#ifndef _COCUCKOO_H_
#define _COCUCKOO_H_

#include <getopt.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sched.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <assert.h>
#include <errno.h>

#include <string>

#include "union_find.h"
#include "MurmurHash3.h"
#include "spinlock.hpp"

#include <time.h>
#include <sys/time.h>
#include <inttypes.h>

#include <cmath>
#include <vector>
#include <array>
#include <set>
#include <queue>
#include <iostream>
#include <algorithm>

#include <chrono>
#include <cstdint>
#include <fstream>
#include <limits>
#include <thread>
#include <utility>
#include <cstdio>
#include <sstream>
#include <atomic>
#include <bitset>
#include <sstream>
#include <random>
#include <map>

using namespace std;


#define TABLE_COUNT 2
#define STASH_SIZE 4

// #define COCUCKOO_DEBUG 1             /* Uncooment this line for thread-local logs */
#define NEIGHBOURS_PRE_ALLOCATE_NUM 10  /* Resource reservation for deletion-enabled CoCuckoo */

/*
 * A switch to enable deletion.
 * When deletion is enabled, the tail latency for write-heavy workloads
 * increases quickly when the load ratio increases. The reason is the
 * time-consuming BFS to split a large subgraph.
 */
// #define DELETE_ENABLE 1              /* Uncomment this line to enable deletion */

/*
 * For the OneEmpty case, this implementation leverages locks to control
 * concurrent writes to the empty bucket. The lock is only required for a
 * key larger than 8 bytes, which cannot be atomically inserted using one
 * instruction. Locks ensure the correctness but at the cost of performance
 * degradation for insert-heavy workloads. One possible optimization is to
 * leverage fine-grained locks, e.g., bucket-grained locks.
 *
 * For further optimizations, there are two lock-free schemes.
 *
 * The first scheme is to insert the first 8-byte of the key atomically (via CAS)
 * and remaining content with non-atomic writes. The CAS for the first 8-byte
 * part of the key ensures the exclusive write access to the empty bucket.
 * If the CAS fails, which indicates another thread is modifying the bucket,
 * redo the insertion. When the bucket is being modified, other threads may
 * observe the intermediate states during assignment, which however does not
 * affect the correctness. Specifically, due to partial keys from intermediate
 * states, search/delete operations are expected to not find the new item,
 * since the assignment is not completed yet. For the kick-out
 * operations in insertion, they may find the resident bucket of a key is
 * not one of the two candidate positions due to the two-step assignment.
 * The threads performing kick-out operations need to spin until the
 * assignment completes.
 *
 * The second scheme is to store pointers (to key-value items) in the table.
 * In 64-bit architecture, the pointer is 8 bytes and can be updated atomically.
 * The actual key-value items are dynamically allocated out of the hash table.
 */
#define ONE_EMPTY_LOCK 1

// Used for key/value larger than 8 bytes
template <size_t N>
struct FixedLengthStr
{
    static const string empty_s;
    char x[N];
    static const size_t size = N;

    FixedLengthStr() : x{'\0'} {}
    FixedLengthStr(const char* s) { memcpy(x, s, N); }
    FixedLengthStr(const string &s) { memcpy(x, s.c_str(), N); }
    FixedLengthStr(const FixedLengthStr &other) { memcpy(x, other.x, N); }
    FixedLengthStr& operator=(const FixedLengthStr &other) {
        memcpy(x, other.x, N);
        return *this;
    }
    bool operator!=(const FixedLengthStr &other) const { return memcmp(x, other.x, N); }
    bool operator==(const FixedLengthStr &other) const { return !(*this!=other); }
    bool operator==(const char* s) const { return !memcmp(x, s, N); }
    bool empty() const { return (*this==empty_s.c_str()); }
    void clear() { memset(x, 0, N); }

};

template <size_t N> const string FixedLengthStr<N>::empty_s(N, '\0');


template <size_t NKEY, size_t NVAL, size_t POWER>
class CoCuckoo {
public:
    vector<fstream> thread_logs;

    typedef uint32_t BucketIndex_t;

    typedef FixedLengthStr<NKEY> Key;
    typedef FixedLengthStr<NVAL> Data;

    CoCuckoo(int t_num) :
      thread_num(t_num), table(TABLE_COUNT, vector<Key>(max_size)), subgraph_num_locks(max_size),
      group_isfull(max_size, -2), pset(max_size), stash(STASH_SIZE), thread_sub_queue_front(thread_num, 0),
      thread_logs(t_num) {

        cout << "Construct CoCuckoo, parameters: key_len = " << key_len << ", val_len = " << val_len
            << ", power = " << power << ", thread_num = " << thread_num << endl;

        // init seeds of hash functions
        random_device dev;
        mt19937 gen(dev());
        uniform_int_distribution<uint32_t> dist(0, UINT32_MAX);
        for (int i = 0; i < 2; i++) {
            seeds[i] = dist(gen);
        }

        node = (Node **)malloc(TABLE_COUNT * sizeof(Node *));
        for (size_t i = 0; i < TABLE_COUNT; i++) {
            node[i] = (Node *)malloc(max_size * sizeof(Node));
            for (size_t j = 0; j < max_size; j++) {
                node_sub(i, j).store(-1);
                node[i][j].neighbours.reserve(NEIGHBOURS_PRE_ALLOCATE_NUM);
            }
        }

        #ifdef COCUCKOO_DEBUG
        for (int i = 0; i < thread_num; i++) {
            stringstream ss;
            ss << "thread-" << i << ".log";
            thread_logs[i].open(ss.str(), fstream::out);
        }
        #endif
    }

    size_t get_capacity() { return max_size * TABLE_COUNT; }

    atomic<int>& node_sub(size_t table, size_t hh) { return node[table][hh].sub; }

    int subgraph_num(size_t table, size_t hh) { return pset.find(node_sub(table, hh)); }

    unsigned long get_timestamp() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return (1000000 * tv.tv_sec + tv.tv_usec);
    }

    void walk_table_for_graph_size() {
        vector<size_t> graph_size(max_size, 0);

        for (size_t i = 0; i < TABLE_COUNT; i++) {
            for (size_t j = 0; j < max_size; j++) {
                int sub = node_sub(i, j);
                if (sub > -1) {
                    sub = subgraph_num(i, j);
                    graph_size[sub]++;
                }
            }
        }

        set<int> graphs;
        fstream fs("graph-size.csv", fstream::out);
        size_t counter = 0, total_nodes = 0;
        for (size_t i = 0; i < graph_size.size(); i++) {
            if (graph_size[i] > 0) {
                fs << i << "," << graph_size[i] << endl;
                total_nodes += graph_size[i];
                counter++;

                graphs.insert(i);
            }
        }

        size_t max_counter = 0, non_max_counter = 0;
        for (set<int>::iterator it = graphs.begin(); it != graphs.end(); it++) {
            if (group_isfull[*it] == 0) {
                max_counter++;
            }
            else if (group_isfull[*it] == -1) {
                non_max_counter++;
            }
            else {
                cerr << "walk_table_for_graph_size error!!! subgraph_num: " << *it << endl;
            }
        }
        cout << "Recorded " << counter << " subgraphs, " << total_nodes << " nodes in total" << endl;
        cout << "maximal graphs counter: " << max_counter << ", non-maximal graph counter: " << non_max_counter << endl;
    }

    void walk_table_for_neighbours_statistics() {
        vector<vector<size_t> > neighbours_num(TABLE_COUNT, vector<size_t> (max_size, 0));

        fstream fs("neighbours_num.csv", fstream::out);
        int non_zero_counter = 0, total_counter = 0, max_neighbour_num = 0;

        map<size_t, size_t> bins;

        for (size_t i = 0; i < TABLE_COUNT; i++) {
            for (size_t j = 0; j < max_size; j++) {
                neighbours_num[i][j] = node[i][j].neighbours.size();

                if (bins.find(neighbours_num[i][j]) == bins.end()) {
                    bins[neighbours_num[i][j]] = 0;
                }
                bins[neighbours_num[i][j]]++;

                if (neighbours_num[i][j] > 0) {
                    non_zero_counter++;
                    total_counter += neighbours_num[i][j];
                    if (neighbours_num[i][j] > max_neighbour_num) {
                        max_neighbour_num = neighbours_num[i][j];
                    }
                }
            }
        }


        cout << "===========================================================" << endl;
        cout << "#non_zero_neighbours: " << non_zero_counter << endl;
        cout << "Total neighbours_num: " << total_counter << endl;
        cout << "Average: " << 1.0 * total_counter / TABLE_COUNT / max_size
            << ", average (non-empty): " << 1.0 * total_counter / non_zero_counter << ", maximal: " << max_neighbour_num << endl;

        size_t total_bin = 0;
        for (map<size_t, size_t>::iterator it = bins.begin(); it != bins.end(); it++) {
            total_bin += it->second;
            cout << "Bin " << it->first << ": " << it->second << "\tproportion: " << it->second * 1.0 / TABLE_COUNT / max_size << endl;
            fs << it->first << "," << it->second << endl;
        }

        cout << "Total bins: " << total_bin << endl;
    }

    void Hash(const Key &k, uint32_t* pa, uint32_t* pb) {
        *pa = defaultHash(k, seeds[0]) % max_size;
        *pb = defaultHash(k, seeds[1]) % max_size;
    }

    int search(const Key &m) {
        uint32_t ha, hb;

        Hash(m, &ha, &hb);

        if (table[0][ha] == m || table[1][hb] == m) return 1;

        for (size_t i = 0; i < stash.size(); i++) {
            if (m == stash[i]) {
                return 1;
            }
        }

        return 0;
    }

    int insert(const Key &m, int task_id, int thread_id) {
        int ha_sub, hb_sub, v_num, ha_sub_dc, hb_sub_dc;
        int isfull_one, isfull_two;
        int temp;
        int g;
        uint32_t ha, hb;
        Hash(m, &ha, &hb);

        if (table[0][ha] == m || table[1][hb] == m) {
            return 0;
        }

        while (1) {

            for (int i = 0; i < STASH_SIZE; i++) {
                if (m == stash[i]) {
                    // find in stash
                    return 0;
                }
            }

            while (1) {
                // Get snapshots for subgraph numbers of two candidate positions.
                // Since the subgraph numbers may change between read and locking,
                // double checking for subgraph number is necessary.

                ha_sub = node_sub(0, ha).load();
                hb_sub = node_sub(1, hb).load();

                if (ha_sub > -1) {
                    ha_sub = subgraph_num(0, ha);
                }
                if (hb_sub > -1) {
                    hb_sub = subgraph_num(1, hb);
                }

                #ifndef DELETE_ENABLE
                // failure prediction: if deletion is disabled, maximal subgraphs would not become non-maximal ones
                if (ha_sub > -1 && hb_sub > -1 && group_isfull[ha_sub] == 0 && group_isfull[hb_sub] == 0 && stash_full) {
                    return 2;
                }
                #endif

                if(ha_sub==hb_sub)
                {
                    // [two_empty, same_non, same_max]
                    if(ha_sub > -1) {
                        // [same_non, same_max]
                        lock(ha_sub);
                    }
                }
                else if(ha_sub<hb_sub)
                {
                    // If two subgraphs are different and non-empty, we first lock
                    // the subgraph with the smaller subgraph number and then the other one.
                    // Locking in order avoids deadlocks in synchronization.
                    if(ha_sub > -1)
                    {
                        // [diff_non_non, diff_non_max, diff_max]
                        lock(ha_sub);
                        lock(hb_sub);
                    }else if(hb_sub > -1) {
                        // [one_empty]
                        #ifdef ONE_EMPTY_LOCK
                        lock(hb_sub);
                        #endif

                        break;
                    }
                }
                else
                {
                    // If two subgraphs are different and non-empty, we first lock
                    // the subgraph with the smaller subgraph number and then the other one.
                    // Locking in order avoids deadlocks in synchronization.
                    if(hb_sub > -1)
                    {
                        // [diff_non_non, diff_non_max, diff_max]
                        lock(hb_sub);
                        lock(ha_sub);
                    }else if(ha_sub > -1) {
                        // [one_empty]
                        #ifdef ONE_EMPTY_LOCK
                        lock(ha_sub);
                        #endif

                        break;
                    }
                }

                ha_sub_dc = node_sub(0, ha);
                hb_sub_dc = node_sub(1, hb);

                if (node_sub(0, ha) > -1) {
                    ha_sub_dc = subgraph_num(0, ha);
                }
                if (node_sub(1, hb) > -1) {
                    hb_sub_dc = subgraph_num(1, hb);
                }

                if (ha_sub_dc == ha_sub && hb_sub_dc == hb_sub) {
                    // lock successfully!
                    break;
                }
                else {
                    // union_set modified by other threads, relock!
                    unlock_2_sub_in_order(ha_sub, hb_sub);
                }
            }

            v_num = judge_v_num(ha_sub, hb_sub);

            // Two_empty/v+2
            if (v_num == 2)
            {
                // Generate a subgraph number.
                temp = find_sub_num(thread_id);

                lock(temp);

                bool ha_sub_succ = false;
                if ((ha_sub_succ = atomic_compare_exchange_strong(&node_sub(0, ha), &ha_sub, temp))
                    && atomic_compare_exchange_strong(&node_sub(1, hb), &hb_sub, temp)) {
                    // Successed
                    table[0][ha] = m;

                    #ifdef DELETE_ENABLE
                    insert_bucket_index(ha, hb);
                    #endif

                    group_isfull[temp] = -1;

                    assert(node_sub(0, ha) == temp && node_sub(1, hb) == temp);
                    unlock(temp);

                    return 0;
                }
                else {
                    // Failed, redo
                    // reset node_sub(0, ha) if necessary
                    if (ha_sub_succ) {
                        node_sub(0, ha) = -1;
                    }

                    unlock(temp);

                    continue;
                }
            }

            // One_empty/v+1
            else if (v_num == 1)
            {
                if (ha_sub == -1)
                {
                    assert(hb_sub != -1);
                    if (!atomic_compare_exchange_strong(&node_sub(0, ha), &ha_sub, hb_sub)) {
                        // Failed , redo
                        #ifdef ONE_EMPTY_LOCK
                        unlock(hb_sub);
                        #endif

                        continue;
                    }
                    //Succeed

                    table[0][ha] = m;

                    #ifdef DELETE_ENABLE
                    insert_bucket_index(ha, hb);
                    #endif

                    #ifdef ONE_EMPTY_LOCK
                    unlock(hb_sub);
                    #endif
                }
                else
                {
                    assert(ha_sub != -1);
                    if (!atomic_compare_exchange_strong(&node_sub(1, hb), &hb_sub, node_sub(0, ha).load())) {
                        #ifdef ONE_EMPTY_LOCK
                        unlock(ha_sub);
                        #endif

                        continue;
                    }

                    table[1][hb] = m;

                    #ifdef DELETE_ENABLE
                    insert_bucket_index(ha, hb);
                    #endif

                    #ifdef ONE_EMPTY_LOCK
                    unlock(ha_sub);
                    #endif
                }

                return 1;
            }

            // Non-empty/v+0
            else
            {
                if (ha_sub == -1 || hb_sub == -1)
                    printf("v+0 error, ha_sub: %d, hb_sub: %d\n", ha_sub, hb_sub);

                isfull_one = subgraph_isfull(0, ha);
                isfull_two = subgraph_isfull(1, hb);


                if ((isfull_one == -2) || (isfull_two == -2))
                {
                    printf("insert: subgraph_isfull error\n");
                    printf("isfull_one: %d, isfull_two: %d, ha: %d, hb: %d, ha_sub: %d, hb_sub: %d\n",
                        isfull_one, isfull_two, ha, hb, ha_sub, hb_sub);
                    printf("table[0][ha]=%s, table[1][hb]=%s\n", table[0][ha].x, table[1][hb].x);
                    unlock_2_sub_in_order(ha_sub, hb_sub);

                    continue;
                }

                // Two maximal subgraphs.
                // If two candidate positions are both in maximal subgraphs, the insertion must be a failure.
                if ((isfull_one == 0) && (isfull_two == 0))
                {
                    unlock_2_sub_in_order(ha_sub, hb_sub);

                    if (!stash_full) {
                        for (int i = 0; i < STASH_SIZE; i++) {
                            if (stash[i].empty()) {
                                stash[i] = m;
                                return 0;
                            }
                        }
                        stash_full = true;
                    }

                    return 2;
                }

                #ifdef DELETE_ENABLE
                insert_bucket_index(ha, hb);
                #endif

                // Two non-maximal subgraphs.
                if ((isfull_one == -1) && (isfull_two == -1))
                {
                    // Same non-maximal subgraphs.
                    if (ha_sub == hb_sub)
                    {
                        g = pset.find(node_sub(0, ha));
                        node_sub(0, ha) = g;
                        node_sub(1, hb) = g;
                        group_isfull[ha_sub] = 0;
                        group_isfull[g] = 0;
                        assert(g == ha_sub || g == hb_sub);

                        if(ha_sub > -1)
                            unlock(ha_sub);

                        kick_out(m, 0, ha, ha_sub, hb_sub, isfull_one, isfull_two, thread_id);

                        return 3;
                    }
                    // Different non-maximal subgraphs.
                    else
                    {
                        if (ha_sub < hb_sub) {
                            kick_out(m, 1, hb, ha_sub, hb_sub, isfull_one, isfull_two, thread_id);
                            node_sub(1, hb) = ha_sub;
                        }
                        else {
                            kick_out(m, 0, ha, ha_sub, hb_sub, isfull_one, isfull_two, thread_id);
                            node_sub(0, ha) = hb_sub;
                        }

                        group_isfull[ha_sub] = -1;
                        group_isfull[hb_sub] = -1;

                        // different non-maximal graphs must be merged with lock
                        pset.merge(ha_sub, hb_sub);

                        unlock_2_sub_in_order(ha_sub, hb_sub);

                        return 4;
                    }
                }


                // One subgraph is maximal and the other is non-maximal.

                // Update graph state with lock.
                group_isfull[ha_sub] = 0;
                group_isfull[hb_sub] = 0;

                unlock_2_sub_in_order(ha_sub, hb_sub);

                // Merge and other processing done without lock.
                pset.merge(ha_sub, hb_sub);

                if (ha_sub < hb_sub) {
                    node_sub(1, hb) = ha_sub;
                }
                else {
                    node_sub(0, ha) = hb_sub;
                }

                if (isfull_one == -1)
                {
                    kick_out(m, 0, ha, ha_sub, hb_sub, isfull_one, isfull_two, thread_id);
                }
                else
                {
                    kick_out(m, 1, hb, ha_sub, hb_sub, isfull_one, isfull_two, thread_id);
                }

                return 5;
            }
        }
    }

    int cocuckoo_delete(const Key &m, int task_id, int thread_id) {
        uint32_t ha, hb;
        int ha_sub, hb_sub, ha_sub_dc, hb_sub_dc;

        Hash(m, &ha, &hb);

        for (int i = 0; i < STASH_SIZE; i++) {
            if (m == stash[i]) {
                // find in stash, delete
                stash[i].clear();
                stash_full = false;

                return 0;
            }
        }


        ha_sub = node_sub(0, ha);
        hb_sub = node_sub(1, hb);

        if (ha_sub > -1) {
            ha_sub = subgraph_num(0, ha);
        }
        if (hb_sub > -1) {
            hb_sub = subgraph_num(1, hb);
        }

        if(ha_sub==hb_sub) {
            if(ha_sub > -1) {
                lock(ha_sub);

                ha_sub_dc = node_sub(0, ha);
                hb_sub_dc = node_sub(1, hb);

                if (node_sub(0, ha) > -1) {
                    ha_sub_dc = subgraph_num(0, ha);
                }
                if (node_sub(1, hb) > -1) {
                    hb_sub_dc = subgraph_num(1, hb);
                }

                if (ha_sub_dc == ha_sub && hb_sub_dc == hb_sub) {
                    // lock successfully!
                    if (table[0][ha] == m) {
                        table[0][ha].clear();
                        delete_bucket_index(ha, hb);

                        if (node[0][ha].neighbours.empty()) {
                            node_sub(0, ha) = -1;
                            // group_isfull[ha_sub] = -2;
                        }
                        else {
                            group_isfull[ha_sub] = -1;
                        }

                        if (node[1][hb].neighbours.empty()) {
                            node_sub(1, hb) = -1;
                        }
                        else if (!node[0][ha].neighbours.empty()) {
                            vector<pair<int, BucketIndex_t>> visited_buckets;
                            vector<vector<bool>> visited(TABLE_COUNT, vector<bool>(max_size, false));
                            bool is_subgraph_full = BFS_traverse_subgraph(1, hb, visited_buckets, visited);

                            #ifdef COCUCKOO_DEBUG
                            thread_logs[thread_id] << "  BFS_traverse subgraph-" << ha_sub << ": visited_size = " << visited_buckets.size() << endl;
                            #endif

                            if (!visited[0][ha]) {
                                // A new subgraph is generated
                                int sub_id = find_sub_num(thread_id);
                                group_isfull[sub_id] = (is_subgraph_full ? 0 : -1);

                                for (size_t i = 0; i < visited_buckets.size(); i++) {
                                    node_sub(visited_buckets[i].first, visited_buckets[i].second) = sub_id;
                                }
                            }
                        }
                    }
                    else if (table[1][hb] == m) {
                        table[1][hb].clear();
                        delete_bucket_index(ha, hb);

                        if (node[1][hb].neighbours.empty()) {
                            node_sub(1, hb) = -1;
                            // group_isfull[hb_sub] = -2;
                        }
                        else {
                            group_isfull[hb_sub] = -1;
                        }

                        if (node[0][ha].neighbours.empty()) {
                            node_sub(0, ha) = -1;
                        }
                        else if (!node[1][hb].neighbours.empty()){
                            vector<pair<int, BucketIndex_t>> visited_buckets;
                            vector<vector<bool>> visited(TABLE_COUNT, vector<bool>(max_size, false));
                            bool is_subgraph_full = BFS_traverse_subgraph(0, ha, visited_buckets, visited);

                            #ifdef COCUCKOO_DEBUG
                            thread_logs[thread_id] << "  BFS_traverse subgraph-" << ha_sub << ": visited_size = " << visited_buckets.size() << endl;
                            #endif

                            if (!visited[1][hb]) {
                                // A new subgraph is generated
                                int sub_id = find_sub_num(thread_id);
                                group_isfull[sub_id] = (is_subgraph_full ? 0 : -1);

                                for (size_t i = 0; i < visited_buckets.size(); i++) {
                                    node_sub(visited_buckets[i].first, visited_buckets[i].second) = sub_id;
                                }
                            }
                        }
                    }
                    // Not found. The reason can be
                    // 1. Item not exists;
                    // 2. Other threads may delete the item;
                    // Finish here.
                }

                unlock(ha_sub);
                return 0;
            }
        }

        return 0;
    }

    CoCuckoo(const CoCuckoo&) = delete;

    ~CoCuckoo() {
        #ifdef COCUCKOO_DEBUG
        for (int i = 0; i < thread_num; i++) {
            thread_logs[i].close();
        }
        #endif

        for (size_t i = 0; i < TABLE_COUNT; i++) {
            free(node[i]);
            node[i] = NULL;
        }
        free(node);
        node = NULL;
    }


private:
    struct Node {
        atomic<int> sub;
        vector<BucketIndex_t> neighbours;
        Node(int v) : sub(v) {}
		void add_neighbour(BucketIndex_t bucket_index) {
			neighbours.push_back(bucket_index);
		}
		void remove_neighbour(BucketIndex_t bucket_index) {
			auto it = find(neighbours.begin(), neighbours.end(), bucket_index);
			if (it != neighbours.end()) {
				neighbours.erase(it);
			}
		}
    };



    size_t key_len = NKEY;
    size_t val_len = NVAL;
    size_t power = POWER;

    int thread_num;

    int max_size = 1 << power;  // number of buckets in one sub_table

    vector<vector<Key> > table;
	Node** node;
    vector<spinlock> subgraph_num_locks;
    vector<int> group_isfull;   // -2: uninitialized, -1: non-maximal, 0: maximal

    UFSet pset;

    bool stash_full = false;
    vector<Key> stash;
    vector<size_t> thread_sub_queue_front;

    uint32_t seeds[2];


    void lock(uint32_t ha) {
        subgraph_num_locks[ha].lock();
    }

    void unlock(uint32_t ha) {
        subgraph_num_locks[ha].unlock();
    }

    uint32_t defaultHash(const Key &k, uint32_t seed) {
        uint32_t h;
        MurmurHash3_x86_32(k.x, k.size, seed, &h);
        return h;
    }


    int judge_v_num(int ha_sub, int hb_sub) {
        if ((ha_sub == -1) && (hb_sub == -1))
        {   // Two_empty/v+2
            return 2;
        }
        if ((ha_sub != -1) && (hb_sub != -1))
        {   // Non empty/v+0
            return 0;
        }
        // One_empty/v+1
        return 1;
    }

    // Check if the subgraph number of specific bucket is maximal.
    // -1: non-maximal, 0: maximal.
    int subgraph_isfull(int t, int h) {
        int subnum;

        subnum = node_sub(t, h);
        assert(subnum > -1);
        int setnum = pset.find(subnum);

        if (group_isfull[setnum] == -2)
            printf("subgraph_isfull, t: %d, h: %d, subnum: %d, setnum: %d\n", t, h, subnum, setnum);

        return group_isfull[setnum];
    }

    void kick_out(const Key &m, int table_num, uint32_t ha, int ha_sub, int hb_sub, int isfull_one, int isfull_two, int thread_id) {
        uint32_t i;
        uint32_t curr_hash;
        uint32_t next_hash;
        Key temp, num;
        uint32_t a, b;
        int kick_count;

        kick_count = 0;

        i = table_num;
        curr_hash = ha;
        num = m;

        while (!table[i][curr_hash].empty())
        {
            kick_count++;
            temp = table[i][curr_hash];
            table[i][curr_hash] = num;

            Hash(temp, &a, &b);
            num = temp;
            if (curr_hash != a)
            {
                next_hash = a;
            }
            else
            {
                next_hash = b;
            }

            i = !i;
            curr_hash = next_hash;

            #if defined(COCUCKOO_DEBUG) && defined(DELETE_ENABLE)
            if (i == table_num && ha == curr_hash)
            {
                vector<pair<int, BucketIndex_t>> visited_buckets;
                vector<vector<bool>> visited(TABLE_COUNT, vector<bool>(max_size, false));
                bool is_subgraph_full = BFS_traverse_subgraph(0, ha, visited_buckets, visited);
                thread_logs[thread_id] << "[check] kick_counter " << kick_count
                    << " BFS_traverse subgraph-" << ha_sub << ": visited_size = " << visited_buckets.size() << endl;

                assert(i != table_num && ha != curr_hash);
            }
            #endif

            #ifdef COCUCKOO_DEBUG
            if (kick_count > 10000)
            {
                #ifdef DELETE_ENABLE
                vector<pair<int, BucketIndex_t>> visited_buckets;
                vector<vector<bool>> visited(TABLE_COUNT, vector<bool>(max_size, false));
                bool is_subgraph_full = BFS_traverse_subgraph(0, ha, visited_buckets, visited);
                thread_logs[thread_id] << "[check] BFS_traverse subgraph-" << ha_sub << ": visited_size = " << visited_buckets.size() << endl;
                #endif

                thread_logs[thread_id] << "kick_count over 10000! table: " << table_num << ", ha: " << ha << ", ha_sub: " << ha_sub
                    << ", hb_sub: " << hb_sub << ", subgraph_isfull: " << subgraph_isfull(table_num, ha)
                    << ", isfull_one: " << isfull_one << ", isfull_two: " << isfull_two
                    << ", actual isfull_one: " << group_isfull[pset.find(ha_sub)]
                    << ", isfull_two: " << group_isfull[pset.find(hb_sub)] << endl;
                throw "kick_count over 10000";
            }
            #endif
        }

        table[i][curr_hash] = num;
    }

    int find_sub_num(int thread_id) {
        int temp = thread_sub_queue_front[thread_id];
        int i = thread_sub_queue_front[thread_id] * thread_num + thread_id;

        thread_sub_queue_front[thread_id] = (thread_sub_queue_front[thread_id] + 1) % (max_size / thread_num);

        if (temp > thread_sub_queue_front[thread_id])
        {
            printf("[Thread-%d]sub_num rewind, old: %d, new: %lu\n",
                thread_id, temp, thread_sub_queue_front[thread_id]);
        }

        return i;
    }

    void unlock_2_sub_in_order(int ha_sub, int hb_sub) {
        if(ha_sub==hb_sub)
        {
            if(ha_sub > -1)
                unlock(ha_sub);
        }
        else if(ha_sub<hb_sub)
        {
            if(ha_sub > -1)
            {
                unlock(ha_sub);
                unlock(hb_sub);
            }else if(hb_sub > -1)
                unlock(hb_sub);
        }
        else
        {
            if(hb_sub > -1)
            {
                unlock(hb_sub);
                unlock(ha_sub);
            }else if(ha_sub > -1)
                unlock(ha_sub);
        }
    }

    void insert_bucket_index(BucketIndex_t ha, BucketIndex_t hb) {
        node[0][ha].add_neighbour(hb);
        node[1][hb].add_neighbour(ha);
    }

    void delete_bucket_index(BucketIndex_t ha, BucketIndex_t hb) {
        node[0][ha].remove_neighbour(hb);
        node[1][hb].remove_neighbour(ha);
    }

    bool BFS_traverse_subgraph(int table_num, BucketIndex_t root, vector<pair<int, BucketIndex_t> > &visited_buckets, vector<vector<bool> > &visited) {
        bool is_subgraph_full = true;

        visited_buckets.clear();
        visited_buckets.push_back(make_pair(table_num, root));

        for (size_t offset = 0; offset != visited_buckets.size(); offset++) {
            int table_id = visited_buckets[offset].first;
            BucketIndex_t bucket_index = visited_buckets[offset].second;

            visited[table_id][bucket_index] = true;

            if (is_subgraph_full && table[table_id][bucket_index].empty()) {
                is_subgraph_full = false;
            }

            vector<BucketIndex_t> &neighbours = node[table_id][bucket_index].neighbours;
            if (neighbours.size() > 0) {
                // Linked nodes exist.
                for (size_t i = 0; i < neighbours.size(); i++) {
                    if (!visited[!table_id][neighbours[i]]) {
                        visited_buckets.push_back(make_pair(!table_id, neighbours[i]));
                    }
                }
            }
        }

        return is_subgraph_full;
    }

};


#endif  // _COCUCKOO_H_
