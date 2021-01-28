#ifndef _UNION_FIND_H_
#define _UNION_FIND_H_

#include <vector>
#include <iostream>
#include <time.h>
#include <sys/time.h>
#include "spinlock.hpp"


struct UFSet{
public:
    UFSet(int N) : id(N, -1) {}

    int find(int idx) {
        if (idx < 0 || idx >= id.size()) {
            std::cerr << "FIND_ERROR, idx: " << idx << ", N: " << id.size() << std::endl;
            throw "FIND_ERROR";
        }

        int parent = id[idx];

        if (parent < 0) {
            return idx;
        }
        else {
            int r = find(parent);
            return id[idx] = r;
        }
    }

    bool connected(int p, int q) {
        return find(p) == find(q);
    }

    int merge(int p, int q) {
        int i = find(p);
        int j = find(q);
        if (i == j) return i;

        if (i < j) {
            return id[j] = i;
        }
        else {
            return id[i] = j;
        }
    }

private:
    std::vector<int> id;
};

#endif // _UNION_FIND_H_
