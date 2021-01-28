#include "../CoCuckoo.h"

#define N_KEY 16
#define N_VAL 16
#define HASH_POWER 20

#define THREAD_NUM 16
#define QUERY_NUM 2000000

typedef CoCuckoo<N_KEY, N_VAL, HASH_POWER> Table;

static size_t success_insert_counter;
static vector<double> thread_exe_timer;
static vector<size_t> thread_success_insert_counter;

static inline double timeval_diff_ns(struct timespec *start,
    struct timespec *end)
{
    double r = (end->tv_sec - start->tv_sec) * 1000000000.0;

    if (end->tv_nsec > start->tv_nsec)
        r += (end->tv_nsec - start->tv_nsec);
    else if (end->tv_nsec < start->tv_nsec)
        r -= (start->tv_nsec - end->tv_nsec);
    return r;
}

static inline double timeval_diff(struct timespec *start,
    struct timespec *end)
{
    return timeval_diff_ns(start, end) / 1000000000.0;
}

void do_insert(Table *cht, int thread_id, size_t begin, size_t count) {
    size_t end = begin + count;
    printf("thread %d: process insert ops...begin:%lu end:%lu\n", thread_id, begin, end);

    struct timespec tv_s, tv_e;
    try {
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tv_s);
        for (size_t i = begin; i < end; i++) {
            if (cht->insert(to_string(i), i, thread_id) != 2) {
                thread_success_insert_counter[thread_id]++;
            };
        }
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tv_e);
        thread_exe_timer[thread_id] = timeval_diff(&tv_s, &tv_e);
    }
    catch (const char* msg) {
        cerr << cht->get_timestamp() <<  ", Thread-" << thread_id << " get EXCEPTION, " << msg << endl;
        return;
    }

    printf("thread %d: run %zu #operations in %.2f sec \n",
        thread_id, count, thread_exe_timer[thread_id]);
    printf("thread %d: tput = %.2f\n\n",
        thread_id, count / thread_exe_timer[thread_id]);
}

void do_search(Table *cht, int thread_id, size_t begin, size_t count) {
    size_t end = begin + count;
    printf("thread %d: process search ops...begin:%lu end:%lu\n", thread_id, begin, end);

    struct timespec tv_s, tv_e;
    try {
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tv_s);
        for (size_t i = begin; i < end; i++) {
            cht->search(to_string(i));
        }
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tv_e);
        thread_exe_timer[thread_id] = timeval_diff(&tv_s, &tv_e);
    }
    catch (const char* msg) {
        cerr << cht->get_timestamp() <<  ", Thread-" << thread_id << " get EXCEPTION, " << msg << endl;
        return;
    }

    printf("thread %d: run %zu #operations in %.2f sec \n",
        thread_id, count, thread_exe_timer[thread_id]);
    printf("thread %d: tput = %.2f\n\n",
        thread_id, count / thread_exe_timer[thread_id]);
}

int main()
{
    int thread_num = THREAD_NUM;

    thread_exe_timer = vector<double>(thread_num, 0);
    thread_success_insert_counter = vector<size_t>(thread_num, 0);

    size_t keys_per_thread;
    size_t begin;

    Table *cht = new Table(thread_num);

    size_t num_queries = QUERY_NUM;
    keys_per_thread = num_queries / thread_num;

    /* testing insertion */
    printf("==================insertions begin=================\n");
    vector<thread> writer_threads;
    begin = 0;
    for (int i = 0; i < thread_num; i++) {
        if(i == thread_num - 1)
            keys_per_thread = num_queries - keys_per_thread * (thread_num - 1);
        writer_threads.emplace_back(do_insert, cht, i, begin, keys_per_thread);
        begin += keys_per_thread;
    }
    for (int i = 0; i < thread_num; i++) {
        writer_threads[i].join();
    }

    for (int i = 0; i < thread_num; i++) {
        success_insert_counter += thread_success_insert_counter[i];
    }
    printf("==================insertions end=================\n");

    printf("capacity: %zu\n", cht->get_capacity());
    printf("Total success_insert_counter: %lu, load ratio: %f\n",
        success_insert_counter, success_insert_counter * 1.0 / (cht->get_capacity()));

    double max_timer = 0, temp;
    for (int i = 0; i < thread_num; i ++) {
        temp = thread_exe_timer[i];
        if (temp > max_timer)
            max_timer = temp;
    }

    printf("execution time: %lf\n", max_timer);
    printf("throuput: %lf\n", num_queries / max_timer);

    /* testing search */
    printf("==================searches begin=================\n");
    vector<thread> reader_threads;
    begin = 0;
    for (int i = 0; i < thread_num; i++) {
        if(i == thread_num - 1)
            keys_per_thread = num_queries - keys_per_thread * (thread_num - 1);
        reader_threads.emplace_back(do_search, cht, i, begin, keys_per_thread);
        begin += keys_per_thread;
    }
    for (int i = 0; i < thread_num; i++) {
        reader_threads[i].join();
    }

    printf("==================searches end=================\n");

    max_timer = 0;
    for (int i = 0; i < thread_num; i ++) {
        temp = thread_exe_timer[i];
        if (temp > max_timer)
            max_timer = temp;
    }

    printf("execution time: %lf\n", max_timer);
    printf("throuput: %lf\n", num_queries / max_timer);

    return 0;
}