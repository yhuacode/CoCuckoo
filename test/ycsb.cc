#include "../CoCuckoo.h"

#define N_KEY 16
#define N_VAL 32
#define HASH_POWER 20

/* type of each query */
enum query_types{
    query_insert=0,
    query_update,
    query_read,
    query_delete,
};

typedef struct __attribute__((__packed__)) {
        char hashed_key[N_KEY];
        enum query_types type;
} query;

typedef CoCuckoo<N_KEY, N_VAL, HASH_POWER> Table;


typedef struct {
    size_t num_insert;
    size_t num_update;
    size_t num_read;
    size_t num_delete;
    double tput;
    double time;
} stat_t;


static size_t key_len;
static size_t val_len;
static size_t num_queries;
static size_t num_records;
static char* loadfile = NULL;
static char* runfile = NULL;
static vector<double> thread_exe_timer;
static size_t success_insert_counter;
static vector<size_t> thread_success_insert_counter;


static void usage(char* binname)
{
    printf("%s -l <load_trace> -r <run_trace> -n <thread_num>\n\n", binname);
    printf("\t-l <load_trace>: an input file for the load phase of YCSB\n");
    printf("\t-r <run_trace>: an input file for the run phase of YCSB\n");
    printf("\t-n <thread_num>: number of threads for the run phase\n");
    printf("\t-h: show usage\n");
}

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


/* init all queries from the ycsb trace file before issuing them */
static query *load_trace(char* filename, int phase)
{
    FILE *input;

    /**
     * Both load_trace and run_trace follow the same format as follows.
     *
     * Field            Type        Number      Description
     * ------------------------------------------------------
     * key_len          size_t      1           size in bytes
     * val_len          size_t      1           size in bytes
     * num_records      size_t      1
     * num_queries      size_t      1
     * queries          query       N           array of query
     *
     * "N" indicates the value of "num_records" or "num_queries" for
     * "run_trace" or "load_trace", respectively.
     */

    input = fopen(filename, "rb");
    if (input == NULL) {
        perror("can not open file");
        perror(filename);
        exit(1);
    }

    int n;
    n = fread(&key_len, sizeof(key_len), 1, input);
    if (n != 1)
       perror("fread error");

    n = fread(&val_len, sizeof(val_len), 1, input);
    if (n != 1)
        perror("fread error");

    n = fread(&num_records, sizeof(num_records), 1, input);
    if (n != 1)
        perror("fread error");

    n = fread(&num_queries, sizeof(num_queries), 1, input);
    if (n != 1)
        perror("fread error");

    printf("trace(%s):\n", filename);
    printf("\tkey_len = %zu\n", key_len);
    printf("\tval_len = %zu\n", val_len);
    printf("\tnum_records = %zu\n", num_records);
    printf("\tnum_queries = %zu\n", num_queries);
    printf("\n");

    size_t cnt = 0;
    if (phase == 0) {
        cnt = num_records;
    } else {
        cnt = num_queries;
    }

    query *queries = (query *)malloc(sizeof(query) * cnt);
    if (queries == NULL) {
        perror("not enough memory to init queries\n");
        exit(-1);
    }

    size_t num_read;
    size_t offset = 0;
    while ( (num_read = fread(queries + sizeof(query) * offset, sizeof(query), cnt-offset, input)) > 0) {
        offset += num_read;
        if (offset >= cnt) break;
    }
    if (offset < cnt) {
        fprintf(stderr, "num_read: %zu\n", offset);
        perror("can not read all queries\n");
        fclose(input);
        exit(-1);
    }

    fclose(input);
    printf("queries_init...done\n");
    return queries;
}


void run_trace(Table *cht, query *queries, int phase)
{
    printf("run_trace...begin\n");

    struct timespec tv_s, tv_e;
    stat_t result;

    result.num_insert = 0;
    result.num_update = 0;
    result.num_read = 0;
    result.num_delete = 0;
    result.tput = 0.0;
    result.time = 0.0;

    size_t cnt = 0;
    if (phase == 0) {
        cnt = num_records;
    }
    else {
        cnt = num_queries;
    }

    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tv_s);
    for (size_t i = 0 ; i < cnt; i++) {
        enum query_types type = queries[i].type;
        char *key = queries[i].hashed_key;

        if (type == query_insert) {
            if (cht->insert(key, i, 0) != 2) {
                success_insert_counter++;
            }
            result.num_insert++;
        } else if (type == query_update) {
            cht->search(key);
            result.num_update++;
        } else if (type == query_read) {
            cht->search(key);
            result.num_read++;
        } else if (type == query_delete) {
            cht->cocuckoo_delete(key, i, 0);
            result.num_delete++;
        }
        else {
            fprintf(stderr, "unknown query type\n");
        }
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tv_e);
    result.time = timeval_diff(&tv_s, &tv_e);

    size_t nops = result.num_insert + result.num_update + result.num_read + result.num_delete;
    result.tput = nops / result.time;

    printf("run %zu #operations in %.2f sec \n",
        nops, result.time);
    printf("#insert = %zu, #update = %zu, #read = %zu, #delete = %zu\n",
        result.num_insert,
        result.num_update,
        result.num_read,
        result.num_delete);
    printf("tput = %.2f\n",  result.tput);
    printf("\n");

    printf("run_trace...done\n");
}


void do_queries(Table *cht, query *queries, int thread_id, size_t begin, size_t count) {

    size_t end = begin + count;
    printf("thread %d: run_trace...begin:%lu end:%lu\n", thread_id, begin, end);

    struct timespec tv_s, tv_e;
    stat_t result;


    result.num_insert = 0;
    result.num_update = 0;
    result.num_read = 0;
    result.num_delete = 0;
    result.tput = 0.0;
    result.time = 0.0;

    try {

        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tv_s);
        for (size_t i = begin; i < end; i++) {
            enum query_types type = queries[i].type;
            char *key = queries[i].hashed_key;
            if (type == query_insert) {
                if (cht->insert(key, i, thread_id) != 2) {
                    thread_success_insert_counter[thread_id]++;
                };
                result.num_insert++;
            } else if (type == query_update) {
                cht->search(key);
                result.num_update++;
            } else if (type == query_read) {
                cht->search(key);
                result.num_read++;
            } else if (type == query_delete) {
                cht->cocuckoo_delete(key, i, thread_id);
                result.num_delete++;
            } else {
                fprintf(stderr, "unknown query type\n");
            }
        }
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tv_e);
        result.time = timeval_diff(&tv_s, &tv_e);
    }
    catch (const char* msg) {
        cerr << cht->get_timestamp() <<  ", Thread-" << thread_id << " get EXCEPTION, " << msg << endl;

        return;
    }

    size_t nops = result.num_insert + result.num_update + result.num_read + result.num_delete;
    result.tput = nops / result.time;

    printf("thread %d: run %zu #operations in %.2f sec \n",
        thread_id, nops, result.time);
    printf("thread %d: #insert = %zu, #update = %zu, #read = %zu, #delete = %zu\n",
        thread_id,
        result.num_insert,
        result.num_update,
        result.num_read,
        result.num_delete);
    printf("thread %d: tput = %.2f\n", thread_id,  result.tput);
    printf("\n");

    printf("thread %d: run_trace...done\n", thread_id);

    thread_exe_timer[thread_id] = result.time;
}



int main(int argc, char **argv)
{
    if (argc <= 1) {
        usage(argv[0]);
        exit(-1);
    }

    int thread_num;

    char ch;
    while ((ch = getopt(argc, argv, "l:r:n:")) != -1) {
        switch (ch) {
            case 'l': loadfile   = optarg; break;
            case 'r': runfile    = optarg; break;
            case 'n':
                thread_num = strtoul(optarg, 0, 10);
                if (thread_num < 0 || thread_num > 16)
                    assert(0);
                break;
            case 'h': usage(argv[0]); exit(0); break;
            default:
                  usage(argv[0]);
                  exit(-1);
        }
    }

    if ( loadfile == NULL || runfile == NULL || thread_num == 0) {
        usage(argv[0]);
        exit(-1);
    }

    thread_exe_timer = vector<double>(thread_num, 0);
    thread_success_insert_counter = vector<size_t>(thread_num, 0);

    size_t keys_per_thread;
    size_t begin;

    Table *cht = new Table(thread_num);

    printf("==================load phase begin=================\n");

    query *loadTraces = load_trace(loadfile, 0);
    run_trace(cht, loadTraces, 0);
    free(loadTraces);

    printf("==================load phase end=================\n");

    printf("==================run phase begin=================\n");

    query *runTraces = load_trace(runfile, 1);

    vector<thread> threads;
    keys_per_thread = num_queries / thread_num;

    begin = 0;
    for (int i = 0; i < thread_num; i++) {
        if(i == thread_num - 1)
            keys_per_thread = num_queries - keys_per_thread * (thread_num - 1);
        threads.emplace_back(do_queries, cht, runTraces, i, begin, keys_per_thread);
        begin += keys_per_thread;
    }
    for (int i = 0; i < thread_num; i++) {
        threads[i].join();
    }

    free(runTraces);
    printf("==================run phase end=================\n");

    printf("key_len: %zu\n", key_len);
    printf("val_len: %zu\n", val_len);
    printf("capacity: %zu\n", cht->get_capacity());

    for (int i = 0; i < thread_num; i++) {
        success_insert_counter += thread_success_insert_counter[i];
    }
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

    return 0;
}
