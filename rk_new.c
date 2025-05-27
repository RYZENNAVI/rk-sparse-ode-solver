/**
 * Programm uses explicit and embedded runge-kutta methods for calculating solutions to TASEP sparse linear ODE systems
 * A example usage of explicit RK4 for a chain of 4 sites(n = 4): ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e4.mat
 * 
 * INPUT:   -sparse matrix file in:
 *          crs format (genmat -> gencrs)               example: test_n_4.crs
 *          diagonal format (genmat -> gendia)          example: test_n_4_new.crs
 *          new crs format (genmat -> gencrs_new)       example: test_n_4_dia.mat
 * 
 *          -initial vector file (geninit)              example: test_n_4.init
 * 
 *          -Butcher tableau in:
 *          MAT format (for all exp. and emb. methods)  example: butcher_e4.mat(RK4), butcher_e5.mat(E5), butcher_b45.mat(RKF45), butcher_b23(BS23)
 *          CRS format (only for RK4 + CRS matrix)      example: butcher_c.crs (RK4)
 * 
 * OUTPUT:
 *          -vector calculated with runge kutta-method after specified start time, # of steps and time period
 * 
 * Basic USAGE:     
 *          -original RK4 (no Butcher tableau required)
 *          example: ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init
 * 
 *          -explicit runge-kutta methods with CRS Butcher tableau(only RK4 + CRS or new CRS format)
 *          examples:
 *          RK4(only CRS): ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_c.crs
 * 
 *          -explicit runge-kutta methods with MAT Butcher tableau(RK4 and E5)
 *          examples:
 *          RK4: ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e4.mat
 *          E5:  ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e5.mat
 * 
 *          -richardson extrapolation (RE RK4 and RE E5)
 *          active RE and tolerance: [-r -l tolerance]
 *          examples:
 *          RE RK4:  ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e4.mat -r -l 0.00001
 *          RE E5:   ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e5.mat -r -l 0.00001

 *          -embedded runge-kutta methods (RKF45 and BS23)
 *          tolerance: [ -l tolerance]
 *          examples:
 *          RKF45: ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_b45.mat -l 0.1
 *          BS23:  ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_b23.mat -l 0.1
 * 
 *          -all sparse matrix storage formats (CRS, new CRS and diagonal)
 *          formats: [-m matrix format]
 *          examples:
 *          CRS:        ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e4.mat
 *          new CRS:    ./rk_new -m test_n_4_new.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e4.mat
 *          diagonal:   ./rk_new -m test_n_4_dia.mat -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e4.mat
 * 
 * Advanced USAGE:   
 *          -multiple threads:
 *          manual mode: Set available CPU_ID in global array cpu_ids manually 
 *          auto mode:   [ -u Number of threads required ]
 *          When auto mode inactive, then the first four threads of the CPU are called by default (IDs: 0,1,2,3).
 *          example: ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e4.mat -u 16
 * 
 *          -Loop tiling (only for explicit methods with CRS and new CRS format):
 *          version 1:   [ -V 1 -s tile size ]
 *          version 2:   [ -V 2 -s tile size ]
 *          example: ./rk_new -m test_n_4.crs -o test_n_4.out -T 5 -t 0 -n 1000 -v test_n_4.init -B butcher_e4.mat -V 1 -s 4
 * 
 *         
 * 
 **/



#include <getopt.h>
#include <string.h> 
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <assert.h> 
#include <omp.h> 
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <sched.h>

#define a_max_embedded 2.5
#define a_min_embedded 0.05
#define a_max_richardson 4
#define a_min_richardson 0.05
int chunk_size = 120;

// test in i7-13700k, which has 24 threads
// P-threads id (with hyperthreading)   : 0-15 
// E-threads id                         : 16-23

// first 4 threads
int cpu_ids[] = {0,1,2,3};
// all P-threads
//int cpu_ids[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
// all E-threads
//int cpu_ids[] = {16,17,18,19,20,21,22,23};
// all threads
//int cpu_ids[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};

unsigned int min(unsigned int a, unsigned int b)
{
    return (a>b)? b:a;
}
unsigned int max(unsigned int a, unsigned int b)
{
    return (a<b)? b:a;
}

void set_thread_affinity(int *cpu_ids, int num_ids) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);

    // Set affinity to allow threads to run on all CPU IDs provided
    for (int i = 0; i < num_ids; i++) {
        CPU_SET(cpu_ids[i], &cpuset);
    }

    pthread_t current_thread = pthread_self();
    if (pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset) != 0) {
        perror("pthread_setaffinity_np failed");
    }
}

void print_affinity(int threads) {
    cpu_set_t cpuset;
    pthread_t current_thread = pthread_self();

    if (pthread_getaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset) == 0) {
        printf("cpu_ids: ");
        for (unsigned int i = 0; i < min(CPU_SETSIZE,threads); ++i) {
            if (CPU_ISSET(i, &cpuset)) {
                printf(" %d", i);
            }
        }
        printf("\n");
    } else {
        perror("Error getting affinity");
    }
}

void print_affinity_details() {
    cpu_set_t cpuset;
    pthread_t current_thread = pthread_self();

    if (pthread_getaffinity_np(current_thread, sizeof(cpu_set_t), &cpuset) == 0) {
        printf("Thread %d running on CPUs:", omp_get_thread_num());
        for (int i = 0; i < CPU_SETSIZE; ++i) {
            if (CPU_ISSET(i, &cpuset)) {
                printf(" %d", i);
            }
        }
        printf("\n");
    } else {
        perror("Error getting affinity");
    }
}

static void usage(const char *progname)
{
  printf("Usage: %s [-m matrix file] [-B Butcher tableau file][-v start vector file] [-o output file] [-t start time] [-T  endtime] [-n steps] [-h]\n", progname);
  printf("[-a Create file of result vectors during calculation of solution]\n");
  printf("  the given number is the amount of time steps between snapshots\n");
}

struct timeval start, stop;

static inline uint64_t log2_u64(uint64_t n)
{
  int result = 0;
  while (n >>= 1) ++result;
  return result;
}

int main(int argc, char *argv[])
{   
    int opt;
    char *matrix = NULL;
    int matrix_found = 0;
    int butcher_crs_matrix_found = 0;
    char *butcher_matrix = NULL;
    int butcher_matrix_found = 0;
    char matrix_type, methode_type;
    int CRS_found = 0;
    int CRS_new_found = 0;
    int Diagonal_found = 0;
    int explicit_found = 0;
    int embedded_found = 0;
    int tol_found = 0;
    int r=0;

    char *init_vector = NULL;
    int iv_found = 0;
    char *output_filename = NULL;
    int of_found = 0;
    double end_time;
    int t_found = 0, tile_size_found = 0, richardson_found = 0, tile_version_found = 0, threads_found = 0;
    uint64_t steps, counter_iter = 0, tile_size = 0, tile_version = 0;
    int threads = sizeof(cpu_ids) / sizeof(cpu_ids[0]);

    int n_found = 0;
    double h, h_new, t0, checksum, counter = 0, tol = 0;
    long double stop_condition;
    int t0_found = 0;
    uint64_t val_size = 0, butcher_crs_val_size = 0, dim = 0, butcher_crs_dim = 0, butcher_stage = 0,
    steps_at_stop __attribute__((unused))= 0, *col = NULL, *rowPtr = NULL, *butcher_crs_col = NULL, *butcher_crs_rowPtr = NULL, 
    ord = 0, *addmax = NULL, **index = NULL;
    double *x = NULL, *x1 = NULL, *x2 = NULL, *swap_helper1 = NULL, * x_half = NULL, * x_bs = NULL, 
    *x_last = NULL, *val = NULL, *butcher_crs_val = NULL, *k1 = NULL, *k2 = NULL, *k3 = NULL, *k4 = NULL, 
    *iter_tiling = NULL, **T_tiling = NULL,**ki = NULL, **ki_half = NULL, **diag = NULL, 
    **A = NULL, *b = NULL, *bs = NULL, *bbs = NULL, *err_absolute = NULL, *y_scale = NULL;
    FILE *F;
    long int beg_of_val = 0, butcher_beg_of_val = 0;
    int display __attribute__((unused))= 0, early_stop = 0;
    char *t0_str = NULL, *end_time_str = NULL, *steps_str = NULL, 
    *tol_str = NULL, *tile_size_str = NULL, *tile_version_str = NULL, *thread_str = NULL;
    double secs = 0.;

    while ((opt = getopt(argc, argv, "m:B:C:l:s:V:ru:v:o:t:T:n:h?")) != -1) 
    {
        switch(opt)
        {
        case 'm':
            matrix = strdup(optarg);
            matrix_found = 1;
            break;
        case 'B': // Butcher tableau for explicit and embedded method
            butcher_matrix = strdup(optarg);
            butcher_matrix_found = 1;
            break;
        case 'l':
            tol_str = strdup(optarg);
            tol = atof(tol_str);
            tol_found = 1;
            break;
        case 's':
            tile_size_str = strdup(optarg);
            tile_size = atoi(tile_size_str);
            tile_size_found = 1;
            break;
        case 'V':
            tile_version_str = strdup(optarg);
            tile_version = atoi(tile_version_str);
            tile_version_found = 1;
            break;
        case 'r':
            richardson_found = 1;
            break;
        case 'u':
            thread_str = strdup(optarg);
            threads = atoi(thread_str);
            threads_found = 1;
            break;
        case 'v':
            init_vector = strdup(optarg);
            iv_found = 1;
            break;
        case 'o':
            output_filename = strdup(optarg);
            of_found = 1;
            break;
        case 'T':
            end_time_str = strdup(optarg);
            end_time = atoi(end_time_str);
            t_found = 1;
            free(end_time_str);
            break;
        case 'n':
            steps_str = strdup(optarg);
            steps = atoi(steps_str);
            n_found = 1;
            free(steps_str);
            break;
        case 't':
            t0_str = strdup(optarg);
            t0 = atoi(t0_str);
            t0_found = 1;
            free(t0_str);
            break;
        default:
            usage(argv[0]);
            exit(-1);
        }
    }

    if(matrix_found){        
        if(fopen(matrix, "r") == NULL){   
            printf("ERROR: Can't open matrix.crs file: %s\n", matrix);
            exit(-1);
        }
    }
    else
    {
        printf("No matrix file found. Aborting.\n");
        exit(-1);
    }
    if(butcher_matrix_found){// .mat
        if(fopen(butcher_matrix, "r") == NULL){   
            printf("ERROR: Can't open .mat file for butcher tableau: %s\n", butcher_matrix);
            exit(-1);
        }
    }
    else if(butcher_crs_matrix_found){// .crs
        if(fopen(butcher_matrix, "r") == NULL){   
            printf("ERROR: Can't open .crs file for butcher tableau: %s\n", butcher_matrix);
            exit(-1);
        }
    }
    else
    {
        printf("No butcher matrix file found. so execute the rk4.\n");
        //exit(-1);
    }

    // crs matrix analyse
    F = fopen(matrix, "r");
    fscanf(F, "\n%c" PRIu64 , &matrix_type);
    fscanf(F, "%" PRIu64 "\n%" PRIu64, &val_size, &dim);
    if(matrix_type == 'c')
    {
        CRS_found = 1;
        printf("Matrix is in CRS format.\n");
    }
    if(matrix_type == 'n')
    {
        CRS_new_found = 1;
        printf("Matrix is in new CRS format.\n");
    }
    if(matrix_type == 'd')
    {
        Diagonal_found = 1;
        r = log(dim)/log(2);
        printf("Matrix is in Diagonal format with %d rows.\n",r+1);
    }
    printf("Matrix has %lu entries and dimension %lu.\n", val_size, dim);
    beg_of_val = ftell(F);
    fclose(F);

    // butcher tableau analyse
    if(butcher_matrix_found) // .mat butcher tableau
    {        
        F = fopen(butcher_matrix, "r");
        //fscanf(F, "%" PRIu64 "\n%c" PRIu64 " %d" PRIu64, &butcher_stage, &methode_type, &ord);
        fscanf(F, "%" PRIu64 "\n%c" PRIu64 , &butcher_stage, &methode_type);

        printf("butcher tableau has dimension %lu.\n", butcher_stage);
        //butcher_beg_of_val = ftell(F);
        //fclose(F);

        //printf("%c\n",methode_type);
        if(methode_type == 'e')
        {
            fscanf(F, "\n%" PRIu64 , &ord);
            butcher_beg_of_val = ftell(F);
            fclose(F);

            explicit_found = 1;
            printf("order of method: %" PRIu64 "\n",ord);
        }
        if(methode_type == 'b')
        {
            fscanf(F, "\n%" PRIu64 , &ord);
            butcher_beg_of_val = ftell(F);
            fclose(F);

            embedded_found = 1;
            printf("order of method: %" PRIu64 "\n",ord);
        }
        if(methode_type == 'c')
        {
            butcher_crs_dim = butcher_stage;
            //printf("crs!!!!!!!!!!!!!!!!!!\n");
            fscanf(F, "\n%" PRIu64 , &ord);
            fscanf(F, "\n%" PRIu64 , &butcher_crs_val_size);

            //fscanf(F, "\n%" PRIu64 "\n%" PRIu64 , &crs_ord, &butcher_crs_val_size);
            butcher_beg_of_val = ftell(F);
            fclose(F);
            explicit_found = 1;
            butcher_crs_matrix_found = 1;
            printf("order of method: %" PRIu64 "\n",ord);
            printf("butcher tableau has %lu entries and dimension %lu.\n", butcher_crs_val_size, butcher_crs_dim);
            printf("butcher tableau is crs type\n");
            printf("butcher tableau has order %lu \n",ord);

        }
        iter_tiling = (double *) malloc(((size_t)dim) * sizeof(double));
        if(iter_tiling == NULL)
        {
            fprintf(stdout, "ERROR: Not enough memory to allocate array iter_tiling.\n");
            exit(-1);
        }
        for(unsigned int i=0;i<dim;i++)
        {
            iter_tiling[i] = 0;
        }

        
        T_tiling = (double **)malloc((size_t) (butcher_stage) * sizeof(double*));
        for(unsigned int i = 0; i < butcher_stage; i++){
            T_tiling[i] = (double*)malloc((dim) * sizeof(double));
        }
        for(unsigned int i=0;i<butcher_stage;i++)
        {
            for(unsigned int j=0;j<dim;j++)
            {
                T_tiling[i][j] = 0;
            }
        }

        if(T_tiling == NULL)
        {
            fprintf(stdout, "ERROR: Not enough memory to allocate array T_tiling.\n");
            exit(-1);
        }


        ki = (double **)malloc((size_t) (butcher_stage) * sizeof(double*));
        for(unsigned int i = 0; i < butcher_stage; i++){
            ki[i] = (double*)malloc((dim) * sizeof(double));
        }

        if(ki == NULL)
        {
            fprintf(stdout, "ERROR: Not enough memory to allocate array ki.\n");
            exit(-1);
        }

        ki_half = (double **)malloc((size_t) (butcher_stage) * sizeof(double*));
        for(unsigned int i = 0; i < butcher_stage; i++){
            ki_half[i] = (double*)malloc((dim) * sizeof(double));
        }

        if(ki_half == NULL)
        {
            fprintf(stdout, "ERROR: Not enough memory to allocate array ki_half.\n");
            exit(-1);
        }


        A = (double **)malloc((size_t) (butcher_stage) * sizeof(double*));
        for(unsigned int i = 0; i < butcher_stage; i++){
            A[i] = (double*)malloc((butcher_stage) * sizeof(double));
        }

        if(A == NULL)
        {
            fprintf(stdout, "ERROR: Not enough memory to allocate matrix A.\n");
            exit(-1);
        }
        b = (double *) malloc(((size_t)butcher_stage + 1) * sizeof(double));
        if(b == NULL)
        {
            fprintf(stdout, "ERROR: Not enough memory to allocate vector b.\n");
            exit(-1);
        }

        err_absolute = (double *) malloc(((size_t)dim+1) * sizeof(double));
        if(err_absolute == NULL)
        {
            fprintf(stdout, "ERROR: Not enough memory to allocate err_absolute array.\n");
            exit(-1);
        }
        y_scale = (double *) malloc(((size_t)dim+1) * sizeof(double));
        if(y_scale == NULL)
        {
            fprintf(stdout, "ERROR: Not enough memory to allocate y_scale array.\n");
            exit(-1);
        }

        if(embedded_found==1)
        {
            bs = (double *) malloc(((size_t)butcher_stage + 1) * sizeof(double));
            if(bs == NULL)
            {
                fprintf(stdout, "ERROR: Not enough memory to allocate vector bs.\n");
                exit(-1);
            }

            bbs = (double *) malloc(((size_t)butcher_stage+1) * sizeof(double));
            if(bbs == NULL)
            {
                fprintf(stdout, "ERROR: Not enough memory to allocate bbs array.\n");
                exit(-1);
            }

            x_bs = (double *) malloc(((size_t)dim+1) * sizeof(double));
            if(x_bs == NULL)
            {
                fprintf(stdout, "ERROR: Not enough memory to allocate x_bs array.\n");
                exit(-1);
            }
        }
        if(butcher_crs_matrix_found == 1)
        {
            butcher_crs_val = (double *) malloc(((size_t)butcher_crs_val_size + 1) * sizeof(double));
            if(butcher_crs_val == NULL)
            {
                fprintf(stdout, "ERROR: Not enough memory to allocate value array.\n");
                exit(-1);
            }

            butcher_crs_col = (uint64_t *) malloc(((size_t)butcher_crs_val_size + 1) * sizeof(uint64_t));
            if(butcher_crs_col == NULL)
            {
                fprintf(stdout, "ERROR: Not enough memory to allocate column array.\n");
                exit(-1);
            }

            butcher_crs_rowPtr = (uint64_t *) malloc(((size_t)butcher_crs_val_size + 1) * sizeof(uint64_t));
            if(butcher_crs_rowPtr == NULL)
            {
                fprintf(stdout, "ERROR: Not enough memory to allocate rowPtr array.\n");
                exit(-1);
            }
        }

    }

    // crs matrix analyse
    val = (double *) malloc(((size_t)val_size + 1) * sizeof(double));
    if(val == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate value array.\n");
        exit(-1);
    }

    col = (uint64_t *) malloc(((size_t)val_size + 1) * sizeof(uint64_t));
    if(col == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate column array.\n");
        exit(-1);
    }

    rowPtr = (uint64_t *) malloc(((size_t)dim + 1) * sizeof(uint64_t));
    if(rowPtr == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate rowPtr array.\n");
        exit(-1);
    }

    addmax = (uint64_t *) malloc((dim+1) * sizeof(long unsigned));
    if(addmax == NULL)
    {
        fprintf(stdout, "ERROR: Not enough Memory to allocate addmax.\n");
        exit(-1);
    }

    index = (long unsigned **)malloc((size_t) (dim+1) * sizeof(long unsigned*));
    for(unsigned int i = 0; i < dim+1; i++){
        index[i] = (long unsigned*)malloc((dim+1) * sizeof(long unsigned));
    }
    if(index == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate index array.\n");
        exit(-1);
    }

    diag = (double **)malloc((size_t) (r+1) * sizeof(double*));
    for(int i = 0; i < r+1; i++){
        diag[i] = (double*)malloc((dim+1) * sizeof(double));
    }
    if(diag == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate diag array.\n");
        exit(-1);
    }

    x = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(x == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate x array.\n");
        exit(-1);
    }

    x_half = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(x_half == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate x_half array.\n");
        exit(-1);
    }

    x_last = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(x_last == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate x_last array.\n");
        exit(-1);
    }

    x1 = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(x1 == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate x1 array.\n");
        exit(-1);
    }

    x2 = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(x2 == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate x2 array.\n");
        exit(-1);
    }

    swap_helper1 = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(swap_helper1 == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate swap_helper1 array.\n");
        exit(-1);
    }

    k1 = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(k1 == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate k1 array.\n");
        exit(-1);
    }

    k2 = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(k2 == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate k2 array.\n");
        exit(-1);
    }

    k3 = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(k3 == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate k3 array.\n");
        exit(-1);
    }

    k4 = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(k4 == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate k4 array.\n");
        exit(-1);
    }

    if(of_found){
        if(fopen(output_filename, "w") == NULL){
            printf("ERROR: Can't open output file %s. Aborting!\n", output_filename);
            exit(-1);
        }
    }
    else
    {
        printf("No output file found. Final vector will be displayed in Terminal\n");
        display = 1;
    }

    if(!t_found){        
        printf("No end time found. Please enter: ");
        scanf("%lf", &end_time);
    }

    if(!t0_found){
        printf("No start time found. Please enter: ");
        scanf("%lf", &t0);
    }

    if(!n_found){
        printf("No step number found. Please enter: ");
        scanf("%" PRIu64, &steps);
    }

    if(iv_found){        
        if(fopen(init_vector, "r") == NULL){   
            printf("ERR_absolute: Can't open start vector file: %s\n", init_vector);
            exit(-1);
        }
    }
    else
    {
        printf("No initial vector found. Aborting.\n");
        exit(-1);
    }

    F = fopen(init_vector, "r");
    if((F = fopen(init_vector, "r")) == NULL)
    {
      fprintf(stdout, "ERROR: Can't open init_vec file '%s'.\n", init_vector);
      exit(-1);
    }

    unsigned int i = 0;
    while(!feof(F)){
        fscanf(F, "%le", &x[i]);
        if(richardson_found)
        {
            x_half[i] = x[i];
            x_last[i] = x[i];
        }
        if(embedded_found)
        {
            x_bs[i] = x[i];
            x_last[i] = x[i];
        }
        i++;
    }
    fclose(F);
    
    if(matrix_type == 'c' || matrix_type == 'n')
    {
        // val, col, row for the Matrix
        F = fopen(matrix, "r");
        fseek(F, beg_of_val, SEEK_SET);
        
        for(uint64_t i = 0; i < val_size; i++)
            fscanf(F, "%lf" ,&val[i]);
        
        for(uint64_t i = 0; i < val_size; i++)
            fscanf(F, "%" PRIu64 ,&col[i]);

        for(uint64_t i = 0; i < dim + 1; i++)
            fscanf(F, "%" PRIu64 ,&rowPtr[i]);

        if(matrix_type == 'n')
        {
            for(uint64_t i = 0; i < dim; i++)
                fscanf(F, "%" PRIu64 ,&addmax[i]);

            for(uint64_t i = 0; i < dim; i++)
            {
                for(uint64_t j = 0; j < addmax[i]; j++)
                {
                    fscanf(F, "%" PRIu64 ,&index[i][j]);
                }

            }
        }
    }
    if(matrix_type == 'd')
    {
        F = fopen(matrix, "r");
        fseek(F, beg_of_val, SEEK_SET);

        int rows = r+1; // 
        int cols = dim+1; // 
        //float diag[rows][cols]; // 

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                fscanf(F, "%le"  ,&diag[i][j]);                
            }
        }    
    }

    // A, b, bs, bbs for .mat butcher tableau
    if(butcher_matrix_found && !butcher_crs_matrix_found)
    {
        F = fopen(butcher_matrix, "r");
        fseek(F, butcher_beg_of_val, SEEK_SET);
        for(uint64_t i = 0; i < butcher_stage; i++)
        {
            for(uint64_t j = 0; j < butcher_stage; j++)
            {
                fscanf(F, "%lf",&A[i][j]);
            }
        }
        for(uint64_t i = 0; i < butcher_stage; i++)
            fscanf(F, "%lf",&b[i]);
        if(embedded_found == 1)
        {
            for(uint64_t i = 0; i < butcher_stage; i++)
                fscanf(F, "%lf",&bs[i]);
        }

        fclose(F);
        printf("butcher mat format output.\n ");
        printf("A:\n ");
        for(unsigned int i=0;i<butcher_stage;i++)
        {
            for(unsigned int j=0;j<butcher_stage;j++)
            {
                printf("%le ",A[i][j]);
            }
            printf("\n ");
        }
        printf("b:\n ");
        for(unsigned int i=0;i<butcher_stage;i++)
        {
            printf("%le ",b[i]);
        }
        if(embedded_found == 1)
        {
            printf("\n ");
            printf("bs:\n ");
            for(unsigned int i=0;i<butcher_stage;i++)
            {
                printf("%le ",bs[i]);
            }
        }
        printf("\n");
    }
    if(threads_found)
    {
        print_affinity(threads);                       
        printf("thread number: %d\n", threads);
    }
    else
    {
        set_thread_affinity(cpu_ids, threads); 
        printf("cpu_ids: ");
        for(unsigned int i=0;i<sizeof(cpu_ids)/sizeof(cpu_ids[0]);i++)
        {
            printf("%d ",cpu_ids[i]);
        }                           
        //print_affinity();
        printf("\nthread number: %d\n", threads);
    }
    //omp_set_num_threads(threads);
    if(tile_version_found)
    {
        if(tile_version == 2)
            printf("the second version of loop-tiling (inspired by version 1))\n");
        else if(tile_version == 1)
            printf("the first version of loop-tiling (inspired by adaptive sparse tiling)\n");
        else 
            {
                tile_version = 1;
                printf("only version 1 or 2, so defaulted execute the version 1\n");
            }
    }
    else
    {
        tile_version = 1;   
        //printf("no input of version, so defaulted execute the version 1\n"); 
    }

    if(tile_size_found)
    {
        if(tile_size > dim) tile_size = dim; // set the max of tile size
        if(tile_version == 2)
        {
            if(tile_size == 1)
                printf("tile size: %" PRIu64 " (no tiling)\n", tile_size);
            else
                printf("tile size: %" PRIu64 "\n", tile_size);
        }
        else if(tile_version == 1)
        {
            if(tile_size == 1)
                printf("panel number: %" PRIu64 " (no tiling)\n", tile_size);
            else
                printf("panel number: %" PRIu64 "\n", tile_size);
        }
    }
    else
    {
        tile_size = 1; 
        if(tile_version == 2)
            printf("tile size: %" PRIu64 " (no loop-tiling)\n", tile_size);  
        else
            printf("panel number: %" PRIu64 " (no loop-tiling)\n", tile_size);  
    }

    // val, col, row for .crs butcher tableau
    if(butcher_crs_matrix_found)
    {
        F = fopen(butcher_matrix, "r");
        fseek(F, butcher_beg_of_val, SEEK_SET);
        
        for(uint64_t i = 0; i < butcher_crs_val_size; i++)
            fscanf(F, "%lf" ,&butcher_crs_val[i]);
        
        for(uint64_t i = 0; i < butcher_crs_val_size; i++)
            fscanf(F, "%" PRIu64 ,&butcher_crs_col[i]);
        // maybe explicit(+2 for b)

        for(uint64_t i = 0; i < butcher_crs_dim + 2; i++)
            fscanf(F, "%" PRIu64 ,&butcher_crs_rowPtr[i]);
        

        fclose(F);
        printf("butcher crs format output.\n ");
        for(unsigned int i=0;i<butcher_crs_val_size;i++)
        {
            printf("%le ",butcher_crs_val[i]);
        }
        printf("\n ");
        for(unsigned int i=0;i<butcher_crs_val_size;i++)
        {
            printf("%" PRIu64 " ",butcher_crs_col[i]);
        }
        printf("\n ");
        // maybe explicit(+2 for b)
        for(unsigned int i=0;i<butcher_crs_dim+2;i++) 
        {
            printf("%" PRIu64 " ",butcher_crs_rowPtr[i]);
        }
        printf("\n");

    }
    // measurement
    h = (end_time - t0)/((long double) steps);
    double h_min = h* 0.1;
    double gesamt = (end_time - t0);
    stop_condition = 0.01 *(long double) h / ((long double) dim);

    //start timer
    gettimeofday(&start, NULL);
    
    printf("Solving ODE\n");
    if(!butcher_crs_matrix_found && !butcher_matrix_found){ // rk4
        printf("classic RK4 method\n");
        for(uint64_t i = 0; i < steps; i++){
            checksum = 0;

            #pragma omp parallel for num_threads(threads)
            for(unsigned int k = 0; k < dim; k++)
            {
                double T = 0;
                k1[k] = 0;
                int addMax = rowPtr[k+1] - rowPtr[k];
                for(int i = 0; i < addMax; i++)
                    T += val[rowPtr[k] + i] * x[col[rowPtr[k] + i]];
                x1[k] = (0.5 * h * T) + x[k];
                k1[k] = h * T;
            }


            #pragma omp parallel for num_threads(threads)
            for(unsigned int k = 0; k < dim; k++)
            {
                double T = 0;
                k2[k] = 0;
                int addMax = rowPtr[k+1] - rowPtr[k];
                for(int i = 0; i < addMax; i++)
                    T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];
                x2[k] = (0.5 * h * T) + x[k];
                k2[k] = h * T;
            }
        
            #pragma omp parallel for num_threads(threads)
            for(unsigned int k = 0; k < dim; k++)
            {
                double T = 0;
                k3[k] = 0;
                int addMax = rowPtr[k+1] - rowPtr[k];
                for(int i = 0; i < addMax; i++)
                    T += val[rowPtr[k] + i] * x2[col[rowPtr[k] + i]];
                x1[k] = (h * T) + x[k];
                k3[k] = h * T;
            }

            #pragma omp parallel for num_threads(threads)
            for(unsigned int k = 0; k < dim; k++)
            {
                double T = 0;
                k4[k] = 0;
                int addMax = rowPtr[k+1] - rowPtr[k];
                for(int i = 0; i < addMax; i++)
                {
                    T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];
                }
                //x2[k] = (h * T) + x[k];
                k4[k] = h * T;
            }

            for(unsigned int k = 0; k < dim; k++){
                x[k] += k1[k]*0.1667 + k2[k] * 0.3333 + k3[k] * 0.3333 + k4[k]*0.1667;
                //x[k] += (k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.;
                checksum = fabs(k1[k]*0.1667 + k2[k] * 0.3333 + k3[k] * 0.3333 + k4[k]*0.1667);
                //checksum = fabs((k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.);
                assert(fabs(checksum - 1.0) > 1e-6);
            }

            t0 += h;
            if(checksum < stop_condition){
            early_stop++;
            steps_at_stop = i;
            if(early_stop > 10)
                break;
            }
        }
        for(unsigned int i=0;i<dim;i++)
        {
            printf("x[%d] = %le\n",i,x[i]);
        }
    }

    else if(explicit_found && butcher_crs_matrix_found) // explicit method with crs butcher tableau
    {
        printf("explicit method with crs butcher tableau\n");
        if(Diagonal_found){ // diagonal not acceptable for loop-tiling
            printf("diagonal format unacceptable! only accept crs or new crs formats!\n");
            exit(-1);
        }                            

        for(uint64_t i = 0; i < steps; i++){
            //checksum = 0;
            double local_checksum = 0;
            double T;
            double iter;

            for(unsigned int cur_stage=0;cur_stage<butcher_crs_dim;cur_stage++)
            {
                int addButcher = 0;
                if(cur_stage!=butcher_crs_dim-1)
                {
                    addButcher = butcher_crs_rowPtr[cur_stage+2]-butcher_crs_rowPtr[cur_stage+1];
                }
                //#pragma omp parallel for
                #pragma omp parallel for private(iter,T) num_threads(threads)
                for(unsigned int k = 0; k < dim; k++)
                {
                    T = 0;
                    iter = 0;
                    ki[cur_stage][k] = 0;
                    int addMax = rowPtr[k+1] - rowPtr[k];
                    if(cur_stage == 0)
                    {
                        for(int i = 0; i < addMax; i++)
                        {T += val[rowPtr[k] + i] * x[col[rowPtr[k] + i]];}
                        ki[cur_stage][k] = h*T;
                    }
                    else
                    {
                        for(int i = 0; i < addMax; i++)
                        {T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];}
                        ki[cur_stage][k] = h*T;
                    }
                    if(cur_stage!=butcher_stage-1)
                    {
                        for(int i=0;i<addButcher;i++)
                        {
                            iter+= butcher_crs_val[butcher_crs_rowPtr[cur_stage+2]-1+i] * ki[butcher_crs_col[butcher_crs_rowPtr[cur_stage+1]+i]][k];
                            // sum of ali * Fi
                        }
                        x2[k] = iter + x[k];
                    }
                }
                // Swap x1 and x2, such that 
                double *swap_helper = x1;
                x1 = x2;
                x2 = swap_helper;

            }

            //int addButcher = butcher_crs_rowPtr[butcher_crs_dim+1]-butcher_crs_rowPtr[butcher_crs_dim];
            int start_index = butcher_crs_rowPtr[butcher_crs_dim];
            
            #pragma omp parallel for reduction(+:local_checksum) num_threads(threads)
            for(unsigned int k = 0; k < dim; k++)
            {
                local_checksum = 0;
                for(unsigned int i=0;i<butcher_crs_dim;i++)
                {
                    x[k] += ki[i][k] * butcher_crs_val[start_index+i]; // sum of Fi * bi
                    local_checksum += fabs(ki[i][k] *butcher_crs_val[start_index+i]);
                }
                assert(fabs(local_checksum - 1.0) > 1e-6);
            }

            t0 += h;
            if(local_checksum < stop_condition){
            early_stop++;
            steps_at_stop = i;
            if(early_stop > 10)
                break;
            }
        }
        for(unsigned int i=0;i<dim;i++)
        {
            printf("x[%d] = %le\n",i,x[i]);
        }

    }

    else if(explicit_found && butcher_matrix_found && richardson_found) // richardson extrapolation
    {
        printf("richardson expolation aktived!\n");
        while(counter < gesamt){ 
            double T = 0;
            double iter = 0;

            for(unsigned int cur_stage=0; cur_stage< butcher_stage; cur_stage++)
            {
                #pragma omp parallel for private(iter,T) num_threads(threads) 
                for(unsigned int k = 0; k < dim; k++)
                {
                    T = 0;
                    iter = 0;
                    ki[cur_stage][k] = 0;
                    if(CRS_found) // CRS
                    {
                        int addMax = rowPtr[k+1] - rowPtr[k];
                        if(cur_stage == 0)
                        #pragma simd
                        for(int i = 0; i < addMax; i++) // loop-unswitching
                        {   T += val[rowPtr[k] + i] * x[col[rowPtr[k] + i]];}
                        else
                        #pragma simd
                        for(int i = 0; i < addMax; i++)
                        {   T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];}
                        ki[cur_stage][k] = h*T;
                    }

                    if(CRS_new_found) // new CRS
                    {
                        if(cur_stage == 0)
                        #pragma simd
                        for(unsigned int i = 0; i < addmax[k]; i++) // loop-unswitching
                        {   T += val[rowPtr[k] + i] * x[index[k][i]];}
                        else
                        #pragma simd
                        for(unsigned int i = 0; i < addmax[k]; i++)
                        {   T += val[rowPtr[k] + i] * x1[index[k][i]];}
                        ki[cur_stage][k] = h*T;
                    }

                    if(Diagonal_found) // Diagonal
                    {
                        #pragma simd
                        for (int d = 0; d <= r; d++) // r+1
                        {
                            int diag_index = (int)diag[d][0]; // The Position of diagonal
                            int x_index = k + diag_index; // the index of the needed element in x

                            if (x_index >= 0 && (unsigned int)x_index < dim) // make sure the index is in range
                            { 
                                if(diag[d][k+1]==0)
                                    continue;
                                if (cur_stage == 0) {
                                    T += diag[d][k + 1] * x[x_index];
                                } else {
                                    T += diag[d][k + 1] * x1[x_index];
                                }
                            }
                        }
                        ki[cur_stage][k] = h * T;
                    }

                    if(cur_stage!= (butcher_stage-1))
                    {
                        for(unsigned int i=0;i<butcher_stage;i++)
                        {
                            iter+= A[cur_stage+1][i] * ki[i][k];
                        }
                        x2[k] = iter + x[k];
                    }
                }

                // Swap x1 and x2, such that 
                double *swap_helper = x1;
                x1 = x2;
                x2 = swap_helper;
            }

            for(unsigned int k = 0; k < dim; k++)
            {
                for(unsigned int i = 0; i < butcher_stage; i++)
                {
                    x[k] += ki[i][k] * b[i]; // sum of Fi * bi
                }
            }

            // 2 times calculate for x_half
            for(int i = 0;i<2;i++)
            {
                for(unsigned int cur_stage=0; cur_stage< butcher_stage; cur_stage++)
                {
                    #pragma omp parallel for private(iter,T) num_threads(threads) 
                    for(unsigned int k = 0; k < dim; k++)
                    {
                        T = 0;
                        iter = 0;
                        ki_half[cur_stage][k] = 0;
                        if(CRS_found) // CRS
                        {
                            int addMax = rowPtr[k+1] - rowPtr[k];
                            for(int i = 0; i < addMax; i++)
                            {
                                if(cur_stage==0)
                                    T += val[rowPtr[k] + i] * x_half[col[rowPtr[k] + i]];
                                else 
                                    T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];
                            }
                            ki_half[cur_stage][k] = (h/2)*T;
                        }

                        if(CRS_new_found) // new CRS
                        {
                            for(unsigned int i = 0; i < addmax[k]; i++)
                            {
                                if(cur_stage==0)
                                    T += val[rowPtr[k] + i] * x_half[index[k][i]];
                                else 
                                    T += val[rowPtr[k] + i] * x1[index[k][i]];
                            }
                            ki_half[cur_stage][k] = (h/2)*T;
                        }

                        if(Diagonal_found) // Diagonal
                        {
                            for (int d = 0; d <= r; d++) 
                            {
                                int diag_index = (int)diag[d][0]; // The Position of diagonal
                                int x_index = k + diag_index; // the index of the needed element in x

                                if (x_index >= 0 && (unsigned int)x_index < dim) // make sure the index is in range
                                { 
                                    if(diag[d][k+1]==0)
                                        continue;
                                    if (cur_stage == 0) {
                                        T += diag[d][k + 1] * x_half[x_index];
                                    } else {
                                        T += diag[d][k + 1] * x1[x_index];
                                    }
                                }
                            }
                            ki_half[cur_stage][k] = (h/2) * T;
                        }

                        if(cur_stage!= (butcher_stage-1))
                        {
                            for(unsigned int i=0;i<butcher_stage;i++)
                            {
                                iter+= A[cur_stage+1][i] * ki_half[i][k];
                            }
                            x2[k] = iter + x_half[k];
                        }
                    }

                    // Swap x1 and x2, such that 
                    double *swap_helper = x1;
                    x1 = x2;
                    x2 = swap_helper;
                }
                for(unsigned int k = 0; k < dim; k++)
                {
                    for(unsigned int i = 0; i < butcher_stage; i++)
                    {
                        x_half[k] += ki_half[i][k] * b[i]; // sum of Fi * h * bsi (= y*i)
                    }
                }
            }

            for(unsigned int k=0;k<dim;k++)
            {
                err_absolute[k] = fabs(x[k] - x_half[k]);// absolute error, u_h - u_2x(h/2)
            }
            // calculate scale vector
            #pragma omp parallel for num_threads(threads)
            for (unsigned int i = 0; i < dim; i++)
            {
                double atol = 1, rtol = 1;

                if(fabs(x[i]) <= fabs(x_half[i]))
                {
                    y_scale[i] = atol + x_half[i] * rtol + 1.0e-50;
                }
                else 
                {
                    y_scale[i] = atol + x_last[i] * rtol + 1.0e-50;
                }
            }

            double err_max = 0; //max relative err
            double temp;

            for(unsigned int i=0;i<dim;i++)
            {
                temp = fabs(err_absolute[i]) / ((pow(2,ord) - 1) * y_scale[i]);
                if (temp > err_max)
                {
                    err_max = temp;
                }
            }
            if(tol == 0)
            {
                tol = 1;// default tol, if no input
            }
            err_max = err_max / tol;

            temp = min(a_max_richardson, max(a_min_richardson, 0.08*pow(1.0 / err_max, 1.0 / (butcher_stage + 1.0))));

            if(temp > 1.5) temp = 1.5;
            if(err_max <= 1.0) 
            {
                #pragma omp parallel for num_threads(threads)
                for(unsigned int k=0;k<dim;k++)
                {
                    x[k] = x_half[k] + (x_half[k] - x[k]) / (pow(2,butcher_stage)-1);
                }
                #pragma omp parallel for num_threads(threads)
                for(unsigned int i=0;i<dim;i++)
                {
                    x_last[i] = x[i];
                }

                counter += h;                
            }
            else
            {
                #pragma omp parallel for num_threads(threads)
                for(unsigned int i=0;i<dim;i++)
                {
                    x[i] = x_last[i];
                }
                counter += h/2;
            }
            h *= temp;
            if(h < h_min) h = h_min;
            //printf("use new step h = %le in %d iteration\n",h, counter_iter);
            counter_iter++;
        }
        for(unsigned int i=0;i<dim;i++)
        {
            printf("x[%d] = %le\n",i,x[i]);
            //printf("%le\n",i,x[i]);
        }
        printf("tolerance = %le\n", tol);
        printf("totally %" PRIu64 " iterations in richardson extrapolation\n",counter_iter);

        if(counter < gesamt)
        {
            printf("stoped early!!!!\n");
            printf("counter: %le, steps: %le\n",counter,gesamt);
        }
    }

    else if(explicit_found && butcher_matrix_found && (tile_size >1) && tile_version == 2) // explicit method with loop-tiling version 2
    {
        printf("explicit method with .mat butcher tableau\n");
        //printf("richardson expolation inaktived!\n");
        printf("loop-tiling version 2\n");
        if(Diagonal_found){ // diagonal not acceptable for loop-tiling
            printf("diagonal format unacceptable! only accept crs or new crs formats!\n");
            exit(-1);
        }                            

        for(uint64_t i = 0; i < steps; i++){
            //checksum = 0;
            double local_checksum = 0;
            double iter = 0;
            unsigned int kk = 0;
            unsigned int tile_id = 0, tile_idi = 0, tile_pos = 0;
            unsigned int k = 0;

            for(unsigned int cur_stage=0; cur_stage< butcher_stage; cur_stage++)
            {
                #pragma omp parallel for private(tile_id, k, iter, tile_idi, tile_pos) num_threads(threads)
                for(kk = 0; kk < dim; kk+=tile_size) {
                    for(tile_id = 0; tile_id < dim; tile_id+=tile_size) {
                        for(k = kk; k < min(dim, kk + tile_size); k++) {
                            iter = 0;
                            tile_pos = 0;
                            if(tile_id == 0) {
                                T_tiling[cur_stage][k] = 0;
                            }
                            if(CRS_found) {
                                int addMax = rowPtr[k+1] - rowPtr[k];
                                for(int i = 0; i < addMax; i++) {
                                    for(tile_idi = tile_id;tile_idi < min(dim, tile_id+tile_size);tile_idi++)
                                    {
                                        tile_pos = tile_idi;
                                        if(tile_idi == col[rowPtr[k] + i]) 
                                            T_tiling[cur_stage][k] += (cur_stage == 0 ? val[rowPtr[k] + i] * x[col[rowPtr[k] + i]] : val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]]);
                                    }
                                }
                                if(tile_pos == dim -1) 
                                    ki[cur_stage][k] = h * T_tiling[cur_stage][k];
                            }

                            if(CRS_new_found) {
                                for(unsigned int i = 0; i < addmax[k]; i++) {
                                    for(tile_idi = tile_id;tile_idi < min(dim, tile_id+tile_size);tile_idi++)
                                    {
                                        tile_pos = tile_idi;
                                        if(tile_idi == col[rowPtr[k] + i]) 
                                            T_tiling[cur_stage][k] += (cur_stage == 0 ? val[rowPtr[k] + i] * x[col[rowPtr[k] + i]] : val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]]);
                                    }
                                }
                                if(tile_pos == dim -1) 
                                    ki[cur_stage][k] = h * T_tiling[cur_stage][k];
                            }

                            if(cur_stage != (butcher_stage-1)) {
                                for(unsigned int i = 0; i < butcher_stage; i++) {
                                    iter += A[cur_stage+1][i] * ki[i][k];
                                }
                                x2[k] = iter + x[k];
                            }
                        }
                    }
                }

                //Swap x1 and x2, such that 
                double *swap_helper = x1;
                x1 = x2;
                x2 = swap_helper;
            
            }

            #pragma omp parallel for num_threads(threads) reduction(+:local_checksum)
            for(unsigned int k = 0; k < dim; k++)
            {
                local_checksum = 0;
                for(unsigned int i = 0; i < butcher_stage; i++)
                {
                    x[k] += ki[i][k] * b[i]; // sum of Fi * bi
                    local_checksum += fabs(ki[i][k] * b[i]);
                }
                assert(fabs(local_checksum - 1.0) > 1e-6);

            }

            t0 += h;
            if(local_checksum < stop_condition){
            early_stop++;
            steps_at_stop = i;
            if(early_stop > 10)
                break;
            }
        }
        for(unsigned int i=0;i<dim;i++)
        {
            printf("x[%d] = %le\n",i,x[i]);
        }
    }

    else if(explicit_found && butcher_matrix_found && (tile_size >1) && tile_version == 1) // explicit method with loop-tiling version 1
    {
        printf("explicit method with .mat butcher tableau\n");
        printf("loop-tiling version 1\n");
        if(Diagonal_found){ // diagonal not acceptable for loop-tiling
            printf("diagonal format unacceptable! only accept crs or new crs formats!\n");
            exit(-1);
        }                            

        // calculate panel size 
        int panel_size = 0;
        int panel_size_last = 0;
        unsigned int remaining_blocks = 0;
        if(dim % tile_size == 0)
        {
            panel_size = dim / tile_size;
            panel_size_last = panel_size;
        }
        // else // has rest
        // {   panel_size = dim / (tile_size - 1);
        //     panel_size_last = dim - (tile_size - 1) * panel_size;
        // }
        else // has rest
        {   
            //Calculate the number of basic rows for each panel
            int basic_panel_size = dim / (tile_size);
            // Calculate the number of remaining panel (rest)
            remaining_blocks = dim % tile_size;
            // The number of rows in the remaining panels is equal to the number of base rows
            panel_size_last = basic_panel_size;
            // Add one row to each of the panels before remaining panels
            panel_size = basic_panel_size + 1;
            printf("remaining_block: %d\n", remaining_blocks);
        }
        printf("panel size is: %d\n", panel_size);
        printf("panel size in last panel is: %d\n", panel_size_last);

        for(uint64_t i = 0; i < steps; i++){
            //checksum = 0;
            double local_checksum = 0;
            double iter = 0;
            unsigned int panel_id;
            unsigned int tile_id = 0, tile_idi = 0, tile_pos = 0, k = 0;
            unsigned int row_id;
            for(unsigned int cur_stage=0; cur_stage< butcher_stage; cur_stage++)
            {
                #pragma omp parallel for private(iter, panel_id, k, tile_id, row_id, tile_idi, tile_pos) num_threads(threads)
                for(panel_id = 0; panel_id < tile_size; panel_id++) // iterate all row panels
                {
                    unsigned int num_row = 0;// calculate the number of rows in this panel
                    //if(panel_id == tile_size - 1)
                    if(panel_id > remaining_blocks - 1)
                        num_row = panel_size_last;
                    else 
                        num_row = panel_size;

                    for(tile_id = 0; tile_id < dim; tile_id+=tile_size)
                    {
                        for(row_id=0;row_id < num_row;row_id++)
                        {
                            //T = 0;
                            tile_pos = 0;
                            iter = 0;
                            //ki[cur_stage][k] = 0;
                            k = 0;
                            if(panel_id > remaining_blocks - 1)
                                k = remaining_blocks * panel_size + (panel_id - remaining_blocks) * panel_size_last + row_id;
                            else
                                k = num_row * panel_id + row_id;
                            if(tile_id == 0)
                                T_tiling[cur_stage][k] = 0;
                            if(CRS_found) {
                                int addMax = rowPtr[k+1] - rowPtr[k];
                                for(int i = 0; i < addMax; i++) {
                                    for(tile_idi = tile_id;tile_idi < min(dim, tile_id+tile_size);tile_idi++)
                                    {
                                        tile_pos = tile_idi;
                                        if(tile_idi == col[rowPtr[k] + i]) 
                                            T_tiling[cur_stage][k] += (cur_stage == 0 ? val[rowPtr[k] + i] * x[col[rowPtr[k] + i]] : val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]]);
                                    }
                                }
                                if(tile_pos == dim -1) 
                                    ki[cur_stage][k] = h * T_tiling[cur_stage][k];
                            }

                            if(CRS_new_found) {
                                for(unsigned int i = 0; i < addmax[k]; i++) {
                                    for(tile_idi = tile_id;tile_idi < min(dim, tile_id+tile_size);tile_idi++)
                                    {
                                        tile_pos = tile_idi;
                                        if(tile_idi == col[rowPtr[k] + i]) 
                                            T_tiling[cur_stage][k] += (cur_stage == 0 ? val[rowPtr[k] + i] * x[col[rowPtr[k] + i]] : val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]]);
                                    }
                                }
                                if(tile_pos == dim -1) 
                                    ki[cur_stage][k] = h * T_tiling[cur_stage][k];
                            }

                            if(cur_stage!= (butcher_stage-1))
                            {
                                for(unsigned int i=0;i<butcher_stage;i++)
                                {
                                    iter+= A[cur_stage+1][i] * ki[i][k];
                                    //iter_tiling[k]+= A[cur_stage+1][i] * ki[i][k];
                                }
                                x2[k] = iter + x[k];
                                //printf("x2[%d]: %le\n",k, x2[k]);
                            }

                            // if(k==(dim-1) && cur_stage == (butcher_stage - 1))
                            // //if(k==(dim-1))
                            // {
                            //     double *swap_helper = x1;
                            //     x1 = x2;
                            //     x2 = swap_helper;
                            // }
                        }
                    }
                }

                //Swap x1 and x2, such that 
                double *swap_helper = x1;
                x1 = x2;
                x2 = swap_helper;
            
            }

            #pragma omp parallel for num_threads(threads) reduction(+:local_checksum)
            for(unsigned int k = 0; k < dim; k++)
            {
                local_checksum = 0;
                for(unsigned int i = 0; i < butcher_stage; i++)
                {
                    x[k] += ki[i][k] * b[i]; // sum of Fi * bi
                    local_checksum += fabs(ki[i][k] * b[i]);
                }
                assert(fabs(local_checksum - 1.0) > 1e-6);

            }

            t0 += h;
            if(local_checksum < stop_condition){
            early_stop++;
            steps_at_stop = i;
            if(early_stop > 10)
                break;
            }
        }
        for(unsigned int i=0;i<dim;i++)
        {
            printf("x[%d] = %le\n",i,x[i]);
            //printf("%le\n",i,x[i]);
        }
    }
    
    else if(explicit_found && butcher_matrix_found && (tile_size <=1)) // explicit method with mat butcher tableau
    {
        printf("explicit method with .mat butcher tableau\n");
        for(uint64_t i = 0; i < steps; i++){
            //checksum = 0;
            double local_checksum = 0;
            double T = 0;
            double iter = 0;

            for(unsigned int cur_stage=0; cur_stage< butcher_stage; cur_stage++)
            {
                #pragma omp parallel for private(iter,T) num_threads(threads) //schedule(static, chunk_size)//schedule(auto)////
                for(unsigned int k = 0; k < dim; k++)
                {
                    T = 0;
                    iter = 0;
                    ki[cur_stage][k] = 0;
                    if(CRS_found) // CRS
                    {
                        int addMax = rowPtr[k+1] - rowPtr[k];
                        if(cur_stage == 0)
                        #pragma simd
                        for(int i = 0; i < addMax; i++) // loop-unswitching
                        {   T += val[rowPtr[k] + i] * x[col[rowPtr[k] + i]];}
                        else
                        #pragma simd
                        for(int i = 0; i < addMax; i++)
                        {   T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];}
                        ki[cur_stage][k] = h*T;
                    }

                    if(CRS_new_found) // new CRS
                    {
                        if(cur_stage == 0)
                        #pragma simd
                        for(unsigned int i = 0; i < addmax[k]; i++) // loop-unswitching
                        {   T += val[rowPtr[k] + i] * x[index[k][i]];}
                        else
                        #pragma simd
                        for(unsigned int i = 0; i < addmax[k]; i++)
                        {   T += val[rowPtr[k] + i] * x1[index[k][i]];}
                        ki[cur_stage][k] = h*T;
                    }

                    if(Diagonal_found) // Diagonal
                    {
                        #pragma simd
                        for (int d = 0; d <= r; d++) // r+1
                        {
                            int diag_index = (int)diag[d][0]; // The Position of diagonal
                            int x_index = k + diag_index; // the index of the needed element in x

                            if (x_index >= 0 && (unsigned int)x_index < dim) // make sure the index is in range
                            { 
                                if(diag[d][k+1]==0)
                                    continue;
                                if (cur_stage == 0) {
                                    T += diag[d][k + 1] * x[x_index];
                                } else {
                                    T += diag[d][k + 1] * x1[x_index];
                                }
                            }
                        }
                        ki[cur_stage][k] = h * T;
                    }

                    if(cur_stage!= (butcher_stage-1))
                    {
                        for(unsigned int i=0;i<butcher_stage;i++)
                        {
                            iter+= A[cur_stage+1][i] * ki[i][k];
                        }
                        x2[k] = iter + x[k];
                        //printf("x2[%d]: %le\n",k, x2[k]);
                    }
                }

                // Swap x1 and x2, such that 
                double *swap_helper = x1;
                x1 = x2;
                x2 = swap_helper;
            
            }

            #pragma omp parallel for num_threads(threads) reduction(+:local_checksum)
            for(unsigned int k = 0; k < dim; k++)
            {
                local_checksum = 0;
                for(unsigned int i = 0; i < butcher_stage; i++)
                {
                    x[k] += ki[i][k] * b[i]; // sum of Fi * bi
                    local_checksum += fabs(ki[i][k] * b[i]);
                }
                assert(fabs(local_checksum - 1.0) > 1e-6);
            }

            t0 += h;
            if(local_checksum < stop_condition){
            early_stop++;
            steps_at_stop = i;
            if(early_stop > 10)
                break;
            }
        }
        for(unsigned int i=0;i<dim;i++)
        {
            printf("x[%d] = %le\n",i,x[i]);
        }
    }

    else if(embedded_found) // embedded method
    {
        printf("embedded method with .mat butcher tableau\n");
        if(counter < gesamt)
        {
            double T = 0;
            double iter = 0;
            do
            {
                for(unsigned int cur_stage=0; cur_stage< butcher_stage; cur_stage++)
                {
                    #pragma omp parallel for private(iter,T) num_threads(threads) 
                    for(unsigned int k = 0; k < dim; k++)
                    {
                        T = 0;
                        iter = 0;
                        ki[cur_stage][k] = 0;
                        if(CRS_found) // CRS
                        {
                            int addMax = rowPtr[k+1] - rowPtr[k];
                            if(cur_stage == 0)
                            #pragma simd
                            for(int i = 0; i < addMax; i++) // loop-unswitching
                            {   T += val[rowPtr[k] + i] * x[col[rowPtr[k] + i]];}
                            else
                            #pragma simd
                            for(int i = 0; i < addMax; i++)
                            {   T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];}
                            ki[cur_stage][k] = h*T;
                        }

                        if(CRS_new_found) // new CRS
                        {
                            if(cur_stage == 0)
                            #pragma simd
                            for(unsigned int i = 0; i < addmax[k]; i++) // loop-unswitching
                            {   T += val[rowPtr[k] + i] * x[index[k][i]];}
                            else
                            #pragma simd
                            for(unsigned int i = 0; i < addmax[k]; i++)
                            {   T += val[rowPtr[k] + i] * x1[index[k][i]];}
                            ki[cur_stage][k] = h*T;
                        }

                        if(Diagonal_found) // Diagonal
                        {
                            #pragma simd
                            for (int d = 0; d <= r; d++) // r+1
                            {
                                int diag_index = (int)diag[d][0]; // The Position of diagonal
                                int x_index = k + diag_index; // the index of the needed element in x

                                if (x_index >= 0 && (unsigned int)x_index < dim) // make sure the index is in range
                                { 
                                    if(diag[d][k+1]==0)
                                        continue;
                                    if (cur_stage == 0) {
                                        T += diag[d][k + 1] * x[x_index];
                                    } else {
                                        T += diag[d][k + 1] * x1[x_index];
                                    }
                                }
                            }
                            ki[cur_stage][k] = h * T;
                        }

                        if(cur_stage!= (butcher_stage-1))
                        {
                            for(unsigned int i=0;i<butcher_stage;i++)
                            {
                                iter+= A[cur_stage+1][i] * ki[i][k];
                            }
                            x2[k] = iter + x[k];
                        }
                    }

                    // Swap x1 and x2, such that 
                    double *swap_helper = x1;
                    x1 = x2;
                    x2 = swap_helper;
                }

                #pragma omp parallel for num_threads(threads)
                for(unsigned int k = 0; k < dim; k++)
                {
                    for(unsigned int i = 0; i < butcher_stage; i++)
                    {
                        x[k] += ki[i][k] * b[i]; // sum of Fi * bi
                        x_bs[k] += ki[i][k] * bs[i]; // sum of Fi * h * bsi (= y*i)
                    }
                    err_absolute[k] = fabs(x[k] - x_bs[k]);// absolute error, yi - y*i
                }

                // calculate scale  
                #pragma omp parallel for num_threads(threads)
                for (unsigned int i = 0; i < dim; i++)
                {
                    double atol = 1, rtol = 1; // use tol as tolerance

                    if(fabs(x[i]) <= fabs(x_last[i]))
                    {
                        y_scale[i] = atol + x_last[i] * rtol + 1.0e-50;
                    }
                    else 
                    {
                        y_scale[i] = atol + x[i] * rtol + 1.0e-50;
                    }
                }

                double err_max = 0; //max relative err
                double temp;

                for(unsigned int i=0;i<dim;i++)
                {
                    temp = fabs(err_absolute[i] / y_scale[i]);
                    if (temp > err_max)
                    {
                        err_max = temp;
                    }
                }
                if(tol == 0)
                {
                    tol = 1;// default tol, if no input
                }
                err_max = err_max / tol;
                if(err_max == 0) 
                    err_max = 1.0e-50;
                
                temp = min(a_max_embedded, max(a_min_embedded, 0.06*pow(1.0 / err_max, 1.0 / (ord))));
                if(err_max <= 1.0) // accept
                {
                    #pragma omp parallel for num_threads(threads)
                    for(unsigned int i=0;i<dim;i++)
                    {
                        x_last[i] = x[i];
                    }
                    counter += h;                

                    h_new = h * temp;
                    if(h_new < h_min) h_new = h_min; 
                    h = h_new;
                }
                else // reject
                {
                    #pragma omp parallel for num_threads(threads)
                    for(unsigned int i=0;i<dim;i++)
                    {
                        x[i] = x_last[i];
                    }
                    counter += h/2;                

                    h_new = h * temp;
                    if(h_new < h_min) h_new = h_min; 
                    h = h_new;
                }
                counter_iter++;
            }while(counter < gesamt);
        }
        if(counter < gesamt)
        {
            printf("stoped early!!!!\n");
            printf("counter: %le, steps: %le\n",counter,gesamt);
        }
        for(unsigned int i=0;i<dim;i++)
        {
            printf("x[%d] = %le\n",i,x[i]);
        }
        printf("tolerance = %le\n", tol);
        printf("totally %" PRIu64 " iterations in embedded method\n",counter_iter);
    }

    //stop timer
    gettimeofday(&stop, NULL);
    secs = (double)(stop.tv_usec - start.tv_usec)/1000000 + (double)(stop.tv_sec - start.tv_sec);

    printf("Time: %f\n", secs);
    
    checksum = 0;
    
    // write the solution into the output document
    F = fopen(output_filename, "w");
    for(unsigned int i = 0; i < dim; i++)
        fprintf(F, "%.20g\n", x[i]);
    
    fclose(F);
    
    // calculate the sum of the init vector
    for(unsigned int i = 0; i < dim; i++)
        checksum += x[i];
    
    printf("Checksum over result vector: %lf\n", checksum);
    
    if(of_found){
        int size = strlen(output_filename);
        size += 15; 
        
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);    
    free(x);
    free(x_half);
    free(x_last);
    free(x1);
    free(x2);
    free(swap_helper1);
    free(val);
    free(col);
    free(rowPtr);
    free(addmax);
    for(unsigned int i=0;i<dim+1;i++)
    {
        free(index[i]);
    }
    free(index);
    for(int i=0;i<r+1;i++)
    {
        free(diag[i]);
    }
    free(diag);

    if(matrix_found)
        free(matrix);
    if(iv_found)
        free(init_vector);
    if(of_found)
        free(output_filename);
    if(tol_found)
        free(tol_str);
    if(tile_size_found)
        free(tile_size_str);
    if(tile_version_found)
        free(tile_version_str);
    if(threads_found)
        free(thread_str);

    if(butcher_matrix_found)
    {
        free(butcher_matrix);
        for(unsigned int i=0;i<butcher_stage;i++)
        {
            free(T_tiling[i]);
            free(ki[i]);
            free(ki_half[i]);
            free(A[i]);
        }
        free(iter_tiling);
        free(T_tiling);
        free(ki);
        free(ki_half);
        free(A);
        free(b);
        free(err_absolute);
        free(y_scale);

        if(embedded_found)
        {
            free(bs);
            free(bbs);
            free(x_bs);
        }
    }
    if(butcher_crs_matrix_found)
    {
        free(butcher_crs_val);
        free(butcher_crs_col);
        free(butcher_crs_rowPtr);
    }
    return EXIT_SUCCESS;
}
