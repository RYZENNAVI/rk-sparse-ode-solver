/**
 * INPUT:   -matrix file in crs format (can be created by genmat -> gencrs)
 *          -start vector file (textfile containing double values delimited by whitespaces)
 * 
 * OUTPUT:
 *          -vector calculated with runge kutta-method after specified start time, # of steps and time period
 *          
 * Programm uses runge kutta-method for calculating solutions to TASEP-Problems
 * */





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


static void usage(const char *progname)
{
  printf("Usage: %s [-m matrix file] [-v start vector file] [-o output file] [-t start time] [-T  endtime] [-n steps] [-h]\n", progname);
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

double * condense(double *x, uint64_t sizeofx)
{
    int c_size = log2_u64(sizeofx);
    double *condensed_x;
    condensed_x = (double *) malloc(((size_t)c_size) * sizeof(double));
    if(condensed_x == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate condensed value array.\n");
        return NULL;
    }

    for(int i = 0; i < c_size; i++)
    {
        condensed_x[i] = 0.;
    }
    
    int bit;
    for(uint64_t k = 0; k < sizeofx; k++)
    {
        for(int i = 0; i < c_size; i++){
            uint64_t temp = k;
            if((temp >> i)&1)
                condensed_x[i] += x[k]; 
        }
    }
    return condensed_x;
}

int main(int argc, char *argv[])
{
    int opt;
    char *matrix = NULL;
    int matrix_found = 0;
    char *init_vector = NULL;
    int iv_found = 0;
    char *output_filename = NULL;
    int of_found = 0;
    double end_time;
    int t_found = 0;
    uint64_t steps;
    int n_found = 0;
    double h, t0, checksum;
    long double stop_condition;
    int t0_found = 0;
    uint64_t val_size = 0, dim = 0, steps_at_stop = 0, *col, *rowPtr;
    double * x, * x1, *val, *k1, *k2, *k3, *k4;
    FILE *F;
    long int beg_of_val;
    int display = 0, early_stop = 0;
    char *t0_str, *end_time_str, *steps_str;
    double secs = 0.;
    int animate = 0, a_found = 0, a_granularity;
    char* animate_temp;
    double* animate_vec;
    int save_choice = 1;
    int frames;

    while ((opt = getopt(argc, argv, "m:v:o:t:T:n:a:h?")) != -1) 
    {
        switch(opt)
        {
        case 'm':
            matrix = strdup(optarg);
            matrix_found = 1;
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
        case 'a':
            animate = 1;
            a_found = 1;
            animate_temp = strdup(optarg);
            a_granularity = atoi(animate_temp);
            free(animate_temp);
            break;
        default:
            usage(argv[0]);
            exit(-1);
        }
    }
    if(a_found){
        frames = (steps / a_granularity);
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

    F = fopen(matrix, "r");
    fscanf(F, "%" PRIu64 "\n%" PRIu64, &val_size, &dim);
    printf("Matrix has %lu entries and dimension %lu.\n", val_size, dim);
    beg_of_val = ftell(F);
    fclose(F);

    int chain_length = log2_u64(dim);

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

    x = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(x == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate x array.\n");
        exit(-1);
    }

    x1 = (double *) malloc(((size_t)dim+1) * sizeof(double));
    if(x1 == NULL)
    {
        fprintf(stdout, "ERROR: Not enough memory to allocate x1 array.\n");
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
            printf("ERROR: Can't open start vector file: %s\n", init_vector);
            exit(-1);
        }
    }
    else
    {
        printf("No initial vector found. Aborting.\n");
        exit(-1);
    }

    if(a_found){
        animate_vec = (double *) malloc((frames * chain_length * sizeof(double)));
        printf("%d\n", frames* chain_length);
        if(animate_vec == NULL){
            printf("ERROR! Not enough space to save animate Vector, continue to save animation data to hard disk?\nCAUTION: This will take subsantially more time[y/n]  ");
            char save_anim_vec;
            scanf( "%c",save_anim_vec);
            if(save_anim_vec == 'y' || save_anim_vec == 'Y')
                save_choice = 0;
            else if(save_anim_vec == 'n' || save_anim_vec == 'N')
                save_choice = 1;
            else
            {
                printf("Selection not clear. Aborting\n");
                exit(-1);
            }
        }
        else{
            for(int i = 0; i < frames * chain_length; i++)
                animate_vec[i] = 0.;
        }
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
        i++;
    }
    fclose(F);
    
    
    F = fopen(matrix, "r");
    fseek(F, beg_of_val, SEEK_SET);
    
    for(uint64_t i = 0; i < val_size; i++)
        fscanf(F, "%lf" ,&val[i]);
    
    for(uint64_t i = 0; i < val_size; i++)
        fscanf(F, "%" PRIu64 ,&col[i]);

    for(uint64_t i = 0; i < dim + 1; i++)
        fscanf(F, "%" PRIu64 ,&rowPtr[i]);

    //for(uint i = 0; i < dim; i++)
    //   printf("x[%d] = %lf\n", i, x[i]);

    fclose(F);
    
    h = (end_time - t0)/((long double) steps);
    stop_condition = 0.01 *(long double) h / ((long double) dim);

    //start timer
    gettimeofday(&start, NULL);

    printf("Solving ODE\n");
    if(!a_found){
        for(uint64_t i = 0; i < steps; i++){
            checksum = 0;

            //#pragma omp parallel for
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

            //#pragma omp parallel for
            for(unsigned int k = 0; k < dim; k++)
            {
                double T = 0;
                k2[k] = 0;
                int addMax = rowPtr[k+1] - rowPtr[k];
                for(int i = 0; i < addMax; i++)
                    T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];
                x1[k] = (0.5 * h * T) + x[k];
                k2[k] = h * T;
            }
        
            //#pragma omp parallel for
            for(unsigned int k = 0; k < dim; k++)
            {
                double T = 0;
                k3[k] = 0;
                int addMax = rowPtr[k+1] - rowPtr[k];
                for(int i = 0; i < addMax; i++)
                    T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];
                x1[k] = (h * T) + x[k];
                k3[k] = h * T;
            }

            //#pragma omp parallel for
            for(unsigned int k = 0; k < dim; k++)
            {
                double T = 0;
                k4[k] = 0;
                int addMax = rowPtr[k+1] - rowPtr[k];
                for(int i = 0; i < addMax; i++)
                    T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];
                x1[k] = (h * T) + x[k];
                k4[k] = h * T;
            }

            for(unsigned int k = 0; k < dim; k++){
            x[k] += (k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.;
            checksum = fabs((k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.);
            }

            t0 += h;
            if(checksum < stop_condition){
            early_stop++;
            steps_at_stop = i;
            if(early_stop > 10)
                break;
            }
        }
        for(int i=0;i<dim;i++)
        {
            printf("x[%d] = %le\n",i,x[i]);
        }

    }
    else{
        int anim_arr_marker = 0;
        for(uint64_t i = 0; i < steps; i++){
            checksum = 0;

            #pragma omp parallel for
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

            #pragma omp parallel for
            for(unsigned int k = 0; k < dim; k++)
            {
                double T = 0;
                k2[k] = 0;
                int addMax = rowPtr[k+1] - rowPtr[k];
                for(int i = 0; i < addMax; i++)
                    T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];
                x1[k] = (0.5 * h * T) + x[k];
                k2[k] = h * T;
            }
        
            #pragma omp parallel for
            for(unsigned int k = 0; k < dim; k++)
            {
                double T = 0;
                k3[k] = 0;
                int addMax = rowPtr[k+1] - rowPtr[k];
                for(int i = 0; i < addMax; i++)
                    T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];
                x1[k] = (h * T) + x[k];
                k3[k] = h * T;
            }

            #pragma omp parallel for
            for(unsigned int k = 0; k < dim; k++)
            {
                double T = 0;
                k4[k] = 0;
                int addMax = rowPtr[k+1] - rowPtr[k];
                for(int i = 0; i < addMax; i++)
                    T += val[rowPtr[k] + i] * x1[col[rowPtr[k] + i]];
                x1[k] = (h * T) + x[k];
                k4[k] = h * T;
            }

            for(unsigned int k = 0; k < dim; k++){
            x[k] += (k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.;
            checksum = fabs((k1[k] + (2. * k2[k]) + (2. * k3[k]) + k4[k]) / 6.);
            }

            if(i%a_granularity == 0){
                double *temp_condensed = condense(x, dim);
                if(save_choice){
                    for(int l = 0; l < chain_length; l++){
                        animate_vec[anim_arr_marker * chain_length + l] = temp_condensed[l];
                    }
                    anim_arr_marker++;
                }
                else{
                    int size = strlen(output_filename);
                    size += 24; 
                    char cavrianceMatrix_ppm_filename[size];
                    snprintf(cavrianceMatrix_ppm_filename, size, "%s_animation_raw_data.vec", output_filename);
                    F = fopen(cavrianceMatrix_ppm_filename, "w");
                    for(int l = 0; l < chain_length; l++)
                        fprintf(F, "%f, ", temp_condensed[l]);
                    fprintf(F, "\n");
                    fclose(F);
                }
                free(temp_condensed);
            }
            t0 += h;
            if(checksum < stop_condition){
            early_stop++;
            steps_at_stop = i;
            if(early_stop > 10)
                break;
            }
        }
        if(save_choice){
            int size = strlen(output_filename);
            size += 24; 
            char cavrianceMatrix_ppm_filename[size];
            snprintf(cavrianceMatrix_ppm_filename, size, "%s_animation_raw_data.vec", output_filename);
            F = fopen(cavrianceMatrix_ppm_filename, "w");
            frames = steps_at_stop / a_granularity;
            for(int l = 1; l <= frames * chain_length; l++){
                fprintf(F, "%f, ", animate_vec[l-1]);
                if(l%chain_length == 0 && l > 5)
                    fprintf(F, "\n");
            }
            double *animate_vec_temp = condense(x, dim);
            for(int l = 0; l < chain_length; l++)
                fprintf(F, "%f, ", animate_vec_temp[l]);
            free(animate_vec_temp);
            fprintf(F, "\n");
            fclose(F);
        }
    }
        //stop timer
    gettimeofday(&stop, NULL);
    secs = (double)(stop.tv_usec - start.tv_usec)/1000000 + (double)(stop.tv_sec - start.tv_sec);

    printf("Time: %f\n", secs);

    if(early_stop){
        printf("Convergence after %" PRIu64 " steps\nConfidience |\u0394x| = %.Lg\n", steps_at_stop, stop_condition);
        printf("Simulation time elapsed: %lf\n", t0);
    }
    
    checksum = 0;
    
    if(display){
        for(unsigned int i = 0; i < dim; i++)
            printf("%lf\n" ,x[i]);
    }
    else{
        /**F = fopen(output_filename, "w");
        for(unsigned int i = 0; i < dim; i++)
            fprintf(F, "%.20g\n", x[i]);
        
        fclose(F);**/
    }
    
    
    for(unsigned int i = 0; i < dim; i++)
        checksum += x[i];
    
    printf("Checksum over result vector: %lf\n", checksum);
    
    
    if(of_found){
        int size = strlen(output_filename);
        size += 15; 
        char expected_values_filename[size];
        
        //sprintf(expected_values_filename, "%s_condensed.txt", output_filename); not used because not safe
        snprintf(expected_values_filename, size, "%s_condensed.vec", output_filename);
        double check = 0.;
        double * y = condense(x, dim);

        if(fopen(expected_values_filename, "w") != NULL){
            F = fopen(expected_values_filename, "w");
            for(unsigned int i = 0; i < log2_u64(dim) - 1; i++)
                fprintf(F, "%lf,", y[i]);

            fprintf(F, "%lf", y[log2_u64(dim) - 1]);
        }
        fclose(F);
        free(y);
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);    
    free(x);
    free(x1);
    free(val);
    free(col);
    free(rowPtr);
    if(a_found)
        free(animate_vec);
    if(matrix_found)
        free(matrix);
    if(iv_found)
        free(init_vector);
    if(of_found)
        free(output_filename);

    return EXIT_SUCCESS;
}
