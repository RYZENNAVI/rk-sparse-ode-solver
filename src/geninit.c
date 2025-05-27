/**
 * The program generates a start vector file in .init format(.init) with normal distribution
 * USAGE:
 *        [-d size] [-n outputfile] [-h]
 *        example for a chain with 4 sites (n = 4): ./geninit -d 4 -n test_n_4.init
**///

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>

#define POW2(BITPOS) (1ULL << (BITPOS))

static void usage(const char *progname)
{
    printf("Usage: %s [-d size] [-n normal distribution] !NOT YET AVALIABLE BEYOND THIS PIONT![-g gaussian distribution] [-r random distribution] [-v valley distribution]\n", progname);
}

int main(int argc, char *argv[])
{
    int opt, select = -1; 
    char *n_filename = 0;//*g_filename, *r_filename, *v_filename;
    int size = 0;
    int log_size;

    while ((opt = getopt(argc, argv, "d:n:g:r:v:h?")) != -1)
    {
        switch(opt)
        {
        case 'd':
        {
            select = 1;
            char* size_str = strdup(optarg);
            log_size = atoi(size_str);
            free(size_str);
            size = POW2(log_size);
            break;
        }
        case 'n':
            select = 0;
            n_filename = strdup(optarg);
            break;
        // case 'g':
        //     select = 1;
        //     g_filename = strdup(optarg);
        //     break;
        // case 'r':
        //     select = 2;
        //     r_filename = strdup(optarg);
        //     break;
        // case 'v':
        //     select = 3;
        //     v_filename = strdup(optarg);
        //     break;
        default:
            usage(argv[0]);
            exit(-1);
        }
    }

    if(select == -1)
    {
        printf("No size specified. Aborting!\n");
        exit(-1);
    }

    switch(select)
    {
    case 0:
    {
        FILE *f;
        long double check = 0;
        if(fopen(n_filename, "w") == NULL){
            printf("Error. Could not open file %s. Aborting!\n", n_filename);
            exit(-1);
        }
        else{
            f = fopen(n_filename, "w");
            double val = 1./(double) size;
            printf("Generating inital vector\nSize:%d\nValue:%f\n", size, val);
            for(int i = 0; i < size; i++){
                fprintf(f, "%.20g\n", val); 
                check += val;
            }
        }
        printf("%s created.\nChecksum: %Lf\n", n_filename, check);
        free(n_filename);
        break;
    }
    default:
        exit(-1);
    }
}
