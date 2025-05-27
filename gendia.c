/**
 * The program converts a sparse matrix file from COO to diagonal format(.mat -> .mat)
 * USAGE:
 *        [-i inputfile] [-o outputfile] [-h]
 *        example for a chain with 4 sites (n = 4): ./gendia -i test_n_4.mat -o test_n_4_dia.mat
**/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

int main(int argc, char *argv[]) {
    char *inputFileName = NULL;
    char *outputFileName = NULL;
    int opt;

    // Process command line arguments
    while ((opt = getopt(argc, argv, "i:o:")) != -1) {
        switch (opt) {
            case 'i':
                inputFileName = optarg;
                break;
            case 'o':
                outputFileName = optarg;
                break;
            default: /* '?' */
                fprintf(stderr, "Usage: %s -i inputfile -o outputfile\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    if (inputFileName == NULL || outputFileName == NULL) {
        fprintf(stderr, "Input and output files must be specified\n");
        exit(EXIT_FAILURE);
    }

    FILE *file = fopen(inputFileName, "r");
    FILE *outfile = fopen(outputFileName, "w");
    if (file == NULL || outfile == NULL) {
        printf("Error opening file.\n");
        if (file) fclose(file);
        if (outfile) fclose(outfile);
        return 1;
    }
    int size, r, row, col, non_zero_count = 0;
    float value;

    // scan the size of matrix
    fscanf(file, "%d\n%d\n", &size, &r);

    // create original matrix
    float **matrix = (float **)malloc(size * sizeof(float *));
    for (int i = 0; i < size; i++) {
        matrix[i] = (float *)malloc(size * sizeof(float));
        for (int j = 0; j < size; j++) {
            matrix[i][j] = 0;
        }
    }
    // create dia_matrix 
    float **dia_matrix = (float **)malloc((r + 1) * sizeof(float *));
    for (int i = 0; i < r + 1; i++) {
        dia_matrix[i] = (float *)malloc((size + 1) * sizeof(float)); // 改为 size + 1
        for (int j = 0; j < size + 1; j++) {
            dia_matrix[i][j] = 0;
        }
    }    
    // number of remaining diagonals is r-1
    int *index_diagonal = (int*)malloc((r-1)* sizeof(int));

    // scan the non-zero elements
    while (fscanf(file, "%d %d %f\n", &row, &col, &value) == 3) {
        //matrix[row][col] = '*';
        matrix[row][col] = value;

        non_zero_count++;
    }

    // store the values of a specific diagonal
    // for example size = 8, r= 3
    int j=0;
    dia_matrix[0][0] = -size/2 ;
    for(int i= (size/2) ;i<size;i++)
    {
        dia_matrix[0][i+1] = matrix[i][j];
        j++;
    }
    // store the values of main diagonal
    for(int i=0;i<size;i++)
    {
        dia_matrix[1][i+1] = matrix[i][i];
    }
    //Find other diagonals with non-zero elements
    int diagonal_index = 0;
    for (int diag = 1; diag < size; diag++) { // from 1, because the main diagonal is processed
        for (int x = 0; x < size - diag; x++) {
            int y = x + diag;
            if (matrix[x][y] != 0) {
                index_diagonal[diagonal_index++] = diag;
                break; // find the non-zero element, break the loop
            }
            
        }
        if (diagonal_index >= r-1) break; // found all diagonals
    }

// Fill the found diagonal elements into dia_matrix
for (int k = 0; k < r - 1; k++) {
    int diag = index_diagonal[k];
    if (diag < size) { // Make sure diag is within a reasonable range
        dia_matrix[k + 2][0] = (float)diag; // Store diagonal index
        for (int i = 0; i < size - diag; i++) {
            dia_matrix[k + 2][i + 1] = matrix[i][i + diag];
        }
    }
}

int count = 0;
    for (int i = 0; i < r+1; i++) {
        for (int j = 1; j < size + 1; j++) {
            if(dia_matrix[i][j]!=0)
            {
                count++;
            }
        }
    }

    // for(int i=0;i<(r-1);i++)
    // {
    //     printf("the index of remaining diagonals: %d\n",index_diagonal[i]);

    // }
    // store original matrix in output
    // for (int i = 0; i < size; i++) {
    //     for (int j = 0; j < size; j++) {
    //         fprintf(outfile, "%le ", matrix[i][j]);
    //     }
    //     fprintf(outfile, "\n");
    // }
    // fprintf(outfile, "\n");

    // store dia_matrix in output
    fprintf(outfile, "%c\n", 'd');
    fprintf(outfile, "%d\n", count);
    fprintf(outfile, "%d\n", size);

    for (int i = 0; i < r+1; i++) {
        for (int j = 0; j < size + 1; j++) {
            fprintf(outfile, "%le ", dia_matrix[i][j]);
        }
        fprintf(outfile, "\n");
    }

    // output test
    // for (int i = 0; i < r+1; i++) {
    //     for (int j = 0; j < size+1; j++) {
    //         fprintf(outfile, "%le ", dia_matrix[i][j]);
    //     }
    //     fprintf(outfile, "\n");
    // }

    // free memory and close file
    for (int i = 0; i < r+1; i++) free(dia_matrix[i]);
    free(dia_matrix);
    for (int i = 0; i < size; i++) free(matrix[i]);
    free(matrix);
    free(index_diagonal);
    fclose(file);
    fclose(outfile);

    return 0;
}
