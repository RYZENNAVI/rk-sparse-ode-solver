/**
 * The program generates the visualization (.mat file) of sparse matrix in MAT format(.mat)
 * USAGE:
 *        [-i inputfile] [-o outputfile] [-h]
 *        example for a chain with 4 sites (n = 4): ./vis_new -i test_n_4.mat -o test_n_4_vis_new.mat
**/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
    char *inputFilename = NULL;
    char *outputFilename = NULL;
    int opt;

    while ((opt = getopt(argc, argv, "i:o:")) != -1) {
        switch (opt) {
            case 'i':
                inputFilename = optarg;
                break;
            case 'o':
                outputFilename = optarg;
                break;
            default:
                fprintf(stderr, "Usage: %s -i inputfile -o outputfile\n", argv[0]);
                return 1;
        }
    }

    if (inputFilename == NULL || outputFilename == NULL) {
        fprintf(stderr, "Input and output files must be specified.\n");
        return 1;
    }

    FILE *file = fopen(inputFilename, "r");
    FILE *outfile = fopen(outputFilename, "w");
    int size, r, row, col, non_zero_count = 0;
    float value;

    // Read matrix size
    fscanf(file, "%d\n%d\n", &size, &r);

    // Create a two-dimensional array to represent input matrix
    char **matrix = (char **)malloc(size * sizeof(char *));
    for (int i = 0; i < size; i++) {
        matrix[i] = (char *)malloc(size * sizeof(char));
        for (int j = 0; j < size; j++) {
            matrix[i][j] = ' ';
        }
    }

    // Read non-zero elements
    while (fscanf(file, "%d %d %f\n", &row, &col, &value) == 3) {
        matrix[row][col] = '*';
        non_zero_count++;
    }

    // Output matrix to file
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            fprintf(outfile, "%c ", matrix[i][j]);
        }
        fprintf(outfile, "\n");
    }

    // free memory and close file
    for (int i = 0; i < size; i++) free(matrix[i]);
    free(matrix);
    fclose(file);
    fclose(outfile);

    return 0;
}
