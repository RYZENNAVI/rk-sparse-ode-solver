#include <stdio.h>
#include <stdlib.h>

int main() {
    FILE *file = fopen("test_n_4.mat", "r");
    FILE *outfile = fopen("test_n_4_vis_new.mat", "w");
    if (file == NULL || outfile == NULL) {
        printf("Error opening file.\n");
        if (file) fclose(file);
        if (outfile) fclose(outfile);
        return 1;
    }

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
