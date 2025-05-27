/**
 * The program generates a sparse matrix file in COO format(.mat)
 * USAGE:
 *        -n n [-i inputfile] [-o outputfile] [-h]
 *        example for a chain with 4 sites (n = 4): ./genmat -n 4 -i zahlen.txt -o test_n_4.mat  
**/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

#define POW2(BITPOS) (1UL << (BITPOS))

static inline int test_bit(uint64_t value, uint64_t i)
{
  return (value & POW2(i)) > 0;
}

static inline void set_bit(uint64_t *value, uint64_t i)
{
  *value |= POW2(i);
} 

static inline void clear_bit(uint64_t *value, uint64_t i)
{
  *value &= ~POW2(i);
}

static void usage(const char *progname)
{
  printf("Usage: %s -n n [-i inputfile] [-o outputfile] [-h]\n", progname);
}

int main(int argc, char *argv[])
{
  int opt;
  uint64_t n = 1;
  int n_found = 0;
  char *input_filename = NULL;
  int if_found = 0;
  char *output_filename = NULL;
  int of_found = 0;
  FILE *f;
  uint64_t i, j, r;
  double *h, *col;

  while ((opt = getopt(argc, argv, "n:i:o:h?")) != -1)
  {
    switch (opt)
    {
    case 'n':
      n = strtoull(optarg, NULL, 10);
      n_found = 1;
      break;
    case 'i':
      input_filename = strdup(optarg);
      if_found = 1;
      break;
    case 'o':
      output_filename = strdup(optarg);
      of_found = 1;
      break;
    default: /* '?' */
      usage(argv[0]);
      exit(-1);
    }
  }

  if (!n_found)
  {
    printf("ERROR: Missing argument: -n n\n");
    usage(argv[0]);
    exit(-1);
  }

  if (n > 63)
  {
    printf("ERROR: n must be less than 64.\n");
    exit(-1);
  }

  h = (double *) malloc((n + 1) * sizeof(double));
  if (h == NULL)
  {
    fprintf(stdout, "ERROR: Not enough memory to allocate h.\n");
    exit(-1);
  }
  
  if (if_found)
  {
    if ((f = fopen(input_filename, "r")) == NULL)
    {
      fprintf(stdout, "ERROR: Can't open input file '%s'.\n", input_filename);
      free(h);
      exit(-1);
    }
  }
  else
  {
    f = stdin;
    printf("Please enter %" PRIu64 " probabilities:\n", n+1);
  }

  for (i = 0; i <= n; i++) // TODO: Check that file contains enough values.
    fscanf(f, "%lf", &h[i]);

  if (if_found)
  {
    fclose(f);
  }

  if (of_found)
  {
    if ((f = fopen(output_filename, "w")) == NULL)
    {
      fprintf(stdout, "ERROR: Can't open output file '%s'.\n", output_filename);
      free(h);
      exit(-1);
    }
  }
  else
  {
    f = stdout;
  }

  r = POW2(n); // number of states = matrix dimension
  fprintf(f, "%" PRIu64 "\n%" PRIu64 "\n\n", r, n);
  
  // column can have at most one more entry than there are bit positions
  col = (double *) malloc((n + 1) * sizeof(double)); 
  if (col == NULL)
  {
    fprintf(stdout, "ERROR: Not enough memory to allocate col.\n");
    fclose(f);
    exit(-1);
  }

  for(j = 0; j < r; j++)
  {
    for(i = 0; i <= n; i++)
      col[i] = 0.0;

    if (!test_bit(j, n-1))
    {
      //z = j;
      //set_bit(&z, n-1);
      col[n] = h[0]; // m(z,j) = h[0];
      col[0] -= h[0]; // m(j,j) -= h[0];
    }
      
    for(i = 0; i <= n - 2; i++)
    {
      if (test_bit(j, n-1-i) && !test_bit(j, n-1-i-1)) // 1 0
      {
        //z = j;
        //clear_bit(&z, n-1-i); //a[i] = 0;
        //set_bit(&z, n-1-i-1); // a[i+1] = 1;
        col[n-1-i] = h[i+1]; // m(z,j) = h[i+1];
        col[0] -= h[i+1]; // m(j,j) -= h[i+1];
      }
    }

    if (test_bit(j, 0))
    {
      //z = j;
      //clear_bit(&z, 0);
      col[1] = h[n]; // m(z,j) = h[n];
      col[0] -= h[n]; // m(j,j) -= h[n];
    }

    for (i = n - 1; i > 0; i--)
    {
      if (col[i] != 0.0)
        fprintf(f, "%" PRIu64 " %" PRIu64 " %f\n", j - POW2(i-1), j, col[i]);
    }
    fprintf(f, "%" PRIu64 " %" PRIu64 " %f\n", j, j, col[0]);
    if (col[n] != 0.0)
      fprintf(f, "%" PRIu64 " %" PRIu64 " %f\n", j + POW2(n-1), j, col[n]);
    fprintf(f, "\n");
  }

  free(h);
  free(col);

  if (of_found)
    fclose(f);

  if (if_found) free(input_filename);
  if (of_found) free(output_filename);

  return EXIT_SUCCESS;
}
