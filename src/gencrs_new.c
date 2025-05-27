/**
 * The program converts a sparse matrix file from COO to new CRS format(.mat -> .crs)
 * USAGE:
 *        [-i inputfile] [-o outputfile] [-h]
 *        example for a chain with 4 sites (n = 4): ./gencrs_new -i test_n_4.mat -o test_n_4_new.crs
**/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>

#define POW2(BITPOS) (1UL << (BITPOS))

static void usage(const char *progname)
{
  printf("Usage: %s [-i inputfile] [-o outputfile] [-h]\n", progname);
}

int main(int argc, char *argv[])
{
  int opt;
  uint64_t n = 0, d = 0, r = 0;
  char *input_filename = NULL;
  int if_found = 0;
  char *output_filename = NULL;
  int of_found = 0;
  FILE *f;
  uint64_t * col, *col_new, * row, * row_ptr, * hist, *addmax, **index;
  double * val, * val_new;

  while ((opt = getopt(argc, argv, "i:o:h?")) != -1)
  {
    switch (opt)
    {
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

  if(!if_found)
  {
    printf("ERROR: Missing argument: -i\n");
    usage(argv[0]);
    exit(-1);
  }

  f = fopen(input_filename, "r");

  fscanf(f, "%lu" "%lu", &d, &r); // read dimension of Matrix
  //fprintf(stdout, "Matrix in file %s has dimension %lu\n", input_filename, d);
  
  //not an accurate calculation, only crude upper bound
  n = POW2(r) * r;

  //allocate Menory
  col = (uint64_t *) malloc((n+1) * sizeof(long unsigned));
  if(col == NULL)
  {
    fprintf(stdout, "ERROR: Not enough memory to allocate col.\n");
    fclose(f);
    exit(-1);
  }

  row = (uint64_t *) malloc((n+1) * sizeof (long unsigned));
  if(row == NULL)
  {
    fprintf(stdout, "ERROR: Not enough memory to allocate row.\n");
    fclose(f);
    exit(-1);
  }

  val = (double *) malloc((n+1) * sizeof(double));
  if(val == NULL)
  {
    fprintf(stdout, "ERROR: Not enough memory to allocate val.\n");
    fclose(f);
    exit(-1);
  }

  row_ptr = (uint64_t *) malloc((d+1) * sizeof(long unsigned));
  if(row_ptr == NULL)
  {
    fprintf(stdout, "ERROR: Not enough Memory to allocate row_ptr.\n");
    fclose(f);
    exit(-1);
  }

  addmax = (uint64_t *) malloc((d+1) * sizeof(long unsigned));
  if(addmax == NULL)
  {
    fprintf(stdout, "ERROR: Not enough Memory to allocate addmax.\n");
    fclose(f);
    exit(-1);
  }

  index = (long unsigned **)malloc((size_t) (d+1) * sizeof(long unsigned*));
  for(unsigned int i = 0; i < d+1; i++){
      index[i] = (long unsigned*)malloc((d+1) * sizeof(long unsigned));
  }
  if(index == NULL)
  {
      fprintf(stdout, "ERROR: Not enough memory to allocate index array.\n");
      exit(-1);
  }

  //Histogram for entries in row
  hist = (uint64_t *) malloc((d) * sizeof(long unsigned));
  if(hist == NULL)
  {
    fprintf(stdout, "ERROR: Not enough Memory to allocate hist.\n");
    fclose(f);
    exit(-1);
  }

  //set hist entries to 0
  for(uint64_t k = 0; k < d; k++)
    hist[k] = 0;
    
  //read triplets into RAM, create histogram
  
  uint64_t ctr = 0;
  while(!feof(f))
  {
    if (fscanf(f, "%lu %lu %le", &row[ctr], &col[ctr], &val[ctr]) != 3)
      break;
    //printf("%lu %lu %e\n",row[ctr],col[ctr],val[ctr]);
    hist[row[ctr]]++;
    ctr++;
  }

  //fprintf(stdout, "\n\n");
  //only for debugging purposes
  //for(int i = 0; i < d; i++)
    //fprintf(stdout, "Number of %i  entries in hist:%" PRIu64 "\n",i, hist[i]);
  
  //fprintf(stdout, "\n\n");
  
  fclose(f);

  //check for errors
  for(uint32_t k = 0; k < d; k++)
  {
    if(hist[k] > d){
      fprintf(stdout, "\nERROR: A problem has occured while calculating hist at pos %u\n", k);
      exit(-1);
    }
  }
  //calculate row_ptr
  row_ptr[0] = 0;
  for(uint32_t k = 0; k < d; k++)
  {
    row_ptr[k+1] = row_ptr[k] + hist[k];
  }

  // calculate addmax 
  for(unsigned int i=0;i<d;i++)
  {
    addmax[i] = hist[i];
  }

  //resize n to the correct size

  n = row_ptr[d];

  if(of_found)
  {
    if((f = fopen(output_filename, "w")) == NULL)
    {
      fprintf(stdout, "ERROR: Can't open output file '%s'.\n", output_filename);
      exit(-1);
    }
  }

  else
  {
    f = stdout;
  }

  val_new = (double *) malloc((n) * sizeof(double));
  if(val_new == NULL)
  {
    fprintf(stdout, "ERROR: Not enough memory to allocate val_new.\n");
    fclose(f);
    exit(-1);
  }

  col_new = (uint64_t *) malloc((n+1) * sizeof(long unsigned));
  if(col_new == NULL)
  {
    fprintf(stdout, "ERROR: Not enough memory to allocate col.\n");
    fclose(f);
    exit(-1);
  }

  ctr = 0;

  for(uint64_t k = 0; k < n; k++)
  {
    val_new[row_ptr[row[k]+1] - hist[row[k]]] = val[k];
    col_new[row_ptr[row[k]+1] - hist[row[k]]] = col[k];
    hist[row[k]]--;
  }

  free(hist);
  free(row);

  fprintf(f, "%c\n", 'n');
  fprintf(f, "%lu\n", n);
  fprintf(f, "%lu\n", d);
  // val 
  for(uint64_t k = 0; k < n; k++)
    fprintf(f, "%f ", val_new[k]);

  fprintf(f, "\n");
  free(val);
  free(val_new);
  // col 
  for(uint64_t k = 0; k < n; k++)
    fprintf(f, "%lu ", col_new[k]);
  
  fprintf(f, "\n");
  free(col);
  // row_ptr 
  for(uint64_t k = 0; k <= d; k++)
    fprintf(f, "%lu ", row_ptr[k]);

  fprintf(f, "\n");

  // addmax 
  for(uint64_t k = 0; k < d; k++)
    fprintf(f, "%lu ", addmax[k]);

  // index 
  // calculate index 
  for(unsigned int i=0;i<d;i++)
  {
    for(unsigned int j=0;j<addmax[i];j++)
    {
      index[i][j] = col_new[j+row_ptr[i]];
    }
  }
  fprintf(f, "\n");
  for(uint64_t i = 0; i < d; i++)
  {
    for(uint64_t j = 0; j < addmax[i]; j++)
    {
      fprintf(f, "%lu ", index[i][j]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n");
  free(col_new);
  free(row_ptr);
  free(addmax);
  for(unsigned int i=0;i<d+1;i++)
  {
      free(index[i]);
  }
  free(index);
  free(input_filename);
  if (of_found)
    free(output_filename);
}

