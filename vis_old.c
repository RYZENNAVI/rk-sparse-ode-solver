/**
 * The program generates the visualization (.ppm file) of sparse matrix in MAT format(.mat)
 * USAGE:
 *        [-r resolution] [-i inputfile] [-o outputfile] [-l] [-h]
 *        example for a chain with 4 sites (n = 4): ./vis -i testXXX.mat -o testXXX.ppm
**/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include <string.h>
#include <inttypes.h>

#include "graphlib.h"

static inline int scale(uint64_t v, uint64_t n, int res)
{
  return (int) ((double) v * (double) res / (double) (n));
}

static void usage(const char *progname)
{
  printf("Usage: %s [-r resolution] [-i inputfile] [-o outputfile] [-l] [-h]\n", progname);
}

int main(int argc, char *argv[])
{
  int opt;
  char *input_filename = NULL;
  int if_found = 0;
  char *output_filename = NULL;
  int of_found = 0;
  int lines_found = 0;
  screen_t *S;
  unsigned int res = 1024;
  uint64_t r, c, n, i, p;
  double d;
  FILE *f;

  while ((opt = getopt(argc, argv, "r:i:o:lh?")) != -1)
  {
    switch (opt)
    {
    case 'r':
      res = atoi(optarg);
      break;
    case 'i':
      input_filename = strdup(optarg);
      if_found = 1;
      break;
    case 'o':
      output_filename = strdup(optarg);
      of_found = 1;
      break;
    case 'l':
      lines_found = 1;
      break;
    default: /* '?' */
      usage(argv[0]);
      return EXIT_FAILURE;
    }
  }

  if (if_found)
  {
    if ((f = fopen(input_filename, "r")) == NULL)
    {
      fprintf(stdout, "ERROR: Can't open input file '%s'.\n", input_filename);
      return EXIT_FAILURE;
    }
  }
  else
  {
    f = stdin;
  }

  fscanf(f, "%" SCNu64, &n);

  if (n < res)
    res = n;
  
  S = create_screen(res, res);
  clear_screen(S);

  if (lines_found)
  {
    const pixel_t lightred = {255, 230, 230};
    const pixel_t lightblue = {230, 230, 255};

    for (i = 0; i < n; i++)
    {
      for (p = 1; p < (n >> 1); p <<= 1)
        set_pixel(S, scale(i, n, res), scale(i - p, n, res),
                  lightred);
      set_pixel(S, scale(i, n, res), scale(i + p, n, res),
                lightred);
    }

    for (i = 0; i < n; i++)
    {
      set_pixel(S, scale(i, n, res), scale(n >> 1, n, res),
                lightblue);
      set_pixel(S, scale(n >> 1, n, res), scale(i, n, res),
                lightblue);
    }
  }

  while (!feof(f))
  {
    fscanf(f, "%" SCNu64, &r);
    fscanf(f, "%" SCNu64, &c);
    fscanf(f, "%lf", &d);

    set_pixel(S, scale(c, n, res), scale(r, n, res),
              black);
  }

  if (if_found)
    fclose(f);

  write_ppm(S, output_filename);

  destroy_screen(S);

  if (if_found) free(input_filename);
  if (of_found) free(output_filename);

  return EXIT_SUCCESS;
}

