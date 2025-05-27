#ifndef GRAPHLIB_H
#define GRAPHLIB_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct
{
  unsigned char r, g, b;
} pixel_t;

extern const pixel_t white;
extern const pixel_t black;
extern const pixel_t red;
extern const pixel_t green;
extern const pixel_t blue;

typedef struct
{
  pixel_t *data;
  int rx;
  int ry;
} screen_t;

screen_t *create_screen(int rx, int ry);
void destroy_screen(screen_t *screen);
void clear_screen(screen_t *screen);
void fill_screen(screen_t *screen, const pixel_t col);
void write_ppm(screen_t *screen, const char *name);
void set_pixel(screen_t *screen, int x, int y, pixel_t col);

#endif /* GRAPHLIB_H */
