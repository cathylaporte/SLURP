#ifndef SNAKE_H
#define SNAKE_H

#include "image.h"

typedef struct snake
{
  double * pts_x;
  double * pts_y;
  double * norm_x;
  double * norm_y;
  unsigned int npts;
  } snake;

void snake_alloc(snake * s, unsigned int npts);
void snake_free(snake * s);
void snake_init(snake * s, const double *pt_data);
void snake_save(const snake * s, double * pt_data);

void snake_create_adaptive(double *img, double * egrad, int ncols, int nrows, const double * init_pts, int nanchors, const int * delta, double alpha, double beta, double lambda1, double lambda2, double band_penalty, int use_band_energy, double * snake_pts, double * snake_energy, double * snake_internal_energy, double * snake_external_energy); 

#endif
