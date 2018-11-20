/* EdgeTrak snake mex-gateway */

#include <mex.h>
#include <stdio.h>
#include "snake.h"

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#if defined(_WIN32) || defined(_WIN64)
#define fmax max
#define fmin min
#pragma warning (disable:4996)
#endif


/*#include <sys/time.h>*/

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  /* will need to deal with options (parameters) at some point */
  double *img, *transpose_init_pts;
  double *init_pts, *snake_pts, *transpose_snake_pts;
  double *snake_energy, *transpose_snake_energy;
  double *internal_snake_energy, *transpose_internal_snake_energy;
  double *external_snake_energy, *transpose_external_snake_energy;
  double * egrad;
  double * transpose_delta;
  int rows, cols, ninit, dim;
  int rows2,cols2;
  int * delta;
  double lambda1, lambda2, alpha, beta, band_penalty;
  int use_band_energy;
  int i;
  int j;
  
  if(nrhs < 8)
    mexErrMsgTxt("Not enough input arguments.\n");
  else if(nrhs > 8)
    mexErrMsgTxt("Too many input arguments.\n");

  /* Validate input matrices */
  /* we handle only doubles */
  if(!mxIsDouble(prhs[0])) 
    mexErrMsgTxt("Input image should be double.\n");
  
  img = mxGetPr(prhs[0]);
  rows = mxGetN(prhs[0]);
  cols = mxGetM(prhs[0]);

  /* get external (Gradient) energy data (faster in Matlab) */
  egrad = mxGetPr(prhs[1]);

  /* Get initial anchor points for snake */
  transpose_init_pts = mxGetPr(prhs[2]);
  ninit = mxGetM(prhs[2]);

  if(ninit < 4)
    mexErrMsgTxt("Not enough points in the initial snake.\n");
  dim = mxGetN(prhs[2]);
  if(dim != 2)
    mexErrMsgTxt("Point array must have 2 rows\n");

  /* convert snake point data to row indexed format */
  init_pts = calloc(ninit*2, sizeof(double));
  for(j = 0; j < ninit; j++)
    {
      init_pts[j*2] = transpose_init_pts[j];
      init_pts[j*2+1] = transpose_init_pts[ninit+j];      
    }

  /* Get delta values for snaxel search regions */
  transpose_delta = mxGetPr(prhs[3]);
  if(mxGetM(prhs[3]) != ninit)
    mexErrMsgTxt("Wrong number of delta values.\n");
  delta = calloc(ninit, sizeof(int));
  for(i = 0; i < ninit; i++)
    {
      delta[i] = (int)(transpose_delta[i]);      
    }
  
  if(nlhs > 4)
    mexErrMsgTxt("Too many output arguments.\n");

  snake_pts = calloc(ninit*2, sizeof(double));
  snake_energy = calloc(ninit, sizeof(double));
  internal_snake_energy = calloc(ninit, sizeof(double));
  external_snake_energy = calloc(ninit, sizeof(double));
  
  /* Do we want to optimize band energy or not?  Not optimizing it reduces computation 
     time */
  band_penalty = (double)(mxGetPr(prhs[4])[0]);
  
  alpha = (double)(mxGetPr(prhs[5])[0]);
  beta = 1-alpha;

  lambda1 = (double)(mxGetPr(prhs[6])[0]);
  lambda2 = 1-lambda1;

  use_band_energy = (int)(mxGetPr(prhs[7])[0]);  
  
  snake_create_adaptive(img, egrad, cols, rows, init_pts, ninit, delta, alpha, beta, lambda1, lambda2, band_penalty, use_band_energy, snake_pts, snake_energy, internal_snake_energy, external_snake_energy);
  
  
  /* Create output matrices */
  
  
  if(nlhs >= 1)
    {      
      plhs[0] = mxCreateDoubleMatrix(ninit, 2, mxREAL); 
      transpose_snake_pts = mxGetPr(plhs[0]);
      
      for(j = 0; j < ninit; j++)
	{
	  transpose_snake_pts[j] = snake_pts[j*2];
	  transpose_snake_pts[ninit+j] = snake_pts[j*2+1];
	}
      
      if(nlhs >= 2)
	{
	  plhs[1] = mxCreateDoubleMatrix(ninit, 1, mxREAL);
	  transpose_snake_energy = mxGetPr(plhs[1]);
	  for(j = 0; j < ninit; j++)
	    transpose_snake_energy[j] = snake_energy[j];	  
		  
	  if(nlhs >= 3)
	    {
	      plhs[2] = mxCreateDoubleMatrix(ninit, 1, mxREAL);
	      transpose_internal_snake_energy = mxGetPr(plhs[2]);
	      for(j = 0; j < ninit; j++)
		transpose_internal_snake_energy[j] = internal_snake_energy[j];
	

	      if(nlhs >= 4)
		{
		  plhs[3] = mxCreateDoubleMatrix(ninit, 1, mxREAL);
		  transpose_external_snake_energy = mxGetPr(plhs[3]);
		  for(j = 0; j < ninit; j++)
		    transpose_external_snake_energy[j] = external_snake_energy[j];
		}
	    }
	}

    }
      
  free(init_pts);
  free(snake_pts);
  free(snake_energy);
  free(internal_snake_energy);
  free(external_snake_energy);
  free(delta);  
}
