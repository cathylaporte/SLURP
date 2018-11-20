#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "image.h"
#include <mex.h>
#ifndef M_PI
#  define M_PI 3.141592653589793238
#endif


#define MAX_STR_LEN 70
#define _USE_MATH_DEFINES

unsigned int sub2ind(unsigned int x, unsigned int y, unsigned int ncols)
{
  return y*ncols+x;
}

void image_alloc(image * img, unsigned int nrows, unsigned int ncols)
{
  img->nrows = nrows;
  img->ncols = ncols;
  img->pixels = calloc(nrows*ncols, sizeof(double));
}

void image_free(image * img)
{
  free(img->pixels);
}

/* Note: this function assumes that the cropping is done correctly, i.e. x1, x2, y1, y2 are
   all within the size of the original image */
void image_crop(const image * img, 
		unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2,
		image * img_cropped)
{
  int i, j;

  /* allocate space for new image */
  unsigned int nrows_cropped = y2-y1 + 1;
  unsigned int ncols_cropped = x2-x1 + 1;
  image_alloc(img_cropped, nrows_cropped, ncols_cropped);
  
  for(i = y1; i <= y2; i++)
    for(j = x1; j <= x2; j++)
      {
	img_cropped->pixels[(i-y1)*ncols_cropped+(j-x1)] =
	  img->pixels[i*img->ncols + j];
      }
}  

void image_read_pgm(const char * filename, image * img)
{
  FILE * img_file;
  char str[MAX_STR_LEN];
  unsigned int nrows, ncols;  
  unsigned int i, j;
  unsigned char * pix_array;
  
  img_file = fopen(filename, "r");
  
  str[0] = '#';
  
  /* read in magic number and throw it away*/
  while(str[0] == '#')
    {
      fscanf(img_file, "%s", str);
    }
  
  /* now read in the number of columns */
  str[0] = '#';
  while(str[0] == '#')
    {
      fscanf(img_file, "%s", str);
    }
  ncols = atoi(str);

  /* now read in the number of rows */
  str[0] = '#';
  while(str[0] == '#')
    {
      fscanf(img_file, "%s", str);
    }
  nrows = atoi(str);
  
  /* now read in the maximum pixel value and throw it away (we assume it doesn't exceed 255) */
  str[0] = '#';
  while(str[0] == '#')
    {
      fscanf(img_file, "%s", str);
    }

  /* make space for the image */
  image_alloc(img, nrows, ncols);

  pix_array = calloc(nrows*ncols, sizeof(unsigned char));
  /* read in image pixels (we assume there are no comments after the max val) */
  fread(pix_array, sizeof(unsigned char), nrows*ncols, img_file);

  /* convert to floating point */
 
  for(i = 0; i < ncols*nrows; i++)
    img->pixels[i] = (double)(pix_array[i]/255);
  
  free(pix_array);

  fclose(img_file);
}

void image_write_pgm(const char * filename, const image * img)
{
  unsigned char * pix_array;
  unsigned int i, j;
  double maximum = -1000;
  double minimum = 1000;
  FILE * img_file;
  img_file = fopen(filename, "w");
  fprintf(img_file, "P5\n");
  fprintf(img_file, "%d %d\n", img->ncols, img->nrows);
  fprintf(img_file, "255\n");

  /* Convert to uint8 */
  pix_array = calloc(img->nrows*img->ncols, sizeof(unsigned char));

  for(i = 0; i < img->ncols*img->nrows; i++)
    {
      if(img->pixels[i] < minimum)
	minimum = img->pixels[i];
      if(img->pixels[i] > maximum)
	maximum = img->pixels[i];
    }

  for(i = 0; i < img->ncols*img->nrows; i++)
    {
      pix_array[i] = (unsigned char)((img->pixels[i]-minimum)*255.0/(maximum-minimum));
    }

  fwrite(pix_array, sizeof(unsigned char), img->nrows*img->ncols, img_file);

  free(pix_array);

  fclose(img_file);
}

void image_derivative_x(const image * img, double sigma, image * img_dx)
{
  unsigned int mask_size = (unsigned int)(6*sigma+1);
  double * mask = calloc(mask_size*mask_size, sizeof(double));
  int i, j, k, l;
  double x, y;

  /* compute derivative of Gaussian mask along x direction */
  for(i = 0; i < mask_size; i++)
    {
      for(j = 0; j < mask_size; j++)
	{
	  x = (double)j - 3 * sigma;
	  y = (double)i - 3 * sigma;
	  mask[sub2ind(j, i, mask_size)] = -(x/(2*M_PI*pow(sigma, 4)))*exp(-(x*x+y*y)/(2*sigma*sigma));
	}
    }

  /* filter the image */
  for(i = (int)(3*sigma); i < img->nrows - (int)(3*sigma); i++)
    for(j = (int)(3*sigma); j < img->ncols - (int)(3*sigma); j++)
      for(k = (int)(-3*sigma); k < (int)(3*sigma+1); k++)
	for(l = (int)(-3*sigma); l < (int)(3*sigma+1); l++)
	    {	      
	      img_dx->pixels[sub2ind(j, i, img->ncols)] += 
		img->pixels[sub2ind(j+l, i+k, img->ncols)]*
		mask[sub2ind(l+(int)(3*sigma), 
			     k+(int)(3*sigma), mask_size)];	      
	    }
  free(mask);

}

void image_derivative_y(const image * img, double sigma, image * img_dy)
{
  unsigned int mask_size = (unsigned int)(6*sigma+1);
  double * mask = calloc(mask_size*mask_size, sizeof(double));
  int i, j, k, l;
  double x, y;

  for(i = 0; i < mask_size; i++)
    {
      for(j = 0; j < mask_size; j++)
	{
	  x = (double)j - 3 * sigma;
	  y = (double)i - 3 * sigma;
	  mask[sub2ind(j, i, mask_size)] = -(y/(2*M_PI*pow(sigma, 4)))*exp(-(x*x+y*y)/(2*sigma*sigma));

	}
    }
  
  /* filter the image */
  for(i = (int)(3*sigma); i < img->nrows - (int)(3*sigma); i++)
      for(j = (int)(3*sigma); j < img->ncols - (int)(3*sigma); j++)
	for(k = (int)(-3*sigma); k < (int)(3*sigma+1); k++)
	  for(l = (int)(-3*sigma); l < (int)(3*sigma+1); l++)
	    {
	      img_dy->pixels[sub2ind(j, i, img->ncols)] += 
		img->pixels[sub2ind(j+l, i+k, img->ncols)]*
		mask[sub2ind(l+(int)(3*sigma), 
			     k+(int)(3*sigma), mask_size)];	      
	    }   

  free(mask);
}

void image_normalized_gradient_energy(const image * img, double sigma, image * egrad)
{
  unsigned int i;
  double maximum = 0;  
  image dx, dy;
  image_alloc(&dx, img->nrows, img->ncols);
  image_alloc(&dy, img->nrows, img->ncols);

  image_derivative_x(img, sigma, &dx);
  image_derivative_y(img, sigma, &dy);


  for(i = 0; i < img->nrows*img->ncols; i++)
    {
      egrad->pixels[i] = sqrt(dx.pixels[i]*dx.pixels[i] + dy.pixels[i]*dy.pixels[i]);
      if(egrad->pixels[i] > maximum)
	maximum = egrad->pixels[i];
    }

  for(i = 0; i < img->nrows*img->ncols; i++)
    {
      egrad->pixels[i] = 1.0 - egrad->pixels[i]/maximum;
    }

  image_free(&dx);
  image_free(&dy);
}


