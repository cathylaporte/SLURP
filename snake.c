/* EdgeTrak snake */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



/*#include <sys/time.h>*/

#include <mex.h>
#include "image.h"
#include "snake.h"
#include "pnpoly.h"
#include "cmath.h"

#if defined(_WIN32) || defined(_WIN64)
#define isfinite _finite
#pragma warning (disable:4996)
#endif


void snake_update_normals(snake * s);
void snake_redistribute_pts(snake * s);
double snake_avg_step_size(const snake * s); 
double internal_energy(double v_x, double v_y, 
		       double prev_x, double prev_y, 
		       double next_x, double next_y,double d,
		       double lambda1, double lambda2);

double fast_band_energy(double v_x, double v_y,
			double prev_x, double prev_y,
			double next_x, double next_y,
			const image * img,
			double d, double band_penalty);

void snake_alloc(snake * s, unsigned int npts)
{
  s->npts = npts;
  s->pts_x = calloc(npts, sizeof(double));
  s->pts_y = calloc(npts, sizeof(double));
  s->norm_x = calloc(npts, sizeof(double));
  s->norm_y = calloc(npts, sizeof(double));
}

void snake_free(snake * s)
{
  free(s->pts_x);
  free(s->pts_y);
  free(s->norm_x);
  free(s->norm_y);
}

void snake_init(snake * s, const double *pt_data)
{
  unsigned int i;
  for(i = 0; i < s->npts; i++)
    {
      s->pts_x[i] = pt_data[2*i];
      s->pts_y[i] = pt_data[2*i+1];
    }

  snake_update_normals(s);
}


void snake_save(const snake * s, double *pt_data)
{
  unsigned int i;
  for(i = 0; i < s->npts; i++)
    {
      pt_data[2*i] = s->pts_x[i];
      pt_data[2*i+1] = s->pts_y[i];
    }
}

void snake_redistribute_pts(snake * s)
{
  double * arclength = calloc(s->npts, sizeof(double));

  double *bx = calloc(s->npts, sizeof(double)); 
  double *cx = calloc(s->npts, sizeof(double));
  double *dx = calloc(s->npts, sizeof(double));
  double *by = calloc(s->npts, sizeof(double)); 
  double *cy = calloc(s->npts, sizeof(double));
  double *dy = calloc(s->npts, sizeof(double));
  double diff_x, diff_y;
  int spl_bool;
  unsigned int i;
  unsigned int first_good = 0;
  unsigned int last_good = s->npts-1;
  unsigned int nknots = s->npts;
  int lastx=0;
  int lasty=0;
  double inc = 0;
  double * new_pts_x;
  double * new_pts_y;

  for(i = 0; s->pts_x[i] < 0 || s->pts_y[i] < 0; i++)
    first_good++;

  for(i = s->npts-1; s->pts_x[i] < 0 || s->pts_y[i] < 0; i--)
    last_good--;
  
  for(i = first_good+1; i < last_good+1; i++)
    {
      diff_x = s->pts_x[i] - s->pts_x[i-1];
      diff_y = s->pts_y[i] - s->pts_y[i-1];
      arclength[i] = arclength[i-1] + sqrt(diff_x*diff_x + diff_y*diff_y);
    }  

  /* increment to use for interpolation */
  inc = arclength[last_good]/(s->npts-1);

  nknots = last_good - first_good + 1;
  spline(nknots,0,0,1,1,&arclength[first_good],&(s->pts_x)[first_good],bx,cx,dx,&spl_bool);
  spline(nknots,0,0,1,1,&arclength[first_good],&(s->pts_y)[first_good],by,cy,dy,&spl_bool);
 
  new_pts_x = calloc(s->npts, sizeof(double));
  new_pts_y = calloc(s->npts, sizeof(double));

  for(i = 0; i < s->npts; i++)
    {
      new_pts_x[i] = seval(nknots,(double)i*inc,&arclength[first_good],&(s->pts_x)[first_good],bx,cx,dx,&lastx);
      new_pts_y[i] = seval(nknots,(double)i*inc,&arclength[first_good],&(s->pts_y)[first_good],by,cy,dy,&lasty);
     
      
    }

  memcpy(s->pts_x, new_pts_x, s->npts*sizeof(double));
  memcpy(s->pts_y, new_pts_y, s->npts*sizeof(double));

  snake_update_normals(s);

  free(arclength);
  free(new_pts_x);
  free(new_pts_y);
  free(bx);
  free(cx);
  free(dx);
  free(by);
  free(cy);
  free(dy);
}

void snake_update_normals(snake * s)
{
  unsigned int i;
  double * tg_x = calloc(s->npts-2, sizeof(double));
  double * tg_y = calloc(s->npts-2, sizeof(double));
  
  double * v0_x = calloc(s->npts-2, sizeof(double));
  double * v0_y = calloc(s->npts-2, sizeof(double));
  
  double * v1_x = calloc(s->npts-2, sizeof(double));
  double * v1_y = calloc(s->npts-2, sizeof(double));
  
  double * u1_x = calloc(s->npts-2, sizeof(double));
  double * u1_y = calloc(s->npts-2, sizeof(double));
  
  double * u0_x = calloc(s->npts-2, sizeof(double));
  double * u0_y = calloc(s->npts-2, sizeof(double));
  
  
  double lenv0;
  double lenv1;
  double len;

  for(i = 0; i < s->npts-2; i++)
    {
      /* compute normalised tangents */
      
      v0_x[i] = s->pts_x[i+1]-s->pts_x[i];
      v0_y[i] = s->pts_y[i+1]-s->pts_y[i];
      
      lenv0 = sqrt(v0_x[i]*v0_x[i] + v0_y[i]*v0_y[i]);   
      
      u0_x[i] = v0_x[i]/lenv0;
      u0_y[i] = v0_y[i]/lenv0;
      
      v1_x[i] = s->pts_x[i+2]-s->pts_x[i+1];
      v1_y[i] = s->pts_y[i+2]-s->pts_y[i+1];
      
      lenv1 = sqrt(v1_x[i]*v1_x[i] + v1_y[i]*v1_y[i]);
      
      u1_x[i] = v1_x[i]/lenv1;
      u1_y[i] = v1_y[i]/lenv1;
      
      tg_x[i] =u0_x[i]+  u1_x[i];
      tg_y[i] =u0_y[i]+  u1_y[i];
      
      len = fabs(tg_y[i]);
      
      if (fabs(tg_x[i])>fabs(tg_y[i])){
        len = fabs(tg_x[i]);
      }
      
      
      s->norm_x[i+1] = -tg_y[i]/len;
      s->norm_y[i+1] = tg_x[i]/len;
      
    }

  /* take care of boundary points */
  s->norm_x[0] = s->norm_x[1];
  s->norm_y[0] = s->norm_y[1];
  s->norm_x[s->npts-1] = s->norm_x[s->npts-2];
  s->norm_y[s->npts-1] = s->norm_y[s->npts-2];

  free(tg_x);
  free(tg_y);
  free(v0_x);
  free(v0_y);
  free(v1_x);
  free(v1_y);
  free(u0_x);
  free(u0_y);
  free(u1_x);
  free(u1_y);
      
  
}

double snake_avg_step_size(const snake * s)
{
  unsigned int i;
  double total_len = 0;
  double diff_x, diff_y;
  for(i = 0; i < s->npts-1; i++)
    {
      diff_x = s->pts_x[i+1]-s->pts_x[i];
      diff_y = s->pts_y[i+1]-s->pts_y[i];
      total_len += sqrt(diff_x*diff_x + diff_y*diff_y);
    }

  return total_len/(s->npts-1);
}

double internal_energy(double v_x, double v_y, 
		       double prev_x, double prev_y, 
		       double next_x, double next_y,double d,
		       double lambda1, double lambda2)
{
  double diff1_x, diff1_y, diff2_x, diff2_y,len1,len2,energy, bending, tensile;
  
  diff1_x = v_x - prev_x;
  diff1_y = v_y - prev_y;
  diff2_x = next_x - v_x;
  diff2_y = next_y - v_y;  

  len1 = sqrt(diff1_x*diff1_x + diff1_y*diff1_y);
  len2 = sqrt(diff2_x*diff2_x + diff2_y*diff2_y);
 
  energy = lambda1*(1.0 - (diff1_x*diff2_x + diff1_y*diff2_y)/(len1*len2)) + lambda2*fabs(len1-d)/d;
    
   if(!isfinite(energy))
    return 0;   
  else
    return energy; 
}

double fast_band_energy(double  v_x, double v_y,
			 double prev_x, double prev_y,
			 double next_x, double next_y,
			 const image * img,
			double d, double band_penalty)
{
  double tg_x, tg_y, len, curve_x, curve_y;
  double norm_x, norm_y;
  double dif = 0;
  int i, j;
  double eband;
  unsigned int count = 0;
  int search_offset;
  double top_vert_x[4];
  double top_vert_y[4];
  double bottom_vert_x[4];
  double bottom_vert_y[4];

  tg_x = next_x - prev_x;
  tg_y = next_y - prev_y;
  len = sqrt(tg_x*tg_x + tg_y*tg_y);
  tg_x /= len;
  tg_y /= len;
  
  norm_x = (d+1)*tg_y;
  norm_y = -(d+1)*tg_x;
/* define band regions (quadrilaterals) based on scaled normals */
  top_vert_x[0] = v_x;
  top_vert_x[1] = v_x + norm_x;
  top_vert_x[2] = next_x + norm_x;
  top_vert_x[3] = next_x;
  top_vert_y[0] = v_y;
  top_vert_y[1] = v_y + norm_y;
  top_vert_y[2] = next_y + norm_y;
  top_vert_y[3] = next_y;

  bottom_vert_x[0] = v_x;
  bottom_vert_x[1] = next_x;
  bottom_vert_x[2] = next_x - norm_x;
  bottom_vert_x[3] = v_x - norm_x;
  bottom_vert_y[0] = v_y;
  bottom_vert_y[1] = next_y;
  bottom_vert_y[2] = next_y - norm_y;
  bottom_vert_y[3] = v_y - norm_y;

/* this could be changed eventually */
  search_offset = (int)(d+1);
  search_offset = 20;

  /* compute normalised difference between pixels in upper and lower band regions by 
filling the quadrilaterals using a fast scan-line method - 
adapted from public domain code by Darel Rex Finley, 2007 */

int  nodes, nodeX[4], pixelX, pixelY, swap ;

/* Loop through the rows of the search region */
 for(pixelY = (int)v_y-search_offset; pixelY < (int)v_y + search_offset+1; pixelY++)
   {
     /* build lists of nodes */
     nodes = 0; j = 3;
     for(i = 0; i < 4; i++)
       {
	 if((int)top_vert_y[i] <= pixelY && (int)(top_vert_y[j]+1) > pixelY ||
	    (int)top_vert_y[j] <= pixelY && (int)(top_vert_y[i]+1) > pixelY) 
	   {
	     
	     if((int)top_vert_y[i] == (int)top_vert_y[j])
	       nodeX[nodes] = (int)top_vert_x[i];
	     else
	       nodeX[nodes] = (int)(top_vert_x[i] + (pixelY - top_vert_y[i])/(top_vert_y[j]-top_vert_y[i])*(top_vert_x[j]-top_vert_x[i]));
	     nodes++;
	   }
	 j = i;
       }

     /* sort the nodes via a simple bubble sort */
     i = 0;
     while(i < nodes -1) 
       {
	 if(nodeX[i] > nodeX[i+1])
	   {
	     swap = nodeX[i];
	     nodeX[i] = nodeX[i+1];
	     nodeX[i+1] = swap;
	     if(i > 0)
	       i--;
	   }
	 else
	   i++;
       }

     /* Gather pixel values between node pairs */
     for(i = 0; i < nodes; i += 2)
       {
	 if(nodeX[i] >= (int)v_x + search_offset)
	   break;
	 if(nodeX[i+1] > (int)v_x - search_offset) 
	   {
	     if(nodeX[i] < (int)v_x - search_offset)
	       nodeX[i] = (int)v_x - search_offset;
	     if(nodeX[i+1] > (int)v_x + search_offset)
	       nodeX[i+1] = (int)v_x + search_offset;
	     for(j = nodeX[i]; j <= nodeX[i+1]; j++)
	       {
		 count++;
		 dif += img->pixels[sub2ind(j, pixelY, img->ncols)];		 
	       }
	   }
       }
   }


/* Loop through the rows of the search region */
 for(pixelY = (int)v_y - search_offset; pixelY < (int)v_y + search_offset; pixelY++)
   {
     /* build lists of nodes */
     nodes = 0; j = 3;
     for(i = 0; i < 4; i++)
       {
	 if((int)bottom_vert_y[i] <= pixelY && (int)(bottom_vert_y[j]+1) > pixelY ||
	    (int)bottom_vert_y[j] <= pixelY && (int)(bottom_vert_y[i]+1) > pixelY) 
	   {
	     if((int)bottom_vert_y[i] == (int)bottom_vert_y[j])
	       nodeX[nodes] = (int)bottom_vert_x[i];
	     else
	     nodeX[nodes] = (int)(bottom_vert_x[i] + (pixelY - bottom_vert_y[i])/(bottom_vert_y[j]-bottom_vert_y[i])*(bottom_vert_x[j]-bottom_vert_x[i]));
	     nodes++;
	   }
	 j = i;
       }

     /* sort the nodes via a simple bubble sort */
     i = 0;
     while(i < nodes -1) 
       {
	 if(nodeX[i] > nodeX[i+1])
	   {
	     swap = nodeX[i];
	     nodeX[i] = nodeX[i+1];
	     nodeX[i+1] = swap;
	     if(i > 0)
	       i--;
	   }
	 else
	   i++;
       }

     /* Gather pixel values between node pairs */
     for(i = 0; i < nodes; i += 2)
       {
	 if(nodeX[i] >= (int)v_x + search_offset)
	   break;
	 if(nodeX[i+1] > (int)v_x - search_offset) 
	   {
	     if(nodeX[i] < (int)v_x - search_offset)
	       nodeX[i] = (int)v_x - search_offset;
	     if(nodeX[i+1] > (int)v_x + search_offset)
	       nodeX[i+1] = (int)v_x + search_offset;
	     for(j = nodeX[i]; j <= nodeX[i+1]; j++)
	       {
		 count++;
		 dif -= img->pixels[sub2ind(j, pixelY, img->ncols)];
	       }
	   }
       }
   }

 int c1, c2, c3, c4;

 if(count == 0) /* the quadrilateral is just a thin strip */
    {
      dif += img->pixels[sub2ind((int)(v_x+norm_x), (int)(v_y+norm_y), img->ncols)];
      dif -= img->pixels[sub2ind((int)(v_x-norm_x), (int)(v_y-norm_y), img->ncols)];

/* Loop through the rows of the search region */
 for(pixelY = (int)v_y-search_offset; pixelY < (int)v_y + search_offset+1; pixelY++)
   {
     /* build lists of nodes */
     nodes = 0; j = 3;
     for(i = 0; i < 4; i++)
       {
	 if((int)top_vert_y[i] <= pixelY && (int)(top_vert_y[j]+1) >= pixelY ||
	    (int)top_vert_y[j] <= pixelY && (int)(top_vert_y[i]+1) >= pixelY) 
	   {
	     if((int)top_vert_y[i] == (int)top_vert_y[j])
	       nodeX[nodes] = (int)(top_vert_x[i]);
	     else
	       nodeX[nodes] = (int)(top_vert_x[i] + (pixelY - top_vert_y[i])/(top_vert_y[j]-top_vert_y[i])*(top_vert_x[j]-top_vert_x[i]));
	     nodes++;
	   }
	 j = i;
       }

     /* sort the nodes via a simple bubble sort */
     i = 0;
     while(i < nodes -1) 
       {
	 if(nodeX[i] > nodeX[i+1])
	   {
	     swap = nodeX[i];
	     nodeX[i] = nodeX[i+1];
	     nodeX[i+1] = swap;
	     if(i > 0)
	       i--;
	   }
	 else
	   i++;
       }

     /* Gather pixel values between node pairs */
     for(i = 0; i < nodes; i += 2)
       {
	 if(nodeX[i] >= (int)v_x + search_offset)
	   break;
	 if(nodeX[i+1] > (int)v_x - search_offset) 
	   {
	     if(nodeX[i] < (int)v_x - search_offset)
	       nodeX[i] = (int)v_x - search_offset;
	     if(nodeX[i+1] > (int)v_x + search_offset)
	       nodeX[i+1] = (int)v_x + search_offset;
	     for(j = nodeX[i]; j <= nodeX[i+1]; j++)
	       {
		 count++;
		 dif += img->pixels[sub2ind(j, pixelY, img->ncols)];		
	       }
	   }
       }
   }

       }
  else
    {
      dif /= (0.5*count);
    }

 
     if(dif < 0)
       eband = band_penalty;
     else
       eband = 1-dif;

  return eband;
}




void snake_create_adaptive(double *img, double * egrad, int ncols, int nrows, const double * init_pts, int ninit, const int * delta, double alpha, double beta, double lambda1, double lambda2, double band_penalty, int use_band_energy, double * snake_pts, double * snake_energy, double * snake_internal_energy, double * snake_external_energy)
{
   int first = true;
  int converge = false;
  int iteration = 0;
  double Eold = 10000;
  double Etotal = 0;
  double Eband, Einternal, minEnergy, Energy;
  double d;
  int minConfig;

  int * nstates;
  int i, j, k, good;
  int yuck = 0;
  int prev, next, cur;  
  unsigned int ncalls = 0;
  double timer_spent = 0;
  int ** search_space ;
  double ** states_x, ** states_y;
  unsigned int ** states_ix;
  double ** EnergyTable;
  unsigned int ** ConfigTable;
  image my_img;
  image Egrad;
  snake my_snake;  

    
  my_img.ncols = ncols;
  my_img.nrows = nrows;
  my_img.pixels = img;

  Egrad.ncols = ncols;
  Egrad.nrows = nrows;
  Egrad.pixels = egrad;

  snake_alloc(&my_snake, ninit);
  snake_init(&my_snake, init_pts);
 

  nstates = calloc(my_snake.npts, sizeof(int));


  for(i = 0; i < my_snake.npts; i++)
    {
      nstates[i] = 2*delta[i]+1;
    }


  search_space = (int **)malloc(my_snake.npts*sizeof(int *));
  for(i = 0; i < my_snake.npts; i++)
    {
      search_space[i] = (int *)calloc(nstates[i], sizeof(int));
    }
  /* Create array of increments for efficient state computation */
  for(i = 0; i < my_snake.npts; i++)
    for(j = -delta[i]; j <= delta[i]; j++)      
      search_space[i][j+delta[i]] = j;

  /* Create state tables for possible x,y positions of each snaxel */
  states_x = (double **)calloc(my_snake.npts, sizeof(double *));
  states_y = (double **)calloc(my_snake.npts, sizeof(double *));
  states_ix = (unsigned int **)calloc(my_snake.npts, sizeof(unsigned int *));
  for(i = 0; i < my_snake.npts; i++)
    {
      states_x[i] = (double *)calloc(nstates[i], sizeof(double));
      states_y[i] = (double *)calloc(nstates[i], sizeof(double));
      states_ix[i] = (unsigned int *)calloc(nstates[i], sizeof(unsigned int));
    }
  /* Create Energy and Configuration tables for dynamic programming */
  EnergyTable = (double **)calloc(my_snake.npts, sizeof(double *));
  ConfigTable = (unsigned int **)calloc(my_snake.npts, sizeof(unsigned int *));
  for(i = 0; i < my_snake.npts-1; i++)
    {
      EnergyTable[i] = (double *)calloc(nstates[i]*nstates[i+1], sizeof(double));
      ConfigTable[i] = (unsigned int *)calloc(nstates[i]*nstates[i+1], sizeof(unsigned int));
    }
  EnergyTable[my_snake.npts-1] = (double *)calloc(nstates[my_snake.npts-1]*nstates[my_snake.npts-1], sizeof(double));
  ConfigTable[my_snake.npts-1] = (unsigned int *)calloc(nstates[my_snake.npts-1]*nstates[my_snake.npts-1], sizeof(unsigned int));


  while(Eold-Etotal > 0) 
      {	
	/* make sure the snake doesn't go beyond image boundaries */
	/* find index of first good point in the snake */
	for(good = 0; my_snake.pts_x[good] <= 20 || my_snake.pts_y[good] <= 20 || my_snake.pts_x[good] >= ncols-20 || my_snake.pts_y[good] >= nrows-20; good++);
	if(good > 0)
	  yuck = 1;		   
	for(i = 0; i < good; i++)
	  {
	    my_snake.pts_x[i] = -1;
	    my_snake.pts_y[i] = -1;
	  }

	/* find index of last good point in the snake*/
	for(good = my_snake.npts-1; my_snake.pts_x[good] <= 20 || my_snake.pts_y[good] <= 20 || my_snake.pts_x[good] >= ncols-20 || my_snake.pts_y[good] >= nrows-20; good--);
	if(good < my_snake.npts-1)
	  yuck = 1;
	for(i = my_snake.npts-1; i > good; i--)
	  {
	    my_snake.pts_x[i] = -1;
	    my_snake.pts_y[i] = -1;
	  }


      /* Redistribute snake points and update normals*/
      snake_redistribute_pts(&my_snake);

    
      d = snake_avg_step_size(&my_snake);
      
      /* Compute possible states from snake normals */
      for(i = 0; i < my_snake.npts; i++)
  	for(j = 0; j < nstates[i]; j++)
  	  {
  	    states_x[i][j] = my_snake.pts_x[i] + search_space[i][j]*my_snake.norm_x[i];
  	    states_y[i][j] = my_snake.pts_y[i] + search_space[i][j]*my_snake.norm_y[i];

  	    /* may want to add code to deal with image boundaries */
	    
  	    states_ix[i][j] = sub2ind((unsigned int)states_x[i][j], (unsigned int)states_y[i][j],
  				      Egrad.ncols);
  	  }
      
      if(!first)
  	Eold = Etotal;
      else
  	first = 0;

      /* Reset Energy and configuration arrays */
      for(i = 0; i < my_snake.npts-1; i++)
  	{
  	  memset(EnergyTable[i], 0, nstates[i]*nstates[i+1]*sizeof(double));
  	  memset(ConfigTable[i], 0, nstates[i]*nstates[i+1]*sizeof(unsigned int));
  	}
      memset(EnergyTable[my_snake.npts-1], 0, nstates[my_snake.npts-1]*nstates[my_snake.npts-1]*sizeof(double));
      memset(ConfigTable[my_snake.npts-1], 0, nstates[my_snake.npts-1]*nstates[my_snake.npts-1]*sizeof(unsigned int));
      
      

      /* ------------- FORWARD PASS -------------*/
      
      /* Do first point in the snake */
      
      for(next = 0; next < nstates[1]; next++)
      	{
      	  for(cur = 0; cur < nstates[0]; cur++)
      	    {
	      if(use_band_energy)
		{
		  
		  Eband = fast_band_energy(states_x[0][cur], states_y[0][cur],
				      states_x[0][cur], states_y[0][cur],
				      states_x[1][next], states_y[1][next],
					   &my_img, d, band_penalty);
		 
		  EnergyTable[0][next*nstates[0]+cur] = beta*Egrad.pixels[states_ix[0][cur]] *Eband; 
		}
	      else
		EnergyTable[0][next*nstates[0]+cur] = beta*Egrad.pixels[states_ix[0][cur]];	
      	    }
      	}
      
      /* Do all the middle points */
      for(i = 1; i < my_snake.npts-1; i++)
      	{
      	  for(next = 0; next < nstates[i+1]; next++)
      	    {
      	      for(cur = 0; cur < nstates[i]; cur++)
      		{
      		  minEnergy = 10000;
      		  minConfig = 0;
      		  for(prev = 0; prev < nstates[i-1]; prev++)
      		    {
		      
      		      Einternal = internal_energy(states_x[i][cur], states_y[i][cur],
      						  states_x[i-1][prev], states_y[i-1][prev],
      						  states_x[i+1][next], states_y[i+1][next],d,
						  lambda1, lambda2);      		 
		      
		      if(use_band_energy)
			{
			  
			  Eband = fast_band_energy(states_x[i][cur], states_y[i][cur],
      					  states_x[i-1][prev], states_y[i-1][prev],
      					  states_x[i+1][next], states_y[i+1][next],
						   &my_img, d, band_penalty); 
			  
			  
			  Energy = EnergyTable[i-1][cur*nstates[i-1]+prev] + 
			    alpha*Einternal +
			    beta*Egrad.pixels[states_ix[i][cur]] *Eband; 
			}
		      else
			{
			  Energy = EnergyTable[i-1][cur*nstates[i-1]+prev] + 
			    alpha*Einternal +
			    beta*Egrad.pixels[states_ix[i][cur]]; 
			}
      		      if(Energy < minEnergy)
      			{
      			  minEnergy = Energy;
      			  minConfig = prev;
      			}
      		    }
      		  
		    EnergyTable[i][next*nstates[i]+cur] = minEnergy;
		    ConfigTable[i][next*nstates[i]+cur] = minConfig;
		    
      		}
      	    }
      	}
      
      /* Do last point in the snake */
      for(cur = 0; cur < nstates[my_snake.npts-1]; cur++)
  	{
  	  minEnergy = 10000;
  	  minConfig = 0;
  	  for(prev = 0; prev < nstates[my_snake.npts-2]; prev++)
  	    {
  	      Einternal = internal_energy(states_x[my_snake.npts-1][cur],
  					  states_y[my_snake.npts-1][cur],
  					  states_x[my_snake.npts-2][prev],
  					  states_y[my_snake.npts-2][prev],
  					  states_x[my_snake.npts-1][cur],
  					  states_y[my_snake.npts-1][cur],d,
					  lambda1, lambda2);
	      
	      
  	      Energy = EnergyTable[my_snake.npts-2][cur*nstates[my_snake.npts-2]+prev] + Einternal;
  	      if(Energy < minEnergy)
  		{
  		  minEnergy = Energy;
  		  minConfig = prev;
  		}
  	    }
	  
  	  EnergyTable[my_snake.npts-1][cur] = minEnergy;
  	  ConfigTable[my_snake.npts-1][cur] = minConfig;
  	}

      /* ------------- BACKWARD PASS ------------*/
           
      minEnergy = 10000;
      minConfig = 0;
      for(i = 0; i < nstates[my_snake.npts-1]; i++)
  	{
  	  if(EnergyTable[my_snake.npts-1][i] < minEnergy)
  	    {
  	      minEnergy = EnergyTable[my_snake.npts-1][i];
  	      minConfig = i;
  	    }
  	}
      
      Etotal = minEnergy;

      my_snake.pts_x[my_snake.npts-1] = states_x[my_snake.npts-1][minConfig];
      my_snake.pts_y[my_snake.npts-1] = states_y[my_snake.npts-1][minConfig];
      cur = minConfig;
      next = 0;
      for(i = my_snake.npts-2; i >= 0; i--)
  	{
  	  prev = ConfigTable[i+1][next*nstates[i+1]+cur];
  	  next = cur;
  	  cur = prev;

  	  my_snake.pts_x[i] = states_x[i][cur];
  	  my_snake.pts_y[i] = states_y[i][cur];
  	}
      
     
      
 }  


  /* Save and clean up */
  snake_save(&my_snake, snake_pts);

  /* Compute energy snaxel-wise (could probably be done inside the DP algorithm) */
  d = snake_avg_step_size(&my_snake);

 
          Eband = fast_band_energy(my_snake.pts_x[0], my_snake.pts_y[0], my_snake.pts_x[0], my_snake.pts_y[0],
				   my_snake.pts_x[1], my_snake.pts_y[1], &my_img, d, band_penalty);
	  
	  

	  snake_internal_energy[0] = 0;
	  snake_external_energy[0] =  Egrad.pixels[sub2ind((int)my_snake.pts_x[0], (int)my_snake.pts_y[0],Egrad.ncols)] * Eband; 

	  snake_energy[0] =  beta*snake_external_energy[0];
	    

  for(i = 1; i < my_snake.npts-1; i++)
  {
   
    Eband = fast_band_energy(my_snake.pts_x[i], my_snake.pts_y[i], 
			    my_snake.pts_x[i-1], my_snake.pts_y[i-1],
			     my_snake.pts_x[i+1], my_snake.pts_y[i+1], &my_img, d, band_penalty);
    
    
	
    snake_internal_energy[i] = internal_energy(my_snake.pts_x[i], my_snake.pts_y[i],
					       my_snake.pts_x[i-1], my_snake.pts_y[i-1],
					       my_snake.pts_x[i+1], my_snake.pts_y[i+1],d,
					       lambda1, lambda2);
    snake_external_energy[i] = Egrad.pixels[sub2ind((int)my_snake.pts_x[i], (int)my_snake.pts_y[i],Egrad.ncols)]* Eband;
    
    snake_energy[i] = alpha * snake_internal_energy[i] + beta * snake_external_energy[i];
    
    
  }
  
  snake_internal_energy[my_snake.npts-1] = internal_energy(my_snake.pts_x[my_snake.npts-1], 
						  my_snake.pts_y[my_snake.npts-1],
						  my_snake.pts_x[my_snake.npts-2], 
						  my_snake.pts_y[my_snake.npts-2],
						  my_snake.pts_x[my_snake.npts-1], 
							   my_snake.pts_y[my_snake.npts-1],d,
							   lambda1, lambda2);
  snake_external_energy[my_snake.npts-1] = 0;
  snake_energy[my_snake.npts-1] = snake_internal_energy[my_snake.npts-1];

  
  
  for(i = 0; i < my_snake.npts; i++)
    {
      free(search_space[i]);
      free(states_x[i]);
      free(states_y[i]);
      free(states_ix[i]);
      free(EnergyTable[i]);
      free(ConfigTable[i]);
    }
  free(states_x);
  free(states_y);
  free(states_ix);
  free(EnergyTable);
  free(ConfigTable);
  free(search_space);
  free(nstates);

  snake_free(&my_snake);
}


