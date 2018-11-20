#ifndef IMAGE_H
#define IMAGE_H

typedef struct image {
  unsigned int nrows;
  unsigned int ncols;
  double * pixels;
} image;

unsigned int sub2ind(unsigned int x, unsigned int y, unsigned int ncols);
void image_alloc(image * img, unsigned int nrows, unsigned int ncols);
void image_free(image * img);
void image_crop(const image * img, 
		unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2,
		image * img_cropped);  
void image_read_pgm(const char * filename, image * img);
void image_write_pgm(const char * filename, const image * img);
void image_derivative_x(const image * img, double sigma, image * img_dx);
void image_derivative_y(const image * img, double sigma, image * img_dy);
void image_normalized_gradient_energy(const image * img, double sigma, image * egrad);

#endif
