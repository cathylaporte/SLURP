#include "pnpoly.h"

/*
License to Use

Copyright (c) 1970-2003, Wm. Randolph Franklin

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

    Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimers.
    Redistributions in binary form must reproduce the above copyright notice in the documentation and/or other materials provided with the distribution.
    The name of W. Randolph Franklin may not be used to endorse or promote products derived from this Software without specific prior written permission. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
*/

int pnpoly(unsigned int nvert, double *vertx, double *verty, double testx, double testy)
{
  unsigned int i, j;
  int c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}

int pnparallelogram(double Px, double Py, double Qx, double Qy, double Rx, double Ry, double test_x, double test_y)
{
  double detPAPQ = (test_x - Px)*(Qy-Py) - (Qx-Px)*(test_y-Py);
  double detPAPR = (test_x - Px)*(Ry-Py) - (Rx-Px)*(test_y-Py);
  double detPQPR = (Qx-Px)*(Ry-Py) - (Rx-Px)*(Qy-Py);

  double r1 = -detPAPQ/detPQPR;
  double r2 = detPAPR/detPQPR;

  return r1 >= 0 && r1 <= 1 && r2 >= 0 && r2 <= 1;
}

