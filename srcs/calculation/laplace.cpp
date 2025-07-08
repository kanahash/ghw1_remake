#include "ghw1.hpp"

// 有限差分による 2D Laplace 演算子: y = del^2 x
void laplace( AXY fi, AXY fo )
{
  int i, j, im, ip, jm, jp, im2, ip2, jm2, jp2;
  int i0 = 0, i1 = nx1, j0 = 0, j1 = ny1;
  double aa, bb, cc, dd;

  if (b_4th) { aa = -60.; bb = +16.; cc = -1.; dd = hysq; } // 4次
  else { aa = -4.; bb = +1.; cc = 0.; dd = hy2; }; // 2次

#pragma omp parallel for private(im,ip,im2,ip2,j,jm,jp,jm2,jp2)
  for (i=i0; i<=i1; i++) {
    nextp(i,i0,i1,im,ip,im2,ip2);

    for (j=0; j<=ny1; j++) {
      nextp(j,j0,j1,jm,jp,jm2,jp2);

      fo[i][j] = aa*fi[i][j] + bb*(fi[ip][j]+fi[i][jp]+fi[im][j]+fi[i][jm])
    + cc*(fi[ip2][j]+fi[im2][j]+fi[i][jp2]+fi[i][jm2]);
      fo[i][j]*= dd;
    }
  }
}
