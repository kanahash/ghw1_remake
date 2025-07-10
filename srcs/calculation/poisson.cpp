#include "ghw1.hpp"

// FFT "ポアソン方程式" ソルバー: del^2y = x の y を k 空間で解く
// cp で定義された一般的な引数を持つ (純粋なポアソン問題だけでなく)。
void poisson( AXY fi, const AXY cp, AXY fo )
{
  int i,j,k;
  int nyp = ny/2+1;
  fftw_complex wwk[nx][nyp];
  double wwx[nx][ny];
  fftw_plan hinpsub, herpsub;

#ifdef _OPENMP
  omp_set_num_threads(npar);
  fftw_plan_with_nthreads(npar);
#endif

  fftw_import_wisdom_from_string(wisdom_sf);
  hinpsub=fftw_plan_dft_r2c_2d(nx,ny,*wwx,*wwk,FFTW_ESTIMATE);
  herpsub=fftw_plan_dft_c2r_2d(nx,ny,*wwk,*wwx,FFTW_ESTIMATE);

  int i0 = 0, i1=nx1;

#pragma omp parallel for private(j)
  for (i=i0; i<=i1; ++i) 
  {
	for (j=0; j<=ny1; ++j) 
		wwx[i][j] = fi[i][j]; 
  }
  fftw_execute(hinpsub);

#pragma omp parallel for private(j,k)
  for (i=i0; i<=i1; ++i) 
  {
	for (j=0; j<=nyp-1; ++j)
      for (k=0; k<=1; ++k) wwk[i][j][k] *= cp[i][j];
  }
  fftw_execute(herpsub);

#pragma omp parallel for private(j)
  for (i=i0; i<=i1; ++i) for (j=0; j<=ny1; ++j) fo[i][j] = wwx[i][j];

  fftw_destroy_plan(hinpsub);
  fftw_destroy_plan(herpsub);
}
