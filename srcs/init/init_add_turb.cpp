#include"ghw1.hpp"

// 擬似乱流浴摂動を初期化する:
void init_add_turb( AXY ne, AXY ni )
{
  int i, j, ik, jk;
  int nk = 32, mmm = int(aspect);
  double zuf, kk11, turb, ppk, lx = ly*aspect, phs[nk+1][nk+1], cbath[nx][ny];
  fftw_complex wwk[nx][ny], ppx[nx][ny];
  fftw_plan plan_ini_back;

  printf("| 乱流浴初期条件...\n");

  // 擬似乱流浴初期化のスペクトル:
  double ccc = 1./(lx*ly);
  for (i=0; i<=nx-1; ++i) for (j=0; j<=ny-1; ++j)
    {
      ik = (i>=nx/2) ? i - nx : i;
      jk = (j>=ny/2) ? j - ny : j;
      kk11 =  TwoPi*TwoPi*(  ( double(ik*ik)/(lx*lx) ) + ( double(jk*jk)/(ly*ly) )  );
      cbath[i][j] =  ccc*kk11/(1. + kk11*kk11*kk11 );
    }

  // 再現性を希望する場合、同じ擬似「ランダム」シードを使用してください！

#ifdef _OPENMP
  fftw_plan_with_nthreads(npar);
#endif
  plan_ini_back=fftw_plan_dft_2d(nx,ny,&wwk[0][0],&ppx[0][0],FFTW_BACKWARD,FFTW_ESTIMATE);

  // srand(time(NULL) + getpid()); // 再現可能な実行のためにはランダムを使用しないでください！

  // (rand 上での omp parallel は再現性がない)
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) 
    {
      zuf = TwoPi*0.05*(rand() % 20);
      wwk[i][j][0] = amp*cbath[i][j]*cos(zuf);
      wwk[i][j][1] = amp*cbath[i][j]*sin(zuf);
    }
  fftw_execute(plan_ini_back);
  fftw_destroy_plan(plan_ini_back);

#pragma omp parallel for private(j)
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j) 
    {
      ne[i][j] += ppx[i][j][0]*delinv;
      ni[i][j] += ppx[i][j][0]*delinv;
    }

  if (vorfree==1.) 
  { 
    poisson(ne,cvfis,ni);
  }

} // init_add_turb 関数の終了
