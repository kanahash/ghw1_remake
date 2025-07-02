#include "ghw1.hpp"

// ghw1.hpp で 'extern' 宣言されたグローバル変数の定義
// 初期化子を持つ変数はここで初期化します。
// 配列もここで定義します。
unsigned int fflag = FFTW_ESTIMATE; // ghw1.hpp から初期化子を移動

int nx, ny, nxh, nyh, nx1, ny1, ict, jct;
int itstp, itmax, npar, ntrace, jnull_old;
int ipade, ihype;
int maxtrace;

double p_amp = 0., p_old = 0., t_amp = 0., t_old = 0.; // ghw1.hpp から初期化子を移動
double t_1, t_2; // ghw1.hpp から初期化子を削除
int jsgn, jsgn_old;

double pkxavg[nxm/2+1];
double pkyavg[nym/2+1];
double pkyavgn[nxm/2+1];
double pkyavgw[nym/2+1];
double pkyavge[nym/2+1];

double hy, hy2, hysq, xyz, dr, amp, wsrc, sigma, xlwl;
double diff, n00, ly, dt, dtt, ddtt, incon;
double chat, diss, neavg, niavg, peavg, vorfree;
double aspect, eno, enwo, delta, delinv, g_n;
double pdist, vorset, clin, vde_old, vde_avg;
double aae, aai,  mue, mui, taue, taui, zze, zzi;

// 静電ポテンシャル pe、ジャイロ平均 pi、ww:
AXY pe, pi, ww;

// 電子およびイオンのジャイロ中心密度 ne と ni、および格納された値:
AXY nn, ne, ni, ne0, ne1, ne2, ni0, ni1, ni2;

// ジャイロ遮蔽イオン密度:
AXY gni;

// ne および ni 連続方程式の右辺を定義する項
AXY fne0, fne1, fne2, fni0, fni1, fni2;

// ラプラシアンを格納する配列 (例: 粘性用)
AXY vis_e, vis_i, hyve, hyvi;

// FFT カーネル
AXY cvort, cpoti, clap, chyv, cvfis;

char *wisdom_sf;

bool printtraces, b_mhw, b_4th;


// 関数の指定:

// ---------------------------------------------------------------------------
// 入力データファイル "ghw1.inp" を読み込み、さらにいくつかのパラメータを初期化する



// ---------------------------------------------------------------------------



// ---------------------------------------------------------------------------



// ---------------------------------------------------------------------------



// ---------------------------------------------------------------------------



// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// 出力量を計算し書き出す:
// 例: エネルギー、輸送、k スペクトル、2D 配列、プロファイルなど
void diagnose( double ttt, double t00, int it, int itmax,
           AXY n_n, AXY n_e, AXY n_i, AXY p_p, AXY n_g )
{
  int i,j, ik, im,ip,jm,jp;
  double enn, enp, enw, eeb, ezf, fne, fnsol, zfx, dum, rey, freq, grow;
  double etran_grad, etran_adia, etran_visc;
  double lap[nxm][nym];
  double avg[nxm], avp[nxm];
  FILE *f, *g, *h, *g1, *g2, g3, *g4, *g5, *g6, *g7, *g8;

#pragma omp parallel for private(j) shared(avg,avp)
  for (i=0; i<=nx1; i++) // ne と phi の帯状平均
    { avg[i] = 0.; avp[i] = 0.;
      for (j=0; j<=ny1; j++) {avg[i] += n_e[i][j]; avp[i] += p_p[i][j];}
      avg[i]/=ny; avp[i]/=ny;
    }
  enn = 0.; enp = 0.; enw = 0.; eeb = 0.; ezf = 0.; fne = 0.;
  etran_grad = 0.; etran_adia= 0.; etran_visc = 0.;

  // 渦度
  // laplace(p_p,lap);
  poisson(p_p,clap,lap);

  // エネルギー量と輸送
#pragma omp parallel for private(j,ip,dum) reduction(+:enn,ezf,eeb,fne,enw,etran_grad,etran_adia,etran_visc)
  for (i=0; i<=nx1; i++) {
    ip = (i==nx1) ? 0   : i+1;

    for (j=0; j<=ny1; j++) {

      // 全熱自由エネルギー:
      enn += n_e[i][j]*n_e[i][j] + taui*n_i[i][j]*n_i[i][j];

      // (主に帯状) 流れのエネルギー:
      dum = hy*(avp[ip]-avp[i]);
      ezf += dum*dum;

      // 乱流エントロフィー:
      enw += lap[i][j]*lap[i][j];

      // 全運動エネルギー:
      dum = DFDY(p_p,i,j);
      eeb += dum*dum;
      dum = DFDX(p_p,i,j);
      eeb += dum*dum;

      // 電子 ExB 粒子輸送:
      fne += - n_e[i][j]*DFDY(p_p,i,j);

      // 勾配駆動 / 転送:
      etran_grad+= - n_e[i][j]*DFDY(pe,i,j);
      etran_grad+= - n_i[i][j]*DFDY(pi,i][j])*taui;

      // 「並列」抵抗散逸シンク / 転送:
      dum = pe[i][j] - n_e[i][j];
      if (b_mhw) dum-= avp[i] - avg[i];
      etran_adia+= -chat*dum*dum;

      // 垂直ハイパー粘性シンク / 転送:
      dum = hyve[i][j]*(pe[i][j] - n_e[i][j]);
      dum-= hyvi[i][j]*(pi[i][j] + taui*n_i[i][j]);
      etran_visc+= diff*dum;

    }
  }

  enn *= .5*xyz; enp *= xyz; enw *= .5*xyz; eeb *= .5*xyz; ezf *= .5*xyz;
  if (enw==0.) enw=1.e-12; if (ezf==0.) ezf=1.e-12; if (eeb==0.) eeb=1.e-12;

  eeb*=delinv*delinv; ezf*=delinv*delinv; enw*=delinv*delinv; enn*=delinv*delinv;
  fne*= xyz*delinv*delinv;

  etran_grad*= xyz; etran_adia*= xyz; etran_visc*= xyz;
  double etran_tot = etran_grad + etran_adia + etran_visc;

  // クラッシュ制御: nan または inf の場合停止
  // if ( (isnan(enn)) || (enn<1.e-16) )
  if (isnan(enn))
    {
      printf("\n|t=%.3f (%.2f): enn=%.3e\n", ttt,(t00+itmax*dt),enn);
      printf("|HGW1 クラッシュしました！  enn=%.3e,   eno=%.3e\n", enn, eno);
#ifdef _OPENMP
      fftw_cleanup_threads();
#endif
      exit(1);
    }
  eno = enn;

  // グローバルエネルギー時系列出力:
  f = fopen("eng.dat","a");
  if (ttt>0.) fprintf(f,"%.8f  %.6e  %.6e  %.6e  %.6e  %.6e\n",
              ttt, enn, eeb, ezf, (enn+eeb), ezf/eeb);
  fclose(f);
  // グローバルエネルギー転送時系列出力:
  f = fopen("eng-transfer.dat","a");
  if (ttt>0.) fprintf(f,"%.8f  %.6e  %.6e  %.6e  %.6e\n",
              ttt, etran_grad, etran_adia, etran_visc, etran_tot);
  fclose(f);
  // グローバル輸送時系列出力:
  f = fopen("fne.dat","a");
  if (ttt>0.) fprintf(f,"%.8f  %.6e\n", ttt, fne);
  fclose(f);

  // (非)線形成長率 "grow" と線形周波数推定 "freq":
  if (incon==5.) {
    // (近似) 線形ドリフト波周波数 freq を計算する

    double dt_amp = 0.;
    freq = 0.; grow = 0.;
    jsgn_old = jsgn;
    if (p_p[nxh][nyh]==fabs(p_p[nxh][nyh])) {jsgn = 1;} else {jsgn = -1;};

    if (jsgn != jsgn_old) {t_old = t_amp; t_amp = ttt;};
    dt_amp = t_amp - t_old;
    if ( (dt_amp !=0.) && (dt_amp != ddtt) && (it>2) )  freq = 0.5*TwoPi/dt_amp;

  }

  if (enwo == 1.e-12) enwo = enw;
  grow = (sqrt(enw)-sqrt(enwo))/(sqrt(enwo)*ddtt);
  enwo = enw;
  f = fopen("gamma.dat","a");
  if (it>1) fprintf(f,"%.3f  %.6e  %.6e\n", ttt, freq, grow);
  fclose(f);

  // プロット出力

  // 帯状流 (t,x) 2D プロット
  h = fopen( "zfx.dat", "a" );
  for (i=0; i<=nx1; i++)  {
    im = (i==0)   ? nx1 : i-1;
    ip = (i==nx1) ? 0   : i+1;
    for (zfx=0., j=0; j<=ny1; j++)  {
      // zfx += p_p[i][j]/(ny); // 帯状ポテンシャル
      zfx += p_p[ip][j]-p_p[im][j]; // 帯状流 <Vy>
    }
    zfx *= .5*hy/(ny); // pot または vor. の場合はコメントアウト
    if (ttt>0.) fprintf(h,"%.3f  %d  %.6e \n", ttt,i,zfx);
  }
  fprintf(h,"\n");
  fclose(h);

  // y=ny/2 での x-カット:
  g = fopen( "cutx.dat", "w" );
  for (i=0; i<=nx1; i++)
    for (j=nyh; j<=nyh; j++)
      fprintf(g,"%f  %.6e   %.6e   %.6e   %.6e   %.6e\n",
          double(i)/hy,n_e[i][j],n_i[i][j], p_p[i][j],lap[i][j],n_g[i][j]);
  fclose( g );

  g = fopen( "cuty.dat", "w" );
  for (j=0; j<=ny1; j++)
      fprintf(g,"%f  %.6e   %.6e   %.6e   %.6e   %.6e\n",
          double(j)/hy,n_e[nxh][j],n_i[nxh][j], p_p[nxh][j],lap[nxh][j],n_g[nxh][j]);
  fclose( g );

  // y で平均された x プロファイル:
  g = fopen( "xprof.dat", "w" );
  double peprof, neprof, niprof, nsprof, wwprof;
  double pexprof[nxm];
  for (i=0; i<=nx1; i++) {
    peprof = 0.; neprof = 0.; niprof = 0.; nsprof = 0.; wwprof=0.; rey = 0.;
    im = (i==0)   ? nx1 : i-1;
    ip = (i==nx1) ? 0   : i+1;
    for (j=0; j<=ny1; j++) {
      jm = (j==0)   ? ny1 : j-1;
      jp = (j==ny1) ? 0   : j+1;
      // レイノルズ応力:
      rey+= .25*hy*hy*(p_p[ip][j]-p_p[im][j])*(p_p[i][jp]-p_p[i][jm]);
      // ポテンシャル、密度、渦度
      peprof += p_p[i][j];
      neprof += n_e[i][j];
      niprof += n_i[i][j];
      wwprof += lap[i][j];
    }
    peprof/=ny; neprof/=ny; niprof/=ny; nsprof/=ny; wwprof/=ny; rey/=ny;
    fprintf(g,"%f  %.6e   %.6e   %.6e   %.6e   %.6e\n",
        double(i)/hy,peprof,neprof,niprof,wwprof,rey);
  }
  fclose( g );

  // 密度、渦度、ポテンシャルの 2D (x,y) プロット
  g1 = fopen( "n2d0.dat", "w" ); g2 = fopen( "w2d0.dat", "w" );
  g3 = fopen( "p2d0.dat", "w" ); g4 = fopen( "i2d0.dat", "w" );
  for (i=0; i<=nx1; i++) {
    for (j=0; j<=ny1; j++) {
      fprintf(g1,"%d  %d  %.6e\n",i,j,n_e[i][j]);
      fprintf(g2,"%d  %d  %.6e\n",i,j,lap[i][j]);
      fprintf(g3,"%d  %d  %.6e\n",i,j,p_p[i][j]);
      fprintf(g4,"%d  %d  %.6e\n",i,j,(p_p[i][j]-n_e[i][j]));
    }
    fprintf(g1,"\n"); fprintf(g2,"\n"); fprintf(g3,"\n");
    fprintf(g4,"\n");
  }
  fclose( g1 ); fclose( g2 ); fclose( g3 ); fclose( g4 );
  rename("n2d0.dat", "n2d.dat"); rename("w2d0.dat", "w2d.dat");
  rename("p2d0.dat", "p2d.dat"); rename("i2d0.dat", "i2d.dat");

  // フーリエ ky スペクトル
  double py[ny];
  // fftw_complex ky[ny/2], sumky[ny/2];
  fftw_complex ky[ny/2+1], sumky[ny/2+1];
#ifdef _OPENMP
  fftw_plan_with_nthreads(npar);
#endif
  fftw_plan hindfty;
  hindfty = fftw_plan_dft_r2c_1d(ny, &py[0], &ky[0], FFTW_ESTIMATE);

  // ポテンシャル phi
  for (j=0; j<=ny/2+1; j++) sumky[j][0] = 0.;
  for (i=0; i<=nx1; i++)
    {
      for (j=0; j<=ny1; j++) py[j] = p_p[i][j];
      fftw_execute(hindfty);
      // for (j=0; j<=ny/2; j++) sumky[j][0] += ky[j][0]*ky[j][0];
      for (j=0; j<=ny/2+1; j++) sumky[j][0] += fabs(ky[j][0]);
    }
  for (j=0; j<=ny/2+1; j++) sumky[j][0] /= double(nx);
  if (incon==2)
    { for (j=1; j<=ny/2; j++) pkyavg[j] = ( (it-1+1)*pkyavg[j] + sumky[j][0] ) / (it+1); }
  g = fopen( "pky_p.dat", "w" );
  for (j=1; j<ny/2; j++)
    fprintf(g,"%.6e  %.6e  %.6e\n",TwoPi*double(j)/ly,sumky[j][0]*ny/hy/hy,pkyavg[j]*ny/hy/hy);
  fclose( g );

  // 密度 ne
  for (j=0; j<=ny/2+1; j++) sumky[j][0] = 0.;
  for (i=0; i<=nx1; i++)
    {
      for (j=0; j<=ny1; j++) py[j] = ne[i][j];
      fftw_execute(hindfty);
      // for (j=0; j<=ny/2; j++) sumky[j][0] += ky[j][0]*ky[j][0];
      for (j=0; j<=ny/2+1; j++) sumky[j][0] += fabs(ky[j][0]);
    }
  for (j=0; j<=ny/2+1; j++) sumky[j][0] /= double(nx);
  if (incon==2)
    { for (j=1; j<=ny/2; j++) pkyavgn[j] = ( (it-1+1)*pkyavgn[j] + sumky[j][0] ) / (it+1); }
  g = fopen( "pky_n.dat", "w" );
  for (j=1; j<ny/2; j++)
    fprintf(g,"%.6e  %.6e  %.6e\n",TwoPi*double(j)/ly,sumky[j][0]*ny/hy/hy,pkyavgn[j]*ny/hy/hy);
  fclose( g );

  // 渦度 w
  for (j=0; j<=ny/2+1; j++) sumky[j][0] = 0.;
  for (i=0; i<=nx1; i++)
    {
      for (j=0; j<=ny1; j++) py[j] = lap[i][j];
      fftw_execute(hindfty);
      for (j=0; j<=ny/2+1; j++) sumky[j][0] += ky[j][0]*ky[j][0];
    }
  for (j=0; j<=ny/2+1; j++) sumky[j][0] /= double(nx);
  if (incon==2)
    { for (j=1; j<=ny/2; j++) pkyavgw[j] = ( (it-1+1)*pkyavgw[j] + sumky[j][0] ) / (it+1); }
  g = fopen( "pky_w.dat", "w" );
  for (j=1; j<ny/2; j++)
    fprintf(g,"%.6e  %.6e  %.6e\n",TwoPi*double(j)/ly,.5*sumky[j][0]*ny/hy/hy,.5*pkyavgw[j]*ny/hy/hy);
  fclose( g );

  // 運動エネルギー E
  for (j=0; j<=ny/2+1; j++) sumky[j][0] = 0.;
  for (i=0; i<=nx1; i++)  {
    im = (i==0)   ? nx1 : i-1;
    ip = (i==nx1) ? 0   : i+1;
    for (j=0; j<=ny1; j++) {
      jm = (j==0)   ? ny1 : j-1;
      jp = (j==ny1) ? 0   : j+1;
      dum = .5*hy*(p_p[i][jp]-p_p[i][jm]);
      py[j] = dum*dum;
      dum = .5*hy*(p_p[ip][j]-p_p[im][j]);
      py[j]+= dum*dum;
    }
    fftw_execute(hindfty);
    for (j=0; j<=ny/2+1; j++) sumky[j][0] += ky[j][0]*ky[j][0];
  }
  for (j=0; j<=ny/2+1; j++) sumky[j][0] /= (nx);
  if (incon==2)
    { for (j=1; j<=ny/2; j++) pkyavge[j] = ( (it-1+1)*pkyavge[j] + sumky[j][0] ) / (it+1); }
  g = fopen( "pky_e.dat", "w" );
  for (j=1; j<ny/2; j++)
    fprintf(g,"%.6e  %.6e  %.6e\n",TwoPi*double(j)/ly,.5*sumky[j][0]*ny/hy/hy,.5*pkyavge[j]*ny/hy/hy);
  fclose( g );
  fftw_destroy_plan(hindfty);

  // フーリエ kx スペクトル
  double px[nx];
  fftw_complex kx[nx/2+1], sumkx[nx/2+1];
  fftw_plan hindftx;
  hindftx = fftw_plan_dft_r2c_1d(nx, &px[0], &kx[0], FFTW_ESTIMATE);
  for (i=0; i<=nx/2; i++) sumkx[i][0] = 0.;
  for (j=0; j<=ny1; j++)
    {
      for (i=0; i<=nx1; i++) px[i] = pe[i][j];
      fftw_execute(hindftx);
      for (i=0; i<=nx/2+1; i++) sumkx[i][0] += fabs(kx[i][0]);
    }
  for (i=0; i<=nx/2+1; i++) sumkx[i][0] /= (ny);
  fftw_destroy_plan(hindftx);
  if (incon==2)
    { for (i=1; i<=nx/2+1; i++) pkxavg[i] = ( (it-1+1)*pkxavg[i] + sumkx[i][0] ) / (it+1); }
  g = fopen( "pkx_p.dat", "w" );
  for (i=1; i<nx/2; i++)
    fprintf(g,"%.6e  %.6e  %.6e\n",TwoPi*double(i)/ly,sumkx[i][0]*nx/hy/hy,pkxavg[i]*nx/hy/hy);
  fclose( g );

#ifdef _OPENMP
  fftw_cleanup_threads();
#endif
  // スペクトル解析の終了

  // タイムスタンプ出力
  printf("|t=%.3f (%.3f): enn=%.3e\n", ttt,(t00+itmax*dt),enn);

}

// -----------------------------------------------------------------------------
// ブラケットの 2次 Arakawa スキーム (Arakawa 1966, JCP 1997 に再録を参照)
// [u,v] = (du/dx)(dv/dy)-(du/dy)(dv/dx)
void arakawa( AXY uuu, AXY vvv, AXY out )
{
  int i00,j00,ip1,jp1,im1,jm1,im2,ip2,jm2,jp2;
  double aj2, afac = .25*hy*hy/3.;
  int i0 = 0, i1 = nx1, j0 = 0, j1 = ny1;

#pragma omp parallel for private(j00,im1,im2,jm1,jm2,ip1,ip2,jp1,jp2,aj2)
  for (i00=0; i00<=nx1; ++i00)  {
    nextp(i00,i0,i1,im1,ip1,im2,ip2);

    for (j00=0; j00<=ny1; ++j00) {
      nextp(j00,j0,j1,jm1,jp1,jm2,jp2);

      // J++ :
      aj2  = (uuu[ip1][j00] - uuu[im1][j00]) * (vvv[i00][jp1] - vvv[i00][jm1]);
      aj2 -= (uuu[i00][jp1] - uuu[i00][jm1]) * (vvv[ip1][j00] - vvv[im1][j00]);

      // J+x :
      aj2 += uuu[ip1][j00] * ( vvv[ip1][jp1] - vvv[ip1][jm1] );
      aj2 -= uuu[im1][j00] * ( vvv[im1][jp1] - vvv[im1][jm1] );
      aj2 -= uuu[i00][jp1] * ( vvv[ip1][jp1] - vvv[im1][jp1] );
      aj2 += uuu[i00][jm1] * ( vvv[ip1][jm1] - vvv[im1][jm1] );

      // Jx+ :
      aj2 += uuu[ip1][jp1] * ( vvv[i00][jp1] - vvv[ip1][j00] );
      aj2 -= uuu[im1][jm1] * ( vvv[im1][j00] - vvv[i00][jm1] );
      aj2 -= uuu[im1][jp1] * ( vvv[i00][jp1] - vvv[im1][j00] );
      aj2 += uuu[ip1][jm1] * ( vvv[ip1][j00] - vvv[i00][jm1] );

      // 2次部分:
      out[i00][j00] = aj2*afac;
    };
  }
} // arakawa 関数の終了


// -----------------------------------------------------------------------------
// ブラケットの 4次 Arakawa スキーム (Arakawa 1966, JCP 1997 に再録を参照)
// [u,v] = (du/dx)(dv/dy)-(du/dy)(dv/dx)
// void arakaw4(double (*uuu)[nym], double (*vvv)[nym], double (*out)[nym])
void arakaw4( AXY uuu, AXY vvv, AXY out )
{
  int i00,j00,ip1,jp1,im1,jm1,im2,ip2,jm2,jp2;
  double aj2, aj4, afac = .25*hy*hy/3.;
  int i0 = 0, i1 = nx1, j0 = 0, j1 = ny1;

#pragma omp parallel for private(j00,im1,im2,jm1,jm2,ip1,ip2,jp1,jp2,aj2,aj4)
  for (i00=0; i00<=nx1; ++i00)  {
    nextp(i00,i0,i1,im1,ip1,im2,ip2);

    for (j00=0; j00<=ny1; ++j00) {
      nextp(j00,j0,j1,jm1,jp1,jm2,jp2);

      // 2次部分:
      // J++ :
      aj2  = (uuu[ip1][j00] - uuu[im1][j00]) * (vvv[i00][jp1] - vvv[i00][jm1]);
      aj2 -= (uuu[i00][jp1] - uuu[i00][jm1]) * (vvv[ip1][j00] - vvv[im1][j00]);

      // J+x :
      aj2 += uuu[ip1][j00] * ( vvv[ip1][jp1] - vvv[ip1][jm1] );
      aj2 -= uuu[im1][j00] * ( vvv[im1][jp1] - vvv[im1][jm1] );
      aj2 -= uuu[i00][jp1] * ( vvv[ip1][jp1] - vvv[im1][jp1] );
      aj2 += uuu[i00][jm1] * ( vvv[ip1][jm1] - vvv[im1][jm1] );

      // Jx+ :
      aj2 += uuu[ip1][jp1] * ( vvv[i00][jp1] - vvv[ip1][j00] );
      aj2 -= uuu[im1][jm1] * ( vvv[im1][j00] - vvv[i00][jm1] );
      aj2 -= uuu[im1][jp1] * ( vvv[i00][jp1] - vvv[im1][j00] );
      aj2 += uuu[ip1][jm1] * ( vvv[ip1][j00] - vvv[i00][jm1] );

      aj2 *= afac;

      // 4次部分:
      // Jxx ;
      aj4  = (uuu[ip1][jp1] - uuu[im1][jm1]) * (vvv[im1][jp1] - vvv[ip1][jm1]);
      aj4 -= (uuu[im1][jp1] - uuu[ip1][jm1]) * (vvv[ip1][jp1] - vvv[im1][jm1]);

      // J4x+ :
      aj4 += uuu[ip1][jp1] * ( vvv[i00][jp2] - vvv[ip2][j00] );
      aj4 -= uuu[im1][jm1] * ( vvv[im2][j00] - vvv[i00][jm2] ); // Arakawa 1966 で "+" が間違っている！
      aj4 -= uuu[im1][jp1] * ( vvv[i00][jp2] - vvv[im2][j00] );
      aj4 += uuu[ip1][jm1] * ( vvv[ip2][j00] - vvv[i00][jm2] );

      // J4+x :
      aj4 += uuu[ip2][j00] * ( vvv[ip1][jp1] - vvv[ip1][jm1] );
      aj4 -= uuu[im2][j00] * ( vvv[im1][jp1] - vvv[im1][jm1] );
      aj4 -= uuu[i00][jp2] * ( vvv[ip1][jp1] - vvv[im1][jp1] );
      aj4 += uuu[i00][jm2] * ( vvv[ip1][jm1] - vvv[im1][jm1] );

      aj4 *= .5*afac;

      // 合計 (4次):
      out[i00][j00] = 2.*aj2 - aj4 ;

      // テストのためにコメント解除: 2次部分のみ:
      // out[i00][j00] = aj2;
    };
  }
} // arakaw4 関数の終了


// ---------------------------------------------------------------------------
// 4次有限差分に対する座標 i,j の 2近傍を定義する:
// NOTE: この関数の定義は ghw1.hpp にインライン関数として移動されました。
// この空のブロックは、プレースホルダーまたは以前の場所を示すために残されています。
/*
inline void nextp (int i, int i0, int i1, int& im, int& ip, int& im2, int& ip2)
{
  im = (i==i0) ? i1 : i-1;  im2 = (i==i0+1) ? i1 : i-2;  if (i==i0) im2 = i1-1;
  ip = (i==i1) ? i0 : i+1;  ip2 = (i==i1-1) ? i0 : i+2;  if (i==i1) ip2 = i0+1;
}
*/

// ---------------------------------------------------------------------------
// 4次 y 導関数:
// NOTE: この関数の定義は ghw1.hpp にインライン関数として移動されました。
/*
inline double DFDY( AXY ain, int i, int j)
{
  int jm, jp, jm2, jp2;
  nextp(j,0,ny1,jm,jp,jm2,jp2);
  if (b_4th) { return ( hy*r12*( 8.*( ain[i][jp]-ain[i][jm] ) + ain[i][jm2] - ain[i][jp2] ) ); }
  else {return ( hy*.5*( ain[i][jp]-ain[i][jm] ) ) ; };
}
*/

// ---------------------------------------------------------------------------
// 4次 x 導関数:
// NOTE: この関数の定義は ghw1.hpp にインライン関数として移動されました。
/*
inline double DFDX( AXY ain, int i, int j)
{
  int im, ip, im2, ip2;
  nextp(i,0,nx1,im,ip,im2,ip2);
  if (b_4th) { return ( hy*r12*( 8.*( ain[ip][j]-ain[im][j] ) + ain[im2][j] - ain[ip2][j] ) ) ; }
  else {  return ( hy*.5*( ain[ip][j]-ain[im][j] ) ) ;  };
}
*/

// ---------------------------------------------------------------------------
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

// ---------------------------------------------------------------------------
// FFT "ポアソン方程式" ソルバー: del^2y = x の y を k 空間で解く
// cp で定義された一般的な引数を持つ (純粋なポアソン問題だけでなく)。
void poisson( AXY fi, AXY cp, AXY fo )
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
  for (i=i0; i<=i1; ++i) for (j=0; j<=ny1; ++j) { wwx[i][j] = fi[i][j]; }
  fftw_execute(hinpsub);

#pragma omp parallel for private(j,k)
  for (i=i0; i<=i1; ++i) for (j=0; j<=nyp-1; ++j)  {
      for (k=0; k<=1; ++k) wwk[i][j][k] *= cp[i][j];
    }
  fftw_execute(herpsub);

#pragma omp parallel for private(j)
  for (i=i0; i<=i1; ++i) for (j=0; j<=ny1; ++j) fo[i][j] = wwx[i][j];

  fftw_destroy_plan(hinpsub);
  fftw_destroy_plan(herpsub);
}

// ---------------------------------------------------------------------------
// NOTE: この関数の定義は ghw1.hpp にインライン関数として移動されました。
/*
inline void f_copy2darray( AXY arrinp, AXY arrout, int iend, int jend )
{
#pragma omp parallel for
  for (int i=0; i<=iend; ++i) for (int j=0; j<=jend; ++j) arrout[i][j] = arrinp[i][j];
}
*/

// ---------------------------------------------------------------------------
// NOTE: この関数の定義は ghw1.hpp にインライン関数として移動されました。
/*
double gamma0( double bkk)
// ジャイロ流体 Gamma0 演算子 Gamma0 = exp(-b)*I0(b) を評価する:
// Olivares et al, Journal of Physics: Conference series 1043(2018) 012003. の cosh 近似に基づく
// 修正ベッセル関数 I0 に基づく。
// (Gamma 演算子で cosh を exp 形式で記述することで、範囲エラーを回避できる！)
{
  return  0.5*(1.+exp(-2.*bkk))*(1.+0.24273*bkk*bkk)
    / ((1.+0.43023*bkk*bkk)*pow((1.+0.25*bkk*bkk),0.25));
}
*/

// ---------------------------------------------------------------------------
// NOTE: この関数の定義は ghw1.hpp にインライン関数として移動されました。
/*
inline double timer_start()
{ // 最初の omp 時間を測定する
#ifdef _OPENMP
  t_1 = omp_get_wtime();
#endif
  return t_1;
}
*/

// ---------------------------------------------------------------------------
// NOTE: この関数の定義は ghw1.hpp にインライン関数として移動されました。
/*
inline double timer_stop(double t_1)
{ // 2番目の omp 時間を測定し、時間差を返す
#ifdef _OPENMP
  t_2 = omp_get_wtime();
#endif
  return (t_2 - t_1);
}
*/

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
