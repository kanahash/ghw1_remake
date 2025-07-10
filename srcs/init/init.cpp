#include "ghw1.hpp"

// FFTWプランとウィズダム文字列を初期化する関数
void init_FFTW(int nx, int ny, unsigned int fflag, int npar,
                    fftw_plan& hinp, fftw_plan& herp, char*& wisdom_sf)
{
  // 並列化されたFFTWを初期化する:
	#ifdef _OPENMP
  	fftw_plan_with_nthreads(npar);
	#endif

  	int nyp = ny / 2 + 1;
  // wwk と wwx は関数内で一時的に使用されるため、ローカル変数として宣言
  	fftw_complex wwk[nx][nyp];
  	double wwx[nx][ny];

  // プランの作成
  	hinp = fftw_plan_dft_r2c_2d(nx, ny, *wwx, *wwk, fflag);
  	herp = fftw_plan_dft_c2r_2d(nx, ny, *wwk, *wwx, fflag);

  // ウィズダム文字列のエクスポート
  	wisdom_sf = fftw_export_wisdom_to_string();
}

// エネルギーを初期化する関数
void init_energy(AXY ne, double xyz, double& eno)
{
    eno = 0.; // エネルギーの初期化
    int i;
    int j;
    
    // OpenMP を使って eno を並列計算
    #pragma omp parallel for private(i,j) reduction (+:eno)
    for (i = 0; i <= nx1; i++)
	{
        for (j = 0; j <= ny1; j++)
            eno += exp(2. * ne[i][j]);
    }
    eno *= .5 * xyz; // 最終的なエネルギー値の計算
}

// 時間ステップを初期化する関数
// 戻り値: シミュレーションを続行すべきか (0以外) または終了すべきか (0)
int init_time_stepping(int itstp, int itmax, double dt, double cf,
                           int& itn, int& itend, double& dtt, double& ddtt,
                           int& itrace, int& jtrace)
{
    printf("| 時間ステップを開始します ...\n");
    itn = 0; // イテレーションカウンターを初期化

    // 終了イテレーション数を計算
    if (itstp != 0)
        itend = itmax / itstp;
    else
        itend = 0; // itstp が 0 の場合、itend も 0 に設定

    // 時間ステップ関連の変数を計算
    dtt = dt * cf;
    ddtt = static_cast<double>(itstp) * dt; // static_cast で明示的に double に変換

    // トレース関連のカウンターを初期化
    itrace = 0;
    jtrace = 0;

    // itmax が 0 の場合、シミュレーションは実行されないため、ここで終了
    if (itmax == 0)
        return 0; // 終了を示す
    return 1; // 続行を示す
}

// 単一ブロブ摂動を初期化する:
void init_add_blob( AXY ne, AXY ni )
{
  int i, j, ixloc;
  double dr = wsrc*wsrc;

  printf("| ブロブ初期条件...\n");
  ixloc = nxh;
#pragma omp parallel for private(i,j)
  for (i=0; i<=nx1; ++i) 
  {
    for (j=0; j<=ny1; ++j) 
    {
      ne[i][j] += amp*( exp(-(i-ixloc)*(i-ixloc)/dr)*exp(-(j-nyh)*(j-nyh)/dr) );
      ni[i][j] = ne[i][j];
    }
  }   
  if (vorfree==1.) { poisson(ne,cvfis,ni); }
} // init_add_blob 関数の終了

// 二重渦摂動を初期化する:
void init_add_dual( AXY ne, AXY ni )
{
  int i, j, ixloc;
  double dr = wsrc*wsrc;
  double elong = 1.0; // 例: elong=1 で円形
  double dry = elong*dr;
  double xxx1, xxx2, vorset;
  printf("| 合体初期条件...\n");

  vorset = 0.01; // ジャイロ密度差から渦度を定義するため
  g_n = 0.; // 背景密度勾配なし
  ixloc = nxh;
#pragma omp parallel for private(i,j,pdist,xxx1,xxx2)
  for (i=0; i<=nx1; ++i)
  {
    for (j=0; j<=ny1; ++j)
    {
      pdist = 4.0*wsrc; // ここで距離を設定、例: 4*wsrc
      xxx1 = (i-ixloc -0.5*pdist);
      xxx2 = (i-ixloc +0.5*pdist);
      ne[i][j] += +amp*( exp(-xxx1*xxx1/dr)*exp(-(j-ny1*1/2)*(j-ny1*1/2)/dry) );
      ne[i][j] += +amp*( exp(-xxx2*xxx2/dr)*exp(-(j-ny1*1/2)*(j-ny1*1/2)/dry) );
    }
  }
  for (i=0; i<=nx1; ++i)
  {
    for (j=0; j<=ny1; ++j)
        ni[i][j] = (1.-vorset)*aai*ne[i][j];
  }

} // init_add_dual 関数の終了

// シアー流摂動を初期化する:
void init_add_flow( AXY ne, AXY ni )
{
  int i, j, ixloc = nxh;
  printf("| シアー流初期条件...\n");

  g_n = 0.; 
  chat = 0.; // 背景密度勾配なし、HW 駆動なし
#pragma omp parallel for private(i,j)
  for (i=0; i<=nx1; ++i)
  {
    for (j=0; j<=ny1; ++j)
    {
      ne[i][j] = amp*cos(TwoPi*1.*double(i-nxh)/double(nx));
      ne[i][j]+= amp*0.001*sin(TwoPi*5.*double(j)/double(ny));
      ni[i][j] = 0.;
    }
  }
} // init_add_flow 関数の終了

// 線形単一 ky モードを初期化する
void init_add_mode( AXY ne, AXY ni )
{
  int i, j;
  printf("| 単一 ky モード初期条件...\n");

  // ここで par[6]=sigma はモード数 (rho_s ky) である。
  // Ly は周期性 (初期化で実行される; m=1, ...) を得るために m*TwoPi*sigma に設定される！
  clin = 0.;
#pragma omp parallel for private(i,j)
  for (i=0; i<=nx1; ++i) for (j=0; j<=ny1; ++j)
   {
      ne[i][j] = amp*sin(sigma*double(j)/hy);
      // ne[i][j]*= sin(TwoPi*0.5*double(i)/double(nx)); // (kx!=0 の場合)
      ni[i][j] = 0.;
    }

  // kx=0 の解析的分散関係を書き出す (cf. Camargo 95)
  // "diagnose" で計算された周波数および成長率と比較するため:
  /*
  double ky, a_omr, a_gam, a_sig, a_lam0, a_lam1, a_aaa, a_bbb;
  f = fopen("formula-gamma.dat", "w");
  for (i=1; i<=1000; ++i) {
    ky = 0.01*double(i);
    a_lam0 = chat*(1.+ky*ky)/(ky*ky);
    a_lam1 = 2.*diff*pow(ky,4.);
    a_sig = chat/ky;
    a_aaa = a_lam0/(sqrt(2.)*2.);
    a_bbb = sqrt(1.+16.*a_sig*a_sig/(a_lam0*a_lam0*a_lam0*a_lam0));
    a_omr = a_aaa*sqrt(a_bbb - 1.);
    a_gam = - .5*(a_lam0+a_lam1) + a_aaa*sqrt(a_bbb + 1.);
    fprintf(f,"%.3e  %.3e  %.3e \n",ky,a_omr,a_gam);
  }
  fclose(f);
  */
} // init_add_mode 関数の終了
