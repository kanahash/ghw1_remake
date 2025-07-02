#include "ghw1.hpp"

void init_parameters(void)
{
  int i,j, ith; char s[80]; FILE *g, *f;
  double rel, para[50];

  for (int i=0;i<=40;i++) printf("_");
  printf("\n| GHW1 開始...\n");

  printf("| 初期化中 ... \n");

  // ghw1.inp 内のすべての "=" 記号の後に出現するすべての数値を順次読み取る
  // 注意: ghw1.inp に追加する場合、ナンバリングが変更される可能性があることに注意してください
  g = fopen( "ghw1.inp", "r" );
  i = 0;
  while (!feof(g)) { fscanf(g, "%s ", s); if (strstr(s,"=")) i++; };
  rewind(g);
  for (j=1; j<=i; j++) {
    while (!strstr(s,"=")) fscanf(g,"%s ",s);
    fscanf(g,"%s ",s); para[j]=atof(s);
  };
  fclose( g );

  chat  = para[1];   // 非断熱 HW カップリングパラメータ (最適な選択: 0.01 < chat < 10)
                     // (Numata et al. の "alpha"、Camargo et al. の C に相当)
  b_mhw = true;
  if (chat<0.)
  {
    chat = fabs(chat); 
    b_mhw = false;
    };  // 通常の Hasegawa-Wakatani

  taui  = para[2];   // イオン温度比 taui = T_i / ( Z_i * T_e )
  xlwl = 1.; 
  if (taui<0.) 
    {
         taui = fabs(taui); xlwl = 0.; 
    };

  delta = para[3];   // ドリフトスケール rho_s / Lperp (Lperp = (d/dx) ln no)
                     // (例: delta = 1/64 = 0.0156)
                     // 時間スケールを t*cs/Lperp に設定する (Camargo et al. と同様)
                     // (Numata et al. は t*omega_ci 規格化を使用するため、delta==1)

  diff  = para[4];   // 垂直ハイパー粘性の強度、グリッド解像度によって非自明にスケールする
                     // ガイドライン: 解像度を2倍にすると、diff を8で割る。

  incon = para[5];   // 初期条件: ブロブ、乱流、渦、シアー流から選択
                     // (乱流はブロブから発達するが、「乱流」浴からの方が速い)

  sigma = para[6];   // 初期ブロブまたは渦の幅

  amp   = para[7];   // 初期摂動振幅 (ブロブ、渦、または「乱流」浴)

  vorfree = para[8]; // 渦度自由ブロブ (特殊な研究用、通常は関連しない)

  ly = para[9];       // y ドメインの長さ (rho_s 単位)
  if (incon==5) 
    ly = 1.*TwoPi/sigma;

  nx = int(para[10]); // x 方向のグリッド点数
  ny = int(para[11]); // y 方向のグリッド点数

  hy = double(ny)/ly; // 有限差分因子、正方形グリッドを仮定
  wsrc = sigma*double(ny)/ly;
  aspect = double(nx)/double(ny); // x-y ボックスのアスペクト比 (浴の初期化のみで使用)
  clin = 1.;          // ExB 非線形性の場合 1、線形化の場合 0
  delinv = 1./delta;

  dt    = para[12];       // 時間ステップ
  itstp = int(para[13]);  // 出力間の内部ステップ
  itmax = int(para[14]);  // 終了までの合計ステップ
  ntrace= int(para[15]);  // 高速トレース出力間のステップ
  maxtrace = int(para[14]/para[15]); // 記録されるトレースステップ数
  if (ntrace==0) 
    maxtrace = 1;


  ipade = int(para[16]);  // ジャイロ演算子に対する Pade 近似またはベッセル関数
  ihype = int(para[17]);  // ハイパー粘性の次数と種類 (p=2 が通常の選択)
  // fflag = FFTW_ESTIMATE;  // FFTW プランの推定または測定 -- グローバル定義に移動
  if (para[18]==1.) 
    fflag = FFTW_MEASURE;
  if (para[18]==2.) 
    fflag = FFTW_PATIENT;
  npar  = int(para[19]);  // 並列 OpenMP スレッド数 (利用可能なコア数による)

  b_4th = true; 
  if (para[20] == 0.) 
    b_4th = false; // 2次または4次 Arakawa スキーム

  // プラズマ種パラメータ
  // (コードに種を追加する場合のみ変更する)
  taue= -1.;              // 固定電子温度 tau_e = T_e / ( Z_e * T_e ) == -1
  mue =  0.;              // 一般に mue = m_e/(Z_e*m_i) は非常に小さい、ゼロに設定可能
  // mue =- 1./(2.*1836);    // 重水素主イオン種に対する電子質量比
  mui = +1.;              // 固定イオン質量比 mui = m_i/(Z_i*m_i) == +1
  aae = -1.;              // 固定電子密度比
  aai = +1.;              // 固定イオン密度比
  zze = -1.;              // 固定電子電荷
  zzi = mui/fabs(mui); if (mui==0.) zzi=0.; // 固定イオン電荷

  // Open_MP の初期化 ...................
#ifdef _OPENMP
  int np, kp, npmax;
  ith = fftw_init_threads();
  if (ith==0)
  {
     printf("| ith=%d : スレッド初期化に失敗しました... 終了します\n",ith);
     exit(1);
    }
  npmax = omp_get_max_threads();
  if (npmax<npar)
    { 
        printf("| npar = %d > %d の最大利用可能スレッド... 終了します\n",npar,npmax); 
        exit(1);
    }
  omp_set_num_threads(npar);
  printf("| 並列スレッド ");
#pragma omp parallel private(kp)
  {
    np = omp_get_num_threads();
    kp = omp_get_thread_num();
    printf("%d/",kp);
  }
  printf(" のうち %d がアクティブ。\n",np);
  double dnp=double(np);
#else
  printf("| シングルコア処理\n");
  double dnp=1.;
#pragma omp barrier
#endif

  // パラメータの初期化 ..................
    printtraces = (ntrace > 0) ? true : false;
    g_n = 1.;    // 背景密度: g_n ~ 逆勾配長 1/Lp (= Numata の kappa)
    if (chat == 0.) g_n = 0.;

    init_output_and_arrays(nx, ny, incon, printtraces);
  
  // FFTW ポアソンソルバーの係数事前計算
    precompute_fftw_coefficients(nx, ny, nxh, nyh, hy, taui, mue, taue, mui, aae, aai, ipade, xlwl);

  // ベッセル出力のテスト: Dorland & Hammett 1993 の Fig.3 と比較する。
  /*
  g = fopen( "bessel.dat", "w" );
  for (int i=0; i<=nxh; ++i)
    {
      ik = (i>=nxh) ? i - nx : i;
      kkqq = (TwoPi*hy*ik*cnx)*(TwoPi*hy*ik*cnx);
      bki = taui*mui*kkqq;
      bke = taue*mue*kkqq;
      gambesi = gamma0(bki);
      gambese = gamma0(bke);

      fprintf(g,"%.8e  %.8e  %.8e  %.8e  %.8e  %.8e\n",sqrt(kkqq),
              1./(1.+.5*taui*mui*kkqq), sqrt( gambesi ), (kkqq*(aai*mui/(1.+taui*mui*kkqq)) + 1.*(aae*mue/(1.+bke))), ((aai/taui)*(1.-gambesi) + 1.*(aae/taue)*(1.-gambese)), gambesi );
    }
  fclose( g );
  // */
} // init_parameters 関数の終了
