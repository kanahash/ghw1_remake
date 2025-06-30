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
    
    // OpenMP を使って eno を並列計算
    #pragma omp parallel for private(i,j) reduction (+:eno)
    for (int i = 0; i <= nx1; i++)
	{
        for (int j = 0; j <= ny1; j++)
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
