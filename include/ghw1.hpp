#ifndef GHW1_HPP // ヘッダーガード: GHW1_H が定義されていなければ、以下の内容を処理する
#define GHW1_HPP // GHW1_H を定義して、二重インクルードを防ぐ

// ライブラリのインクルード: 標準的なLinuxディストリビューションで利用可能 (テスト済み: OpenSuse)
#ifdef _OPENMP // OpenMP が有効な場合
#include <omp.h> // OpenMP ライブラリをインクルード
#endif

// maxtraceの最大値を仮に MAX_TRACE_STEPS として定義
#define MAX_TRACE_K_MAX_SQRT 8 // kt_maxsqrt
#define MAX_TRACE_STEPS 1000 // maxtrace の上限値に合わせて設定

// MAX_TRACE_K_MAX_SQRT と MAX_TRACE_STEPS は、collect_time_trace_data の定義で使用したものを再利用
#define MAX_TRACE_K_MAX_SQRT 8
#define MAX_TRACE_STEPS 1000

// trace_n と trace_p はグローバル変数であると仮定（または main で宣言し、ポインタで渡す）
static double trace_n_global[MAX_TRACE_K_MAX_SQRT * MAX_TRACE_K_MAX_SQRT][MAX_TRACE_STEPS];
static double trace_p_global[MAX_TRACE_K_MAX_SQRT * MAX_TRACE_K_MAX_SQRT][MAX_TRACE_STEPS];


#include <iostream> // C++の入出力ストリーム (例: cout, cin)
#include <cstdio>   // C++スタイルの C標準入出力ヘッダー (stdio.h に対応)
#include <cstdlib>  // C++スタイルの C標準ユーティリティヘッダー (stdlib.h に対応)
#include <cmath>    // C++スタイルの C標準数学関数ヘッダー (math.h に対応)
#include <cstring>  // C++スタイルの C文字列関数ヘッダー (string.h に対応)
#include <fftw3.h>  // FFTW (Fastest Fourier Transform in the West) ライブラリ
#include <stdbool.h>

// 静的定数の定義:
const static double TwoPi = 2. * M_PI, r12 = 1. / 12.; // 円周率の2倍と1/12

// nxm, nym: 配列を初期化するための最大 nx および ny サイズ:
const static int nxm = 512 + 2, nym = 512 + 2; // nx + 2 = 2^n + 2 はFFTに理想的

// 時間トレース記録点の数 (例: 各x軸とy軸あたり8点):
const static int kt_maxsqrt = 8;

// グローバル変数として宣言された ffalg。その定義（および初期化）は .cpp ファイルにある必要があります。
extern unsigned int fflag;

typedef double AXY[nxm][nym]; // 2D配列 AXY を double 型の nxm x nym 配列として定義

// 'extern' で宣言されたグローバル変数。これらは1つの .cpp ファイルで定義されます。
extern int nx, ny, nxh, nyh, nx1, ny1, ict, jct;
extern int itstp, itmax, npar, ntrace, jnull_old;
extern int ipade, ihype;
extern int maxtrace;

extern double p_amp, p_old, t_amp, t_old;
extern double t_1, t_2;
extern int jsgn, jsgn_old;

extern double pkxavg[nxm / 2 + 1], pkyavg[nym / 2 + 1];
extern double pkyavgn[nxm / 2 + 1], pkyavgw[nym / 2 + 1], pkyavge[nym / 2 + 1];

extern double hy, hy2, hysq, xyz, dr, amp, wsrc, sigma, xlwl;
extern double diff, n00, ly, dt, dtt, ddtt, incon;
extern double chat, diss, neavg, niavg, peavg, vorfree;
extern double aspect, eno, enwo, delta, delinv, g_n;
extern double pdist, vorset, clin, vde_old, vde_avg;
extern double aae, aai, mue, mui, taue, taui, zze, zzi;

// 静電ポテンシャル pe, ジャイロ平均された pi, ww:
extern AXY pe, pi, ww;

// 電子およびイオンのジャイロ中心密度 ne および ni, および保存された値:
extern AXY nn, ne, ni, ne0, ne1, ne2, ni0, ni1, ni2;

// ジャイロ遮蔽されたイオン密度:
extern AXY gni;

// ne および ni 連続の式を定義する項
extern AXY fne0, fne1, fne2, fni0, fni1, fni2;

// ラプラシアンを保存するための配列 (例: 粘性用)
extern AXY vis_e, vis_i, hyve, hyvi;

// FFT カーネル
extern AXY cvort, cpoti, clap, chyv, cvfis;

extern char *wisdom_sf;

extern bool printtraces, b_mhw, b_4th;

// 関数プロトタイプ
void init_parameters(void);
void init_add_blob(AXY ne, AXY ni);
void init_add_turb(AXY ne, AXY ni);
void init_add_dual(AXY ne, AXY ni);
void init_add_flow(AXY ne, AXY ni);
void init_add_mode(AXY ne, AXY ni);

void arakawa(AXY uuu, AXY vvv, AXY www);
void arakaw4(AXY uuu, AXY vvv, AXY www);
void poisson(AXY pe, AXY cp, AXY pz);
void laplace(AXY pp, AXY ww);
void diagnose(double ttt, double t00, int it, int itmax,
              AXY nn, AXY ne, AXY ni, AXY pp, AXY gni);

// グローバル変数として定義されていると仮定
extern double trace_n[/* k_max の最大値に対応するサイズ */][/* maxtrace の最大値に対応するサイズ */];
extern double trace_p[/* k_max の最大値に対応するサイズ */][/* maxtrace の最大値に対応するサイズ */];
extern const int kt_maxsqrt; // ghw1.hpp で定義されていると仮定

// インライン関数定義 (定義はここにあり、プロトタイプだけではありません)
inline double DFDY(AXY ain, int i, int j)
{
    int jm, jp, jm2, jp2;
    nextp(j, 0, ny1, jm, jp, jm2, jp2);
    if (b_4th)
    {
        return (hy * r12 * (8. * (ain[i][jp] - ain[i][jm]) + ain[i][jm2] - ain[i][jp2]));
    }
    else
    {
        return (hy * .5 * (ain[i][jp] - ain[i][jm]));
    };
}

inline double DFDX(AXY ain, int i, int j)
{
    int im, ip, im2, ip2;
    nextp(i, 0, nx1, im, ip, im2, ip2);
    if (b_4th)
    {
        return (hy * r12 * (8. * (ain[ip][j] - ain[im][j]) + ain[im2][j] - ain[ip2][j]));
    }
    else
    {
        return (hy * .5 * (ain[ip][j] - ain[im][j]));
    };
}

inline void nextp(int i, int i0, int i1, int &im, int &ip, int &im2, int &ip2)
{
    im = (i == i0) ? i1 : i - 1;
    im2 = (i == i0 + 1) ? i1 : i - 2;
    if (i == i0)
        im2 = i1 - 1;
    ip = (i == i1) ? i0 : i + 1;
    ip2 = (i == i1 - 1) ? i0 : i + 2;
    if (i == i1)
        ip2 = i0 + 1;
}

inline void f_copy2darray(AXY arrinp, AXY arrout, int iend, int jend)
{
#pragma omp parallel for
    for (int i = 0; i <= iend; ++i)
        for (int j = 0; j <= jend; ++j)
            arrout[i][j] = arrinp[i][j];
}

inline double gamma0(double bkk)
// ジャイロ流体 Gamma0 演算子を評価する: Gamma0 = exp(-b)*I0(b):
// Olivares et al, Journal of Physics: Conference series 1043(2018) 012003 の
// 修正ベッセル関数 I0 と cosh 近似に基づく。
// (ガンマ演算子で cosh を exp 形式で記述することで、範囲エラーを回避する！)
{
    return 0.5 * (1. + exp(-2. * bkk)) * (1. + 0.24273 * bkk * bkk) / ((1. + 0.43023 * bkk * bkk) * std::pow((1. + 0.25 * bkk * bkk), 0.25));
}

inline double timer_start()
{ // 最初のOpenMP時間を計測する
#ifdef _OPENMP
    t_1 = omp_get_wtime();
#endif
    return t_1;
}

inline double timer_stop(double t_1)
{ // 2番目のOpenMP時間を計測し、時間差を返す
#ifdef _OPENMP
    t_2 = omp_get_wtime();
#endif
    return (t_2 - t_1);
}

typedef struct {
    // ... 既存のパラメータ ...

    // 新しく追加する派生パラメータ
    // bool printtraces; // 必要であれば追加
    int nxh;
    int nyh;
    int nx1;
    int ny1;
    double hy2;
    double hysq;
    double g_n;
    double xyz;

    // 平均値配列は、構造体に入れるのではなく、別途動的に確保してポインタを渡すのが一般的
    // double* pkxavg;
    // double* pkyavg;
    // ...など
    // ただし、これらの配列はシミュレーションのメインループで継続的に使われるため、
    // ここで初期化するよりも、メインのシミュレーションデータ構造の一部として管理する方が適切かもしれません。

} HWParameters;

//init
//init.potential.cpp
void init_FFTW(int nx, int ny, unsigned int fflag, int npar,
                    fftw_plan& hinp, fftw_plan& herp, char*& wisdom_sf);
void init_energy(AXY ne, double xyz, double& eno);
int init_time_stepping(int itstp, int itmax, double dt, double cf,
                           int& itn, int& itend, double& dtt, double& ddtt,
                           int& itrace, int& jtrace);
void initialize_output_files(int incon_param);
void initialize_average_arrays(int nx_param, int ny_param, double* pkxavg, double* pkyavg, double* pkyavgn, double* pkyavgw, double* pkyavge);

//ghw1.cpp
//main.cpp

//setting
//set_init_density.cpp
void set_init_density_perturbation(int incon, AXY ne, AXY ni);
//set_time.cpp
void init_history_densities(AXY ne, AXY ni, 
                                AXY& ne0, AXY& ne1, AXY& ne2,
                                AXY& ni0, AXY& ni1, AXY& ni2,
                                AXY& fne1, AXY& fne2, AXY& fni1, AXY& fni2) ;
//set_trace.cpp
void set_trace_output_locations(int nx1, int ny1, int kt_maxsqrt, 
                               int iout[], int jout[], int& k_max);

//calculation
//init_potential.cpp
void init_potentials(AXY ni, AXY ne, AXY& pe, AXY& gni, AXY& ww, 
                                const AXY cpoti, const AXY cvort,
                                double aae, double aai);
//density_time_evolution_step.cpp
void update_densities_one_time_step(double ttt, double dt, AXY ne, AXY ni, AXY pe, AXY pi, AXY ane, AXY ani,
                                     AXY hyve, AXY hyvi, AXY vis_e, AXY vis_i,
                                     fftw_plan hinp, fftw_plan herp, double wwx[][nym], double wwk[][nym/2+1][2],
                                     const AXY chyv_coeff, bool use_4th_order, int ihype_mode, int nx_max, int ny_max, int ny_padded);

//gyro_averate_calculate.cpp
void gyro_averaged_ion_density_calculate(AXY ni_in, AXY gni_out, fftw_plan hinp, fftw_plan herp,
                                         double wwx[][nym], double wwk[][nym/2+1][2],
                                         const AXY cpoti_coeff, double taui_param,
                                         int nx_max, int ny_max, int ny_padded);

//solve_polarization_equation.cpp
void solve_polarization_equation(AXY ne_in, AXY gni_in, AXY pe_out,
                                 fftw_plan hinp, fftw_plan herp,
                                 double wwx[][nym], double wwk[][nym/2+1][2],
                                 const AXY cvort_coeff, double zze_param, double zzi_param,
                                 int nx_max, int ny_max, int ny_padded);

//gyro_shielded_potential_calculate.cpp
void calculate_gyro_shielded_potential(AXY pe_in, AXY pi_out,
                                       fftw_plan hinp, fftw_plan herp,
                                       double wwx[][nym], double wwk[][nym/2+1][2],
                                       const AXY cpoti_coeff, double taui_param,
                                       int nx_max, int ny_max, int ny_padded);

//FFT.cpp
void perform_time_trace_analysis(int maxtrace, int k_max, double t00, double dt, int itmax, int ntrace);

//derived_parameters.cpp
void derived_parameters_calculate(HWParameters* params);

//utils
//utils_time.cpp
void run_time_loop();
//utils_history.cpp
void update_history_densities(AXY ne, AXY ne0, AXY ne1, AXY ne2,
                              AXY ni, AXY ni0, AXY ni1, AXY ni2,
                              AXY fne0, AXY fne1, AXY fne2,
                              AXY fni0, AXY fni1, AXY fni2,
                              int nx_max, int ny_max);
//utils_data.cpp
void collect_time_trace_data(double ttt_in, AXY ne_in, AXY pe_in,
                             int iout_loc[], int jout_loc[], int k_max_loc,
                             int& current_itrace, int& current_jtrace, int ntrace_interval,
                             double& td_out_total, double& t_start_measure,
                             bool printtraces_flag);

int	main(void);

#endif // GHW1_H
