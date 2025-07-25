#ifndef GHW1_HPP // ヘッダーガード: GHW1_HPP が定義されていなければ、以下の内容を処理する
#define GHW1_HPP // GHW1_HPP を定義して、二重インクルードを防ぐ

// 1. 必要なライブラリのインクルードとマクロ定義を**最初に**配置
//    FFTW や OpenMP の型が他の場所で使われるため、これらを最初に読み込む必要があります。
#ifdef __cplusplus // C++コンパイラの場合のみ extern "C" を適用
extern "C" {
#endif
#ifdef _OPENMP // OpenMP が有効な場合
#include <omp.h> // OpenMP ライブラリをインクルード
#endif
#ifdef __cplusplus
}
#endif

#include <iostream>  // C++の入出力ストリーム (例: cout, cin)
#include <cstdio>    // C++スタイルの C標準入出力ヘッダー (stdio.h に対応)
#include <cstdlib>   // C++スタイルの C標準ユーティリティヘッダー (stdlib.h に対応)
#include <cmath>     // C++スタイルの C標準数学関数ヘッダー (math.h に対応)
#include <cstring>   // C++スタイルの C文字列関数ヘッダー (string.h に対応)
#include <fftw3.h>   // ★FFTWの型 (fftw_plan, fftw_complex など) の定義のため、AXY より前に配置★
#include <stdbool.h> // bool 型のサポート (C++では不要な場合もあるが互換性のため)

// maxtraceの最大値を仮に MAX_TRACE_STEPS として定義
#define MAX_TRACE_K_MAX_SQRT 8   // kt_maxsqrt
#define MAX_TRACE_STEPS 1000     // maxtrace の上限値に合わせて設定
#define KMAX_VAL            100  // 例：k_max の最大値
#define MAXTRACE_VAL        500  // 例：maxtrace の最大値

// 2. 静的定数や typedef の定義を次に配置
//    AXY型はnxm, nymに依存するため、それらの定義の後に置きます。
const static int nxm = 512 + 2, nym = 512 + 2; // nxm, nym の定義
typedef double AXY[nxm][nym]; // ★AXY の定義を、それを使用する関数プロトタイプより前に移動★

const static double TwoPi = 2. * M_PI, r12 = 1. / 12.; // 円周率の2倍と1/12

// 時間トレース記録点の数 (例: 各x軸とy軸あたり8点):
const static int kt_maxsqrt = 8; // kt_maxsqrt は const static で定義されているので、extern宣言は不要

// 3. 'extern' で宣言されたグローバル変数
//    これらの変数は、上記で定義された型 (AXYなど) や定数 (nxmなど) を使用するため、これらの定義の後に記述します。
extern unsigned int fflag;
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

extern double trace_n_global[MAX_TRACE_K_MAX_SQRT * MAX_TRACE_K_MAX_SQRT][MAX_TRACE_STEPS];
extern double trace_p_global[MAX_TRACE_K_MAX_SQRT * MAX_TRACE_K_MAX_SQRT][MAX_TRACE_STEPS];

// 4. ★インライン関数の定義ブロックを配置★
//    これらの定義は、それが使用する全ての型（AXYなど）が定義された後に来ます。
//    また、呼び出し側のコード（main.cppなど）よりも物理的に前に来ます。
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

inline void f_copy2darray(AXY arrinp, AXY arrout, int iend, int jend)
{
#pragma omp parallel for
    for (int i = 0; i <= iend; ++i)
        for (int j = 0; j <= jend; ++j)
            arrout[i][j] = arrinp[i][j];
}

inline double gamma0(double bkk)
{
    return 0.5 * (1. + exp(-2. * bkk)) * (1. + 0.24273 * bkk * bkk) / ((1. + 0.43023 * bkk * bkk) * std::pow((1. + 0.25 * bkk * bkk), 0.25));
}

inline double timer_start()
{
#ifdef _OPENMP
    t_1 = omp_get_wtime();
#endif
    return t_1;
}

inline double timer_stop(double t_1)
{
#ifdef _OPENMP
    t_2 = omp_get_wtime();
#endif
    return (t_2 - t_1);
}

// 5. 通常の関数プロトタイプ
//    これらの関数は、上記で定義された型やインライン関数を使用できます。
// init.cpp
void init_parameters(void);
void init_add_blob(AXY ne, AXY ni);
void init_add_turb(AXY ne, AXY ni);
void init_add_dual(AXY ne, AXY ni);
void init_add_flow(AXY ne, AXY ni);
void init_add_mode(AXY ne, AXY ni);
void init_FFTW(int nx, int ny, unsigned int fflag, int npar,
               fftw_plan& hinp, fftw_plan& herp, char*& wisdom_sf);
void init_energy(AXY ne, double xyz, double& eno);
int init_time_stepping(int itstp, int itmax, double dt, double cf,
                       int& itn, int& itend, double& dtt, double& ddtt,
                       int& itrace, int& jtrace);
//init_output_and_arrays.cpp
void init_output_and_arrays(int nx, int ny, int incon, bool printtraces);
//setting
//set_init_density.cpp
void set_init_density_perturbation(int incon, AXY ne, AXY ni);
//set_time.cpp
void init_history_density(AXY ne, AXY ni, 
                          AXY& ne0, AXY& ne1, AXY& ne2,
                          AXY& ni0, AXY& ni1, AXY& ni2,
                          AXY& fne1, AXY& fne2, AXY& fni1, AXY& fni2);
//set_trace.cpp
void set_trace_output_locations(int nx1, int ny1, int kt_maxsqrt, 
                                 int iout[], int jout[], int& k_max);

//calculation
//init_potential.cpp
void init_potentials(AXY ni, AXY ne, AXY& pe, AXY& gni, AXY& ww, 
                     AXY cpoti, AXY cvort,
                     double aae, double aai);
//density_time_evolution_step.cpp
void update_densities_one_time_step(double ttt, double dt, AXY ne, AXY ni, AXY pe, AXY pi, AXY ane, AXY ani,
                                     AXY hyve, AXY hyvi, AXY vis_e, AXY vis_i,
                                     fftw_plan hinp, fftw_plan herp, double wwx[][nym], double wwk[][nym/2+1][2],
                                     const AXY chyv_coeff, bool use_4th_order, int ihype_mode, int nx_max, int ny_max, int ny_padded);

//gyro_average_calculate.cpp
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

//precompute_fftw_coefficients.cpp
void precompute_fftw_coefficients(int nx, int ny, int nxh, int nyh, double hy, double taui, double mue, double taue, double mui, double aae, double aai, int ipade, bool xlwl);

//FFT.cpp
void perform_time_trace_analysis(int maxtrace, int k_max, double t00, double dt, int itmax, int ntrace);

//energy_and_transport_calculate.cpp
void calculate_energies_and_transport(double& enn, double& enp, double& enw, double& eeb, double& ezf, double& fne,
                                      double& etran_grad, double& etran_adia, double& etran_visc,
                                      AXY n_e, AXY n_i, AXY p_p, double lap[][nym],
                                      double avg[], double avp[], double hy, double taui,
                                      AXY pe, AXY pi, bool b_mhw, double chat, double diff,
                                      AXY hyve, AXY hyvi, int nx1, int ny1);

//linear_growth_and_frequency_calculate.cpp
void calculate_linear_growth_and_frequency(double ttt, double dt, int it, int itmax, int incon,
                                            double& freq, double& grow,
                                            AXY p_p, double& t_amp, double& t_old, int& jsgn, int& jsgn_old,
                                            double& ddtt, double& enw, double& enwo, int nxh, int nyh);

//calculate_and_write_ky_spectrum.cpp
void calculate_and_write_ky_spectrum(double AXY_data[][nym],
                                      const char *filename, fftw_plan hindfty, double *py_array,
                                      fftw_complex *ky_array, double *sumky_array, double *avgky_array,
                                      int nx1, int ny1, int ny, int nx,
                                      double hy, double TwoPi, double ly,
                                      int incon, int it);

//calculate_and_write_energy_ky_spectrum.cpp
void calculate_and_write_energy_ky_spectrum(const double p_p[][514],
                                            const char *filename, fftw_plan hindfty, double *py_array,
                                            fftw_complex *ky_array, double *sumky_array, double *pkyavge_array,
                                            int nx1, int ny1, int ny, int nx, double hy, double TwoPi, double ly,
                                            int incon, int it);

//calculate_and_write_kx_spectrum.cpp
void calculate_and_write_kx_spectrum(const AXY pe_data, const char *filename,
                                      int nx1, int ny1, int nx, int ny, double hy, double TwoPi, double ly,
                                      int incon, int it, double *pkxavg_array);

//diagnose.cpp
void diagnose( double ttt, double t00, int it, int itmax,
               AXY n_n, AXY n_e, AXY n_i, AXY p_p, AXY n_g );

//arakawa.cpp
void arakawa( AXY uuu, AXY vvv, AXY out );

//arakaw4.cpp
void arakaw4( AXY uuu, AXY vvv, AXY out );

//laplace.cpp
void laplace( AXY fi, AXY fo );

//poisson.cpp
void poisson( AXY fi, const AXY cp, AXY fo );

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
void write_global_outputs(double ttt, double enn, double eeb, double ezf,
                          double etran_grad, double etran_adia, double etran_visc, double etran_tot,
                          double fne);
void write_zonal_flow_data(double ttt, AXY p_p, int nx1, int ny1, double hy);
void write_x_cut_data(AXY n_e, AXY n_i, AXY p_p, double lap[][nym], AXY n_g, int nx1, int nyh, double hy);
void write_y_cut_data(AXY n_e, AXY n_i, AXY p_p, double lap[][nym], AXY n_g, int ny1, int nxh, double hy);
void write_x_profile_data(AXY p_p, AXY n_e, AXY n_i, double lap[][nym],
                          int nx1, int ny1, double hy, int ny);
void write_2d_plot_data(AXY n_e, AXY p_p, double lap[][nym], int nx1, int ny1);

int main(void); // main関数のプロトタイプ宣言

#endif // GHW1_HPP
