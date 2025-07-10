// srcs/main.cpp

#include "ghw1.hpp"
#include <cstdio>   // For printf, fprintf, fopen, fclose, fscanf, feof, stderr
#include <cstdlib>  // For exit, EXIT_FAILURE, atof, atoi
#include <cmath>    // For exp, std::exp
#include <cstddef>  // For NULL

// =========================================================
// ghw1.hpp で extern 宣言されているグローバル変数の定義
// これらの変数が main.cpp の中で一度だけ定義されていることを確認してください。
// =========================================================

unsigned int fflag;

int nx, ny, nxh, nyh, nx1, ny1, ict, jct;
int itstp, itmax, npar, ntrace, jnull_old;
int ipade, ihype;
int maxtrace; 

double p_amp, p_old, t_amp, t_old;
double t_1, t_2;
int jsgn, jsgn_old;

// ★★★ itn と itend をここに追加します！★★★
int itn;
int itend;


double pkxavg[nxm / 2 + 1];
double pkyavg[nym / 2 + 1];
double pkyavgn[nxm / 2 + 1], pkyavgw[nym / 2 + 1], pkyavge[nym / 2 + 1];

double hy, hy2, hysq, xyz, dr, amp, wsrc, sigma, xlwl;
double diff, n00, ly, dt, dtt, ddtt, incon;
double chat, diss, neavg, niavg, peavg, vorfree;
double aspect, eno, enwo, delta, delinv, g_n;
double pdist, vorset, clin, vde_old, vde_avg;
double aae, aai, mue, mui, taue, taui, zze, zzi;

// AXY型の配列定義
AXY pe, pi, ww;
AXY nn, ne, ni, ne0, ne1, ne2, ni0, ni1, ni2;
AXY gni;
AXY fne0, fne1, fne2, fni0, fni1, fni2;
AXY vis_e, vis_i, hyve, hyvi;
AXY cvort, cpoti, clap, chyv, cvfis; // chyv は AXY 型として定義されていることを確認

char *wisdom_sf;

bool printtraces, b_mhw, b_4th;

// trace_n_global と trace_p_global の定義
double trace_n_global[MAX_TRACE_K_MAX_SQRT * MAX_TRACE_K_MAX_SQRT][MAX_TRACE_STEPS];
double trace_p_global[MAX_TRACE_K_MAX_SQRT * MAX_TRACE_K_MAX_SQRT][MAX_TRACE_STEPS];

// wwx_global と wwk_global の定義
AXY wwx_global;
double wwk_global[nxm][nym / 2 + 1][2];

// nyp の定義
int nyp;

// FFTW plan 変数の定義
fftw_plan hinp_main;
fftw_plan herp_main;

// 時間計測変数の定義
double td_init;
double td_update;
double td_pois1;
double td_pois2;
double td_bnd;
double td_mem;
double td_end;
double td_out;
double td_pol;
double td_tot;

// ttt と t00 もグローバル変数として定義
double ttt; 
double t00;


// =========================================================

int main(void)
{
    // main 関数内のローカル変数宣言 (グローバル変数と重複しないように注意)
    double c0 = 18.0 / 11.0;
    double c1 = 9.0 / 11.0;
    double c2 = 2.0 / 11.0;
    double cf = 6.0 / 11.0;

    // td_変数は上記でグローバル変数として定義済みなので、ここでは宣言しません。
    // 初期化は main 関数内で引き続き行います。
    td_init = 0.0;
    td_update = 0.0;
    td_pois1 = 0.0;
    td_pois2 = 0.0;
    td_bnd = 0.0;
    td_mem = 0.0;
    td_end = 0.0;
    td_out = 0.0;
    td_pol = 0.0;
    td_tot = 0.0;

    double bndys[nxm]; // nxm は ghw1.hpp で const static int として定義済み
    double ne_za[nxm], ni_za[nxm];
    char s[80], str[80];
    int nx_i;
    double dx0;
    int rs_check;

    // ループ変数や一時変数 (it, itn, itend, is をここに追加)
    int i, j, k, l, m, n, ik, jk, it, is, itn, itend; 

    // init_time_stepping 関数に渡すためのローカル変数 (int& で参照渡しされるため)
    int itrace, jtrace;

    // トレース出力位置関連の配列と変数
    int iout[kt_maxsqrt * kt_maxsqrt], jout[kt_maxsqrt * kt_maxsqrt]; // kt_maxsqrt は ghw1.hpp で const static int として定義済み
    int k_max; // set_trace_output_locations で値が設定される

    // ttt, t00 はグローバル変数として定義済みなので、ここでは宣言しません。
    double phase, zuf, nue, nui, nuzf; // main()でしか使われないローカル変数
    double neprof, pbar, nbar, dve, gre, gri;

    AXY ane, ani, dummy; // AXY は ghw1.hpp で typedef されています

    // ファイルポインタ。FILE型は <cstdio> で定義されます。
    FILE *f, *g, *h, *g1, *g2, *g3, *g4, *g5, *ft;

    t_1 = timer_start(); // timer_start は ghw1.hpp で inline 関数として定義済み

    // 入力ファイルを読み込み、パラメータを初期化する:
    init_parameters(); // init_parameters は ghw1.hpp でプロトタイプ宣言済み

    // 並列化されたFFTWを初期化する:
    init_FFTW(nx, ny, fflag, npar, hinp_main, herp_main, wisdom_sf); // init_FFTW は ghw1.hpp でプロトタイプ宣言済み。wisdom_sf はグローバル変数

    // オプション: ディリクレ境界条件を定義する (選択した場合):
    for (i = 0; i <= nx1; i++)
    {
        nx_i = static_cast<int>(static_cast<double>(nx1 - i));
        dx0 = 2.0 * 2.0;
        bndys[i] = 1.0 - std::exp(-nx_i * nx_i / dx0) - std::exp(-static_cast<double>((i - 0) * (i - 0)) / dx0); // std::exp を使用
    }

    // 初期密度摂動を設定する (ブロブ、渦、乱流、流れなど)
    set_init_density_perturbation(incon, ne, ni); // set_init_density_perturbation は ghw1.hpp でプロトタイプ宣言済み

    // 初期ポテンシャルを計算する (関数名修正)
    init_potentials(ni, ne, pe, gni, ww, cpoti, cvort, aae, aai); // init_potentials は ghw1.hpp でプロtoタイプ宣言済み

    // 人為的な「以前の」時間値を設定する (多段ソルバー用) (関数名修正)
    init_history_densities(ne, ni, ne0, ne1, ne2, ni0, ni1, ni2, fne1, fne2, fni1, fni2); // init_history_densities は ghw1.hpp でプロトタイプ宣言済み

    // エネルギーを初期化する:
    init_energy(ne, xyz, eno); // init_energy は ghw1.hpp でプロトタイプ宣言済み

    if (incon == 2)
    { // 以前の実行の保存された終了出力からリスタート
        printf("| データセットを継続します... \n");
        rs_check = 1;
        g = fopen("restart.dat", "r");
        if (g == NULL)
        {
            printf("\n ファイル 'restart.dat' が見つかりません。終了します... \n");
            exit(EXIT_FAILURE); // EXIT_FAILURE は <cstdlib> で定義
        }
        int read_count;

        read_count = fscanf(g, "%s ", str);
        if (read_count != 1 && !feof(g)) { fprintf(stderr, "Error reading t00 from restart.dat\n"); exit(EXIT_FAILURE); }
        t00 = atof(str); // atof は <cstdlib> で定義

        while (!feof(g)) // feof は <cstdio> で定義
        {
            read_count = fscanf(g, "%s ", str); // fscanf は <cstdio> で定義
            if (read_count != 1) { break; } // EOFまたはエラーの場合はループを抜ける
            i = atoi(str); // atoi は <cstdlib> で定義

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading j for i=%d from restart.dat\n", i); exit(EXIT_FAILURE); }
            j = atoi(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading eno for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            eno = atof(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading pe for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            pe[i][j] = atof(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading ne for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            ne[i][j] = atof(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading ne1 for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            ne1[i][j] = atof(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading ne2 for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            ne2[i][j] = atof(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading ni for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            ni[i][j] = atof(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading ni1 for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            ni1[i][j] = atof(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading ni2 for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            ni2[i][j] = atof(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading fne1 for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            fne1[i][j] = atof(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading fne2 for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            fne2[i][j] = atof(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading fni1 for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            fni1[i][j] = atof(str);

            read_count = fscanf(g, "%s ", str);
            if (read_count != 1) { fprintf(stderr, "Error reading fni2 for (%d,%d) from restart.dat\n", i, j); exit(EXIT_FAILURE); }
            fni2[i][j] = atof(str);
        }
        fclose(g); // fclose は <cstdio> で定義
        poisson(pe, cpoti, pi); // poisson は ghw1.hpp でプロトタイプ宣言済み
        for (i = 0; i <= nx1; i++)
        {
            for (j = 0; j <= ny1; j++)
            {
                ne0[i][j] = ne[i][j];
                ni0[i][j] = ni[i][j];
            }
        }
    }
    // トレース出力位置を設定する (関数名修正)
    set_trace_output_locations(nx1, ny1, kt_maxsqrt, iout, jout, k_max); // set_trace_output_locations は ghw1.hpp でプロトタイプ宣言済み

    // 時間ステップを初期化する:
    // itn, itend, it, is は main 関数内のローカル変数として宣言済み
    if (init_time_stepping(itstp, itmax, dt, cf, itn, itend, dtt, ddtt, itrace, jtrace) == 0) // init_time_stepping は ghw1.hpp でプロトタイプ宣言済み
        return (0); // シミュレーションを終了

    // 初期プロファイルの出力制御
    diagnose(0.0, 0.0, it, itmax, nn, ne, ni, pe, gni); // diagnose は ghw1.hpp でプロトタイプ宣言済み
    td_init = timer_stop(t_1); // timer_stop は ghw1.hpp で inline 関数として定義済み

    // メインの時間ループを run_time_loop 関数で実行
    run_time_loop(); // run_time_loop は ghw1.hpp でプロトタイプ宣言済み

    // 時間ステップ終了 --------------------------------------------------------
    t_1 = timer_start();
    // 完全なリスタートファイルを書き込む
    g = fopen("restart.dat", "w");
    if (g == NULL)
    {
        fprintf(stderr, "Error: Could not open restart.dat for writing.\n");
        exit(EXIT_FAILURE);
    }
    fprintf(g, "%5e \n", ttt); // ttt は run_time_loop内で最終値が設定される
    for (i = 0; i <= nx1; i++)
        for (j = 0; j <= ny1; j++)
            fprintf(g,
                    "%d %d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e "
                    "%.6e %.6e %.6e %.6e\n",
                    i, j, eno, pe[i][j], ne[i][j],
                    ne1[i][j],
                    ne2[i][j],
                    ni[i][j],
                    ni1[i][j],
                    ni2[i][j],
                    fne1[i][j],
                    fne2[i][j],
                    fni1[i][j],
                    fni2[i][j]);
    fclose(g);

    // 時間トレース / 周波数スペクトル解析:
    if (printtraces)
        perform_time_trace_analysis(maxtrace, k_max, t00, dt, itmax, ntrace); // perform_time_trace_analysis は ghw1.hpp でプロトタイプ宣言済み

    td_out += timer_stop(t_1);

    td_tot = td_init + td_update + td_pois1 + td_pois2 + td_pol + td_bnd + td_out + td_mem;
    printf("|\n| 総実行時間:  %.2f s  =  %.2f min \n", td_tot, td_tot / 60.0);
    td_tot = 100.0 / td_tot;

    // 総時間ボトルネック分析出力:
    printf("|\n| 絶対実行時間: init: %.2f upd: %.2f poi1: %.2f poi2: %.2f pol: %.2f bnd: %.2f out: %.2f mem: %.2f\n",
           td_init, td_update, td_pois1, td_pois2, td_pol, td_bnd, td_out, td_mem);

    // 割合ボトルネック分析出力:
    printf("| 相対実行時間: init: %.2f upd: %.2f poi1: %.2f poi2: %.2f pol: %.2f bnd: %.2f out: %.2f mem: %.2f\n|\n",
           td_init * td_tot, td_update * td_tot,
           td_pois1 * td_tot,
           td_pois2 * td_tot,
           td_pol * td_tot,
           td_bnd * td_tot,
           td_out * td_tot,
           td_mem * td_tot);

#ifdef _OPENMP
    fftw_cleanup_threads(); // fftw_cleanup_threads は <fftw3.h> で定義
#endif
    printf("| GHW1 終了。\n\n");

    return 0;
}
