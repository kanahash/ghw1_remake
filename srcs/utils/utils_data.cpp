#include "ghw1.hpp"

void collect_time_trace_data(double ttt_in, AXY ne_in, AXY pe_in,
                             int iout_loc[], int jout_loc[], int k_max_loc,
                             int& current_itrace, int& current_jtrace, int ntrace_interval,
                             double& td_out_total, double& t_start_measure,
                             bool printtraces_flag)
{
    if (printtraces_flag)
	{ // printtraces フラグに基づいて実行
        t_start_measure = timer_start(); // td_out の測定を開始

        current_itrace++;
        if (current_itrace >= ntrace_interval)
		{ // ntrace が ntrace_interval に対応
            // ファイル出力はコメントアウトされたまま、データを配列に記録するロジック
            current_itrace = 0; // カウンターをリセット

            // 記録されるトレースステップ数が MAX_TRACE_STEPS を超えないようにチェック
            if (current_jtrace < MAX_TRACE_STEPS)
			{
                for (int k = 0; k < k_max_loc; k++)
				{
                    trace_n_global[k][current_jtrace] = ne_in[iout_loc[k]][jout_loc[k]];
                    trace_p_global[k][current_jtrace] = pe_in[iout_loc[k]][jout_loc[k]];
                }
                current_jtrace++; // 記録ステップ数をインクリメント
            }
			else // 配列の容量を超えた場合の警告またはエラー処理
                printf("Warning: Maximum trace steps reached. Data might be overwritten or truncated.\n");
        }
        td_out_total += timer_stop(t_start_measure); // td_out の測定を停止
    }
}

// グローバルなエネルギー、輸送量、時系列データをファイルに書き出す関数
void write_global_outputs(double ttt, double enn, double eeb, double ezf,
                          double etran_grad, double etran_adia, double etran_visc, double etran_tot,
                          double fne)
{
    FILE *f;

    // グローバルエネルギー時系列出力:
    f = fopen("eng.dat", "a");
    if (ttt > 0.) {
        fprintf(f, "%.8f  %.6e  %.6e  %.6e  %.6e  %.6e\n",
                ttt, enn, eeb, ezf, (enn + eeb), ezf / eeb);
    }
    fclose(f);

    // グローバルエネルギー転送時系列出力:
    f = fopen("eng-transfer.dat", "a");
    if (ttt > 0.) {
        fprintf(f, "%.8f  %.6e  %.6e  %.6e  %.6e\n",
                ttt, etran_grad, etran_adia, etran_visc, etran_tot);
    }
    fclose(f);

    // グローバル輸送時系列出力:
    f = fopen("fne.dat", "a");
    if (ttt > 0.) {
        fprintf(f, "%.8f  %.6e\n", ttt, fne);
    }
    fclose(f);
}

// 帯状流 (t,x) 2D プロットデータをファイルに書き出す関数
void write_zonal_flow_data(double ttt, AXY p_p, int nx1, int ny1, double hy)
{
    FILE *h = fopen("zfx.dat", "a");
    double zfx;
    int im, ip;

    for (int i = 0; i <= nx1; i++)
    {
        im = (i == 0) ? nx1 : i - 1;
        ip = (i == nx1) ? 0 : i + 1;
        for (zfx = 0., int j = 0; j <= ny1; j++)
        {
            // zfx += p_p[i][j]/(ny); // 帯状ポテンシャル
            zfx += p_p[ip][j] - p_p[im][j]; // 帯状流 <Vy>
        }
        zfx *= .5 * hy / ny; // pot または vor. の場合はコメントアウト
        if (ttt > 0.) {
            fprintf(h, "%.3f  %d  %.6e \n", ttt, i, zfx);
        }
    }
    fprintf(h, "\n");
    fclose(h);
}

// y=ny/2 での x-カット データをファイルに書き出す関数
void write_x_cut_data(AXY n_e, AXY n_i, AXY p_p, double lap[][nym], AXY n_g, int nx1, int nyh, double hy)
{
    FILE *g = fopen("cutx.dat", "w");
    for (int i = 0; i <= nx1; i++) {
        for (int j = nyh; j <= nyh; j++) {
            fprintf(g, "%f  %.6e   %.6e   %.6e   %.6e   %.6e\n",
                    double(i) / hy, n_e[i][j], n_i[i][j], p_p[i][j], lap[i][j], n_g[i][j]);
        }
    }
    fclose(g);
}

// x=nx/2 での y-カット データをファイルに書き出す関数
void write_y_cut_data(AXY n_e, AXY n_i, AXY p_p, double lap[][nym], AXY n_g, int ny1, int nxh, double hy)
{
    FILE *g = fopen("cuty.dat", "w");
    for (int j = 0; j <= ny1; j++) {
        fprintf(g, "%f  %.6e   %.6e   %.6e   %.6e   %.6e\n",
                double(j) / hy, n_e[nxh][j], n_i[nxh][j], p_p[nxh][j], lap[nxh][j], n_g[nxh][j]);
    }
    fclose(g);
}

// y で平均された x プロファイルデータをファイルに書き出す関数
void write_x_profile_data(AXY p_p, AXY n_e, AXY n_i, double lap[][nym],
                          int nx1, int ny1, double hy, int ny)
{
    FILE *g = fopen("xprof.dat", "w");
    double peprof, neprof, niprof, nsprof, wwprof;
    double rey;
    int im, ip, jm, jp;

    for (int i = 0; i <= nx1; i++) {
        peprof = 0.; neprof = 0.; niprof = 0.; nsprof = 0.; wwprof = 0.; rey = 0.;
        im = (i == 0) ? nx1 : i - 1;
        ip = (i == nx1) ? 0 : i + 1;
        for (int j = 0; j <= ny1; j++) {
            jm = (j == 0) ? ny1 : j - 1;
            jp = (j == ny1) ? 0 : j + 1;
            // レイノルズ応力:
            rey += .25 * hy * hy * (p_p[ip][j] - p_p[im][j]) * (p_p[i][jp] - p_p[i][jm]);
            // ポテンシャル、密度、渦度
            peprof += p_p[i][j];
            neprof += n_e[i][j];
            niprof += n_i[i][j];
            wwprof += lap[i][j];
        }
        peprof /= ny; neprof /= ny; niprof /= ny; nsprof /= ny; wwprof /= ny; rey /= ny;
        fprintf(g, "%f  %.6e   %.6e   %.6e   %.6e   %.6e\n",
                double(i) / hy, peprof, neprof, niprof, wwprof, rey);
    }
    fclose(g);
}

// 密度、渦度、ポテンシャルの 2D (x,y) プロットデータをファイルに書き出す関数
void write_2d_plot_data(AXY n_e, AXY p_p, double lap[][nym], int nx1, int ny1)
{
    FILE *g1 = fopen("n2d0.dat", "w");
    FILE *g2 = fopen("w2d0.dat", "w");
    FILE *g3 = fopen("p2d0.dat", "w");
    FILE *g4 = fopen("i2d0.dat", "w");

    for (int i = 0; i <= nx1; i++) {
        for (int j = 0; j <= ny1; j++) {
            fprintf(g1, "%d  %d  %.6e\n", i, j, n_e[i][j]);
            fprintf(g2, "%d  %d  %.6e\n", i, j, lap[i][j]);
            fprintf(g3, "%d  %d  %.6e\n", i, j, p_p[i][j]);
            fprintf(g4, "%d  %d  %.6e\n", i, j, (p_p[i][j] - n_e[i][j]));
        }
        fprintf(g1, "\n");
        fprintf(g2, "\n");
        fprintf(g3, "\n");
        fprintf(g4, "\n");
    }
    fclose(g1);
    fclose(g2);
    fclose(g3);
    fclose(g4);

    rename("n2d0.dat", "n2d.dat");
    rename("w2d0.dat", "w2d.dat");
    rename("p2d0.dat", "p2d.dat");
    rename("i2d0.dat", "i2d.dat");
}

