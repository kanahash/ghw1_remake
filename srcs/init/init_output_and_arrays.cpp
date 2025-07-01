#include "ghw1.hpp"

// 出力ファイルと配列を初期化する
void init_output_and_arrays(int nx, int ny, int incon, bool printtraces)
{
    int i, j; // ループ変数

    // グローバル変数の初期化（元のコードにあったもの）
    nxh = nx / 2;
    nyh = ny / 2;
    nx1 = nx - 1;
    ny1 = ny - 1;

    hy2 = hy * hy;
    hysq = hy2 / 12.;
    // g_nとchatの関連付けはinit_parametersに残すか、別途処理が必要になる可能性あり
    // g_n = 1.;
    // if (chat == 0.) g_n = 0.;

    xyz = 1. / double(nx * ny);

    // 出力ファイルの初期化
    if (incon != 2)
    {
        FILE *f;
        f = fopen("eng.dat", "w");
        fclose(f);
        f = fopen("eng-transfer.dat", "w");
        fclose(f);
        f = fopen("fne.dat", "w");
        fclose(f);
        f = fopen("gamma.dat", "w");
        fclose(f);
        f = fopen("traces.dat", "w");
        fclose(f);
        f = fopen("time.dat", "w");
        fclose(f);
        f = fopen("zfx.dat", "w");
        fclose(f);
    }

    // 平均配列の初期化
    for (i = 1; i <= nx / 2; i++)
        pkxavg[i] = 0;
    for (j = 1; j <= ny / 2; j++)
    {
        pkyavg[j] = 0;
        pkyavgn[j] = 0;
        pkyavgw[j] = 0;
        pkyavge[j] = 0;
    }
}