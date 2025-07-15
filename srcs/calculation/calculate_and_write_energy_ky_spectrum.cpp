#include "ghw1.hpp"

void calculate_and_write_energy_ky_spectrum(const double p_p[][514], // const double (*p_p)[514] と互換性のある形式
        const char *filename, fftw_plan hindfty, double *py_array, // fftw_plan のまま
        fftw_complex *ky_array, double *sumky_array, double *pkyavge_array,
        int nx1, int ny1, int ny, int nx, double hy, double TwoPi, double ly,
        int incon, int it)
{
    int     im;
    int     ip;
    int     jm;
    int     jp;
    FILE    *g;

    // sumkyの初期化
    for (int j = 0; j <= ny / 2 + 1; j++)
        sumky_array[j] = 0.;
    // X方向のループとFFT、パワースペクトル合計
    for (int i = 0; i <= nx1; i++)
    {
        im = (i == 0) ? nx1 : i - 1;
        ip = (i == nx1) ? 0 : i + 1;
        for (int j = 0; j <= ny1; j++)
        {
            jm = (j == 0) ? ny1 : j - 1;
            jp = (j == ny1) ? 0 : j + 1;
            double dum_y = 0.5 * hy * (p_p[i][jp] - p_p[i][jm]); // y方向勾配
            py_array[j] = dum_y * dum_y;                                 // y方向勾配の二乗
            double dum_x = 0.5 * hy * (p_p[ip][j] - p_p[im][j]); // x方向勾配
            py_array[j] += dum_x * dum_x;                                // x方向勾配の二乗を加算
        }
        fftw_execute(hindfty);
        for (int j = 0; j <= ny / 2 + 1; j++)
            sumky_array[j] += ky_array[j][0] * ky_array[j][0]; // パワースペクトルを加算
    }
    // X方向平均
    for (int j = 0; j <= ny / 2 + 1; j++)
        sumky_array[j] /= (double)nx;
    // 時間平均
    if (incon == 2)
    {
        for (int j = 1; j <= ny / 2; j++)
        {
            pkyavge_array[j] = ((it - 1 + 1) * pkyavge_array[j]
                            + sumky_array[j]) / (it + 1);
        }
    }
    // ファイル出力
    g = fopen(filename, "w");
    if (g == NULL)
    {
        perror("Error opening file for Energy Ky spectrum");
        return ; // エラー処理
    }
    for (int j = 1; j < ny / 2; j++)
    {
        fprintf(g, "%.6e  %.6e  %.6e\n", TwoPi * (double)j / ly, 0.5
            * sumky_array[j] * ny / hy / hy, 0.5 * pkyavge_array[j] * ny / hy
            / hy);
    }
    fclose(g);
}
