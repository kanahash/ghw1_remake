#include "ghw1.hpp"

void	calculate_and_write_ky_spectrum(double AXY_data[nxm][nym],
		const char *filename, fftw_plan hindfty, double *py_array,
		fftw_complex *ky_array, double *sumky_array, double *avgky_array,
		int nx1, int ny1, int ny, int nx, int incon, int it)
{
	FILE	*g;

	// sumkyの初期化
	for (int j = 0; j <= ny / 2 + 1; j++)
	{
		sumky_array[j] = 0.;
	}
	// X方向のループとFFT、スペクトル合計
	for (int i = 0; i <= nx1; i++)
	{
		for (int j = 0; j <= ny1; j++)
		{
			py_array[j] = AXY_data[i][j];
		}
		fftw_execute(hindfty);
		// ここでスペクトルの種類を選択できるよう引数で制御することも可能
		// 現在のコードはfabs(ky[j][0])を使用
		for (int j = 0; j <= ny / 2 + 1; j++)
		{
			sumky_array[j] += fabs(ky_array[j][0]); // またはky_array[j][0]
			*ky_array[j][0]
		}
	}
	// X方向平均
	for (int j = 0; j <= ny / 2 + 1; j++)
	{
		sumky_array[j] /= (double)nx;
	}
	// 時間平均
	if (incon == 2)
	{
		for (int j = 1; j <= ny / 2; j++)
		{
			avgky_array[j] = ((it - 1 + 1) * avgky_array[j] + sumky_array[j])
				/ (it + 1);
		}
	}
	// ファイル出力
	g = fopen(filename, "w");
	if (g == NULL)
	{
		perror("Error opening file for Ky spectrum");
		return ; // エラー処理
	}
	for (int j = 1; j < ny / 2; j++)
	{
		fprintf(g, "%.6e  %.6e  %.6e\n", TwoPi * (double)j / ly, sumky_array[j]
			* ny / hy / hy, avgky_array[j] * ny / hy / hy);
	}
	fclose(g);
}
