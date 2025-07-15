#include "ghw1.hpp"

void	calculate_and_write_kx_spectrum(const double pe_data[][514], const char *filename,
		int nx1, int ny1, int nx, int ny, double hy, double TwoPi, double ly,
		int incon, int it, double *pkxavg_array)
{
	double			px[nx];
	fftw_complex	kx[nx / 2 + 1], sumkx[nx / 2 + 1];
	fftw_plan		hindftx;
	FILE			*g;

	// FFTW のためのローカル変数
	// FFTW Plan の作成 (関数内で完結させるため、ここで作成・破棄)
	// ただし、この計画作成は一度行えばよいので、もし他のKxスペクトル計算があるなら、
	// diagnose 関数の先頭でまとめて作成し、引数で渡す方が効率的です。
	// 今回は、このブロックのみを切り出すという前提でここに含めます。
	hindftx = fftw_plan_dft_r2c_1d(nx, &px[0], &kx[0], FFTW_ESTIMATE);
	// sumkx の初期化
	for (int i = 0; i <= nx / 2; i++)
		sumkx[i][0] = 0.;
	// Y 方向のループと FFT、スペクトル合計
	for (int j = 0; j <= ny1; j++)
	{
		for (int i = 0; i <= nx1; i++)
			px[i] = pe_data[i][j]; // pe_data を使用
		fftw_execute(hindftx);
		for (int i = 0; i <= nx / 2 + 1; i++)
			sumkx[i][0] += fabs(kx[i][0]);
	}
	// Y 方向平均
	for (int i = 0; i <= nx / 2 + 1; i++)
		sumkx[i][0] /= (double)ny;
	// FFTW Plan の破棄
	fftw_destroy_plan(hindftx);
	// 時間平均
	if (incon == 2)
	{
		for (int i = 1; i <= nx / 2 + 1; i++)
		{
			pkxavg_array[i] = ((it - 1 + 1) * pkxavg_array[i] + sumkx[i][0])
				/ (it + 1);
		}
	}
	// ファイル出力
	g = fopen(filename, "w");
	if (g == NULL)
	{
		perror("Error opening file for Kx spectrum");
		return ; // エラー処理
	}
	for (int i = 1; i < nx / 2; i++)
	{
		fprintf(g, "%.6e  %.6e  %.6e\n", TwoPi * (double)i / ly, sumkx[i][0]
			* nx / hy / hy, pkxavg_array[i] * nx / hy / hy);
	}
	fclose(g);
}
