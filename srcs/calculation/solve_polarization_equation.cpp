#include "ghw1.hpp"

void solve_polarization_equation(AXY ne_in, AXY gni_in, AXY pe_out,
                                 fftw_plan hinp, fftw_plan herp,
                                 double wwx[][nym], double wwk[][nym/2+1][2],
                                 const AXY cvort_coeff, double zze_param, double zzi_param,
                                 int nx_max, int ny_max, int ny_padded)
{
    int i, j, k; // ループ変数

	#pragma omp parallel for private(j)
    	for (i = 0; i <= nx_max; ++i) 
		{
        	for (j = 0; j <= ny_max; ++j) // 入力データの構築
            	wwx[i][j] = - zze_param * ne_in[i][j] - zzi_param * gni_in[i][j];
    	}
    	fftw_execute_dft_r2c(hinp, wwx[0], (fftw_complex*)wwk[0]); // 順方向FFT

	#pragma omp parallel for private(j,k)
    	for (i = 0; i <= nx_max; ++i)
		{
        	for (j = 0; j <= ny_padded - 1; ++j)
			{
            	for (k = 0; k <= 1; k++)
                	wwk[i][j][k] *= cvort_coeff[i][j]; // 周波数空間で係数を乗算
        	}
    	}
    	fftw_execute_dft_c2r(herp, (fftw_complex*)wwk[0], wwx[0]); // 逆方向FFT

    	// 結果を正規化してコピー
    	double norm_factor = 1.0 / ( (double)(nx_max + 1) * (ny_max + 1) );
	#pragma omp parallel for private(j)
    	for (i = 0; i <= nx_max; ++i)
		{
        	for (j = 0; j <= ny_max; ++j)
            	pe_out[i][j] = wwx[i][j] * norm_factor;
    	}
}
