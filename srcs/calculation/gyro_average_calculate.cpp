#include "ghw1.hpp"

void gyro_averaged_ion_density_calculate(AXY ni_in, AXY gni_out, fftw_plan hinp, fftw_plan herp,
                                         double wwx[][nym], double wwk[][nym/2+1][2],
                                         const AXY cpoti_coeff, double taui_param,
                                         int nx_max, int ny_max, int ny_padded) 
{
    int i, j, k; // ループ変数

    if (taui_param > 0.)
	{
        // データコピー
		#pragma omp parallel for private(j)
        	for (i = 0; i <= nx_max; ++i)
			{
            	for (j = 0; j <= ny_max; ++j)
                	wwx[i][j] = ni_in[i][j];
        	}
        fftw_execute_dft_r2c(hinp, wwx[0], (fftw_complex*)wwk[0]); // スレッド化された FFTW 呼び出し

        // 周波数空間での操作
		#pragma omp parallel for private(j,k)
        	for (i = 0; i <= nx_max; ++i)
			{
            	for (j = 0; j <= ny_padded - 1; ++j)
				{
                	for (k = 0; k <= 1; k++)
                    	wwk[i][j][k] *= cpoti_coeff[i][j];
            	}
        	}
        	fftw_execute_dft_c2r(herp, (fftw_complex*)wwk[0], wwx[0]); // スレッド化された FFTW 呼び出し

        // 結果を正規化してコピー
        	double norm_factor = 1.0 / ( (double)(nx_max + 1) * (ny_max + 1) );
		#pragma omp parallel for private(j)
        	for (i = 0; i <= nx_max; ++i)
			{
            	for (j = 0; j <= ny_max; ++j)
                	gni_out[i][j] = wwx[i][j] * norm_factor;
        	}
    }
	else // taui == 0 の場合のショートカット: ni を gni に直接コピー
        f_copy2darray(ni_in, gni_out, nx_max, ny_max);
}
