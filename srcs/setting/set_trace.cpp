#include "ghw1.hpp"

// トレース出力位置を設定する関数
void set_trace_output_locations(int nx1, int ny1, int kt_maxsqrt, 
                               int iout[], int jout[], int& k_max) 
{
    int kk = 0;
    // k_max を計算 (配列の総サイズ)
    k_max = kt_maxsqrt * kt_maxsqrt; 
    
    // トレース点のインデックスを計算して格納
    for (int ki = 0; ki < kt_maxsqrt; ++ki)
	{
        for (int kj = 0; kj < kt_maxsqrt; ++kj)
		{
            iout[kk] = nx1 * ki / kt_maxsqrt;
            jout[kk] = ny1 * kj / kt_maxsqrt;
            kk++;
        }
    }
}
