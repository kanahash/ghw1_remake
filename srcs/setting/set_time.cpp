#include "ghw1.hpp"

// 履歴密度を初期設定する関数 (多段ソルバー用)
void init_history_density(AXY ne, AXY ni, 
                                AXY& ne0, AXY& ne1, AXY& ne2,
                                AXY& ni0, AXY& ni1, AXY& ni2,
                                AXY& fne1, AXY& fne2, AXY& fni1, AXY& fni2)
{
    int i;
    int j;

    #pragma omp parallel for private(i,j) // OpenMP による並列化
    for (i = 0; i <= nx1; ++i)
	{
        for (j = 0; j <= ny1; ++j)
		{
            // 現在の密度を履歴にコピー
            ne0[i][j] = ne[i][j];
            ne1[i][j] = ne[i][j];
            ne2[i][j] = ne[i][j];

            ni0[i][j] = ni[i][j];
            ni1[i][j] = ni[i][j];
            ni2[i][j] = ni[i][j];

            // フラックス項をゼロクリア
            fne1[i][j] = 0.;
            fne2[i][j] = 0.;
            fni1[i][j] = 0.;
            fni2[i][j] = 0.;
        }
    }
}
