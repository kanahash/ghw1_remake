#include "ghw1.hpp"

void update_history_densities(AXY ne, AXY ne0, AXY ne1, AXY ne2,
                              AXY ni, AXY ni0, AXY ni1, AXY ni2,
                              AXY fne0, AXY fne1, AXY fne2,
                              AXY fni0, AXY fni1, AXY fni2,
                              int nx_max, int ny_max)
{
    int i, j; // ループ変数

#pragma omp parallel for private(j) collapse(2)
    for (i = 0; i <= nx_max; ++i)
	{
        for (j = 0; j <= ny_max; ++j)
		{
            // 電子の密度履歴を更新
            ne2[i][j] = ne1[i][j];
            ne1[i][j] = ne0[i][j];
            ne0[i][j] = ne[i][j];

            // イオンの密度履歴を更新
            ni2[i][j] = ni1[i][j];
            ni1[i][j] = ni0[i][j];
            ni0[i][j] = ni[i][j];

            // 電子の連続方程式の右辺の履歴を更新
            fne2[i][j] = fne1[i][j];
            fne1[i][j] = fne0[i][j];

            // イオンの連続方程式の右辺の履歴を更新
            fni2[i][j] = fni1[i][j];
            fni1[i][j] = fni0[i][j];
        }
    }
}
