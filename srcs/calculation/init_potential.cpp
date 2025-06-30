#include "ghw1.hpp"

// 初期ポテンシャルを計算する関数
void init_potentials(AXY ni, AXY ne, AXY& pe, AXY& gni, AXY& ww, 
                                const AXY cpoti, const AXY cvort,
                                double aae, double aai)
{
    // ni からジャイロ遮蔽イオン密度 gni を計算する
    poisson(ni, cpoti, gni);

    // ww (渦度に対応) を計算する
    for (int i = 0; i <= nx1; ++i)
	{
        for (int j = 0; j <= ny1; j++)
            ww[i][j] = -aae * ne[i][j] - aai * gni[i][j];
    }

    // ww から初期静電ポテンシャル pe を計算する
    poisson(ww, cvort, pe);
}
