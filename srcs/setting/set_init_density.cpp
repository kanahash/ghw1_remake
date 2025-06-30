#include "ghw1.hpp"

// 初期密度摂動を設定する関数
void set_init_density_perturbation(int incon, AXY ne, AXY ni)
{
    // まず、全ての密度配列を0で初期化します
    for (int i = 0; i <= nx1; ++i)
	{
        for (int j = 0; j <= ny1; ++j)
		{
            ne[i][j] = 0.;
            ni[i][j] = 0.;
        }
    }

    // incon の値に基づいて、適切な初期化関数を呼び出します
    if (incon == 0.)
        init_add_blob(ne, ni); // ガウス型密度ブロブの初期化

    else if (incon == 1.)
        init_add_turb(ne, ni); // 擬似乱流浴の初期化

    // incon==2 はリスタートを意味するため、ここでは処理しません（main関数で別途処理されるため）
    else if (incon == 3.)
        init_add_dual(ne, ni); // 二重渦合体初期化

    else if (incon == 4.)
        init_add_flow(ne, ni); // シアー流初期化

    else if (incon == 5.)
        init_add_mode(ne, ni); // 「線形」単一 ky モード初期化
		
    // その他の incon 値に対するエラーハンドリングやデフォルト動作を追加することも可能です
}