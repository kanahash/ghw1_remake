#include "ghw1.hpp"

// メインの時間ループを簡潔にする関数
void run_time_loop()
{
    for (it = itn / itstp; it < itend; ++it)
	{ // 全体の時間ループ (出力ステップ)
        for (is = 0; is < itstp; ++is)
		{    // 内部時間ループ (出力なし)
            ttt = t00 + (it + 0) * ddtt + (is + 1) * dt;

            // ジャイロ流体密度 ne と ni の時間ステップ更新:
            // すべての複雑な計算を単一の関数呼び出しで処理
            update_densities_one_time_step(ttt, dt, ne, ni, pe, pi, ane, ani,
                                             hyve, hyvi, vis_e, vis_i,
                                             hinp, herp, wwx, wwk,
                                             chyv, b_4th, ihype, nx1, ny1, nyp);

            // 必要に応じて他の診断や更新処理を追加
        }
        // 出力ステップでの処理 (例えば diagnose 関数など)
        // diagnose(ttt, t00, it, itmax, nn, ne, ni, pe, gni);
    }
}
