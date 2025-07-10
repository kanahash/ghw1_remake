#include "ghw1.hpp"

// ExB 対流項を計算する関数
static void exb_convection_calculate(AXY density_e, AXY potential_e, AXY result_e,
                              AXY density_i, AXY potential_i, AXY result_i, bool use_4th_order) 
{
    if (use_4th_order)
	{
        arakaw4(density_e, potential_e, result_e);
        arakaw4(density_i, potential_i, result_i);
    }
	else
	{
        arakawa(density_e, potential_e, result_e);
        arakawa(density_i, potential_i, result_i); // 2次の方が大きな dt に対して安定
    }
}

// ハイパー粘性項を計算する関数
static void hyper_viscosity_calculate(AXY density, AXY hyv_out, AXY vis_tmp,
                               fftw_plan hinp, fftw_plan herp, double wwx[][nym], double wwk[][nym/2+1][2],
                               const AXY chyv, int ihype, int nx1, int ny1, int nyp) 
{
    int i;
    int j;
    int k_idx;

    if (ihype == -1)
	{ // del^p を (p/2) ラプラス演算子で評価する
        laplace(density, hyv_out);
    }
	else if (ihype == 0)
	{ // del^p を (p/2) ラプラス演算子で評価する
        laplace(density, vis_tmp);
        laplace(vis_tmp, hyv_out);
    }
	else if (ihype > 0)
	{ // del^p の FFT を k^p として評価する: 推奨 (p==2)
        // データコピー
        #pragma omp parallel for private(j)
        for (i = 0; i <= nx1; ++i)
		{
            for (j = 0; j <= ny1; ++j)
                wwx[i][j] = density[i][j];
        }
        fftw_execute(hinp); // 順方向FFTを実行

        // 周波数空間での操作
        #pragma omp parallel for private(j, k_idx) // k_idx: FFTW配列の複素数成分インデックス (実部/虚部)
        for (i = 0; i <= nx1; ++i)
		{
            for (j = 0; j <= nyp - 1; ++j)
			{
                for (k_idx = 0; k_idx <= 1; k_idx++)
                    wwk[i][j][k_idx] *= chyv[i][j]; // 周波数空間で係数を乗算
            }
        }
        fftw_execute(herp); // 逆方向FFTを実行

        // 結果をコピー
        #pragma omp parallel for private(j)
        for (i = 0; i <= nx1; ++i)
		{
            for (j = 0; j <= ny1; ++j)
                hyv_out[i][j] = wwx[i][j];
        }
    }
}

// --- 主要な時間ステップ更新関数 ---
void update_densities_one_time_step(double ttt, double dt, AXY ne, AXY ni, AXY pe, AXY pi, AXY ane, AXY ani,
                                     AXY hyve, AXY hyvi, AXY vis_e, AXY vis_i,
                                     fftw_plan hinp, fftw_plan herp, double wwx[][nym], double wwk[][nym/2+1][2],
                                     const AXY chyv_coeff, bool use_4th_order, int ihype_mode, int nx_max, int ny_max, int ny_padded)
{
    timer_start(); // td_update の測定を開始 (ここで時間を計測)

    // ExB 対流項の計算
    exb_convection_calculate(ne, pe, ane, ni, pi, ani, use_4th_order);

    // 電子とイオンのハイパー粘性項を個別に計算
    hyper_viscosity_calculate(ne, hyve, vis_e, hinp, herp, wwx, wwk, chyv_coeff, ihype_mode, nx_max, ny_max, ny_padded);
    hyper_viscosity_calculate(ni, hyvi, vis_i, hinp, herp, wwx, wwk, chyv_coeff, ihype_mode, nx_max, ny_max, ny_padded);

    // ここに時間更新スキーム (例えばルンゲクッタ法など) の残りの部分を追加
    // 通常、各項 (ane, hyveなど) を組み合わせて時間微分項 (fne0, fni0など) を構築し、
    // それらを使って密度の配列 (ne, ni) を更新します。
    // 以下はRK4のようなスキームの一部になる可能性のある概念的なコードです。

    // 仮の時間微分項 fne0 と fni0 (コメントアウトを解除して使用する場合)
    AXY fne0, fni0; // これらは引数で渡されるか、グローバルである必要があります

    // 密度の更新ロジックの例 (単純なオイラー法の場合)
    // for (int i = 0; i <= nx_max; ++i) {
    //     for (int j = 0; j <= ny_max; ++j) {
    //         // 密度連続の式に基づいて fne0 と fni0 を計算
    //         fne0[i][j] = ane[i][j] - hyve[i][j]; // + 他の物理項
    //         fni0[i][j] = ani[i][j] - hyvi[i][j]; // + 他の物理項

    //         ne[i][j] += dt * fne0[i][j];
    //         ni[i][j] += dt * fni0[i][j];
    //     }
    // }

    // ルンゲクッタ法などのより複雑なスキームでは、
    // k1, k2, k3, k4 などの計算がここに入ります。
    // 例えば、ne0, ne1, ne2 などが時間スキームの中間ステップに使われることがあります。

    // timer_stop() は呼び出し元で処理することも考えられる
}
