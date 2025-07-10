#include "ghw1.hpp" // ghw1.hpp をインクルードすることで、extern 宣言されたグローバル変数にアクセスできます

// === ここに、ghw1.hpp に extern 宣言が抜けている（または認識されない）グローバル変数を追加 ===
// もし ghw1.hpp に既に含まれていれば、このブロックは不要です。
// ただし、最後の手段として、ここに直接 extern 宣言を記述することでコンパイルを通すことがあります。

// 時間ステップ制御関連
extern int itn, itend, itstp, itmax;
extern double t00, ddtt, dt;

// FFTW関連のプランとデータ配列
extern fftw_plan hinp_main, herp_main; // main.cppで定義されているFFTWプラン
extern AXY wwx_global;
extern double wwk_global[nxm][nym / 2 + 1][2]; // nxm, nym は ghw1.hpp で const static int として定義済み
extern int nyp;

// 時間計測関連
extern double t_1, t_2;
extern double td_pois1, td_bnd, td_pol, td_pois2, td_mem, td_out;
extern double td_init, td_update, td_tot, ttt; // ttt もグローバル変数として利用

// 物理量配列 (AXY型)
extern AXY ne, ni, pe, pi, nn, gni;
extern AXY ne0, ne1, ne2, ni0, ni1, ni2; // 履歴密度
extern AXY fne0, fne1, fne2, fni0, fni1, fni2; // 連続の式項
extern AXY hyve, hyvi, vis_e, vis_i; // 粘性関連
extern AXY cvort, cpoti; // FFT カーネル関連
// extern AXY chyv; // <--- この行もghw1.hppにAXY型としてexternされているので不要

// その他、グローバルで使われる変数
extern bool printtraces, b_4th;
extern int ihype;
extern int nx1, ny1;
extern double taui, zze, zzi;
// extern double chyv; // <--- ★★★ この行を削除します！★★★

// ===========================================================================


void run_time_loop() {
    // run_time_loop() 関数内でローカルに宣言すべき変数
    int it, is; // ループカウンタ (it は global it をシャドウしないよう注意)
    // ttt はグローバル変数として定義済みなので、ここでは宣言しません。

    // AXY 型のローカル一時配列 (main.cpp から引き継がれる)
    AXY ane, ani;

    // トレース出力位置関連の配列と変数 (main.cpp から引き継がれる)
    // kt_maxsqrt は ghw1.hpp で const static int として定義されているので、直接使えます。
    int iout[kt_maxsqrt * kt_maxsqrt], jout[kt_maxsqrt * kt_maxsqrt];
    int k_max; // set_trace_output_locations で値が設定される

    // トレース関連のローカルカウンタ (main.cpp から引き継がれる)
    int itrace_local, jtrace_local;


    // itn, itstp, itend, t00, ddtt, dt などは ghw1.hpp で extern 宣言されており、
    // main.cpp で定義・初期化されているグローバル変数なので、ここで直接使えます。
    for (it = itn / itstp; it < itend; ++it) // itn, itend がグローバルとして認識されることを期待
    {
        for (is = 0; is < itstp; ++is) // itstp がグローバルとして認識されることを期待
        {
            // グローバル変数 ttt に値を設定
            ttt = t00 + static_cast<double>(it) * ddtt + static_cast<double>(is + 1) * dt; // t00, ddtt, dt がグローバルとして認識されることを期待

            // update_densities_one_time_step の呼び出し
            // 全ての引数はグローバル変数としてアクセスされます
            update_densities_one_time_step(ttt, dt, ne, ni, pe, pi, ane, ani,
                                            hyve, hyvi, vis_e, vis_i, hinp_main, herp_main,
                                            wwx_global, wwk_global, chyv, b_4th, ihype, nx1, ny1, nyp); // 全ての引数がグローバルとして認識されることを期待

            // 時間計測変数はグローバルなので直接加算
            // td_update += timer_stop(t_1); // この行は update_densities_one_time_step 内で td_update の加算が行われている場合、または td_update が参照渡しされるべき
            
            // ni から gni へのジャイロ平均:
            t_1 = timer_start();
            gyro_averaged_ion_density_calculate(ni, gni, hinp_main, herp_main,
                                                wwx_global, wwk_global, cpoti, taui, nx1, ny1, nyp);
            td_pois1 += timer_stop(t_1); // td_pois1 がグローバルとして認識されることを期待

            // 境界条件:
            t_1 = timer_start();
            // 境界条件のコードは省略されていますが、もしあればここに含まれます
            // 例:
            // if (vorfree == 2.0) { /* 境界条件のループ */ }
            td_bnd += timer_stop(t_1); // td_bnd がグローバルとして認識されることを期待

            // delta-f ジャイロ流体分極方程式を解く:
            t_1 = timer_start();
            solve_polarization_equation(ne, gni, pe, hinp_main, herp_main,
                                        wwx_global, wwk_global, cvort, zze, zzi, nx1, ny1, nyp);
            td_pol += timer_stop(t_1); // td_pol がグローバルとして認識されることを期待

            // ジャイロ遮蔽ポテンシャル pi:
            t_1 = timer_start();
            calculate_gyro_shielded_potential(pe, pi, hinp_main, herp_main,
                                            wwx_global, wwk_global, cpoti, taui, nx1, ny1, nyp);
            td_pois2 += timer_stop(t_1); // td_pois2 がグローバルとして認識されることを期待

            // 多段記憶:
            t_1 = timer_start();
            update_history_densities(ne, ne0, ne1, ne2, ni, ni0, ni1, ni2, fne0,
                                     fne1, fne2, fni0, fni1, fni2, nx1, ny1);
            td_mem += timer_stop(t_1); // td_mem がグローバルとして認識されることを期待

            // ローカルなトレースの高解像度ファイル出力:
            if (printtraces) // printtraces がグローバルとして認識されることを期待
            {
                collect_time_trace_data(ttt, ne, pe, iout, jout, k_max,
                                        itrace_local, jtrace_local, ntrace, // ntrace がグローバルとして認識されることを期待
                                        td_out, t_1, printtraces); // td_out がグローバルとして認識されることを期待
            }
        }
        // 低時間分解能診断
        t_1 = timer_start();
        diagnose(ttt, t00, it, itmax, nn, ne, ni, pe, gni); // t00, itmax, nn がグローバルとして認識されることを期待
        td_out += timer_stop(t_1);
    }
}
