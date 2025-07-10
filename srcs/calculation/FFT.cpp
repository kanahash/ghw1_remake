#include "ghw1.hpp" // カスタムヘッダーファイルをインクルード

// ① FFTWの初期化と窓関数の計算
// k_max を引数に追加
static void initialize_fft_for_time_traces(int maxtrace, int k_max, fftw_plan& hindftt, double* pt, fftw_complex* kt, double* func_window)
{
    printf("|\n| %d 個の周波数スペクトル合計を計算中 ... \n", k_max); // k_max は引数として渡される
    double dx0 = static_cast<double>(maxtrace) / 10.0; // 明示的なキャスト

#ifdef _OPENMP
    fftw_plan_with_nthreads(npar); // npar は ghw1.hpp で extern 宣言されているグローバル変数
#endif

    hindftt = fftw_plan_dft_r2c_1d(maxtrace, pt, kt, FFTW_ESTIMATE);

    for (int i = 0; i < maxtrace; i++)
    {
        func_window[i] = 1.0;
        // 明示的なキャストと浮動小数点数リテラルを使用
        func_window[i] -= exp(-(static_cast<double>(i) / dx0) * (static_cast<double>(i) / dx0)) +
                          exp(-(static_cast<double>((maxtrace - i)) / dx0) * (static_cast<double>(maxtrace - i) / dx0));
    }
}

// ② トレースデータのファイル出力
static void write_time_trace_data(double t00, int ntrace, int maxtrace, const double* func_window)
{
    FILE *ft = fopen("traces.dat", "w");
    if (ft == NULL)
    {
        perror("Error opening traces.dat"); // エラーメッセージを出力
        return; // 関数を終了
    }

    for (int i = 0; i < maxtrace; i++)
    {
        fprintf(ft, "%.6e  %.6e  %.6e  %.6e \n",
                (t00 + static_cast<double>(i * ntrace) * dt), // dt は ghw1.hpp で extern 宣言されているグローバル変数
                trace_n_global[1][i] * func_window[i],       // trace_n_global に修正
                trace_p_global[1][i] * func_window[i],       // trace_p_global に修正
                func_window[i]);
    }
    fclose(ft); // ファイルを閉じる
}

// ③ フーリエ周波数スペクトルの計算とファイル出力
static void calculate_and_write_frequency_spectra(int maxtrace, int k_max, double t00, int itmax, int ntrace, fftw_plan hindftt, double* pt, fftw_complex* kt, double* kt_sum, const double* func_window)
{
    FILE *ft;

    // ne のフーリエ周波数スペクトル
    for (int i = 0; i <= maxtrace / 2; i++) kt_sum[i] = 0.0;

    for (int k = 0; k < k_max; k++)
    {
        for (int i = 0; i < maxtrace; i++) // ループ条件を maxtrace に修正 (配列のインデックス)
        {
            pt[i] = trace_n_global[k][i] * func_window[i]; // trace_n_global に修正
        }
        fftw_execute(hindftt); // inner loop is for data prep, execute after loop

        for (int i = 0; i <= maxtrace / 2; i++)
        {
            kt_sum[i] += fabs(kt[i][0]); // 実部の絶対値の合計。必要に応じて複素数の振幅 (sqrt(re*re + im*im)) に変更
        }
    }

    ft = fopen("ptt_n.dat", "w");
    if (ft == NULL)
    {
        perror("Error opening ptt_n.dat");
        return;
    }
    for (int i = 1; i < maxtrace / 2; i++)
        fprintf(ft, "%.6e  %.6e\n", static_cast<double>(i) * TwoPi / (static_cast<double>(itmax) * dt), kt_sum[i] / static_cast<double>(k_max * maxtrace)); // キャスト追加, TwoPi, dt はグローバル
    fclose(ft);

    // phi のフーリエ周波数スペクトル
    for (int i = 0; i <= maxtrace / 2; i++) kt_sum[i] = 0.0;

    for (int k = 0; k < k_max; k++)
    {
        for (int i = 0; i < maxtrace; i++) // ループ条件を maxtrace に修正
        {
            pt[i] = trace_p_global[k][i] * func_window[i]; // trace_p_global に修正
        }
        fftw_execute(hindftt); // inner loop is for data prep, execute after loop

        for (int i = 0; i <= maxtrace / 2; i++)
        {
            kt_sum[i] += fabs(kt[i][0]); // 実部の絶対値の合計
        }
    }

    ft = fopen("ptt_p.dat", "w");
    if (ft == NULL)
    {
        perror("Error opening ptt_p.dat");
        return;
    }
    for (int i = 1; i < maxtrace / 2; i++)
        fprintf(ft, "%.6e  %.6e\n", static_cast<double>(i) * TwoPi / (static_cast<double>(itmax) * dt), kt_sum[i] / static_cast<double>(k_max * maxtrace)); // キャスト追加
    fclose(ft);

    fftw_destroy_plan(hindftt);
}

// 全体をまとめた関数
void perform_time_trace_analysis(int maxtrace, int k_max, double t00, double dt, int itmax, int ntrace)
{
    // 配列宣言。可変長配列 (VLA) はC++標準ではないため、C++11以降の std::vector を使うのが現代的。
    // しかし、このコードが古いCスタイルに基づいているため、VLAがコンパイラ拡張として許可されている場合が多い。
    // 安全のため、VLAではなく固定サイズ配列や動的確保に切り替えることも検討
    double pt[maxtrace], func_window[maxtrace];
    fftw_complex kt[maxtrace / 2 + 1];
    double kt_sum[maxtrace / 2 + 1];

    fftw_plan hindftt;

    // k_max を引数として渡す
    initialize_fft_for_time_traces(maxtrace, k_max, hindftt, pt, kt, func_window);
    // t00, ntrace, maxtrace, dt はグローバル変数または引数として渡す
    write_time_trace_data(t00, ntrace, maxtrace, func_window);
    // k_max, t00, itmax, ntrace, dt はグローバル変数または引数として渡す
    calculate_and_write_frequency_spectra(maxtrace, k_max, t00, itmax, ntrace, hindftt, pt, kt, kt_sum, func_window);
}
