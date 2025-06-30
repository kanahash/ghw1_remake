#include "ghw1.hpp"

void collect_time_trace_data(double ttt_in, AXY ne_in, AXY pe_in,
                             int iout_loc[], int jout_loc[], int k_max_loc,
                             int& current_itrace, int& current_jtrace, int ntrace_interval,
                             double& td_out_total, double& t_start_measure,
                             bool printtraces_flag)
{
    if (printtraces_flag)
	{ // printtraces フラグに基づいて実行
        t_start_measure = timer_start(); // td_out の測定を開始

        current_itrace++;
        if (current_itrace >= ntrace_interval)
		{ // ntrace が ntrace_interval に対応
            // ファイル出力はコメントアウトされたまま、データを配列に記録するロジック
            current_itrace = 0; // カウンターをリセット

            // 記録されるトレースステップ数が MAX_TRACE_STEPS を超えないようにチェック
            if (current_jtrace < MAX_TRACE_STEPS)
			{
                for (int k = 0; k < k_max_loc; k++)
				{
                    trace_n_global[k][current_jtrace] = ne_in[iout_loc[k]][jout_loc[k]];
                    trace_p_global[k][current_jtrace] = pe_in[iout_loc[k]][jout_loc[k]];
                }
                current_jtrace++; // 記録ステップ数をインクリメント
            }
			else // 配列の容量を超えた場合の警告またはエラー処理
                printf("Warning: Maximum trace steps reached. Data might be overwritten or truncated.\n");
        }
        td_out_total += timer_stop(t_start_measure); // td_out の測定を停止
    }
}
