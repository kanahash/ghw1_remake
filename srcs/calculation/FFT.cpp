#include "ghw1.hpp"

// ① FFTWの初期化と窓関数の計算
static void initialize_fft_for_time_traces(int maxtrace, fftw_plan& hindftt, double* pt, fftw_complex* kt, double* func_window)
{
    printf("|\n| %d 個の周波数スペクトル合計を計算中 ... \n", k_max); // k_max はグローバルまたは引数として渡す
    double dx0 = double(maxtrace)/10.;
#ifdef _OPENMP
    fftw_plan_with_nthreads(npar); // npar はグローバル変数または引数として渡す
#endif
    hindftt = fftw_plan_dft_r2c_1d(maxtrace, pt, kt, FFTW_ESTIMATE);

    for (int i=0; i<maxtrace; i++)
    {
        func_window[i] = 1.;
        func_window[i] -= exp(-(double(i)/dx0)*(double(i)/dx0)) + exp(-(double((maxtrace-i)/dx0)*(double(maxtrace-i)/dx0)));
    }
}

// ② トレースデータのファイル出力
static void write_time_trace_data(double t00, int ntrace, int maxtrace, const double* func_window)
{
    FILE *ft = fopen("traces.dat","w");
    if (ft == NULL)
	{
        perror("Error opening traces.dat");
        return;
    }
    for (int i=0; i<maxtrace; i++)
    {
        fprintf(ft,"%.6e  %.6e  %.6e  %.6e \n",
                (t00 + double(i*ntrace)*dt), trace_n[1][i]*func_window[i], trace_p[1][i]*func_window[i], func_window[i]); // dt はグローバル変数または引数として渡す
    }
    fclose(ft);
}

// ③ フーリエ周波数スペクトルの計算とファイル出力
static void calculate_and_write_frequency_spectra(int maxtrace, int k_max, double t00, int itmax, int ntrace, fftw_plan hindftt, double* pt, fftw_complex* kt, double* kt_sum, const double* func_window)
{
    FILE *ft;

    // ne のフーリエ周波数スペクトル
    for (int i=0; i<=maxtrace/2; i++) kt_sum[i] = 0.;
    for (int k=0; k<k_max; k++)
    {
        for (int i=0; i<=maxtrace; i++)
		{
			pt[i] = trace_n[k][i]*func_window[i];
        	fftw_execute(hindftt);
		}
        for (int i=0; i<=maxtrace/2; i++)
			kt_sum[i] += fabs(kt[i][0]);
    }
    ft = fopen( "ptt_n.dat", "w" );
    if (ft == NULL)
	{
        perror("Error opening ptt_n.dat");
        return;
    }
    for (int i=1; i<maxtrace/2; i++)
        fprintf(ft,"%.6e  %.6e\n",double(i)*TwoPi/(double(itmax)*dt),kt_sum[i]/double(k_max*maxtrace));
    fclose( ft );

    // phi のフーリエ周波数スペクトル
    for (int i=0; i<=maxtrace/2; i++) kt_sum[i] = 0.;
    for (int k=0; k<k_max; k++)
    {
        for (int i=0; i<=maxtrace; i++)
		{
			pt[i] = trace_p[k][i]*func_window[i];
        	fftw_execute(hindftt);
		}
        for (int i=0; i<=maxtrace/2; i++)
			kt_sum[i] += fabs(kt[i][0]);
    }
    ft = fopen( "ptt_p.dat", "w" );
    if (ft == NULL)
	{
        perror("Error opening ptt_p.dat");
        return;
    }
    for (int i=1; i<maxtrace/2; i++)
        fprintf(ft,"%.6e  %.6e\n",double(i)*TwoPi/(double(itmax)*dt),kt_sum[i]/double(k_max*maxtrace));
    fclose( ft );

    fftw_destroy_plan(hindftt);
}

// 全体をまとめた関数
void perform_time_trace_analysis(int maxtrace, int k_max, double t00, double dt, int itmax, int ntrace)
{
    double pt[maxtrace], func_window[maxtrace];
    fftw_complex kt[maxtrace/2+1]; // maxtrace/2+1 に変更、FFTW の r2c 出力配列のサイズに合わせる
    double kt_sum[maxtrace/2+1]; // maxtrace/2+1 に変更

    fftw_plan hindftt;

    initialize_fft_for_time_traces(maxtrace, hindftt, pt, kt, func_window);
    write_time_trace_data(t00, ntrace, maxtrace, func_window);
    calculate_and_write_frequency_spectra(maxtrace, k_max, t00, itmax, ntrace, hindftt, pt, kt, kt_sum, func_window);
}
