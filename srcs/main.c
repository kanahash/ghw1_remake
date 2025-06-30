#include "ghw1.hpp"

int	main(void)
{
	double	c0 = 18. / 11., c1;
	double	td_init = 0., td_update = 0., td_pois1;
	double	td_bnd = 0., td_mem;
	double	bndys[nxm];
	double	ne_za[nxm], ni_za[nxm];
	char	s[80], str[80];
	int		nx_i;
	double	dx0;
	int		rs_check;

	c0 = 18. / 11., c1 = 9. / 11., c2;
	td_init = 0., td_update = 0., td_pois1 = 0., td_pois2;
	td_bnd = 0., td_mem = 0., td_end;
	c0 = 18. / 11., c1 = 9. / 11., c2 = 2. / 11., cf;
	td_init = 0., td_update = 0., td_pois1 = 0., td_pois2 = 0., td_out;
	td_bnd = 0., td_mem = 0., td_end = 0., td_pol;
	int i, j, k, l, m, n, ik, jk, it, is, itn, itend;
	int iout[kt_maxsqrt * kt_maxsqrt], jout[kt_maxsqrt * kt_maxsqrt], k_max;
	int im, ip, jm, jp, im1, ip1, im2, ip2, jm2, jp2, icc, ixloc, itrace,
		jtrace;
	c0 = 18. / 11., c1 = 9. / 11., c2 = 2. / 11., cf = 6. / 11.;
	double ttt, t00, phase, zuf, nue, nui, nuzf;
	double neprof, pbar, nbar, dve, gre, gri;
	// プログラムの各部の時間消費を測定するための変数
	td_init = 0., td_update = 0., td_pois1 = 0., td_pois2 = 0., td_out = 0.;
	td_bnd = 0., td_mem = 0., td_end = 0., td_pol = 0., td_tot;
	AXY ane, ani, dummy;
	char const *fn, *fw;
	FILE *f, *g, *h, *g1, *g2, *g3, *g4, *g5, *ft;
	t_1 = timer_start(); // td_init の測定を開始 (コードの各部分の相対的なパフォーマンスを分析するため)
	// 入力ファイルを読み込み、パラメータを初期化する:
	init_parameters();
	// 並列化されたFFTWを初期化する:
	init_FFTW(nx, ny, fflag, npar, hinp_main, herp_main, wisdom_sf_main);
	// オプション: ディリクレ境界条件を定義する (選択した場合):
	for (i = 0; i <= nx1; i++)
	{
		nx_i = double(nx1 - i);
		dx0 = 2. * 2.;
		bndys[i] = 1. - exp(-nx_i * nx_i / dx0) - exp(-double((i - 0) * (i - 0))
				/ dx0);
	}
	// 初期密度摂動を設定する (ブロブ、渦、乱流、流れなど)
	set_init_density_perturbation(incon, ne, ni);
	// 初期ポテンシャルを計算する
	init_potential(ni, ne, pe, gni, ww, cpoti, cvort, aae, aai);
	// 人為的な「以前の」時間値を設定する (多段ソルバー用)
	init_history_density(ne, ni, ne0, ne1, ne2, ni0, ni1, ni2, fne1, fne2, fni1,
		fni2);
	// エネルギーを初期化する:
	init_energy(ne, xyz, eno);
	if (incon == 2.)
	{ // 以前の実行の保存された終了出力からリスタート
		printf("| データセットを継続します... \n");
		rs_check = 1;
		g = fopen("restart.dat", "r");
		if (g == NULL)
		{
			printf("\n ファイル 'restart.dat' が見つかりません。終了します... \n");
			exit(1);
		}
		fscanf(g, "%s ", str);
		t00 = atof(str);
		while (!feof(g))
		{
			rs_check *= fscanf(g, "%s ", str);
			i = atoi(str);
			rs_check *= fscanf(g, "%s ", str);
			j = atoi(str);
			rs_check *= fscanf(g, "%s ", str);
			eno = atof(str);
			rs_check *= fscanf(g, "%s ", str);
			pe[i][j] = atof(str);
			rs_check *= fscanf(g, "%s ", str);
			ne[i][j] = atof(str);
			rs_check *= fscanf(g, "%s ", str);
			ne1[i][j] = atof(str);
			rs_check *= fscanf(g, "%s ", str);
			ne2[i][j] = atof(str);
			rs_check *= fscanf(g, "%s ", str);
			ni[i][j] = atof(str);
			rs_check *= fscanf(g, "%s ", str);
			ni1[i][j] = atof(str);
			rs_check *= fscanf(g, "%s ", str);
			ni2[i][j] = atof(str);
			rs_check *= fscanf(g, "%s ", str);
			fne1[i][j] = atof(str);
			rs_check *= fscanf(g, "%s ", str);
			fne2[i][j] = atof(str);
			rs_check *= fscanf(g, "%s ", str);
			fni1[i][j] = atof(str);
			rs_check *= fscanf(g, "%s ", str);
			fni2[i][j] = atof(str);
		}
		fclose(g);
		poisson(pe, cpoti, pi);
		for (i = 0; i <= nx1; i++)
		{
			for (j = 0; j <= ny1; j++)
			{
				ne0[i][j] = ne[i][j];
				ni0[i][j] = ni[i][j];
			}
		}
	}
	// トレース出力位置を設定する
	set_trace_output_location(nx1, ny1, kt_maxsqrt, iout, jout, k_max);
	// 時間ステップを初期化する:
	if (init_time_stepping(itstp, itmax, dt, cf, itn, itend, dtt, ddtt, itrace,
			jtrace) == 0)
		return (0); // シミュレーションを終了
	// 初期プロファイルの出力制御
	diagnose(0., 0., it, itmax, nn, ne, ni, pe, gni);
	td_init = timer_stop(t_1); // td_init の測定を停止
	// 時間ステップ  --------------------------------------------------------------
	// 3次 Karniadakis 時間スキーム; 対流項には Arakawa ブラケット
	for (it = itn / itstp; it < itend; ++it)
	{
		// 全体の時間ループ (出力ステップ)
		for (is = 0; is < itstp; ++is)
		{
			// 内部時間ループ (出力なし)
			ttt = t00 + (it + 0) * ddtt + (is + 1) * dt;
			// ジャイロ流体密度 ne と ni の時間ステップ更新:
			// td_update の測定は update_densities_one_time_step 内で行われる
			update_densities_one_time_step(ttt, dt, ne, ni, pe, pi, ane, ani,
				hyve, hyvi, vis_e, vis_i, fftw_plan_hinp, fftw_plan_herp,
				wwx_global, wwk_global, chyv, b_4th, ihype, nx1, ny1, nyp);
			td_update += timer_stop(t_1); // td_update の測定を停止
			// テスト: 帯状密度の減衰:
			/*
			nuzf = 0.0; // 0.02
			nuzf *= dtt/double(ny);
		#pragma omp parallel for private(j, ne_za, ni_za)
			for (i=0; i<=nx1; ++i) {
			ne_za[i] = 0.; ni_za[i] = 0.;
			// ne および ni の帯状平均:
			for (j=0; j<=ny1; ++j) { ne_za[i] += ne[i][j]; ni_za[i]
				+= ni[i][j]; };
			// 帯状密度の割合で ne と ni を削減する:
			for (j=0; j<=ny1; ++j) { ne[i][j] -= nuzf*ne_za[i]; ni[i][j]
				-= nuzf*ni_za[i]; };
			}
			*/
			// ni から gni へのジャイロ平均: -------------------------------
			t_1 = timer_start(); // td_pois1 の測定を開始
			gyro_averaged_ion_density_calculate(ni, gni, hinp_main, herp_main,
				wwx_global, wwk_global, cpoti, taui, nx1, ny1, nyp);
			td_pois1 += timer_stop(t_1); // td_pois1 の測定を停止
			// 境界条件 (TIFF コードの局所境界条件と一致するチャネル形状):
			t_1 = timer_start(); // td_bnd の測定を開始
			if (vorfree == 2.)
			{
#pragma omp parallel for private(j)
				for (i = 0; i <= nx1; i++)
				{
					for (j = 0; j <= ny1; j++)
					{
						ne[i][j] *= bndys[i];
						ni[i][j] *= bndys[i];
						gni[i][j] *= bndys[i];
					}
				}
			}
			td_bnd += timer_stop(t_1); // td_bnd の測定を停止
			// delta-f ジャイロ流体分極方程式を解く: ----------------------------
			// k空間でのFFTにより渦度から電位 phi を計算する
			t_1 = timer_start(); // td_pol の測定を開始
			solve_polarization_equation(ne, gni, pe, hinp_main, herp_main,
				wwx_global, wwk_global, cvort, zze, zzi, nx1, ny1, nyp);
			td_pol += timer_stop(t_1); // td_pol の測定を停止
			// ジャイロ遮蔽ポテンシャル pi: ---------------------------------------------
			t_1 = timer_start(); // td_pois2 の測定を開始
			gyro_shielded_potential_calculate(pe, pi, hinp_main, herp_main,
				wwx_global, wwk_global, cpoti, taui, nx1, ny1, nyp);
			td_pois2 += timer_stop(t_1); // td_pois2 の測定を停止
			// 多段記憶: 2ステップ前の値を記憶する -------------------
			t_1 = timer_start();
			// td_mem の測定を開始
			update_history_densities(ne, ne0, ne1, ne2, ni, ni0, ni1, ni2, fne0,
				fne1, fne2, fni0, fni1, fni2, nx1, ny1); // nx1,
			// ny1 はそれぞれ nx_max, ny_max に対応
			td_mem += timer_stop(t_1);
			// td_mem の測定を停止
			// ローカルなトレースの高解像度ファイル出力 (オプション、頻繁な I/O は速度を低下させる)
			// 時間トレース解析用 (例: ある x,y での pe(t) から周波数スペクトルを計算する)
			// /*
			if (printtraces)
			{
				// `t_1 = timer_start();` と `td_out += timer_stop(t_1);` は
				// `collect_time_trace_data` の内部に移動されるか、
				// `td_out` を参照渡しで関数に渡し、関数内で加算するように調整されます。
				// `t_1` は `collect_time_trace_data` 関数の引数としても渡すことが可能。
				collect_time_trace_data(ttt,
										ne,
										pe,
										iout,
										jout,
										k_max,
										itrace_local,
										jtrace_local,
										ntrace,
										td_out,
										t_1,
										// ここで t_1 を計測開始のダミーとして渡し、td_out に加算
										printtraces);
				// グローバルな printtraces フラグを使用
			}
			// */
		}
		// ... 内部時間ループの終了 ..............................................
		// 低時間分解能診断: エネルギーとスナップショット出力
		t_1 = timer_start(); // td_out の測定を開始
		diagnose(ttt, t00, it, itmax, nn, ne, ni, pe, gni);
		td_out += timer_stop(t_1); // td_out の測定を停止
	}
	// ... 外部時間ループの終了 ..........................................
	// 時間ステップ終了 --------------------------------------------------------
	t_1 = timer_start(); // td_out の測定を開始
	// 完全なリスタートファイルを書き込む
	g = fopen("restart.dat", "w");
	fprintf(g, "%5e  \n", ttt);
	for (i = 0; i <= nx1; i++)
		for (j = 0; j <= ny1; j++)
			fprintf(g,
					"%d  %d  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e  %.6e 
						% .6e % .6e % .6e %
						.6e\n ", i, j, eno, pe[i][j], ne[i][j],
						ne1[i][j],
					ne2[i][j],
					ni[i][j],
					ni1[i][j],
					ni2[i][j],
					fne1[i][j],
					fne2[i][j],
					fni1[i][j],
					fni2[i][j]);
	fclose(g);
	// 時間トレース / 周波数スペクトル解析:
	if (printtraces)
		perform_time_trace_analysis(maxtrace, k_max, t00, dt, itmax, ntrace);
	td_out += timer_stop(t_1); // td_out の測定を停止
	// td_out += t_2 - t_1;
	td_tot = td_init + td_update + td_pois1 + td_pois2 + td_pol + td_bnd
		+ td_out + td_mem;
	printf("|\n| 総実行時間:  %.2f s  =  %.2f min \n", td_tot, td_tot / 60.);
	td_tot = 100. / td_tot;
	// 総時間ボトルネック分析出力:
	printf("|\n| 絶対実行時間: init: %.2f upd: %.2f poi1: %.2f poi2: %.2f pol:
				% .2f bnd
			: % .2f out
			: % .2f mem
			: %
					.2f\n ", td_init, td_update, td_pois1,
					td_pois2,
				td_pol,
				td_bnd,
				td_out,
				td_mem);
	// 割合ボトルネック分析出力:
	printf("| 相対実行時間: init: %.2f upd: %.2f poi1: %.2f poi2: %.2f pol: %.2f bnd:
				% .2f out
			: % .2f mem
			: % .2f\n |\n ", td_init * td_tot, td_update * td_tot,
						td_pois1 *
						td_tot,
				td_pois2 * td_tot,
				td_pol * td_tot,
				td_bnd * td_tot,
				td_out * td_tot,
				td_mem * td_tot);
#ifdef _OPENMP
	fftw_cleanup_threads();
#endif
	printf("| GHW1 終了。\n\n");
}
