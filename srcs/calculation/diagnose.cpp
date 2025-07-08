#include "ghw1.hpp"

// 出力量を計算し書き出す:
// 例: エネルギー、輸送、k スペクトル、2D 配列、プロファイルなど
void diagnose( double ttt, double t00, int it, int itmax,
           AXY n_n, AXY n_e, AXY n_i, AXY p_p, AXY n_g )
{
  int i,j, ik, im,ip,jm,jp;
  double enn, enp, enw, eeb, ezf, fne, fnsol, zfx, dum, rey, freq, grow;
  double etran_grad, etran_adia, etran_visc;
  double lap[nxm][nym];
  double avg[nxm], avp[nxm];
  FILE *f, *g, *h, *g1, *g2, g3, *g4, *g5, *g6, *g7, *g8;

#pragma omp parallel for private(j) shared(avg,avp)
  for (i=0; i<=nx1; i++) // ne と phi の帯状平均
    { 
		avg[i] = 0.; avp[i] = 0.;
      	for (j=0; j<=ny1; j++)
			avg[i] += n_e[i][j]; avp[i] += p_p[i][j];
      	avg[i]/=ny; avp[i]/=ny;
    }
  enn = 0.; enp = 0.; enw = 0.; eeb = 0.; ezf = 0.; fne = 0.;
  etran_grad = 0.; etran_adia= 0.; etran_visc = 0.;

  // 渦度
  // laplace(p_p,lap);
  poisson(p_p,clap,lap);

   // エネルギー量と輸送
  calculate_energies_and_transport(enn, enp, enw, eeb, ezf, fne,
                                   etran_grad, etran_adia, etran_visc,
                                   n_e, n_i, p_p, lap,
                                   avg, avp, hy, taui,
                                   pe, pi, b_mhw, chat, diff,
                                   hyve, hyvi, nx1, ny1);

  enn *= .5*xyz; enp *= xyz; enw *= .5*xyz; eeb *= .5*xyz; ezf *= .5*xyz;
  if (enw==0.) enw=1.e-12; if (ezf==0.) ezf=1.e-12; if (eeb==0.) eeb=1.e-12;

  eeb*=delinv*delinv; ezf*=delinv*delinv; enw*=delinv*delinv; enn*=delinv*delinv;
  fne*= xyz*delinv*delinv;

  etran_grad*= xyz; etran_adia*= xyz; etran_visc*= xyz;
  double etran_tot = etran_grad + etran_adia + etran_visc;

  // クラッシュ制御: nan または inf の場合停止
  // if ( (isnan(enn)) || (enn<1.e-16) )
  if (isnan(enn))
    {
      printf("\n|t=%.3f (%.2f): enn=%.3e\n", ttt,(t00+itmax*dt),enn);
      printf("|HGW1 クラッシュしました！  enn=%.3e,   eno=%.3e\n", enn, eno);
#ifdef _OPENMP
      fftw_cleanup_threads();
#endif
      exit(1);
    }
  eno = enn;

   // グローバルエネルギー時系列出力と輸送時系列出力
  write_global_outputs(ttt, enn, eeb, ezf, etran_grad, etran_adia, etran_visc, etran_tot, fne);

  // (非)線形成長率 "grow" と線形周波数推定 "freq":
  calculate_linear_growth_and_frequency(ttt, dt, it, itmax, incon,
                                        freq, grow,
                                        p_p, t_amp, t_old, jsgn, jsgn_old,
                                        ddtt, enw, enwo, nxh, nyh);

  // プロット出力

  // 帯状流 (t,x) 2D プロット
  write_zonal_flow_data(ttt, p_p, nx1, ny1, hy);

  // y=ny/2 での x-カット:
  write_x_cut_data(n_e, n_i, p_p, lap, n_g, nx1, nyh, hy);

  // x=nx/2 での y-カット:
  write_y_cut_data(n_e, n_i, p_p, lap, n_g, ny1, nxh, hy);

  // y で平均された x プロファイル:
  write_x_profile_data(p_p, n_e, n_i, lap, nx1, ny1, hy, ny);

  // 密度、渦度、ポテンシャルの 2D (x,y) プロット
  write_2d_plot_data(n_e, p_p, lap, nx1, ny1);

  // フーリエ ky スペクトル
  double py[ny];
  // fftw_complex ky[ny/2], sumky[ny/2];
  fftw_complex ky[ny/2+1], sumky[ny/2+1];
#ifdef _OPENMP
  fftw_plan_with_nthreads(npar);
#endif
  fftw_plan hindfty;
  hindfty = fftw_plan_dft_r2c_1d(ny, &py[0], &ky[0], FFTW_ESTIMATE);

// ポテンシャル phi
  calculate_and_write_ky_spectrum(p_p, "pky_p.dat",
                                  hindfty, py, ky,
                                  sumky[0], pkyavg, // sumkyとpkyavgは適切な次元で渡す
                                  nx1, ny1, ny, nx, incon, it);
  
  // 密度 ne
  calculate_and_write_ky_spectrum(n_e, "pky_n.dat",
                                  hindfty, py, ky,
                                  sumky[0], pkyavgn,
                                  nx1, ny1, ny, nx, incon, it);

 // 渦度 w
  // 渦度はky[j][0]*ky[j][0]を使うので、別途引数を追加するか、
  // 関数を2種類作るか、内部で条件分岐させる必要があります。
  // ここでは簡単な例としてfabsを仮定しています。
  calculate_and_write_ky_spectrum(lap, "pky_w.dat",
                                  hindfty, py, ky,
                                  sumky[0], pkyavgw,
                                  nx1, ny1, ny, nx, incon, it);

  // 運動エネルギー E
  calculate_and_write_energy_ky_spectrum(p_p, "pky_e.dat",
                                         hindfty, py, ky,
                                         sumky[0], pkyavge, // sumkyとpkyavgeは適切な次元で渡す
                                         nx1, ny1, ny, nx,
                                         hy, TwoPi, ly,
                                         incon, it);

  // 他のKyスペクトル計算がすべて終わった後で、一度だけ hindfty を破棄します
  fftw_destroy_plan(hindfty);

 // フーリエ kx スペクトル
  calculate_and_write_kx_spectrum(pe, "pkx_p.dat",
                                  nx1, ny1, nx, ny,
                                  hy, TwoPi, ly,
                                  incon, it,
                                  pkxavg); // pkxavg は適切な次元で渡す

#ifdef _OPENMP
  fftw_cleanup_threads();
#endif
  // スペクトル解析の終了

  // タイムスタンプ出力
  printf("|t=%.3f (%.3f): enn=%.3e\n", ttt,(t00+itmax*dt),enn);

}
