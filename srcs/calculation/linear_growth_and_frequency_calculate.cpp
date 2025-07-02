#include "ghw1.hpp"

// 線形ドリフト波の成長率と周波数を計算し、ファイルに出力する関数
void calculate_linear_growth_and_frequency(double ttt, double dt, int it, int itmax, int incon,
                                           double& freq, double& grow,
                                           AXY p_p, double& t_amp, double& t_old, int& jsgn, int& jsgn_old,
                                           double& ddtt, double& enw, double& enwo, int nxh, int nyh)
{
    // (近似) 線形ドリフト波周波数 freq を計算する
    if (incon == 5.)
    {
        double dt_amp = 0.;
        freq = 0.;
        grow = 0.;
        jsgn_old = jsgn;
        if (p_p[nxh][nyh] == fabs(p_p[nxh][nyh])){jsgn = 1;} else {jsgn = -1;};

        if (jsgn != jsgn_old) {t_old = t_amp; t_amp = ttt;};
        dt_amp = t_amp - t_old;
        if ((dt_amp != 0.) && (dt_amp != ddtt) && (it > 2)) freq = 0.5 * TwoPi / dt_amp;
    }

    if (enwo == 1.e-12) enwo = enw;
    grow = (sqrt(enw) - sqrt(enwo)) / (sqrt(enwo) * ddtt);
    enwo = enw;

    FILE *f = fopen("gamma.dat", "a");
    if (it > 1) fprintf(f, "%.3f  %.6e  %.6e\n", ttt, freq, grow);
    fclose(f);
}
