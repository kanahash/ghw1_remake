#include "ghw1.hpp"

// エネルギーと輸送量を計算する関数
void calculate_energies_and_transport(double& enn, double& enp, double& enw, double& eeb, double& ezf, double& fne,
                                      double& etran_grad, double& etran_adia, double& etran_visc,
                                      AXY n_e, AXY n_i, AXY p_p, double lap[][nym],
                                      double avg[], double avp[], double hy, double taui,
                                      AXY pe, AXY pi, bool b_mhw, double chat, double diff,
                                      AXY hyve, AXY hyvi, int nx1, int ny1)
{
    int i, j, ip;
    double dum;

    // 変数を初期化
    enn = 0.; enp = 0.; enw = 0.; eeb = 0.; ezf = 0.; fne = 0.;
    etran_grad = 0.; etran_adia= 0.; etran_visc = 0.;

#pragma omp parallel for private(j,ip,dum) reduction(+:enn,ezf,eeb,fne,enw,etran_grad,etran_adia,etran_visc)
    for (i = 0; i <= nx1; i++)
    {
        ip = (i == nx1) ? 0 : i + 1;

        for (j = 0; j <= ny1; j++)
        {
            // 全熱自由エネルギー:
            enn += n_e[i][j] * n_e[i][j] + taui * n_i[i][j] * n_i[i][j];

            // (主に帯状) 流れのエネルギー:
            dum = hy * (avp[ip] - avp[i]);
            ezf += dum * dum;

            // 乱流エントロフィー:
            enw += lap[i][j] * lap[i][j];

            // 全運動エネルギー:
            dum = DFDY(p_p, i, j);
            eeb += dum * dum;
            dum = DFDX(p_p, i, j);
            eeb += dum * dum;

            // 電子 ExB 粒子輸送:
            fne += -n_e[i][j] * DFDY(p_p, i, j);

            // 勾配駆動 / 転送:
            etran_grad += -n_e[i][j] * DFDY(pe, i, j);
            etran_grad += -n_i[i][j] * DFDY(pi, i, j) * taui;

            // 「並列」抵抗散逸シンク / 転送:
            dum = pe[i][j] - n_e[i][j];
            if (b_mhw) dum -= avp[i] - avg[i];
            etran_adia += -chat * dum * dum;

            // 垂直ハイパー粘性シンク / 転送:
            dum = hyve[i][j] * (pe[i][j] - n_e[i][j]);
            dum -= hyvi[i][j] * (pi[i][j] + taui * n_i[i][j]);
            etran_visc += diff * dum;
        }
    }
}
