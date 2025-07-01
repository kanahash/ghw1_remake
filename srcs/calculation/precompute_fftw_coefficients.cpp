#include "ghw1.hpp"

// FFTW ポアソンソルバーの係数を事前計算する
void precompute_fftw_coefficients(int nx, int ny, int nxh, int nyh, double hy, double taui, double mue, double taue, double mui, double aae, double aai, int ipade, bool xlwl)
{
    int ik, jk;
    double kx, ky; // kx, kyは使用されていないようですが、元のコードにあったため残します
    double dxy = 1. / (hy * hy); // この変数は使われていませんが、元のコードにあったため残します
    double cnx = 1. / double(nx);
    double cny = 1. / double(ny);
    double cxy = cnx * cny;
    double kkqq, cinv, bki, bke, taui0;
    double gambese, gambesi; // gambese, gambesiはコメントアウトされた部分でしか使われていませんが、元のコードにあったため残します
    double G0i, G0e, G1i, G1e;

    taui0 = (taui == 0.) ? 1.e-9 : taui;
    double mue0 = (mue == 0.) ? 1.e-9 : mue;

    // for (int i=0; i<=nx1; ++i) for (int j=0; j<=ny1; ++j) { // 元のコメント
    for (int i = 0; i <= nx - 1; ++i)
    { 
        // nx1の代わりにnx-1を使用
        for (int j = 0; j <= nyh; ++j)
        { 
            // nyh+1-1はnyhと等しい
            ik = (i >= nxh) ? i - nx : i;
            jk = (j >= nyh) ? j - ny : j;

            // k_perp^2:
            kkqq = (TwoPi * hy * ik * cnx) * (TwoPi * hy * ik * cnx) + (TwoPi * hy * jk * cny) * (TwoPi * hy * jk * cny);

            bki = taui0 * mui * kkqq;
            bke = taue * mue0 * kkqq;

            // Gamma 演算子 (ジャイロ平均とジャイロ遮蔽):
            if (ipade == 0)
            { 
                // Pade 近似
                G0i = 1. / (1. + xlwl * bki);
                G1i = 1. / (1. + 0.5 * bki);
                G0e = 1. / (1. + xlwl * bke);
                G1e = 1. / (1. + 0.5 * bke);
                cinv = kkqq * (aai * mui * G0i + aae * mue * G0e);
            }
            else if (ipade == 1)
            { 
                // ベッセル関数
                // gamma0 関数がどこかで定義されている必要があります
                G0i = gamma0(bki);
                G1i = sqrt(G0i);
                G0e = gamma0(bke);
                G1e = sqrt(G0e);
                cinv = (aai / taui0) * (1. - G0i) + (aae / taue) * (1. - G0e);
            }

            // グローバル配列に値を設定（これらの配列は関数外で定義されていると仮定します）
            cvort[i][j] = -cxy / cinv;
            cpoti[i][j] = +cxy * G1i;
            cvfis[i][j] = +cxy / G1i;

            // ラプラス
            clap[i][j] = -cxy * kkqq;

            // ハイパー粘性
            if (ihype == 1) 
                chyv[i][j] = +cxy * kkqq; // Navier-Stokes 型粘性
            else if (ihype == 2) 
                chyv[i][j] = cxy * kkqq * kkqq; // k^4 「通常の」ハイパー粘性
            else if (ihype == 4) 
                chyv[i][j] = cxy * kkqq * kkqq * kkqq * kkqq; // k^8 強化
            else if (ihype == 8) 
                chyv[i][j] = cxy * kkqq * kkqq * kkqq * kkqq * kkqq * kkqq; // k^12 (元のコメントではk^8ですが、計算式はk^12)
        }
    }
    cvort[0][0] = 0.;
}
