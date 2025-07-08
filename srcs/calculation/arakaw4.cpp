#include "ghw1.hpp"

// 2次精度のアラカワ・ヤコビアンの各項の合計を計算する関数
// スケールファクター (.25 * hy * hy / 3.) はこの関数内で乗じる
static double calculate_arakawa_jacobian_2nd_order(
    AXY uuu, AXY vvv,
    int i00, int j00,
    int im1, int ip1, int jm1, int jp1,
    double hy
)
{
    double aj2_sum = 0.0;

    // J++ :
    aj2_sum  = (uuu[ip1][j00] - uuu[im1][j00]) * (vvv[i00][jp1] - vvv[i00][jm1]);
    aj2_sum -= (uuu[i00][jp1] - uuu[i00][jm1]) * (vvv[ip1][j00] - vvv[im1][j00]);

    // J+x :
    aj2_sum += uuu[ip1][j00] * ( vvv[ip1][jp1] - vvv[ip1][jm1] );
    aj2_sum -= uuu[im1][j00] * ( vvv[im1][jp1] - vvv[im1][jm1] );
    aj2_sum -= uuu[i00][jp1] * ( vvv[ip1][jp1] - vvv[im1][jp1] );
    aj2_sum += uuu[i00][jm1] * ( vvv[ip1][jm1] - vvv[im1][jm1] );

    // Jx+ :
    aj2_sum += uuu[ip1][jp1] * ( vvv[i00][jp1] - vvv[ip1][j00] );
    aj2_sum -= uuu[im1][jm1] * ( vvv[im1][j00] - vvv[i00][jm1] );
    aj2_sum -= uuu[im1][jp1] * ( vvv[i00][jp1] - vvv[im1][j00] );
    aj2_sum += uuu[ip1][jm1] * ( vvv[ip1][j00] - vvv[i00][jm1] );

    // スケールファクターを乗じる
    return aj2_sum * (.25 * hy * hy / 3.);
}

// 4次精度のアラカワ・ヤコビアンの各項の合計を計算する関数
// スケールファクター (.5 * afac) はこの関数内で乗じる
static double calculate_arakawa_jacobian_4th_order(
    AXY uuu, AXY vvv,
    int i00, int j00,
    int im1, int ip1, int jm1, int jp1,
    int im2, int ip2, int jm2, int jp2, // 4次項に必要な追加のインデックス
    double hy // スケールファクターの一部
)
{
    double aj4_sum = 0.0;
    double afac = .25 * hy * hy / 3.; // 2次部分と同じafacをここで再計算

    // Jxx ;
    aj4_sum  = (uuu[ip1][jp1] - uuu[im1][jm1]) * (vvv[im1][jp1] - vvv[ip1][jm1]);
    aj4_sum -= (uuu[im1][jp1] - uuu[ip1][jm1]) * (vvv[ip1][jp1] - vvv[im1][jm1]);

    // J4x+ :
    aj4_sum += uuu[ip1][jp1] * ( vvv[i00][jp2] - vvv[ip2][j00] );
    aj4_sum -= uuu[im1][jm1] * ( vvv[im2][j00] - vvv[i00][jm2] );
    // "Arakawa 1966 で "+" が間違っている！" のコメントに注意。
    // 元の論文との比較で修正が必要な場合があるかもしれません。
    aj4_sum -= uuu[im1][jp1] * ( vvv[i00][jp2] - vvv[im2][j00] );
    aj4_sum += uuu[ip1][jm1] * ( vvv[ip2][j00] - vvv[i00][jm2] );

    // J4+x :
    aj4_sum += uuu[ip2][j00] * ( vvv[ip1][jp1] - vvv[ip1][jm1] );
    aj4_sum -= uuu[im2][j00] * ( vvv[im1][jp1] - vvv[im1][jm1] );
    aj4_sum -= uuu[i00][jp2] * ( vvv[ip1][jp1] - vvv[im1][jp1] );
    aj4_sum += uuu[i00][jm2] * ( vvv[ip1][jm1] - vvv[im1][jm1] );

    // スケールファクターを乗じる
    return aj4_sum * (.5 * afac);
}

// ブラケットの 4次 Arakawa スキーム (Arakawa 1966, JCP 1997 に再録を参照)
// [u,v] = (du/dx)(dv/dy)-(du/dy)(dv/dx)
// void arakaw4(double (*uuu)[nym], double (*vvv)[nym], double (*out)[nym])
void arakaw4( AXY uuu, AXY vvv, AXY out )
{
  int i00,j00,ip1,jp1,im1,jm1,im2,ip2,jm2,jp2;
  double aj2, aj4, afac = .25*hy*hy/3.;
  int i0 = 0, i1 = nx1, j0 = 0, j1 = ny1;

#pragma omp parallel for private(j00,im1,im2,jm1,jm2,ip1,ip2,jp1,jp2,aj2,aj4)
  for (i00=0; i00<=nx1; ++i00)  
  {
    nextp(i00,i0,i1,im1,ip1,im2,ip2);

    for (j00=0; j00<=ny1; ++j00) 
	{
      nextp(j00,j0,j1,jm1,jp1,jm2,jp2);

      // 2次部分の計算を関数に委譲
      double aj2 = calculate_arakawa_jacobian_2nd_order(
                       uuu, vvv,
                       i00, j00,
                       im1, ip1, jm1, jp1,
                       hy
                   );

      // 4次部分の計算を関数に委譲
      double aj4 = calculate_arakawa_jacobian_4th_order(
                       uuu, vvv,
                       i00, j00,
                       im1, ip1, jm1, jp1,
                       im2, ip2, jm2, jp2, // im2, ip2, jm2, jp2 を渡す
                       hy
                   );

      // 合計 (4次):
      out[i00][j00] = 2.*aj2 - aj4 ;

      // テストのためにコメント解除: 2次部分のみ:
      // out[i00][j00] = aj2;
    };
  }
} // arakaw4 関数の終了
