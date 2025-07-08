#include "ghw1.hpp"

static double calculate_arakawa_jacobian_at_point(
    AXY uuu, AXY vvv,
    int i00, int j00,
    int im1, int ip1, int jm1, int jp1,
    double hy
)
{
    double aj2_term = 0.0; // 計算途中のArakawa Jacobianの値を保持

    // --- Arakawa Jacobian の各項の計算 ---

    // J++ (J_plus_plus) 項: 中央差分近似に基づいたヤコビアン
    // (u_i+1,j - u_i-1,j) * (v_i,j+1 - v_i,j-1)
    // - (u_i,j+1 - u_i,j-1) * (v_i+1,j - v_i-1,j)
    aj2_term  = (uuu[ip1][j00] - uuu[im1][j00]) * (vvv[i00][jp1] - vvv[i00][jm1]);
    aj2_term -= (uuu[i00][jp1] - uuu[i00][jm1]) * (vvv[ip1][j00] - vvv[im1][j00]);

    // J+x (J_plus_cross) 項:
    // u_{i+1,j} * (v_{i+1,j+1} - v_{i+1,j-1})
    // - u_{i-1,j} * (v_{i-1,j+1} - v_{i-1,j-1})
    // - u_{i,j+1} * (v_{i+1,j+1} - v_{i-1,j+1})
    // + u_{i,j-1} * (v_{i+1,j-1} - v_{i-1,j-1})
    aj2_term += uuu[ip1][j00] * (vvv[ip1][jp1] - vvv[ip1][jm1]);
    aj2_term -= uuu[im1][j00] * (vvv[im1][jp1] - vvv[im1][jm1]);
    aj2_term -= uuu[i00][jp1] * (vvv[ip1][jp1] - vvv[im1][jp1]);
    aj2_term += uuu[i00][jm1] * (vvv[ip1][jm1] - vvv[im1][jm1]);

    // Jx+ (J_cross_plus) 項:
    // u_{i+1,j+1} * (v_{i,j+1} - v_{i+1,j})
    // - u_{i-1,j-1} * (v_{i-1,j} - v_{i,j-1})
    // - u_{i-1,j+1} * (v_{i,j+1} - v_{i-1,j})
    // + u_{i+1,j-1} * (v_{i+1,j} - v_{i,j-1})
    aj2_term += uuu[ip1][jp1] * (vvv[i00][jp1] - vvv[ip1][j00]);
    aj2_term -= uuu[im1][jm1] * (vvv[im1][j00] - vvv[i00][jm1]);
    aj2_term -= uuu[im1][jp1] * (vvv[i00][jp1] - vvv[im1][j00]);
    aj2_term += uuu[ip1][jm1] * (vvv[ip1][j00] - vvv[i00][jm1]);

    // 最終的なスケールファクター (.25 * hy * hy / 3.) を適用して返す
    return aj2_term * (.25 * hy * hy / 3.);
}

// ブラケットの 2次 Arakawa スキーム (Arakawa 1966, JCP 1997 に再録を参照)
// [u,v] = (du/dx)(dv/dy)-(du/dy)(dv/dx)
void arakawa( AXY uuu, AXY vvv, AXY out )
{
  int i00,j00,ip1,jp1,im1,jm1,im2,ip2,jm2,jp2;
  double aj2, afac = .25*hy*hy/3.;
  int i0 = 0, i1 = nx1, j0 = 0, j1 = ny1;

#pragma omp parallel for private(j00,im1,im2,jm1,jm2,ip1,ip2,jp1,jp2,aj2)
  for (i00=0; i00<=nx1; ++i00)  
  {
    nextp(i00,i0,i1,im1,ip1,im2,ip2);

    for (j00=0; j00<=ny1; ++j00) 
	{
      nextp(j00,j0,j1,jm1,jp1,jm2,jp2);

     out[i00][j00] = calculate_arakawa_jacobian_at_point(
                          uuu, vvv,
                          i00, j00,
                          im1, ip1, jm1, jp1,
                          hy
                      );
    };
  }
} // arakawa 関数の終了
