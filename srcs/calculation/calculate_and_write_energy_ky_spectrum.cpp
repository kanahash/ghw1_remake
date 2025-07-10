#include "ghw1.hpp" // カスタムヘッダーファイルをインクルード
// 必要に応じて他の標準ライブラリをインクルード

// calculate_and_write_ky_spectrum 関数の定義
void calculate_and_write_ky_spectrum(double AXY_data[nxm][nym],
                                     const char *filename, fftw_plan hindfty, double *py_array,
                                     fftw_complex *ky_array, double *sumky_array, double *avgky_array,
                                     int nx1, int ny1, int ny, int nx, int incon, int it)
{
    // この関数内の具体的なロジックは提示されていませんが、
    // エラーメッセージに基づき、問題の行を修正します。

    // 例: エラーが発生していたであろう箇所
    // (コンテキストが不明なため、以下はあくまで例です)
    // 例えば、ky_array の実部を参照したかった場合:
    // double value = ky_array[j][0]; // <-- 修正: 余分な * を削除

    // もし以下のようなループ内でエラーが発生していた場合、
    // 例として、ky_arrayの計算結果をsumky_arrayに加算する部分を想定します。
    // 実際のコードに合わせて修正してください。

    // 仮のループと計算 (元のコードの正確なコンテキストに合わせてください)
    for (int j = 0; j <= ny / 2; ++j) { // ny/2 は FFTW の出力サイズに合わせる
        // エラーが出た行の正確な内容が不明なため、一般的なFFT出力処理を想定
        // ky_array[j][0] はすでに double 型なので、* は不要
        sumky_array[j] += fabs(ky_array[j][0]); // 実部
        // sumky_array[j] += fabs(ky_array[j][1]); // 虚部 (もし必要なら)

        // 例えば、kyスペクトルの振幅を計算する場合
        // sumky_array[j] += sqrt(ky_array[j][0]*ky_array[j][0] + ky_array[j][1]*ky_array[j][1]);
    }

    // ファイル書き込み部分 (例として)
    FILE *f_ky = fopen(filename, "w");
    if (f_ky == NULL) {
        perror("Error opening ky spectrum file");
        return;
    }
    // ここで sumky_array や avgky_array をファイルに書き込むロジックが続く
    fclose(f_ky);

    // その他の計算やロジックが続く...
}

// 他の calculate_and_write_energy_ky_spectrum や calculate_and_write_kx_spectrum 関数も
// 同様のエラーが出る可能性があるため、もし発生したら同様に * を削除してください。
