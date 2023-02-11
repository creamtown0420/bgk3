#include <stdio.h>
#include <functional>
#include <cmath>
#include <errno.h>
#include <iomanip>
#include <omp.h>
#include <time.h>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include "boost/multi_array.hpp"
#include <fstream>
#include <vector>
#include "Eigen/Core"
#include <iostream>
#include <fstream>
// #include <Eigen/CXX11/Tensor>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
class MeshData // 汎用性の高いメッシュ生成機能を持つクラス
{
private:
    static const int step = 2;                  // 時間データが2層
    std::array<std::vector<double>, step> data; // 二次元配列dataの定義 data[step][]
    int nx;

public:
    std::vector<double> &operator[](int nt) //  operator[] という関数は, オブジェクト変数名[n] という書き方をされた時に, どういう値が得られるか，を定義するもの
    {
        return data[nt % step]; // dataはvectorの配列なので、data[nt]とすると、二次元配列でも配列の一番最初の行から数えていく．nt%stepは1か0
    }
    MeshData(int N) : nx(N) // コンストラクタ nxは引数Nをもち，nxをNに初期化
    {
        for (int nt = 0; nt < step; nt++)
            data[nt].resize(N); // ２次元配列data[nt][]うしろの[]をresizeで配列の要素数を無理やりNにしてる
    }
    int size() // サイズを聞かれたら教える
    {
        return nx;
    }
};
class Solver
{
    MeshData U, V; // MeshからU,Vをつくる　MeshData(int N) : nx(N)というコンストラクタがあるので，MeshData U(N),V(N)ということができる
    double c = 1.0, cfl;
    double dt, dx;
    std::vector<double> X; // Xは初期値用の定義か？ちがいそう
    std::ofstream ofs;

public:
    Solver(double xmax, int nx, const double dt, std::string init, std::string filename)
        : U(nx + 1), V(nx + 1), dt(dt) //  MeshData(int N)のint Nにnx+1をいれてUとVのメッシュを作る　dt(dt)はdtをdtに初期化
    {
        if (nx % 2 != 0)
            throw(std::runtime_error("nxが偶数じゃないぞタコ"));
        boost::filesystem::path folder = getenv("HOME");
        boost::filesystem::current_path(folder / "tes_test");
        dx = xmax / nx;     // 0,1,2,...,nx
        X.resize(U.size()); // XのサイズをUのサイズにする
        for (int i = 0; i < X.size(); i++)
            X[i] = -xmax / 2 + dx * i; // 格子Xをきんとうに分割して作る
        cfl = c * dt / dx;
        if (cfl > 1.0)
            throw(std::runtime_error("cflでかすぎぢゃボケ"));
        // 色んな初期値で遊ぶ
        if (init == "kaku")
        { // カクカクした初期値の例
            for (int i = 0; i < X.size(); i++)
            {
                U[0][i] = (std::abs(X[i]) <= 1.0) ? 1.0 : 0.0;
                V[0][i] = (std::abs(X[i]) <= 1.0) ? 0.0 : 0.0;
            }
        }
        else
        { // ぬるぬるした初期値の例
          // for (int i = 0; i < X.size(); i++)
          // {
          //     U[0][i] = std::exp(-5 * X[i] * X[i]);
          //     V[0][i] = 0.0;
          // }
        }
        ofs.open(filename);
    }
    void Step(int nt)
    {
        // 境界条件を設定
        U[nt + 1][0] = 0.;
        V[nt + 1][X.size()] = 0.;
        // 解きさらせ
        for (int i = 1; i <= X.size(); i++)
            U[nt + 1][i] = U[nt][i] - cfl * (V[nt][i] - V[nt][i - 1]);
        for (int i = 0; i < X.size(); i++)
            V[nt + 1][i] = V[nt][i] - cfl * (U[nt][i + 1] - U[nt][i]);
    }
    void Write(int nt)
    {
        for (int i = 0; i < X.size(); i++)
            ofs << X[i] << " " << U[nt][i] << std::endl;
        ofs << std::endl;
    }
};
int main()
{

    auto tens = JKLM();
    tens.Jmax = 16, tens.Kmax = 16, tens.Lmax = 16, tens.Mmax = 16;
    boost::multi_array<double, 4> Kjklm(boost::extents[tens.Jmax][tens.Kmax][tens.Lmax][tens.Mmax]);

    // ファイルをバイナリ形式で開く
    ifstream binaryFile("test.dat", ios::binary);

    // 読み取り可能な状態であるか確認
    if (binaryFile.good())
    {
        // Kjklmにファイルから読み取った値を代入
        binaryFile.read((char *)Kjklm.data(), sizeof(double) * tens.Jmax * tens.Kmax * tens.Lmax * tens.Mmax);

        // 結果を表示
        //     for (int i = 0; i < tens.Jmax; i++)
        //     {
        //         for (int j = 0; j < tens.Kmax; j++)
        //         {
        //             for (int k = 0; k < tens.Lmax; k++)
        //             {
        //                 for (int l = 0; l < tens.Mmax; l++)
        //                 {
        //                     cout << "Kjklm(" << i << "," << j << "," << k << "," << l << ") = " << Kjklm[i][j][k][l] << endl;
        //                 }
        //             }
        //         }
        //     }
    }
    else
    {
        // ファイルを開けなかった場合
        cout << "Error opening file" << endl;
    }
}

// ファイルを閉じる
binaryFile.close();

return 0;
}