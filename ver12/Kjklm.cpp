#include <iostream>
#include <stdio.h>
#include <functional>
#include <cmath>
#include <errno.h>
#include <iomanip>
#include <omp.h>
#include <time.h>
#include <numbers>
// #include <boost/multiprecision/cpp_dec_float.hpp>
// #include <boost/math/constants/constants.hpp> // 16を任意の精度に変えれる
#include <fstream>
#include <vector>
#include <boost/math/special_functions/ellint_rf.hpp>
#include "3d_tanh_sinh.hpp"
#include "gauss.hpp"
using namespace std;
#define Jmax 4
#define Kmax 4
#define Lmax 4
#define Mmax 4
namespace std
{ // これをいれとくと関数を先に定義したことになるらしい
    float comp_ellint_1f(float k);
    double comp_ellint_1(double k);
    long double comp_ellint_1l(long double k);
}
namespace std
{ // c++のばーじょんがc++17以上じゃないとエラーになる　今はc++98
    float comp_ellint_2f(float k);
    double comp_ellint_2(double k);
    long double comp_ellint_2l(long double k);
}
namespace boost
{
    namespace math
    {

        template <class T1, class T2, class T3>
            calculated - result - type ellint_rf(T1 x, T2 y, T3 z)

                                      template <class T1, class T2, class T3, class Policy>
                                      calculated - result - type ellint_rf(T1 x, T2 y, T3 z, const Policy &)

    }
} // namespaces

double M_PI = 3.14159265358979;
double F1_(double zeta_bar, double theta_bar, double psi, double zeta_j, double theta_k)
{

    return (zeta_bar) * (zeta_bar) + (zeta_j) * (zeta_j)             //
           - 2 * (zeta_bar) * (zeta_j)*cos(theta_bar) * cos(theta_k) //
           - 2 * (zeta_bar) * (zeta_j)*sin(theta_bar) * sin(theta_k) * cos(psi);
}
double F2_(double zeta_bar, double theta_bar, double psi, double zeta_j, double theta_k)
{
    return (zeta_bar) * (zeta_bar) * (zeta_j) * (zeta_j)   //
           - (zeta_bar) * (zeta_bar) * (zeta_j) * (zeta_j) //
                 * (cos(theta_k) * cos(theta_bar) + sin(theta_k) * sin(theta_bar) * cos(psi));
}
double f_for_psi_(double zeta_bar, double theta_bar, double zeta_l, double theta_m, double alpha, double beta, double gamma, double delta)
{
    return ((zeta_bar) - (alpha)) * ((zeta_bar) - (beta)) * ((theta_bar) - (gamma)) * ((theta_bar) - (delta)) / (((zeta_l) - (alpha)) * ((zeta_l) - (beta)) * ((theta_m) - (gamma)) * ((theta_m) - (delta)));
}
double L1(double zeta_bar, double theta_bar, double psi, double zeta_j, double theta_k, double zeta_l, double theta_m, double alpha, double beta, double gamma, double delta)
{
    return 1 / (sqrt(2) * M_PI) * (zeta_bar) * (zeta_bar)*sin(theta_bar) * cos(psi)                                 //
           / sqrt(F1_(zeta_bar, theta_bar, psi, zeta_j, theta_k))                                                   //
                                                                                                                    //                                                                        //
           * exp(-(zeta_bar) * (zeta_bar)                                                                           //
                 + F2_(zeta_bar, theta_bar, psi, zeta_j, theta_k) / F1_(zeta_bar, theta_bar, psi, zeta_j, theta_k)) //
           * f_for_psi_(zeta_bar, theta_bar, zeta_l, theta_m, alpha, beta, gamma, delta);
}
double L2(double zeta_bar, double theta_bar, double psi, double zeta_j, double theta_k, double zeta_l, double theta_m, double alpha, double beta, double gamma, double delta)

{
    return 1 / (2 * sqrt(2) * M_PI) * (zeta_bar) * (zeta_bar)*sin(theta_bar) //
           * sqrt(F1_(zeta_bar, theta_bar, psi, zeta_j, theta_k)) * cos(psi) //
           * exp(-(zeta_bar) * (zeta_bar))                                   //
           * f_for_psi_(zeta_bar, theta_bar, zeta_l, theta_m, alpha, beta, gamma, delta);
}
double L2_ellip(double zeta_bar, double theta_bar, double zeta_j, double theta_k, double zeta_l, double theta_m, double alpha, double beta, double gamma, double delta)

{
    double y = 4 * (zeta_bar) / (zeta_j);
    return 1 / (2 * sqrt(2) * M_PI) * (zeta_bar) * (zeta_bar)*sin(theta_bar)             //
           * sqrt((zeta_bar) * (zeta_bar) + (zeta_j) * (zeta_j)                          //
                  - 2 * (zeta_bar) * (zeta_j)*cos(theta_bar) * cos(theta_k)              //
                  - 2 * (zeta_bar) * (zeta_j)*sin(theta_bar) * sin(theta_k))             //
           * exp(-(zeta_bar) * (zeta_bar))                                               //
           * f_for_psi_(zeta_bar, theta_bar, zeta_l, theta_m, alpha, beta, gamma, delta) //
           * (((y - 2) * comp_ellint_2(y) + 2 * (1 - y) * comp_ellint_1(y)) / (3 * y));
}

int main()
{
    clock_t start = clock(); // 実行時間

    int ier = 0;
    int key = 3;
    int info = 0;
    double xa;
    double xb;
    double ya;
    double yb;
    double za;
    double zb;

    //==========物理的な仮定の変数=========
    double Cap_zeta = 100; // 無限遠のzeta
    double Cap_theta = 100;
    // //=========格子点に関する変数=========
    // int Kmax = 16; // 16分割する
    // int Jmax = 16; // 16分割する
    // int Lmax = 16; // 16分割する
    // int Mmax = 16; // 16分割する

    //================変数=============
    double zeta[Jmax];
    double theta[Kmax];
    static double Kjklm[Jmax][Kmax][Lmax][Mmax]; // セグメント違反が出るのでstaticで

    for (int j = 0; j < Jmax; j++)
    {
        zeta[j] = Cap_zeta * j / Jmax;
    }
    for (int k = 0; k < Kmax; k++)
    {
        theta[k] = M_PI * k / Kmax;
    }

    //=================================================================
    //============forループ内で使う変数の定義============================
    int j, k, l, m;
    double eps = 1e-8;
    double s = 0;
    int count = 0;
    double psi_a = 0, psi_b = 2 * M_PI;
    int a_l, b_l, c_m, d_m;                                     // 積分区間のインデックス
    int a_l_p, a_l_m, b_l_p, b_l_m, c_m_p, c_m_m, d_m_p, d_m_m; // 偶偶用の積分区間のインデックス(p,mはそれぞれplus,minus,(pp,mm,pm,mp))
    double alpha, beta, gamma, delta;                           // psi用の変数
    double zeta_l, theta_m, zeta_j, theta_k;                    // zeta[l]とか，特異点
    double zeta_a_l, zeta_b_l, theta_d_m, theta_c_m;            // 積分範囲lが始まり，mがおわり
    double zeta_a, zeta_b, theta_d, theta_c;                    // 積分範囲lが始まり，mがおわり
    double zeta_a_l_m, zeta_b_l_m, theta_d_m_m, theta_c_m_m;    // 積分範囲lが始まり，mがおわり

    double val_kiki_1, val_kiki_2;
    double val_gugu_pp_1, val_gugu_pp_2, val_gugu_pm_1, val_gugu_pm_2, val_gugu_mp_1, val_gugu_mp_2, val_gugu_mm_1, val_gugu_mm_2;
    double val_kigu_0p_1, val_kigu_0p_2, val_kigu_0m_1, val_kigu_0m_2, val_guki_p0_1, val_guki_p0_2, val_guki_m0_1, val_guki_m0_2;
    double val1, val2, val3, val4; // 特異点を避ける積分範囲の値
    //================================================================
    // #pragma omp parallel for private(j, k, l, m, zeta_j, theta_k, zeta_l, theta_m, zeta_a_l, theta_c_m, a_l, b_l, c_m, d_m, alpha, beta, gamma, delta, s, ier) // epsをprivateにすると何故かバグる
    // #pragma omp parallel for private(k, l, m) // firstprivate(zeta_j, theta_k, zeta_l, theta_m, zeta_a_l, theta_c_m, alpha, beta, gamma, delta) // epsをprivateにすると何故かバグる
    // #pragma omp parallel for private(k, l, m, count) // firstprivate(zeta_j, theta_k, zeta_l, theta_m, zeta_a_l, theta_c_m, alpha, beta, gamma, delta) // epsをprivateにすると何故かバグる
    for (j = 0; j < Jmax; j++)
    {
        for (k = 0; k < Kmax; k++)
        {
            for (l = 0; l < Lmax; l++)
            {
                for (m = 0; m < Mmax; m++)
                {
                    printf("j= %d,k= %d, l= %d, m= %d\n", j, k, l, m);
                    // printf("singurality"); // ここからprintできない　わけわからん なぜかL1を変えると値も出る

                    zeta_j = zeta[j], theta_k = theta[k], zeta_l = zeta[l], theta_m = theta[m];
                    // printf("singurality");

                    if (j == 0 || k == 0 || k == Kmax)
                    {
                        Kjklm[j][k][l][m] = 0.0;
                    }
                    else
                    {
                        if (l % 2 == 1 & m % 2 == 1)
                        {

                            a_l = l - 1, b_l = l + 1, c_m = m - 1, d_m = m + 1;
                            // zeta_a_l = zeta[a_l];
                            // #pragma omp critical                                                                                                   // 積分範囲の格子の位置
                            alpha = zeta[l - 1], beta = zeta[l + 1], gamma = theta[m - 1], delta = theta[m + 1]; // psiの0にする点
                                                                                                                 //===========奇数奇数ようの関数
                            auto L1_kiki = [&](double zeta_bar, double theta_bar, double psi) -> double
                            {
                                return L1(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                            };
                            auto L2_kiki = [&](double zeta_bar, double theta_bar) -> double
                            {
                                return L2_ellip(zeta_bar, theta_bar, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                            };

                            if (0 <= a_l & a_l <= Lmax & 0 <= b_l & b_l <= Lmax & 0 <= c_m & c_m <= Mmax & 0 <= d_m & d_m <= Mmax)
                            {
                                // printf("singurality");

                                zeta_a_l = zeta[l - 1], zeta_b_l = zeta[l + 1], theta_c_m = theta[m - 1], theta_d_m = theta[m + 1];

                                if (zeta[a_l] <= zeta_j & zeta_j <= zeta[b_l] & theta[c_m] <= theta_k & theta_k <= theta[d_m])
                                {
                                    // 特異点があるとき
                                    printf("singurality\n");

                                    val1 = dde3d(L1_kiki, zeta_a_l, zeta_j, theta_c_m, theta_k, psi_a, psi_b, eps); //                 //左下

                                    // // 右下
                                    val2 = dde3d(L1_kiki, zeta_j, zeta_b_l, theta_c_m, theta_k, psi_a, psi_b, eps); //                 //左下

                                    val3 = dde3d(L1_kiki, zeta_a_l, zeta_j, theta_k, theta_d_m, psi_a, psi_b, eps); //                 //左下
                                    //  // 左上

                                    val4 = dde3d(L1_kiki, zeta_j, zeta_b_l, theta_k, theta_d_m, psi_a, psi_b, eps); //                 //左下
                                    //  // 右上
                                    val_kiki_1 = val1 + val2 + val3 + val4;
                                    // strfromf128(buf, sizeof(buf), "%.40g", val_kiki_1);
                                    // printf("val_kiki_2=%s\n", buf);
                                }

                                else
                                {
                                    printf("no singurality\n");
                                    val_kiki_1 = quad3(L1_kiki, zeta_a_l, zeta_b_l, theta_c_m, theta_d_m, psi_a, psi_b, eps); //                 //左下
                                }
                                printf("L2\n");

                                val_kiki_2 = quad2(L2_kiki, zeta_a_l, zeta_b_l, theta_c_m, theta_d_m, eps); //                 //左下

                                // printf("val_kiki_2=%lf\n", val_kiki_2);
                            }
                            else
                            {
                                val_kiki_1 = 0.0, val_kiki_2 = 0.0;
                            }
                            Kjklm[j][k][l][m] = val_kiki_1 + val_kiki_2;
                        }
                        //==================偶数偶数========================================
                        else if (l % 2 == 0 & m % 2 == 0)
                        {
                            a_l_p = l, b_l_p = l + 2, c_m_p = m, d_m_p = m + 2; //// こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です
                            a_l_m = l - 2, b_l_m = l, c_m_m = m - 2, d_m_m = m;
                            //==============pp==================ppの領域が範囲外に飛び出してなかったら                                                                   // //こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です
                            if (0 <= a_l_p & a_l_p <= Lmax & 0 <= b_l_p & b_l_p <= Lmax & 0 <= c_m_p & c_m_p <= Mmax & 0 <= d_m_p & d_m_p <= Mmax) // pp
                            {
                                alpha = zeta[l + 1], beta = zeta[l + 2], gamma = theta[m + 1], delta = theta[m + 2];
                                // 関数定義
                                auto L1_gugu_pp = [&](double zeta_bar, double theta_bar, double psi) -> double
                                {
                                    return L1(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                auto L2_gugu_pp = [&](double zeta_bar, double theta_bar) -> double
                                {
                                    return L2_ellip(zeta_bar, theta_bar, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };

                                //===================特異点===============================
                                if (zeta[a_l_p] <= zeta_j & zeta_j <= zeta[b_l_p] & theta[c_m_p] <= theta_k & theta_k <= theta[d_m_p])
                                { // 特異点があるとき;
                                    // 積分範囲
                                    zeta_a = zeta[a_l_p], zeta_b = zeta[b_l_p], theta_c = theta[c_m_p], theta_d = theta[d_m_p];

                                    printf("singurality\n");
                                    val1 = dde3d(L1_gugu_pp, zeta_a, zeta_j, theta_c, theta_k, psi_a, psi_b, eps); // 左下;

                                    val2 = dde3d(L1_gugu_pp, zeta_j, zeta_b, theta_c, theta_k, psi_a, psi_b, eps); //; 右下

                                    val3 = dde3d(L1_gugu_pp, zeta_a, zeta_j, theta_k, theta_d, psi_a, psi_b, eps); //; 左上

                                    val4 = dde3d(L1_gugu_pp, zeta_j, zeta_b, theta_k, theta_d, psi_a, psi_b, eps); //; 右上

                                    val_gugu_pp_1 = val1 + val2 + val3 + val4;
                                    // singurality();
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    val_gugu_pp_1 = quad3(L1_gugu_pp, zeta_a, zeta_b, theta_c, theta_d, psi_a, psi_b, eps);
                                }
                                printf("L2\n");
                                val_gugu_pp_2 = quad2(L2_gugu_pp, zeta_a, zeta_b, theta_c, theta_d, eps);
                            }
                            //==============pp==================ppの領域が範囲外に飛び出してたら
                            else
                            {
                                printf("hanigai_pp\n");
                                val_gugu_pp_1 = 0.0;
                                val_gugu_pp_2 = 0.0;
                            }
                            Kjklm[j][k][l][m] = val_gugu_pp_1 + val_gugu_pp_2;
                            //=================pm===================pmの領域が範囲内だったら
                            if (0 <= a_l_p & a_l_p <= Lmax & 0 <= b_l_p & b_l_p <= Lmax & 0 <= c_m_m & c_m_m <= Mmax & 0 <= d_m_m & d_m_m <= Mmax)
                            {
                                alpha = zeta[l + 1], beta = zeta[l + 2], gamma = theta[m - 1], delta = theta[m - 2];

                                auto L1_gugu_pm = [&](double zeta_bar, double theta_bar, double psi) -> double
                                {
                                    return L1(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                auto L2_gugu_pm = [&](double zeta_bar, double theta_bar) -> double
                                {
                                    return L2_ellip(zeta_bar, theta_bar, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                //=================特異点=======================
                                if (zeta[a_l_p] <= zeta_j & zeta_j <= zeta[b_l_p] & theta[c_m_m] <= theta_k & theta_k <= theta[d_m_m])
                                {
                                    zeta_a = zeta[a_l_p], zeta_b = zeta[b_l_p], theta_c = theta[c_m_m], theta_d = theta[d_m_m];

                                    printf("singurality\n");
                                    dde3d(L1_gugu_pm, zeta_a, zeta_j, theta_c, theta_k, psi_a, psi_b, eps); // 左下;

                                    dde3d(L1_gugu_pm, zeta_j, zeta_b, theta_c, theta_k, psi_a, psi_b, eps); //; 右下

                                    dde3d(L1_gugu_pm, zeta_a, zeta_j, theta_k, theta_d, psi_a, psi_b, eps); //; 左上

                                    dde3d(L1_gugu_pm, zeta_j, zeta_b, theta_k, theta_d, psi_a, psi_b, eps); //; 右上

                                    val_gugu_pm_1 = val1 + val2 + val3 + val4;
                                    // singurality();
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    quad3(L1_gugu_pm, zeta_a, zeta_b, theta_c, theta_d, psi_a, psi_b, eps);
                                    val_gugu_pm_1 = s;
                                }
                                printf("L2\n");
                                quad2(L2_gugu_pm, zeta_a, zeta_b, theta_c, theta_d, eps);
                                val_gugu_pm_2 = s;

                            } // ppの領域が範囲外に飛び出してたら
                            else
                            {
                                printf("hanigai_pm\n");
                                val_gugu_pm_1 = 0.0;
                                val_gugu_pm_2 = 0.0;
                            }
                            //!-- -- -- -- -- -- -- -- -- -- -- -mp-- -- -- -- -- -- -- -- -- -- -l_m, m_pになる要チェック
                            //
                            if (0 <= a_l_m & a_l_m <= Lmax & 0 <= b_l_m & b_l_m <= Lmax & 0 <= c_m_p & c_m_p <= Mmax & 0 <= d_m_p & d_m_p <= Mmax) // mp
                            {
                                alpha = zeta[l - 1], beta = zeta[l - 2], gamma = theta[m + 1], delta = theta[m + 2];
                                auto L1_gugu_mp = [&](double zeta_bar, double theta_bar, double psi) -> double
                                {
                                    return L1(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                auto L2_gugu_mp = [&](double zeta_bar, double theta_bar) -> double
                                {
                                    return L2_ellip(zeta_bar, theta_bar, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                if (zeta[a_l_m] <= zeta_j & zeta_j <= zeta[b_l_m] & theta[c_m_p] <= theta_k & theta_k <= theta[d_m_p])
                                {
                                    zeta_a = zeta[a_l_m], zeta_b = zeta[b_l_m], theta_c = theta[c_m_p], theta_d = theta[d_m_p];

                                    printf("singurality\n");
                                    val1 = dde3d(L1_gugu_mp, zeta_a, zeta_j, theta_c, theta_k, psi_a, psi_b, eps); // 左下;

                                    val2 = dde3d(L1_gugu_mp, zeta_j, zeta_b, theta_c, theta_k, psi_a, psi_b, eps); //; 右下

                                    val3 = dde3d(L1_gugu_mp, zeta_a, zeta_j, theta_k, theta_d, psi_a, psi_b, eps); //; 左上

                                    val4 = dde3d(L1_gugu_mp, zeta_j, zeta_b, theta_k, theta_d, psi_a, psi_b, eps); //; 右上

                                    val_gugu_mp_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    val_gugu_mp_1 = quad3(L1_gugu_mp, zeta_a, zeta_b, theta_c, theta_d, psi_a, psi_b, eps);
                                }
                                printf("L2\n");
                                val_gugu_mp_2 = quad2(L2_gugu_mp, zeta_a, zeta_b, theta_c, theta_d, eps);

                            } // ppの領域が範囲外に飛び出してたら
                            else
                            {
                                printf("hanigai_mp\n");
                                val_gugu_mp_1 = 0.0;
                                val_gugu_mp_2 = 0.0;
                            }
                            //!-- -- -- -- -- -- -- -- -- -- -- -mm-- -- -- -- -- -- -- -- -- -- -l_m, m_mになる要チェック
                            //
                            if (0 <= a_l_m & a_l_m <= Lmax & 0 <= b_l_m & b_l_m <= Lmax & 0 <= c_m_m & c_m_m <= Mmax & 0 <= d_m_m & d_m_m <= Mmax) // mp
                            {
                                alpha = zeta[l - 1], beta = zeta[l - 2], gamma = theta[m - 1], delta = theta[m - 2];
                                auto L1_gugu_mm = [&](double zeta_bar, double theta_bar, double psi) -> double
                                {
                                    return L1(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                auto L2_gugu_mm = [&](double zeta_bar, double theta_bar) -> double
                                {
                                    return L2_ellip(zeta_bar, theta_bar, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                if (zeta[a_l_m] <= zeta_j & zeta_j <= zeta[b_l_m] & theta[c_m_m] <= theta_k & theta_k <= theta[d_m_m])
                                {
                                    zeta_a = zeta[a_l_m], zeta_b = zeta[b_l_m], theta_c = theta[c_m_m], theta_d = theta[d_m_m];

                                    printf("singurality\n");
                                    val1 = dde3d(L1_gugu_mm, zeta_a, zeta_j, theta_c, theta_k, psi_a, psi_b, eps); // 左下;

                                    val2 = dde3d(L1_gugu_mm, zeta_j, zeta_b, theta_c, theta_k, psi_a, psi_b, eps); //; 右下

                                    val3 = dde3d(L1_gugu_mm, zeta_a, zeta_j, theta_k, theta_d, psi_a, psi_b, eps); //; 左上

                                    val4 = dde3d(L1_gugu_mm, zeta_j, zeta_b, theta_k, theta_d, psi_a, psi_b, eps); //; 右上

                                    val_gugu_mm_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    val_gugu_mm_1 = quad3(L1_gugu_mm, zeta_a, zeta_b, theta_c, theta_d, psi_a, psi_b, eps);
                                }
                                printf("L2\n");
                                val_gugu_mm_2 = quad2(L2_gugu_mm, zeta_a, zeta_b, theta_c, theta_d, eps);

                            } // ppの領域が範囲外に飛び出してたら
                            else
                            {
                                printf("hanigai_mm\n");
                                val_gugu_mm_1 = 0.0;
                                val_gugu_mm_2 = 0.0;
                            }
                            Kjklm[j][k][l][m] = val_gugu_pp_1 + val_gugu_pp_2 + val_gugu_pm_1 + val_gugu_pm_2 //
                                                + val_gugu_mp_1 + val_gugu_mp_2 + val_gugu_mm_1 + val_gugu_mm_2;
                        }

                        //==================奇数偶数====================
                        else if (l % 2 == 1 & m % 2 == 0)
                        {
                            a_l = l - 1, b_l = l + 1, c_m_p = m, d_m_p = m + 2, c_m_m = m - 2, d_m_m = m;
                            //=============0p=========================
                            if (0 <= a_l & a_l <= Lmax & 0 <= b_l & b_l <= Lmax & 0 <= c_m_p & c_m_p <= Mmax & 0 <= d_m_p & d_m_p <= Mmax)
                            {
                                alpha = zeta[l - 1], beta = zeta[l + 1], gamma = theta[m + 1], delta = theta[m + 2];

                                auto L1_kigu_0p = [&](double zeta_bar, double theta_bar, double psi) -> double
                                {
                                    return L1(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                auto L2_kigu_0p = [&](double zeta_bar, double theta_bar) -> double
                                {
                                    return L2_ellip(zeta_bar, theta_bar, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                if (zeta[a_l] <= zeta_j & zeta_j <= zeta[b_l] & theta[c_m_p] <= theta_k & theta_k <= theta[d_m_p])
                                {
                                    printf("singrality");
                                    val1 = dde3d(L1_kigu_0p, zeta_a, zeta_j, theta_c, theta_k, psi_a, psi_b, eps); // 左下;

                                    val2 = dde3d(L1_kigu_0p, zeta_j, zeta_b, theta_c, theta_k, psi_a, psi_b, eps); //; 右下

                                    val3 = dde3d(L1_kigu_0p, zeta_a, zeta_j, theta_k, theta_d, psi_a, psi_b, eps); //; 左上

                                    val4 = dde3d(L1_kigu_0p, zeta_j, zeta_b, theta_k, theta_d, psi_a, psi_b, eps); //; 右上

                                    val_kigu_0p_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    val_kigu_0p_1 = quad3(L1_kigu_0p, zeta_a, zeta_b, theta_c, theta_d, psi_a, psi_b, eps);
                                }
                                printf("L2\n");
                                val_kigu_0p_2 = quad2(L2_kigu_0p, zeta_a, zeta_b, theta_c, theta_d, eps);
                            }
                            else
                            {
                                printf("hanigai_0p\n");
                                val_kigu_0p_1 = 0.0;
                                val_kigu_0p_2 = 0.0;
                            }

                            //================0m=======================
                            if (0 <= a_l & a_l <= Lmax & 0 <= b_l & b_l <= Lmax & 0 <= c_m_m & c_m_m <= Mmax & 0 <= d_m_m & d_m_m <= Mmax)
                            {
                                alpha = zeta[l - 1], beta = zeta[l + 1], gamma = theta[m - 1], delta = theta[m - 2];

                                auto L1_kigu_0m = [&](double zeta_bar, double theta_bar, double psi) -> double
                                {
                                    return L1(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                auto L2_kigu_0m = [&](double zeta_bar, double theta_bar) -> double
                                {
                                    return L2_ellip(zeta_bar, theta_bar, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                if (zeta[a_l] <= zeta_j & zeta_j <= zeta[b_l] & theta[c_m_p] <= theta_k & theta_k <= theta[d_m_p])
                                {
                                    printf("singrality");
                                    val1 = dde3d(L1_kigu_0m, zeta_a, zeta_j, theta_c, theta_k, psi_a, psi_b, eps); // 左下;

                                    val2 = dde3d(L1_kigu_0m, zeta_j, zeta_b, theta_c, theta_k, psi_a, psi_b, eps); //; 右下

                                    val3 = dde3d(L1_kigu_0m, zeta_a, zeta_j, theta_k, theta_d, psi_a, psi_b, eps); //; 左上

                                    val4 = dde3d(L1_kigu_0m, zeta_j, zeta_b, theta_k, theta_d, psi_a, psi_b, eps); //; 右上

                                    val_kigu_0m_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    val_kigu_0m_1 = quad3(L1_kigu_0m, zeta_a, zeta_b, theta_c, theta_d, psi_a, psi_b, eps);
                                }
                                printf("L2\n");
                                val_kigu_0m_2 = quad2(L2_kigu_0m, zeta_a, zeta_b, theta_c, theta_d, eps);
                            }
                            else
                            {
                                printf("hanigai_0m\n");
                                val_kigu_0m_1 = 0.0;
                                val_kigu_0m_2 = 0.0;
                            }
                            Kjklm[j][k][l][m] = val_kigu_0p_1 + val_kigu_0p_2 + val_kigu_0m_1 + val_kigu_0m_2;
                        }

                        //===============偶数奇数===========
                        else
                        {
                            a_l_p = l, b_l_p = l + 2, a_l_m = l - 2, b_l_p = l, c_m = m - 1, d_m = m + 1;
                            //================p0===============
                            if (0 <= a_l_p & a_l_p <= Lmax & 0 <= b_l_p & b_l_p <= Lmax & 0 <= c_m & c_m <= Mmax & 0 <= d_m & d_m <= Mmax)
                            {
                                alpha = zeta[l + 1], beta = zeta[l + 2], gamma = theta[m + 1], delta = theta[m - 1];

                                auto L1_guki_p0 = [&](double zeta_bar, double theta_bar, double psi) -> double
                                {
                                    return L1(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                auto L2_guki_p0 = [&](double zeta_bar, double theta_bar) -> double
                                {
                                    return L2_ellip(zeta_bar, theta_bar, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                //======特異点があるとき=========
                                if (zeta[a_l_p] <= zeta_j & zeta_j <= zeta[b_l_p] & theta[c_m] <= theta_k & theta_k <= theta[d_m])
                                {
                                    printf("singrality");
                                    val1 = dde3d(L1_guki_p0, zeta_a, zeta_j, theta_c, theta_k, psi_a, psi_b, eps); // 左下;

                                    val2 = dde3d(L1_guki_p0, zeta_j, zeta_b, theta_c, theta_k, psi_a, psi_b, eps); //; 右下

                                    val3 = dde3d(L1_guki_p0, zeta_a, zeta_j, theta_k, theta_d, psi_a, psi_b, eps); //; 左上

                                    val4 = dde3d(L1_guki_p0, zeta_j, zeta_b, theta_k, theta_d, psi_a, psi_b, eps); //; 右上

                                    val_guki_p0_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    val_guki_p0_1 = quad3(L1_guki_p0, zeta_a, zeta_b, theta_c, theta_d, psi_a, psi_b, eps);
                                }
                                printf("L2\n");
                                val_guki_p0_2 = quad2(L2_guki_p0, zeta_a, zeta_b, theta_c, theta_d, eps);
                            }
                            else
                            {
                                printf("hanigai_p0\n");
                                val_guki_p0_1 = 0.0;
                                val_guki_p0_2 = 0.0;
                            }
                            //================m0================
                            if (0 <= a_l_m & a_l_m <= Lmax & 0 <= b_l_m & b_l_m <= Lmax & 0 <= c_m & c_m <= Mmax & 0 <= d_m & d_m <= Mmax)
                            {
                                alpha = zeta[l - 1], beta = zeta[l - 2], gamma = theta[m + 1], delta = theta[m - 1];

                                auto L1_guki_m0 = [&](double zeta_bar, double theta_bar, double psi) -> double
                                {
                                    return L1(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                auto L2_guki_m0 = [&](double zeta_bar, double theta_bar) -> double
                                {
                                    return L2_ellip(zeta_bar, theta_bar, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta);
                                };
                                //===========特異点があるとき===========
                                if (zeta[a_l_m] <= zeta_j & zeta_j <= zeta[b_l_m] & theta[c_m] <= theta_k & theta_k <= theta[d_m])
                                {
                                    printf("singrality");
                                    val1 = dde3d(L1_guki_m0, zeta_a, zeta_j, theta_c, theta_k, psi_a, psi_b, eps); // 左下;

                                    val2 = dde3d(L1_guki_m0, zeta_j, zeta_b, theta_c, theta_k, psi_a, psi_b, eps); //; 右下

                                    val3 = dde3d(L1_guki_m0, zeta_a, zeta_j, theta_k, theta_d, psi_a, psi_b, eps); //; 左上

                                    val4 = dde3d(L1_guki_m0, zeta_j, zeta_b, theta_k, theta_d, psi_a, psi_b, eps); //; 右上

                                    val_guki_m0_1 = val1 + val2 + val3 + val4;
                                }
                                else
                                {
                                    printf("no singurality\n");
                                    val_guki_m0_1 = quad3(L1_guki_m0, zeta_a, zeta_b, theta_c, theta_d, psi_a, psi_b, eps);
                                }
                                printf("L2\n");
                                val_guki_m0_2 = (L2_guki_m0, zeta_a, zeta_b, theta_c, theta_d, eps);
                            }
                            else
                            {
                                printf("hanigai_m0\n");
                                val_guki_m0_1 = 0.0;
                                val_guki_m0_2 = 0.0;
                            }
                            Kjklm[j][k][l][m] = val_guki_p0_1 + val_guki_p0_2 + val_guki_m0_1 + val_guki_m0_2;
                        }
                    }
                    printf("s=%lf\n", Kjklm[j][k][l][m]);
                }
            }
        }
    }
    FILE *file;
    file = fopen("test.dat", "wb");
    fwrite(Kjklm, sizeof(Kjklm), 1, file);
    fclose(file);
    return 0;
}