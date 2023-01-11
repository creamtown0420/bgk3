#define __STDC_WANT_IEC_60559_TYPES_EXT__
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
extern void dqag_sk3d_(_Float128 (*)(_Float128 *, _Float128 *, _Float128 *), _Float128 *xa, _Float128 *xb, _Float128 *ya, _Float128 *yb, _Float128 *za, _Float128 *zb, _Float128 *eps, _Float128 *s, int *ier, int *key);
extern void dde3d_(_Float128 (*)(_Float128 *, _Float128 *, _Float128 *), _Float128 *xa, _Float128 *xb, _Float128 *ya, _Float128 *yb, _Float128 *za, _Float128 *zb, _Float128 *eps, _Float128 *s, int *info);

int main()
{
    int ier = 0;
    int key = 3;
    int info = 0;
    _Float128 xa;
    _Float128 xb;
    _Float128 ya;
    _Float128 yb;
    _Float128 za;
    _Float128 zb;
    _Float128 eps = 1e-20;
    _Float128 s = 0;
    char buf[1024];

    //==========物理的な仮定の変数=========
    _Float128 Cap_zeta2 = 100; // 無限遠のzeta
    _Float128 Cap_zetarho = 100;
    //=========格子点に関する変数=========
    int Kmax = 16; // 16分割する
    int Jmax = 16; // 16分割する
    int Lmax = 16; // 16分割する
    int Mmax = 16; // 16分割する

    //================変数=============
    _Float128 zeta2[Jmax];
    _Float128 zetarho[Kmax];
    static _Float128 Kjklm[16][16][16][16]; // セグメント違反が出るのでstaticで

    // _Float128 ****Kjklm;
    // Kjklm = (_Float128 ****)malloc(sizeof(int ***) * Jmax);
    // for (int e1 = 0; e1 < Jmax; e1 = e1 + 1)
    // {
    //     Kjklm[e1] = (_Float128 ***)malloc(sizeof(int **) * Kmax); // int*型のスペースを要素数（[][☆][]）の分だけ確保する。
    //     for (int e2 = 0; e2 < Kmax; e2 = e2 + 1)
    //     {
    //         Kjklm[e1][e2] = (_Float128 **)malloc(sizeof(int *) * Lmax); // int型のスペースを要素数（[][][☆]）の分だけ確保する。
    //         for (int e3 = 0; e3 < Lmax; e3 = e3 + 1)
    //         {
    //             Kjklm[e1][e2][e3] = (_Float128 *)malloc(sizeof(int) * Mmax);
    //         }
    //     }
    // }

    //_Float128 Kjklm[Jmax][Kmax][Lmax][Mmax];
    //_Float128 Kjklm[][][][];
    //================

    for (int j = 0; j < Jmax; j++)
    {
        zeta2[j] = Cap_zeta2 * j / Jmax;
    }
    for (int k = 0; k < Kmax; k++)
    {
        zetarho[k] = Cap_zetarho * k / Kmax;
    }
    //===========関数=============
    _Float128 F1_(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi, _Float128 * zeta2_j, _Float128 * zetarho_k)
    {
        // clang-format off
        return (*zeta2_bar) * (*zeta2_bar) + (*zeta2_j) * (*zeta2_j) + (*zetarho_bar) * (*zetarho_bar) + (*zetarho_k) * (*zetarho_k) \
        - 2 * ((*zeta2_bar) * (*zeta2_bar) + (*zetarho_bar) * (*zetarho_k)*cos(*psi));
        // clang-format on
    }
    _Float128 F2_(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi, _Float128 * zeta2_j, _Float128 * zetarho_k)
    {
        // clang-format off
        return ((*zeta2_bar) * (*zeta2_bar)+ (*zetarho_bar) * (*zetarho_bar)) *( (*zeta2_j) * (*zeta2_j)  + (*zetarho_k) * (*zetarho_k)) \
        -  pow(((*zeta2_bar) * (*zeta2_bar) + (*zetarho_bar) * (*zetarho_k)*cos(*psi)),2.0);
        // clang-format on
    }
    _Float128 f_for_psi_(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * zeta2_l, _Float128 * zetarho_m, _Float128 * alpha, _Float128 * beta, _Float128 * gamma, _Float128 * delta)
    {
        // clang-format off
        return ((*zeta2_bar) - (*alpha)) * ((*zeta2_bar) - (*beta)) * ((*zetarho_bar) - (*gamma)) * ((*zetarho_bar) - (*delta)) \
               / ((*zeta2_l) - (*alpha)) * ((*zeta2_l) - (*beta)) * ((*zetarho_m) - (*gamma)) * ((*zetarho_m) - (*delta));
        // clang-format on
    }
    _Float128 L1(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi, _Float128 * zeta2_j, _Float128 * zetarho_k, _Float128 * zeta2_l, _Float128 * zetarho_m, _Float128 * alpha, _Float128 * beta, _Float128 * gamma, _Float128 * delta)
    // clang-format off
    {
        return 1 / sqrt(2) / M_PI * (*zetarho_bar) / sqrt(F1_(zeta2_bar, zetarho_bar, psi, zeta2_j, zetarho_k)) \
        * exp(-(*zeta2_bar)*(*zeta2_bar)-(*zetarho_bar)*(*zetarho_bar)+F1_(zeta2_bar, zetarho_bar, psi, zeta2_j, zetarho_k)/F2_(zeta2_bar, zetarho_bar, psi, zeta2_j, zetarho_k))\
        * f_for_psi_(zeta2_bar,zetarho_bar,zeta2_l,zetarho_m,alpha,beta,gamma,delta);
        // clang-format on
    }
    _Float128 L2(_Float128 * zeta2_bar, _Float128 * zetarho_bar, _Float128 * psi, _Float128 * zeta2_j, _Float128 * zetarho_k, _Float128 * zeta2_l, _Float128 * zetarho_m, _Float128 * alpha, _Float128 * beta, _Float128 * gamma, _Float128 * delta)
    // clang-format off
    {
        return 1 / 2/sqrt(2) / M_PI * (*zetarho_bar) * sqrt(F1_(zeta2_bar, zetarho_bar, psi, zeta2_j, zetarho_k)) \
        * exp(-(*zeta2_bar)*(*zeta2_bar)-(*zetarho_bar)*(*zetarho_bar))\
        * f_for_psi_(zeta2_bar,zetarho_bar,zeta2_l,zetarho_m,alpha,beta,gamma,delta);
        // clang-format on
    }
    // fsim_(&i, &r);
    // fs_(&func_, &x, &y);
    int j, k, l, m;
    // xa = 0;
    // xb = 1;
    // ya = 0;
    // yb = 1;
    // za = 0;
    // zb = 1;
    // s = 0;
    // _Float128 s1;
    // _Float128 s2;
    // dqag_sk3d_(&func_, &xa, &xb, &ya, &yb, &za, &zb, &eps, &s, &ier, &key);
    // s = Kjklm[0][1][0][0];
#pragma omp parallel for private(j, k, l, m, xa, xb, ya, yb, za, zb, s, ier) // epsをprivateにすると何故かバグる
    for (j = 0; j < Jmax; j++)
    {
        for (k = 0; k < Kmax; k++)
        {
            for (l = 0; l < Lmax; l++)
            {
                for (m = 0; m < Mmax; m++)
                {
                    printf("j= %d,k= %d, l= %d, m= %d\n", j, k, l, m);

                    // dqag_sk3d_(&func_, &xa, &xb, &ya, &yb, &za, &zb, &eps, &s, &ier, &key);
                    //   s1 = s;

                    // dde3d_(&func_, &xa, &xb, &ya, &yb, &za, &zb, &eps, &s, &info);
                    // printf("j= %d,k= %d, l= %d, m= %d", j, k, l, m);

                    Kjklm[j][k][l][m] = s;
                    // printf("j= %d,k= %d, l= %d, m= %d", j, k, l, m);

                    // strfromf128(buf, sizeof(buf), "%.40g", s);
                    strfromf128(buf, sizeof(buf), "%.40g", Kjklm[j][k][l][m]);
                    printf("s=%s\n", buf);
                }
            }
        }
    }
    for (int j = 0; j < 3; j++)
    {
        for (int k = 0; k < 3; k++)
        {
            // strfromf128(buf, sizeof(buf), "%.40g", Kjklm[j][k][l][m]);
            // printf("i=%d j=%d K=%s\n", j, k, buf);
        }
    }

    return 0;
}