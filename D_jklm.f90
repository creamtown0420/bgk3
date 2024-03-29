module DE_integration
   implicit none
   integer, private, parameter::kmin = 3
   integer, private, parameter::kmax = 14
   ! max order of DE integration
   real(16), private, parameter::hr = 6.0_16    ! [-1,1]
   real(16), private, parameter::c0 = 0.01_16  ! safe value for eps
   !real(16), parameter::pi = acos(-1.0_16)
   real(16), private, parameter::pi = 3.141592653589793238462643383279502884_16
   real(16), private, parameter::pi2 = pi/2.0_16
   interface de3d
      module procedure &
         dde3d
   end interface de3d
contains
   subroutine dde3d(f, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, xa, xb, ya, yb, za, zb, eps, s, info)
      ! Integrate for 3d
      implicit none
      real(16), intent(in)::xa, xb, ya, yb, za, zb, eps, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta
      real(16), intent(out)::s
      integer, intent(out)::info
      interface
         function f(x, y, z, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)
            implicit none
            real(16), intent(in)::x, y, z, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta
            real(16)::f
         end function f
      end interface
      integer::j, nc, l
      real(16)::h, s0, xt, wt, t, ft, as, shk, mba, pba, err, seps
      info = 0
      seps = c0*sqrt(eps)
      h = hr*0.5_16

      mba = (xb - xa)*0.5_16
      pba = (xb + xa)*0.5_16
      call dde3dsyz(f, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, ya, yb, za, zb, pba, eps, ft, info)
      s0 = ft*h*pi2*mba
      nc = 1
      do l = 2, kmax
         s = 0.0_16
         nc = 2*nc
         h = h*0.5_16
         do j = 1, nc
            t = dble(2*j - nc - 1)*h
            shk = pi2*sinh(t)

            xt = tanh(shk)
            wt = pi2*cosh(t)/(cosh(shk)*cosh(shk))
          call dde3dsyz(f, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, ya, yb, za, zb, mba*xt + pba, eps, ft, info)
            s = s + ft*wt
         end do
         s = s0*0.5_16 + s*h*mba

         as = abs(s); err = abs(s - s0)
         if (as .ge. 1.0_16) err = err/as

         if (err .le. seps .and. l .ge. kmin) exit
         s0 = s
      end do

      if (l .eq. kmax + 1) info = 1

      return
   end subroutine dde3d

   subroutine dde3dsyz(f, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, ya, yb, za, zb, xc, eps, s, info)
      ! Integrate for 3d about y,z
      implicit none
      real(16), intent(in)::ya, yb, za, zb, xc, eps, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta
      real(16), intent(out)::s
      integer, intent(inout)::info
      interface
         function f(x, y, z, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)
            implicit none
            real(16), intent(in)::x, y, z, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta
            real(16)::f
         end function f
      end interface
      integer::j, nc, l
      real(16)::h, s0, xt, wt, t, ft, as, shk, mba, pba, err, seps
      seps = c0*sqrt(eps)
      h = hr*0.5_16

      ft = 0.0_16; s = 0.0_16
      mba = (yb - ya)*0.5_16
      pba = (yb + ya)*0.5_16
      call dde3dsz(f, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, za, zb, xc, pba, eps, ft, info)
      s0 = ft*h*pi2*mba

      nc = 1
      do l = 2, kmax
         s = 0.0_16
         nc = 2*nc
         h = h*0.5_16
         do j = 1, nc
            t = dble(2*j - nc - 1)*h
            shk = pi2*sinh(t)

            xt = tanh(shk)
            wt = pi2*cosh(t)/(cosh(shk)*cosh(shk))
            call dde3dsz(f, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, za, zb, xc, mba*xt + pba, eps, ft, info)
            s = s + ft*wt
            !print *, "y no sekibun l=", l, "j=", j
         end do
         s = s0*0.5_16 + s*h*mba
         as = abs(s); err = abs(s - s0)
         if (as .ge. 1.0_16) err = err/as
         if (err .le. seps .and. l .ge. kmin) exit

         s0 = s
      end do

      return
   end subroutine dde3dsyz

   subroutine dde3dsz(f, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, za, zb, xc, yc, eps, s, info)
      ! Integrate for 3d about z
      implicit none
      real(16), intent(in)::za, zb, xc, yc, eps, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta
      real(16), intent(out)::s
      integer, intent(inout)::info
      interface
         function f(x, y, z, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)
            implicit none
            real(16), intent(in)::x, y, z, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta
            real(16)::f
         end function f
      end interface
      integer::j, nc, l
      real(16)::h, s0, xt, wt, t, as, shk, mba, pba, err, seps
      seps = c0*sqrt(eps)
      h = hr*0.5_16

      mba = (zb - za)*0.5_16
      pba = (zb + za)*0.5_16
      s0 = f(xc, yc, pba, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)*h*pi2*mba

      s = 0.0_16
      nc = 1
      do l = 2, kmax
         s = 0.0_16
         nc = 2*nc
         h = h*0.5_16
         do j = 1, nc
            t = dble(2*j - nc - 1)*h
            shk = pi2*sinh(t)

            xt = tanh(shk)
            wt = pi2*cosh(t)/(cosh(shk)*cosh(shk))
            s = s + f(xc, yc, mba*xt + pba, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)*wt
            !print *, "z no sekibun l=", l, "j=", j
         end do

         s = s0*0.5_16 + s*h*mba
         as = abs(s); err = abs(s - s0)
         if (as .ge. 1.0_16) err = err/as
         if (err .le. seps .and. l .ge. kmin) exit
         s0 = s
      end do

      if (l .eq. kmax + 1) info = 1

      return
   end subroutine dde3dsz

end module DE_integration
program main
   use DE_integration
   implicit none
   real(16), parameter::pi2 = 1.5707963267948966_16
   !real(16), parameter::pi = acos(-1.0q0)
   real(16), parameter::pi = 3.141592653589793238462643383279502884_16
   !real(16), parameter::pi2 = pi/2.0q0

   real(16):: eps = 1d-24, s = 0q0
   real(16):: alpha, beta, gamma, delta
   real(16):: zeta_l, theta_m, zeta_j, theta_k
   real(16):: val_kiki_1, val_kiki_2
   real(16)::val_gugu_pp_1, val_gugu_pp_2, val_gugu_pm_1, val_gugu_pm_2, val_gugu_mp_1, val_gugu_mp_2, val_gugu_mm_1, val_gugu_mm_2
   real(16)::val_kigu_0p_1, val_kigu_0p_2, val_kigu_0m_1, val_kigu_0m_2, val_guki_p0_1, val_guki_p0_2, val_guki_m0_1, val_guki_m0_2
   real(16):: val1, val2, val3, val4!特異点を避ける積分範囲の値
   !eps = 1d-24; s = 0.0_16
   integer::info
   real(16) :: Cap_theta = 100.0_16
   real(16) :: Cap_zeta = 100.0_16
   !  integer :: display = 10                !if mod(n,display) = 0, then show the situation on the console
   !  integer :: proc                        !For cputime measurement

   !------------------------------------------------------------
   !@@@ Physical and numerical parameters ----------------------
   integer, parameter :: Kmax = 100  !Nz (maximum value of theta)
   integer, parameter :: Jmax = 100 !Nz (maximum value of zeta)
   integer, parameter :: Lmax = 100  !Nz (maximum value of zeta)
   integer, parameter :: Mmax = 100  !Nz (maximum value of zeta)
   !------------------------------------------------------------
   !Physical quantities----------------------------------------
   real(16), allocatable :: zeta(:), theta(:), K_jklm(:, :, :, :)
   !for integral index
   integer::i, j, k, l, m
   integer:: a_l, b_l, c_m, d_m, a_l_p, a_l_m, b_l_p, b_l_m, c_m_p, c_m_m, d_m_p, d_m_m
   allocate (zeta(0:Jmax), theta(0:Kmax))
   allocate (K_jklm(0:Jmax, 0:Kmax, 0:Lmax, 0:Mmax))
   !$omp parallel do
   do i = 0, Jmax
      zeta(i) = Cap_zeta*(i*1q0)/(Jmax*1q0)
   end do
   !$omp parallel do
   do k = 0, Kmax
      theta(k) = pi*(k*1q0)/(Kmax*1q0)
   end do
   !$omp parallel do
   do j = 0, Jmax!j=0では常にL1,L2は0
      do k = 0, Kmax
         do l = 0, Lmax
            do m = 0, Mmax
               zeta_l = zeta(l); theta_m = theta(m); zeta_j = zeta(j); theta_k = theta(k)

               if (j .eq. 0 .or. k .eq. 0 .or. k .eq. Kmax) then ! j=0のときは計算もっと楽
                  K_jklm(j, k, l, m) = 0.0q0

                  !----------------奇数奇数---------------
               else
                  if (mod(l, 2) .eq. 1 .and. mod(m, 2) .eq. 1) then!奇数奇数
                     a_l = l - 1; b_l = l + 1; c_m = m - 1; d_m = m + 1
                     alpha = zeta(l - 1); beta = zeta(l + 1); gamma = theta(m - 1); delta = theta(m + 1)
                     if (0 .le. a_l .and. a_l .le. Lmax .and. 0 .le. b_l .and. b_l .le. Lmax &
                         .and. 0 .le. c_m .and. c_m .le. Mmax .and. 0 .le. d_m .and. d_m .le. Mmax) then

                        if (zeta(a_l) .le. zeta_j .and. zeta_j .le. zeta(b_l) &
                            .and. theta(c_m) .le. theta_k .and. theta_k .le. theta(d_m)) then!特異点があるとき
                           print *, "singurality"
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l), zeta_j, theta(c_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!左下
                           val1 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l), theta(c_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!右下
                           val2 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l), zeta_j, theta_k, theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)!左上
                           val3 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l), theta_k, theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)!右上
                           val4 = s
                           s = 0.0q0
                           val_kiki_1 = val1 + val2 + val3 + val4

                        else
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l), zeta(b_l), theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_kiki_1 = s
                           print *, "no singurality"
                        end if
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, zeta(a_l), zeta(b_l), &
                                   theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_kiki_2 = s
                        s = 0.0q0
                     else
                        val_kiki_1 = 0.0q0
                        val_kiki_2 = 0.0q0
                     end if
                     K_jklm(j, k, l, m) = val_kiki_1 + val_kiki_2

                     !----------------偶数偶数---------------
                  else if (mod(l, 2) .eq. 0 .and. mod(m, 2) .eq. 0) then!偶数偶数
                     a_l_p = l; b_l_p = l + 2; c_m_p = m; d_m_p = m + 2; ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です
                     a_l_m = l - 2; b_l_m = l; c_m_m = m - 2; d_m_m = m  ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です

                     !-----------------------pp---------------------l_p,m_pになる要チェック
                     if (0 .le. a_l_p .and. a_l_p .le. Lmax .and. 0 .le. b_l_p .and. b_l_p .le. Lmax .and. &!pp
                         0 .le. c_m_p .and. c_m_p .le. Mmax .and. 0 .le. d_m_p .and. d_m_p .le. Mmax) then
                        alpha = zeta(l + 1); beta = zeta(l + 2); gamma = theta(m + 1); delta = theta(m + 2)

                        if (zeta(a_l_p) .le. zeta_j .and. zeta_j .le. zeta(b_l_p) &
                            .and. theta(c_m_p) .le. theta_k .and. theta_k .le. theta(d_m_p)) then!特異点があるとき
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_p), zeta_j, theta(c_m_p), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!左下
                           val1 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_p), theta(c_m_p), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!右下
                           val2 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_p), zeta_j, theta_k, theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)!左上
                           val3 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_p), theta_k, theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)!右上
                           val4 = s
                           s = 0.0q0
                           val_gugu_pp_1 = val1 + val2 + val3 + val4
                        else
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_p), zeta(b_l_p), theta(c_m_p), theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_gugu_pp_1 = s
                           s = 0.0q0
                        end if
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l_p), zeta(b_l_p), theta(c_m_p), theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_gugu_pp_2 = s
                        s = 0.0q0
                     else
                        val_gugu_pp_1 = 0.0q0
                        val_gugu_pp_2 = 0.0q0
                     end if
                     !-----------------------pm---------------------l_p,m_mになる要チェック
                     if (0 .le. a_l_p .and. a_l_p .le. Lmax .and. 0 .le. b_l_p .and. b_l_p .le. Lmax .and. &!pm
                         0 .le. c_m_m .and. c_m_m .le. Mmax .and. 0 .le. d_m_m .and. d_m_m .le. Mmax) then
                        alpha = zeta(l + 1); beta = zeta(l + 2); gamma = theta(m - 1); delta = theta(m - 2)

                        if (zeta(a_l_p) .le. zeta_j .and. zeta_j .le. zeta(b_l_p) &
                            .and. theta(c_m_m) .le. theta_k .and. theta_k .le. theta(d_m_m)) then!特異点があるとき
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_p), zeta_j, theta(c_m_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!左下
                           val1 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_p), theta(c_m_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!右下
                           val2 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_p), zeta_j, theta_k, theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)!左上
                           val3 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_p), theta_k, theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)!右上
                           val4 = s
                           s = 0.0q0
                           val_gugu_pm_1 = val1 + val2 + val3 + val4
                        else
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_p), zeta(b_l_p), theta(c_m_m), theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_gugu_pm_1 = s
                           s = 0.0q0
                        end if
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l_p), zeta(b_l_p), theta(c_m_m), theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_gugu_pm_2 = s
                        s = 0.0q0
                     else
                        val_gugu_pm_1 = 0.0q0
                        val_gugu_pm_2 = 0.0q0
                     end if
                     !-----------------------mp---------------------l_m,m_pになる要チェック
                     if (0 .le. a_l_m .and. a_l_m .le. Lmax .and. 0 .le. b_l_m .and. b_l_m .le. Lmax .and. &!mp
                         0 .le. c_m_p .and. c_m_p .le. Mmax .and. 0 .le. d_m_p .and. d_m_p .le. Mmax) then
                        alpha = zeta(l - 1); beta = zeta(l - 2); gamma = theta(m + 1); delta = theta(m + 2)

                        if (zeta(a_l_m) .le. zeta_j .and. zeta_j .le. zeta(b_l_m) &
                            .and. theta(c_m_p) .le. theta_k .and. theta_k .le. theta(d_m_p)) then!特異点があるとき
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_m), zeta_j, theta(c_m_p), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!左下
                           val1 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_m), theta(c_m_p), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!右下
                           val2 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_m), zeta_j, theta_k, theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)!左上
                           val3 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_m), theta_k, theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)!右上
                           val4 = s
                           s = 0.0q0
                           val_gugu_mp_1 = val1 + val2 + val3 + val4
                        else
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_m), zeta(b_l_m), theta(c_m_p), theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_gugu_mp_1 = s
                           s = 0.0q0
                        end if
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l_m), zeta(b_l_m), theta(c_m_p), theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_gugu_mp_2 = s
                        s = 0.0q0
                     else
                        val_gugu_mp_1 = 0.0q0
                        val_gugu_mp_2 = 0.0q0
                     end if
                     !-----------------------mm---------------------l_m,m_mになる要チェック
                     if (0 .le. a_l_m .and. a_l_m .le. Lmax .and. 0 .le. b_l_m .and. b_l_m .le. Lmax .and. &
                         0 .le. c_m_m .and. c_m_m .le. Mmax .and. 0 .le. d_m_m .and. d_m_m .le. Mmax) then
                        alpha = zeta(l - 1); beta = zeta(l - 2); gamma = theta(m - 1); delta = theta(m - 2)

                        if (zeta(a_l_m) .le. zeta_j .and. zeta_j .le. zeta(b_l_m) &
                            .and. theta(c_m_m) .le. theta_k .and. theta_k .le. theta(d_m_m)) then!特異点があるとき
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_m), zeta_j, theta(c_m_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!左下
                           val1 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_m), theta(c_m_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!右下
                           val2 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_m), zeta_j, theta_k, theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)!左上
                           val3 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_m), theta_k, theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)!右上
                           val4 = s
                           s = 0.0q0
                           val_gugu_mm_1 = val1 + val2 + val3 + val4
                        else
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_m), zeta(b_l_m), theta(c_m_m), theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_gugu_mm_1 = s
                           s = 0.0q0
                        end if
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l_m), zeta(b_l_m), theta(c_m_m), theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_gugu_mm_2 = s
                        s = 0.0q0
                     else
                        val_gugu_mm_1 = 0.0q0
                        val_gugu_mm_2 = 0.0q0
                     end if
                     K_jklm(j, k, l, m) = val_gugu_pp_1 + val_gugu_pp_2 + val_gugu_pm_1 + val_gugu_pm_2 &
                                          + val_gugu_mp_1 + val_gugu_mp_2 + val_gugu_mm_1 + val_gugu_mm_2
!-----------------奇数偶数------------------!-----------------奇数偶数------------------!-----------------奇数偶数------------------
                     !-----------------奇数偶数------------------
                  else if (mod(l, 2) .eq. 1 .and. mod(m, 2) .eq. 0) then
                     a_l = l - 1; b_l = l + 1; c_m_p = m; d_m_p = m + 2; ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です
                     c_m_m = m - 2; d_m_m = m  ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です

                     !-----------------------0p---------------------m_pになる要チェック
                     if (0 .le. a_l .and. a_l .le. Lmax .and. 0 .le. b_l .and. b_l .le. Lmax .and. &!pp
                         0 .le. c_m_p .and. c_m_p .le. Mmax .and. 0 .le. d_m_p .and. d_m_p .le. Mmax) then
                        alpha = zeta(l + 1); beta = zeta(l - 1); gamma = theta(m + 1); delta = theta(m + 2)

                        if (zeta(a_l) .le. zeta_j .and. zeta_j .le. zeta(b_l) &
                            .and. theta(c_m_p) .le. theta_k .and. theta_k .le. theta(d_m_p)) then!特異点があるとき
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l), zeta_j, theta(c_m_p), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!左下
                           val1 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l), theta(c_m_p), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!右下
                           val2 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l), zeta_j, theta_k, theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)!左上
                           val3 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l), theta_k, theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)!右上
                           val4 = s
                           s = 0.0q0
                           val_kigu_0p_1 = val1 + val2 + val3 + val4
                        else
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l), zeta(b_l), theta(c_m_p), theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_kigu_0p_1 = s
                           s = 0.0q0
                        end if
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l), zeta(b_l), theta(c_m_p), theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_kigu_0p_2 = s
                        s = 0.0q0
                     else
                        val_kigu_0p_1 = 0.0q0
                        val_kigu_0p_2 = 0.0q0
                     end if
                     !-----------------------0m---------------------m_mになる要チェック
                     if (0 .le. a_l .and. a_l .le. Lmax .and. 0 .le. b_l .and. b_l .le. Lmax .and. &!pm
                         0 .le. c_m_m .and. c_m_m .le. Mmax .and. 0 .le. d_m_m .and. d_m_m .le. Mmax) then
                        alpha = zeta(l + 1); beta = zeta(l - 1); gamma = theta(m - 1); delta = theta(m - 2)

                        if (zeta(a_l) .le. zeta_j .and. zeta_j .le. zeta(b_l) &
                            .and. theta(c_m_m) .le. theta_k .and. theta_k .le. theta(d_m_m)) then!特異点があるとき
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l), zeta_j, theta(c_m_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!左下
                           val1 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l), theta(c_m_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!右下
                           val2 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l), zeta_j, theta_k, theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)!左上
                           val3 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l), theta_k, theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)!右上
                           val4 = s
                           s = 0.0q0
                           val_kigu_0m_1 = val1 + val2 + val3 + val4
                        else
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l), zeta(b_l), theta(c_m_m), theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_kigu_0m_1 = s
                           s = 0.0q0
                        end if
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l), zeta(b_l), theta(c_m_m), theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_kigu_0m_2 = s
                        s = 0.0q0
                     else
                        val_kigu_0m_1 = 0.0q0
                        val_kigu_0m_2 = 0.0q0
                     end if
                     K_jklm(j, k, l, m) = val_kigu_0p_1 + val_kigu_0p_2 + val_kigu_0m_1 + val_kigu_0m_2

                     !------------------偶数奇数-------------------!------------------偶数奇数-------------------!------------------偶数奇数-------------------!------------------偶数奇数-------------------!------------------偶数奇数-------------------
                  else if (mod(l, 2) .eq. 0 .and. mod(m, 2) .eq. 1) then
                     a_l_p = l + 1; b_l_p = l + 2; a_l_m = l - 1; b_l_p = l - 2 ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です
                     c_m = m - 1; d_m = m + 1  ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です

                     !-----------------------p0---------------------l_pになる要チェック
                     if (0 .le. a_l_p .and. a_l_p .le. Lmax .and. 0 .le. b_l_p .and. b_l_p .le. Lmax .and. &!pp
                         0 .le. c_m .and. c_m .le. Mmax .and. 0 .le. d_m .and. d_m .le. Mmax) then
                        alpha = zeta(l + 1); beta = zeta(l + 2); gamma = theta(m + 1); delta = theta(m - 1)

                        if (zeta(a_l_p) .le. zeta_j .and. zeta_j .le. zeta(b_l_p) &
                            .and. theta(c_m) .le. theta_k .and. theta_k .le. theta(d_m)) then!特異点があるとき
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_p), zeta_j, theta(c_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!左下
                           val1 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_p), theta(c_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!右下
                           val2 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_p), zeta_j, theta_k, theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)!左上
                           val3 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_p), theta_k, theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)!右上
                           val4 = s
                           s = 0.0q0
                           val_guki_p0_1 = val1 + val2 + val3 + val4
                        else
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_p), zeta(b_l_p), theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_guki_p0_1 = s
                           s = 0.0q0
                        end if
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l_p), zeta(b_l_p), theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_guki_p0_2 = s
                        s = 0.0q0
                     else
                        val_guki_p0_1 = 0.0q0
                        val_guki_p0_2 = 0.0q0
                     end if
                     !-----------------------m0---------------------l_mになる要チェック
                     if (0 .le. a_l_m .and. a_l_m .le. Lmax .and. 0 .le. b_l_m .and. b_l_m .le. Lmax .and. &!pm
                         0 .le. c_m .and. c_m .le. Mmax .and. 0 .le. d_m .and. d_m .le. Mmax) then
                        alpha = zeta(l - 1); beta = zeta(l - 2); gamma = theta(m - 1); delta = theta(m + 1)

                        if (zeta(a_l_m) .le. zeta_j .and. zeta_j .le. zeta(b_l_m) &
                            .and. theta(c_m) .le. theta_k .and. theta_k .le. theta(d_m)) then!特異点があるとき
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_m), zeta_j, theta(c_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!左下
                           val1 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_m), theta(c_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!右下
                           val2 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_m), zeta_j, theta_k, theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)!左上
                           val3 = s
                           s = 0.0q0
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta_j, zeta(b_l_m), theta_k, theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)!右上
                           val4 = s
                           s = 0.0q0
                           val_guki_m0_1 = val1 + val2 + val3 + val4
                        else
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_m), zeta(b_l_m), theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_guki_m0_1 = s
                           s = 0.0q0
                        end if
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l_m), zeta(b_l_m), theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_guki_m0_2 = s
                        s = 0.0q0
                     else
                        val_guki_m0_1 = 0.0q0
                        val_guki_m0_2 = 0.0q0
                     end if
                     K_jklm(j, k, l, m) = val_guki_p0_1 + val_guki_p0_2 + val_guki_m0_1 + val_guki_m0_2

                  end if

               end if
               print *, "j", j, "k", k, "l", l, "m", m, K_jklm(j, k, l, m)
            end do
         end do
      end do
   end do

   !if (l .eq. 0 .and. m .eq. 0) then j=0のときは計算もっと楽

   ! else if (l .eq. 0 .and. m .eq. Mmax) then
   !    continue
   ! else if (l .eq. Lmax .and. m .eq. 0) then
   !    continue
   ! else if (l .eq. Lmax .and. m .eq. Mmax) then
   !    continu
   !else

   !print *, "j=", j, "k=", k, "l=", l, "m=", m, "K_jklm=", K_jklm(j, k, l, m)
   ! if (mod(l, 2) .eq. 1 .and. mod(m, 2) .eq. 1) then
   !    a_l = l - 1; b_l = l + 1; c_m = m - 1; d_m = m + 1
   !    zeta_l = zeta(j); theta_m = theta(k)
   !    alpha = zeta(j - 1); beta = zeta(j + 1); gamma = theta(k - 1); delta = theta(k + 1)
   !    if (0 .le. a_l .and. a_l .le. Lmax .and. 0 .le. b_l .and. b_l .le. Lmax &
   !        .and. 0 .le. c_m .and. c_m .le. Mmax .and. 0 .le. d_m .and. d_m .le. Mmax) then
   !       call dde3d(L1, zeta_l, theta_m, alpha, beta, gamma, delta, zeta(a_l), zeta(b_l), &
   !                  theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)
   !       val_kiki_1 = s
   !       s = 0.0q0
   !       call dde3d(L2, zeta_l, theta_m, alpha, beta, gamma, delta, zeta(a_l), zeta(b_l), &
   !                  theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)

   !       val_kiki_2 = s
   !       s = 0.0q0
   !    else
   !       val_kiki_1 = 0.0q0
   !       val_kiki_2 = 0.0q0
   !    end if
   !    K_jklm(j, k, l, m) = val_kiki_1 + val_kiki_2
   ! j = 1; k = 1; l = 1; m = 5
   ! a_l = l - 1; b_l = l + 1; c_m = m - 1; d_m = m + 1
   ! zeta_l = zeta(j); theta_m = theta(k)
   ! alpha = zeta(j - 1); beta = zeta(j + 1); gamma = theta(k - 1); delta = theta(k + 1)
   ! call dde3d(L1, zeta_l, theta_m, alpha, beta, gamma, delta, zeta(a_l), zeta(b_l), &
   !            theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)
   ! val_kiki_1 = s
   ! s = 0.0q0
   ! call dde3d(L2, zeta_l, theta_m, alpha, beta, gamma, delta, zeta(a_l), zeta(b_l), &
   !            theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)
   ! val_kiki_2 = s
   ! s = 0.0q0
   ! K_jklm(j, k, l, m) = val_kiki_1 + val_kiki_2
   ! print *, K_jklm(j, k, l, m)
   ! call dde3d(L1, zeta(1), theta(1), alpha, beta, gamma, delta, 0.0_16, 1.0_16, 0.0_16, pi/100.0_16, 0.0_16, 2*pi, eps, s, info)
   ! print *, s! 9.596087731032261608563067133300942E-0010
   ! print *, L1(1.0q0, 1.0q0, 1.0q0, zeta(1), theta(1), alpha, beta, gamma, delta)! -8.022711679671201503181953271830736E-0005
!https://www.wolframalpha.com/input?i=Integrate%5BIntegrate%5BIntegrate%5Bx%5E2%2B1-2*x*cos%28y%29*cos%28%CF%80%2F100%29-2*x*sin%28y%29*sin%28%CF%80%2F100%29*cos%28z%29%2C%7Bz%2C0%2C2*%CF%80%7D%5D%2C%7By%2C0%2C%CF%80%2F100%7D%5D%2C%7Bx%2C0%2C1%7D%5D&lang=ja
   ! -1.961747754100960317717098015371201E-0009
contains
   function F1(zeta_bar, theta_bar, psi, zeta_j, theta_k)!計算できる
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_j, theta_k
      real(16)::F1

      F1 = zeta_bar**2 + zeta_j**2 - 2*zeta_bar*zeta_j*cos(theta_bar)*cos(theta_k) &
           - 2*zeta_bar*zeta_j*sin(theta_bar)*sin(theta_k)*cos(psi)
   end function
   function F2(zeta_bar, theta_bar, psi, zeta_j, theta_k)
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_j, theta_k
      real(16)::F2
      F2 = zeta_bar**2*zeta_j**2 &
           - zeta_bar**2*zeta_j**2 &
           *(cos(theta_k)*cos(theta_bar) + sin(theta_bar)*sin(theta_k)*cos(psi))**2
   end function
   function f_for_psi(zeta_bar, theta_bar, zeta_l, theta_m, alpha, beta, gamma, delta)
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, alpha, beta, gamma, delta, zeta_l, theta_m
      real(16) :: f_for_psi
      f_for_psi = (zeta_bar - alpha)*(zeta_bar - beta)*(theta_bar - gamma)*(theta_bar - delta) &
                  /(zeta_l - alpha)*(zeta_l - beta)*(theta_m - gamma)*(theta_m - delta)
   end function
   function L1(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)!計算できない
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_j, zeta_l, theta_m, theta_k, alpha, beta, gamma, delta
      real(16)::L1

      L1 = zeta_bar*zeta_bar*sin(theta_bar)/sqrt(2.0q0)/pi*cos(psi) &
           /sqrt(F1(zeta_bar, theta_bar, psi, zeta_j, theta_k)) &
           *exp(-zeta_bar**2 + F2(zeta_bar, theta_bar, psi, zeta_j, theta_k) &
                /F1(zeta_bar, theta_bar, psi, zeta_j, theta_k)) &
           *f_for_psi(zeta_bar, theta_bar, zeta_l, theta_m, alpha, beta, gamma, delta)
   end function
   function L1_check(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)!計算できない
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_j, zeta_l, theta_m, theta_k, alpha, beta, gamma, delta
      real(16)::L1_check

      L1_check = 1/sqrt(F1(zeta_bar, theta_bar, psi, zeta_j, theta_k))
   end function!時間かかるだけでちゃんと機能してる
   function L2(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)!計算できない
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta
      real(16)::L2

      L2 = zeta_bar*zeta_bar*sin(theta_bar)/2/sqrt(2.0q0)/pi*cos(psi) &
           *sqrt(F1(zeta_bar, theta_bar, psi, zeta_j, theta_k)) &
           *exp(-zeta_bar**2) &
           *f_for_psi(zeta_bar, theta_bar, zeta_l, theta_m, alpha, beta, gamma, delta)
   end function
   !つねに0
   !https://www.wolframalpha.com/input?i=Integrate%5BIntegrate%5BIntegrate%5Bx%5E2*cos%28z%29*sin%28y%29*exp%28-x%5E2%29*%28x%29*%28x-2%29*%28y%29*%28y-%CF%80%2F50%29%2F%28sqrt%282%29*%CF%80*%281-2%29*%28%CF%80%2F100%29*%28%CF%80%2F100-%CF%80%2F50%29%29%2C%7Bx%2C0%2C2%7D%5D%2C%7By%2C0%2C%CF%80%2F50%7D%5D%2C%7Bz%2C0%2C2*%CF%80%7D%5D&lang=ja
   function L1_j0(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta
      real(16)::L1_j0

      L1_j0 = zeta_bar*sin(theta_bar)/sqrt(2.0q0)/pi*cos(psi) &
              *exp(-zeta_bar**2) &
              *f_for_psi(zeta_bar, theta_bar, zeta_l, theta_m, alpha, beta, gamma, delta)
   end function
   ! function f_check(zeta_bar, theta_bar, psi)
   !    implicit none
   !    real(16), intent(in)::zeta_bar, theta_bar, psi
   !    real(16)::zeta, theta
   !    real(16)::f_check
   !    !zeta = 1.0_16
   !    !theta = pi/100.0_16

   !    f_check = 1 &!((zeta_bar - 0.0_16)*(zeta_bar - 2.0_16)*(theta_bar - 0.0_16)*(theta_bar - pi/50.0_16)) &
   !              /((zeta - 0.0_16)*(zeta - 2.0_16)*(theta - 0.0_16)*(theta - pi/50.0_16))!これが収束しないのはおかしい
   ! end function
end program
