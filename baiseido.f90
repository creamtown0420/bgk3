module DE_integration
   !
   ! Double Exponential (tanh-sinh) Formula
   !       for 1d, 2d, 3d.
   !+------------------------------------+
   ! For 1 dimension integration
   !   call dde1d(f,a,b,eps,s,info)
   !   call cde1d(f,a,b,eps,s,info)
   !   call dde1d_hinf(f,a,eps,s,info)
   !   call cde1d_hinf(f,a,eps,s,info)
   !   call dde1d_ehinf(f,a,eps,s,info)
   !   call cde1d_ehinf(f,a,eps,s,info)
   !   call dde1d_inf(f,eps,s,info)
   !   call cde1d_inf(f,eps,s,info)
   !
   !   call zde1d(f,a,b,eps,s,info)
   !   call zde1d_hinf(f,a,theta,eps,s,info)
   !   call zde1d_inf(f,a,theta,eps,s,info)
   !
   !    d/cde1d,      x:[a,b]
   !       [in]  (dp/cp)f, (dp)a,b, (dp)eps
   !       [out] (dp/cp)s, (int)info
   !
   !    d/cde1d_hinf, x:[a,+infty]
   !       [in]  (dp/cp)f, (dp)a, (dp)eps
   !       [out] (dp/cp)s, (int)info
   !
   !    d/cde1d_inf,  x:[-infty,+infty]
   !       [in]  (dp/cp)f, (dp)eps
   !       [out] (dp/cp)s, (int)info
   !
   !    zde1d,        x:[a,b]
   !       [in]  (cp)f, (cp)a,b, (dp)eps
   !       [out] (cp)s, (int)info
   !
   !    zde1d_hinf,   x:[a,a+infty*e^{i*theta}]
   !       [in]  (cp)f, (cp)a,   (dp)eps, (dp)theta
   !       [out] (cp)s, (int)info
   !
   !    zde1d_inf,    x:[a-infty*e^{i*theta},a+infty*e^{i*theta}]
   !       [in]  (cp)f, (cp)a,   (dp)eps
   !       [out] (cp)s, (int)info
   !
   !+------------------------------------+
   ! For 2 dimension integration
   !   call dde2d(f,xa,xb,ya,yb,eps,s,info)
   !   call cde2d(f,xa,xb,ya,yb,eps,s,info)
   !   call dde2d_inf(f,eps,s,info)
   !   call cde2d_inf(f,eps,s,info)
   !
   !   call dde_polar(f,eps,s,info)
   !   call cde_polar(f,eps,s,info)
   !
   !   call zde2d(f,xa,xb,ya,yb,eps,s,info)
   !   call zde2d_inf(f,xa,thx,ya,thy,eps,s,info)
   !
   !    d/cde2d,      x:[xa,xb] y:[ya,yb]
   !       [in]  (dp/cp)f, (dp)xa,xb,ya,yb, (dp)eps
   !       [out] (dp/cp)s, (int)info
   !
   !    d/cde2d_inf,  x,y:[-infty,+infty]
   !       [in]  (dp/cp)f, (dp)eps
   !       [out] (dp/cp)s, (int)info
   !
   !    d/cde_polar,  r:[0,+infty], theta:[0,2*pi]
   !       [in]  (dp/cp)f, (dp)eps
   !       [out] (dp/cp)s, (int)info
   !
   !    zde2d,        x:[xa,xb] y:[ya,yb]
   !       [in]  (cp)f, (cp)xa,xb,ya,yb, (dp)eps
   !       [out] (cp)s, (int)info
   !
   !    zde2d_inf,  x:[xa-infty*e^{i*thx},xa+infty*e^{i*thx}]
   !                y:[ya-infty*e^{i*thy},ya+infty*e^{i*thy}]
   !       [in]  (cp)f, (cp)xa,ya, (dp)thx,thy, (dp)eps
   !       [out] (cp)s, (int)info
   !
   !+------------------------------------+
   ! For 3 dimension integration
   !   call dde3d(f,xa,xb,ya,yb,za,zb,eps,s,info)
   !   call cde3d(f,xa,xb,ya,yb,za,zb,eps,s,info)
   !   call dde3d_inf(f,eps,s,info)
   !   call cde3d_inf(f,eps,s,info)
   !
   !   call dde_spherical(f,eps,s,info)
   !   call cde_spherical(f,eps,s,info)
   !
   !   call zde3d(f,xa,xb,ya,yb,za,zb,eps,s,info)
   !   call zde3d_inf(f,xa,thx,ya,thy,za,thz,eps,s,info)
   !
   !    d/cde3d,      x:[xa,xb] y:[ya,yb] z:[za,zb]
   !       [in]  (dp/cp)f, (dp)xa,xb,ya,yb,za,zb (dp)eps
   !       [out] (dp/cp)s, (int)info
   !
   !    d/cde3d_inf,  x,y,z:[-infty,+infty]
   !       [in]  (dp/cp)f, (dp)eps
   !       [out] (dp/cp)s, (int)info
   !
   !    d/cde_spherical, r:[0,+infty], theta:[0,pi], phi:[0,2*pi]
   !       [in]  (dp/cp)f, (dp)eps
   !       [out] (dp/cp)s, (int)info
   !
   !    zde3d,        x:[xa,xb] y:[ya,yb] z:[za,zb]
   !       [in]  (cp)f, (cp)xa,xb,ya,yb,za,zb (dp)eps
   !       [out] (cp)s, (int)info
   !
   !    zde3d_inf,  x:[xa-infty*e^{i*thx},xa+infty*e^{i*thx}]
   !                y:[ya-infty*e^{i*thy},ya+infty*e^{i*thy}]
   !                z:[za-infty*e^{i*thz},za+infty*e^{i*thz}]
   !       [in]  (cp)f, (cp)xa,ya,za (dp)thx,thy,thz (dp)eps
   !       [out] (cp)s, (int)info
   !
   !+-------------------------------------+
   ! (dp)  real(16)
   ! (cp)  complex(kind(0.0_16))
   ! (int) integer
   !
   ! Output parameters
   !  s    : integration result
   !  info : 0 --> relative error converge less than 'eps'.
   !       : 1 --> didn't converge in somewhere.
   !
   ! http://slpr.sakura.ne.jp/qp
   !        2016/12/17 (yyyy/mm/dd) by sikino
   !
   implicit none
   integer, private, parameter::kmin = 3
   integer, private, parameter::kmax = 14
   ! max order of DE integration
   real(16), private, parameter::hr = 6.0_16    ! [-1,1]
   real(16), private, parameter::hrhi = 7.6_16 ! [ 0,\infty]
   real(16), private, parameter::hrai = 7.6_16 ! [-\infty,\infty]
   real(16), private, parameter::c0 = 0.01_16  ! safe value for eps
   real(16), private, parameter::pi2 = 1.5707963267948966_16
   interface de3d
      module procedure &
         dde3d
   end interface de3d
contains
   subroutine dde3d(f, xa, xb, ya, yb, za, zb, eps, s, info)
      ! Integrate for 3d
      implicit none
      real(16), intent(in)::xa, xb, ya, yb, za, zb, eps
      real(16), intent(out)::s
      integer, intent(out)::info
      interface
         function f(x, y, z)
            implicit none
            real(16), intent(in)::x, y, z
            real(16)::f
         end function f
      end interface
      integer::j, nc, l
      real(16)::h, s0, xt, wt, t, ft, as, shk, mba, pba, err, seps
      info = 0
      seps = 0.001_16!c0*sqrt(eps)
      h = hr*0.5_16

      mba = (xb - xa)*0.5_16
      pba = (xb + xa)*0.5_16
      call dde3dsyz(f, ya, yb, za, zb, pba, eps, ft, info)
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
            call dde3dsyz(f, ya, yb, za, zb, mba*xt + pba, eps, ft, info)
            s = s + ft*wt
            print *, "s", s
!              s -8.524989881935718142271927179999396E-0021
!  s  1.680295746484329267527464109076638E-0010
!  s -1.823983615679643026042224841586826E-0034
!  s  5.131505265615405285741562101248373E-0018
!  s -4.010222938975687728945849556479320E-0011
!  s -3.997488451027711625848268378315298E-0011
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

   subroutine dde3dsyz(f, ya, yb, za, zb, xc, eps, s, info)
      ! Integrate for 3d about y,z
      implicit none
      real(16), intent(in)::ya, yb, za, zb, xc, eps
      real(16), intent(out)::s
      integer, intent(inout)::info
      interface
         function f(x, y, z)
            implicit none
            real(16), intent(in)::x, y, z
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
      call dde3dsz(f, za, zb, xc, pba, eps, ft, info)
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
            call dde3dsz(f, za, zb, xc, mba*xt + pba, eps, ft, info)
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

   subroutine dde3dsz(f, za, zb, xc, yc, eps, s, info)
      ! Integrate for 3d about z
      implicit none
      real(16), intent(in)::za, zb, xc, yc, eps
      real(16), intent(out)::s
      integer, intent(inout)::info
      interface
         function f(x, y, z)
            implicit none
            real(16), intent(in)::x, y, z
            real(16)::f
         end function f
      end interface
      integer::j, nc, l
      real(16)::h, s0, xt, wt, t, as, shk, mba, pba, err, seps
      seps = c0*sqrt(eps)
      h = hr*0.5_16

      mba = (zb - za)*0.5_16
      pba = (zb + za)*0.5_16
      s0 = f(xc, yc, pba)*h*pi2*mba

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
            s = s + f(xc, yc, mba*xt + pba)*wt
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
   real(16), parameter::pi = 3.141592653589793238462643383279502884_16
   real(16):: zeta, theta_zeta, eps, s
   integer::info!q0
   zeta = 1.0_16
   theta_zeta = pi/100.0_16
   eps = 1d-24; s = 0.0_16
!xa=0.0_16; xb=1.0_16; ya=0.0_16; yb=2.0_16; za=2.0_16; zb=3.0_16
   call dde3d(L1, 0.0_16, 1.0_16, 0.0_16, pi/100.0_16, 0.0_16, 2*pi, eps, s, info)
   !call dde3d(f3, 0.0_16, 1.0_16, 0.0_16, 2.0_16, 2.0_16, 3.0_16, eps, s, info)
   print *, s
   print *, L1(1.0q0, 1.0q0, 1.0q0)! -8.022711679671201503181953271830738E-0005
contains
   function zeta_hiku_zeta_bar(zeta_bar_real, theta_zeta_bar_real, psi_real)!計算できる
      implicit none
      real(16), intent(in)::zeta_bar_real, theta_zeta_bar_real, psi_real
      real(16)::zeta, theta_zeta
      real(16)::zeta_hiku_zeta_bar
      zeta = 1.0_16
      theta_zeta = pi/100.0_16
      zeta_hiku_zeta_bar = zeta_bar_real**2 + zeta**2 - 2*zeta_bar_real*zeta*cos(theta_zeta_bar_real)*cos(theta_zeta) &
                           - 2*zeta_bar_real*zeta*sin(theta_zeta_bar_real)*sin(theta_zeta)*cos(psi_real)
      !   if (zeta_bar_real**2 + zeta**2 - 2*zeta_bar_real*zeta*cos(theta_zeta_bar_real)*cos(theta_zeta) &
      !       - 2*zeta_bar_real*zeta*sin(theta_zeta_bar_real)*sin(theta_zeta)*cos(psi_real) == 0) then
      !      print *, zeta_bar_real**2 + zeta**2 - 2*zeta_bar_real*zeta*cos(theta_zeta_bar_real)*cos(theta_zeta) &
      !         - 2*zeta_bar_real*zeta*sin(theta_zeta_bar_real)*sin(theta_zeta)*cos(psi_real)
      !   end if
   end function zeta_hiku_zeta_bar
   function zeta_gaiseki_zeta_bar(zeta_bar_real, theta_zeta_bar_real, psi_real)
      implicit none
      real(16), intent(in)::zeta_bar_real, theta_zeta_bar_real, psi_real
      real(16):: zeta, theta_zeta
      real(16)::zeta_gaiseki_zeta_bar
      zeta = 1.0_16
      theta_zeta = pi/100.0_16
      zeta_gaiseki_zeta_bar = zeta_bar_real**2*zeta**2 &
                              - zeta_bar_real**2*zeta**2 &
                             *(cos(theta_zeta)*cos(theta_zeta_bar_real) + sin(theta_zeta_bar_real)*sin(theta_zeta)*cos(psi_real))**2
   end function
   function f_for_psi(zeta_bar_real, theta_zeta_bar_real, alpha, beta, gamma, delta)
      implicit none
      real(16), intent(in)::zeta_bar_real, theta_zeta_bar_real, alpha, beta, gamma, delta
      real(16)::zeta, theta_zeta
      real(16) :: f_for_psi
      zeta = 1.0_16
      theta_zeta = pi/100.0_16
      f_for_psi = (zeta_bar_real - alpha)*(zeta_bar_real - beta)*(theta_zeta_bar_real - gamma)*(theta_zeta_bar_real - delta) &
                  /(zeta - alpha)*(zeta - beta)*(theta_zeta - gamma)*(theta_zeta - delta)
   end function
   function L1(zeta_bar_real, theta_zeta_bar_real, psi_real)!計算できない
      implicit none
      real(16), intent(in)::zeta_bar_real, theta_zeta_bar_real, psi_real
      real(16)::zeta, theta_zeta
      real(16)::L1
      zeta = 1.0_16
      theta_zeta = pi/100.0_16

      L1 = zeta_bar_real*zeta_bar_real*sin(theta_zeta_bar_real)/sqrt(2.0q0)/pi*cos(psi_real) &
           /sqrt(zeta_hiku_zeta_bar(zeta_bar_real, theta_zeta_bar_real, psi_real)) &
           *exp(-zeta_bar_real**2 + zeta_gaiseki_zeta_bar(zeta_bar_real, theta_zeta_bar_real, psi_real) &
                /zeta_hiku_zeta_bar(zeta_bar_real, theta_zeta_bar_real, psi_real)) &
           *f_for_psi(zeta_bar_real, theta_zeta_bar_real, 0.0_16, 2.0_16, 0.0_16, pi/50.0_16)
   end function
   function f_check(zeta_bar_real, theta_zeta_bar_real, psi_real)
      implicit none
      real(16), intent(in)::zeta_bar_real, theta_zeta_bar_real, psi_real
      real(16)::zeta, theta_zeta
      real(16)::f_check
      !zeta = 1.0_16
      !theta_zeta = pi/100.0_16

      f_check = 1 &!((zeta_bar_real - 0.0_16)*(zeta_bar_real - 2.0_16)*(theta_zeta_bar_real - 0.0_16)*(theta_zeta_bar_real - pi/50.0_16)) &
                /((zeta - 0.0_16)*(zeta - 2.0_16)*(theta_zeta - 0.0_16)*(theta_zeta - pi/50.0_16))!これが収束しないのはおかしい
   end function
   function f3(x, y, z)
      implicit none
      real(16), intent(in)::x, y, z
      real(16)::f3

      f3 = sin(x*z)*exp(-y)

      return
   end function f3
end program
