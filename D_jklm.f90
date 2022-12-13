module DE_integration
   implicit none
   integer, private, parameter::kmin = 3
   integer, private, parameter::kmax = 14
   ! max order of DE integration
   real(16), private, parameter::hr = 6.0_16    ! [-1,1]
   real(16), private, parameter::c0 = 0.01_16  ! safe value for eps
   real(16), private, parameter::pi2 = 1.5707963267948966_16
   interface de3d
      module procedure &
         dde3d
   end interface de3d
contains
   subroutine dde3d(f, zeta_f, theta_f, alpha, beta, gamma, delta, xa, xb, ya, yb, za, zb, eps, s, info)
      ! Integrate for 3d
      implicit none
      real(16), intent(in)::xa, xb, ya, yb, za, zb, eps, zeta_f, theta_f, alpha, beta, gamma, delta
      real(16), intent(out)::s
      integer, intent(out)::info
      interface
         function f(x, y, z, zeta_f, theta_f, alpha, beta, gamma, delta)
            implicit none
            real(16), intent(in)::x, y, z, zeta_f, theta_f, alpha, beta, gamma, delta
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
      call dde3dsyz(f, zeta_f, theta_f, alpha, beta, gamma, delta, ya, yb, za, zb, pba, eps, ft, info)
      !ここまでは一致する
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
            call dde3dsyz(f, zeta_f, theta_f, alpha, beta, gamma, delta, ya, yb, za, zb, mba*xt + pba, eps, ft, info)
            s = s + ft*wt
            print *, "s", s
!              s -8.524989881935718142271927179999396E-0021
!  s  1.680295746484329267527464109076638E-0010
!  s -1.823983615679643026042224841586826E-0034
!  s  5.131505265615405285741562101248373E-0018
!  s -4.010222938975687728945849556479320E-0011
!  s -3.997488451027711625848268378315298E-0011
!  s -3.563154088907465803585834999578505E-0046
!  s -2.139003248256175712322741044310125E-0026
!  s -3.365281868124206528948217424300352E-0019
!  s  1.293539499958389745697891835628484E-0015
!  s -2.315916291784274698013901791667295E-0012
!  s  2.961524890664837224247027005729159E-0010
!  s  3.061924162148629907113268374645792E-0010
!  s  3.061926148049446586054652861473502E-0010
!  s -6.764203695057187019374457142749371E-0054
!  s -8.964960314536543881738848949584472E-0040
!  s -4.649889460849318695213279303507867E-0030
!  s -2.406333927520666921421127871963126E-0023
!  s -1.171424133454018315127314920096786E-0018
!  s -1.094028869718857775291285907333043E-0018
!  s  1.101560146085826693553389853205558E-0016
!  s  7.262594039840567518078057011095000E-0015
!  s -2.689930992411956222432874245983268E-0013
!  s -1.254407430637054587655196077117998E-0011
!  s  6.757860902121509039841568965448283E-0013
!  s  3.471716789683692563237147589222514E-0010
!  s  3.966912287566047771008145254440072E-0010
!  s  3.980936578170030003617259452215925E-0010
!  s  3.981005417259366469079012315913096E-0010
!  s  3.981005444056825716781853586965566E-0010
!  s -2.455489559932856043105707940051298E-0058
!  s -7.449773976756763305019173490303596E-0050
!  s -7.980390581775544141418837312217674E-0043
!  s -5.377884086624719286016008596231901E-0037
!  s -3.684843481168565161717269232645409E-0032
!  s -3.828341798405611073146571758887075E-0028
!  s -8.415693903341279282915920200713214E-0025
!  s -5.168628589869872572310307700468117E-0022
!  s -1.117314500227204253435396940053238E-0019
!  s -2.021489059990860644681284576401732E-0019
!  s -2.276870457241373481694708597224579E-0018
!  s -1.404538712029818576753263331680402E-0018
!  s  2.439439202271895276037653369109326E-0017
!  s  4.355851880707621905152400644415032E-0016
!  s  3.839951004366279867294662027286201E-0015
!  s  1.410016026466141491523838075991095E-0014
!  s -4.522187319769483467234209911826642E-0014
!  s -9.143061296353893329455926156016769E-0013
!  s -6.480974127307439802389728954410076E-0012
!  s -3.086042052518945926137001993141370E-0011
!  s -7.247299288131502605853784270831846E-0011
!  s  7.451942367910098689478504631000217E-0011
!  s  4.480053586011999509870794898799065E-0010
!  s  7.087234654738819832780494986821974E-0010
!  s  8.046776669155055224069814954750882E-0010
!  s  8.279946786680598189171540597027747E-0010
!  s  8.319366492941674421062616508166730E-0010
!  s  8.323848474783108153526734991338054E-0010
!  s  8.324166395723652315833576718092423E-0010
!  s  8.324179140109959207730632233434306E-0010
!  s  8.324179396133674689601114058634245E-0010
!  s  8.324179398364215219310576915999175E-0010
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

   subroutine dde3dsyz(f, zeta_f, theta_f, alpha, beta, gamma, delta, ya, yb, za, zb, xc, eps, s, info)
      ! Integrate for 3d about y,z
      implicit none
      real(16), intent(in)::ya, yb, za, zb, xc, eps, zeta_f, theta_f, alpha, beta, gamma, delta
      real(16), intent(out)::s
      integer, intent(inout)::info
      interface
         function f(x, y, z, zeta_f, theta_f, alpha, beta, gamma, delta)
            implicit none
            real(16), intent(in)::x, y, z, zeta_f, theta_f, alpha, beta, gamma, delta
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
      call dde3dsz(f, zeta_f, theta_f, alpha, beta, gamma, delta, za, zb, xc, pba, eps, ft, info)
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
            call dde3dsz(f, zeta_f, theta_f, alpha, beta, gamma, delta, za, zb, xc, mba*xt + pba, eps, ft, info)
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

   subroutine dde3dsz(f, zeta_f, theta_f, alpha, beta, gamma, delta, za, zb, xc, yc, eps, s, info)
      ! Integrate for 3d about z
      implicit none
      real(16), intent(in)::za, zb, xc, yc, eps, zeta_f, theta_f, alpha, beta, gamma, delta
      real(16), intent(out)::s
      integer, intent(inout)::info
      interface
         function f(x, y, z, zeta_f, theta_f, alpha, beta, gamma, delta)
            implicit none
            real(16), intent(in)::x, y, z, zeta_f, theta_f, alpha, beta, gamma, delta
            real(16)::f
         end function f
      end interface
      integer::j, nc, l
      real(16)::h, s0, xt, wt, t, as, shk, mba, pba, err, seps
      seps = c0*sqrt(eps)
      h = hr*0.5_16

      mba = (zb - za)*0.5_16
      pba = (zb + za)*0.5_16
      s0 = f(xc, yc, pba, zeta_f, theta_f, alpha, beta, gamma, delta)*h*pi2*mba

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
            s = s + f(xc, yc, mba*xt + pba, zeta_f, theta_f, alpha, beta, gamma, delta)*wt
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

   real(16):: eps = 1d-24, s = 0q0
   real(16):: alpha, beta, gamma, delta
   !eps = 1d-24; s = 0.0_16
   integer::info
   real(16) :: Cap_theta = 100.0_16
   real(16) :: Cap_zeta = 100.0_16
   !  integer :: display = 10                !if mod(n,display) = 0, then show the situation on the console
   !  integer :: proc                        !For cputime measurement

   !------------------------------------------------------------
   !@@@ Physical and numerical parameters ----------------------
   integer, parameter :: Kmax = 100   !Nz (maximum value of theta)
   integer, parameter :: Jmax = 100  !Nz (maximum value of zeta)
   integer, parameter :: Lmax = 100  !Nz (maximum value of zeta)
   integer, parameter :: Mmax = 100  !Nz (maximum value of zeta)
   !------------------------------------------------------------
   !Physical quantities----------------------------------------
   real(16), allocatable :: zeta(:), theta(:)
   !for integral index
   integer::i, j, k, l
   allocate (zeta(0:Jmax), theta(0:Kmax))
   do i = 0, Jmax
      zeta(i) = Cap_zeta*(i*1q0)/(Jmax*1q0)
   end do
   do k = 0, Kmax
      theta(k) = pi*(k*1q0)/(Kmax*1q0)
   end do
   print *, "zeta", zeta(1), "theta", theta(1)
   alpha = 0.0q0; beta = 2.0q0; gamma = 0.0q0; delta = pi/50q0
   !call dde3d(L1,zeta(1),theta(1),zeta(0),zeta(2),theta(0),theta(2), 0.0_16, 1.0_16, 0.0_16, pi/100.0_16, 0.0_16, 2*pi, eps, s, info)
   !!9.033969976609827287099592928924005E-0006
   call dde3d(L1, zeta(1), theta(1), alpha, beta, gamma, delta, 0.0_16, 1.0_16, 0.0_16, pi/100.0_16, 0.0_16, 2*pi, eps, s, info)
   print *, s! 9.596087731032261608563067133300942E-0010
   print *, L1(1.0q0, 1.0q0, 1.0q0, zeta(1), theta(1), alpha, beta, gamma, delta)! -8.022711679671201503181953271830736E-0005

   ! -1.961747754100960317717098015371201E-0009
contains
   function F1(zeta_bar, theta_bar, psi, zeta_f, theta_f)!計算できる
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_f, theta_f
      real(16)::F1

      F1 = zeta_bar**2 + zeta_f**2 - 2*zeta_bar*zeta_f*cos(theta_bar)*cos(theta_f) &
           - 2*zeta_bar*zeta_f*sin(theta_bar)*sin(theta_f)*cos(psi)
   end function
   function F2(zeta_bar, theta_bar, psi, zeta_f, theta_f)
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_f, theta_f
      real(16)::F2
      F2 = zeta_bar**2*zeta_f**2 &
           - zeta_bar**2*zeta_f**2 &
           *(cos(theta_f)*cos(theta_bar) + sin(theta_bar)*sin(theta_f)*cos(psi))**2
   end function
   function f_for_psi(zeta_bar, theta_bar, zeta_f, theta_f, alpha, beta, gamma, delta)
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, alpha, beta, gamma, delta, zeta_f, theta_f
      real(16) :: f_for_psi
      f_for_psi = (zeta_bar - alpha)*(zeta_bar - beta)*(theta_bar - gamma)*(theta_bar - delta) &
                  /(zeta_f - alpha)*(zeta_f - beta)*(theta_f - gamma)*(theta_f - delta)
   end function
   function L1(zeta_bar, theta_bar, psi, zeta_f, theta_f, alpha, beta, gamma, delta)!計算できない
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_f, theta_f, alpha, beta, gamma, delta
      real(16)::L1

      L1 = zeta_bar*zeta_bar*sin(theta_bar)/sqrt(2.0q0)/pi*cos(psi) &
           /sqrt(F1(zeta_bar, theta_bar, psi, zeta_f, theta_f)) &
           *exp(-zeta_bar**2 + F2(zeta_bar, theta_bar, psi, zeta_f, theta_f) &
                /F1(zeta_bar, theta_bar, psi, zeta_f, theta_f)) &
           *f_for_psi(zeta_bar, theta_bar, zeta_f, theta_f, alpha, beta, gamma, delta)
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
