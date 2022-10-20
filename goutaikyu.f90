! call integral(配列の大きさ,配列,刻み幅,積分結果,積分手法)
! という感じです。
! 実際のプログラムでは、以下のように使用してください。配列yは複素数解列(ただし実軸上の積分)でもokです。
!https://slpr.sakura.ne.jp/qp/nc-integral/
! 1次元

! サブルーチンintegralの使い方
! 積分則        呼び方
! 短冊近似(長方形近似)
!  call integral(size(w,1),w,h,s,"box")
! 配列yの大きさに指定はない。
! 台形近似
!  call integral(size(y,1),y,h,s,"trapezoid")
! 配列yの大きさに指定はない。
! シンプソン則
!  call integral(size(y,1),y,h,s,"simpson")
! 配列yの大きさが3,5,7,…個、2n+1個でないといけない
! シンプソン3/8則
!  call integral(size(y,1),y,h,s,"simpson38")
! 配列yの大きさが4,7,10,…個、3n+1個でないといけない
! ブール則
!  call integral(size(y,1),y,h,s,"boole")
! 配列yの大きさが5,9,13,…個、4n+1個でないといけない
! 2次元配列の場合は以下のように指定してください。

!  call integral(size(w,1),size(w,2),w,hx,hy,s,"simpson")
! 3次元の場合はこう指定してください。

!  call integral(size(w,1),size(w,2),size(w,3),w,hx,hy,hz,s,"simpson")
module integral_mod
   !developer --> sikino
   !date --> 2015/04/07
   implicit none
   interface integral
      module procedure &
         dintegral, &
         dintegral2d, &
         dintegral3d, &
         cintegral, &
         cintegral2d, &
         cintegral3d
   end interface integral
contains
   subroutine dintegral(N, y, h, s, method)
      integer, intent(in)::N
      double precision, intent(in)::h, y(1:N)
      character(*), intent(in)::method
      double precision, intent(out)::s
      integer::i
      double precision::y0, y1, y2, y3

      s = 0.d0; y0 = 0.d0; y1 = 0.d0; y2 = 0.d0; y3 = 0.d0

      if (trim(method) .eq. "box") then
         s = h*sum(y(1:N - 1))
      elseif (trim(method) .eq. "trapezoid") then
         y1 = y(1) + y(N)
         do i = 2, N - 1
            y2 = y2 + y(i)
         end do
         s = (y1 + 2.d0*y2)*h*0.5d0
      elseif (trim(method) .eq. "simpson") then
         if (mod(N, 2) .ne. 1) then
            write (6, *) "=====cannot calculation with simpson"
            write (6, *) "=====program stop"
            stop
         end if

         y1 = y(1) + y(N)
         do i = 2, N - 1, 2
            y2 = y2 + y(i)
         end do
         do i = 3, N - 2, 2
            y3 = y3 + y(i)
         end do

         s = (y1 + 4.d0*y2 + 2.d0*y3)*h/3.d0

      elseif (trim(method) .eq. "simpson38") then

         if (mod(N, 3) .ne. 1) then
            write (6, *) "=====cannot calculation with simpson38"
            write (6, *) "=====program stop"
            stop
         end if

         y0 = y(1) + y(N)
         do i = 2, N - 2, 3
            y1 = y1 + y(i)
         end do
         do i = 3, N - 1, 3
            y2 = y2 + y(i)
         end do
         do i = 4, N - 3, 3
            y3 = y3 + y(i)
         end do
         s = (y0 + 3.d0*(y1 + y2) + 2.d0*y3)*3.d0*h/8.d0

      elseif (trim(method) .eq. "boole") then
         if (mod(N, 4) .ne. 1) then
            write (6, *) "=====cannot calculation with boole"
            write (6, *) "=====program stop"
            stop
         end if

         y0 = y(1) + y(N)
         do i = 5, N - 4, 4
            y1 = y1 + y(i)
         end do
         do i = 2, N - 1, 2
            y2 = y2 + y(i)
         end do
         do i = 3, N - 2, 4
            y3 = y3 + y(i)
         end do

         s = (14.d0*y0 + 28.d0*y1 + 64.d0*y2 + 24.d0*y3)*h/45.d0
      else
         write (6, *) "=====cannot calculation in integral"
         write (6, *) "=====program stop"
         stop
      end if

      return
   end subroutine dintegral

   subroutine dintegral2d(Nx, Ny, z, hx, hy, s, method)
      implicit none
      integer, intent(in)::Nx, Ny
      double precision, intent(in)::hx, hy, z(1:Nx, 1:Ny)
      character(*), intent(in)::method
      double precision, intent(out)::s
      integer::i
      double precision::ty(1:Ny), r(1:Nx)

      s = 0.d0
      ty(1:Ny) = 0.d0
      r(1:Nx) = 0.d0
      do i = 1, Nx
         ty(1:Ny) = z(i, 1:Ny)
         call integral(Ny, ty, hy, s, method)
         r(i) = s
      end do
      call integral(Nx, r, hx, s, method)

      return
   end subroutine dintegral2d

   subroutine dintegral3d(Nx, Ny, Nz, w, hx, hy, hz, s, method)
      implicit none
      integer, intent(in)::Nx, Ny, Nz
      double precision, intent(in)::hx, hy, hz, w(1:Nx, 1:Ny, 1:Nz)
      character(*), intent(in)::method
      double precision, intent(out)::s
      integer::i
      double precision::tyz(1:Ny, 1:Nz), r(1:Nx)

      s = 0.d0
      tyz(1:Ny, 1:Nz) = 0.d0
      r(1:Nx) = 0.d0
      do i = 1, Nx
         tyz(1:Ny, 1:Nz) = w(i, 1:Ny, 1:Nz)
         call integral(Ny, Nz, tyz, hy, hz, s, method)
         r(i) = s
      end do
      call integral(Nx, r, hx, s, method)

      return
   end subroutine dintegral3d

   subroutine cintegral(N, y, h, s, method)
      integer, intent(in)::N
      complex(kind(0d0)), intent(in)::y(1:N)
      double precision, intent(in)::h
      character(*), intent(in)::method
      complex(kind(0d0)), intent(out)::s

      double precision::res, ims

      s = dcmplx(0d0, 0d0); res = 0.d0; ims = 0.d0

      call integral(N, dble(y), h, res, trim(method))
      call integral(N, dimag(y), h, ims, trim(method))

      s = dcmplx(res, ims)

      return
   end subroutine cintegral

   subroutine cintegral2d(Nx, Ny, z, hx, hy, s, method)
      integer, intent(in)::Nx, Ny
      complex(kind(0d0)), intent(in)::z(1:Nx, 1:Ny)
      double precision, intent(in)::hx, hy
      character(*), intent(in)::method
      complex(kind(0d0)), intent(out)::s

      double precision::res, ims

      s = dcmplx(0d0, 0d0); res = 0.d0; ims = 0.d0

      call integral(Nx, Ny, dble(z), hx, hy, res, trim(method))
      call integral(Nx, Ny, dimag(z), hx, hy, ims, trim(method))

      s = dcmplx(res, ims)

      return
   end subroutine cintegral2d

   subroutine cintegral3d(Nx, Ny, Nz, w, hx, hy, hz, s, method)
      integer, intent(in)::Nx, Ny, Nz
      complex(kind(0d0)), intent(in)::w(1:Nx, 1:Ny, 1:Nz)
      double precision, intent(in)::hx, hy, hz
      character(*), intent(in)::method
      complex(kind(0d0)), intent(out)::s

      double precision::res, ims

      s = dcmplx(0d0, 0d0); res = 0.d0; ims = 0.d0

      call integral(Nx, Ny, Nz, dble(w), hx, hy, hz, res, trim(method))
      call integral(Nx, Ny, Nz, dimag(w), hx, hy, hz, ims, trim(method))

      s = dcmplx(res, ims)

      return
   end subroutine cintegral3d
end module integral_mod
!**********************************************************
!**********************************************************
!**********************************************************
!**********************************************************
program main
   use integral_mod
   implicit none

   !***************************************************************************
   !*** Definition ************************************************************
   !***************************************************************************

   !------------------------------------------------------------
   !Dummy variables --------------------------------------------
   integer :: i, j, k
   !Auxiliary variables for general part -----------------------
   ! integer :: cputime(1:99, 1:2) = 0       !#1 -> proc, #2 -> 1 is start time, 2 is end time
   !  real(8) :: cputime_sum(1:99) = 0d0     !total cpu time used by proc
   !  integer :: CountPerSec, CountMax       !For cpu time measurement
   real(8) :: pi = 4d0*datan(1d0)        !pi
   real(8) :: spi = dsqrt(4d0*datan(1d0)) !sqrt(pi)
   real(8) :: large = 1d99                !some very large value
   !  real(8) :: small = 1d-16               !some very small value
   real(8) :: Cap_r = 1d2
   real(8) :: Cap_zeta = 1d2
   !  integer :: display = 10                !if mod(n,display) = 0, then show the situation on the console
   !  integer :: proc                        !For cputime measurement

   !------------------------------------------------------------
   !@@@ Physical and numerical parameters ----------------------
   real(8) :: Kn = 1  !仮決め   !Knudsen number*sqrt(pi)/2
   integer :: Imax = 100    !Imax (maximum value of r)
   integer, parameter :: Kmax = 100   !Nz (maximum value of theta_zeta)
   integer, parameter :: Jmax = 100  !Nz (maximum value of zeta)
   integer, parameter :: Khat = Kmax/2 ! ほんとはKmax/2で定義したい
   !------------------------------------------------------------
   !Physical quantities----------------------------------------
   real(8), allocatable :: r(:), dr(:)            !(0:Imax)
   real(8), allocatable :: zeta(:), zeta_bar(:), dzeta(:)            !(0:Imax)
   real(8), allocatable :: thete_zeta(:), theta_bar(:) dthete_zeta(:)            !(0:Kmax)
   real(8), allocatable :: Cap_lambda(:), lambda(:)            !(0:Kmax)
   real(8), allocatable :: r_i(:), theta_zeta_k(:)
   real(8), allocatable :: g(:, :, :)                 !(0:Imax,-Nz:Nz,)
   real(8), allocatable :: uP(:)                  !(0:Imax)
   !  real(8)              :: MP, MP_old = 0d0       !flow rate
   !  real(8)              :: err, threshold = 1d-10 !if "err" < threshold, computation finishs
   !  real(8)              :: gamma1 = 1d0           !for BGK
   !  real(8)              :: k0     = -1.01619d0    !for BGK

!main processing

   call allocating()
   call initialize()
   call lattice_r()
   call lattice_zeta()
   call lattice_thete_zeta()
   do i = 1, 100 !100では結果がすぐ出るのに101ではでない
      call COMPUTE_bc(1)
      call COMPUTE_vdf()
      call COMPUTE_mac()
      call COMPUTE_bc(Imax)
   end do
   call OUTPUT_mac()
contains
   !===========================================================
   !===========================================================
   !===========================================================
   !=== Set up ================================================
   !===========================================================
   !===========================================================
   !===========================================================

   !Allocate variables
   subroutine allocating()
      implicit none
      allocate (r(0:Imax), dr(1:Imax))
      allocate (zeta(0:Jmax), zeta_bar(0:Jmax), dzeta(1:Jmax))
      allocate (thete_zeta(0:Kmax), theta_zeta_bar(0:Kmax), dthete_zeta(1:Kmax))
      allocate (g(0:Imax, 0:Jmax, 0:Kmax))
      allocate (uP(0:Imax))
      allocate (Cap_lambda(1:Imax - 1), lambda(1:Kmax))
      allocate (r_i(0:Khat), theta_zeta_k(0:Imax))!Khat:Kmaxじゃないの？
   end subroutine

   !Set initial setting
   subroutine initialize()
      implicit none
      g = large
      uP = 0d0
   end subroutine

   !===========================================================
   !===========================================================
   !===========================================================
   !=== Lattice system=========================================
   !===========================================================
   !===========================================================
   !===========================================================

   !lattice for r
   subroutine lattice_r()
      implicit none
      do i = 0, Imax
         r(i) = Cap_r*dble(i)/dble(Imax)
      end do
      do i = 1, Imax
         dr(i) = r(i) - r(i - 1)
      end do
      do i = 1, Imax - 1
         Cap_lambda(i) = dr(i + 1)/dr(i)
      end do
   end subroutine
   !lattice for zeta
   subroutine lattice_zeta()
      implicit none
      do i = 0, Jmax
         zeta(i) = Cap_zeta*dble(i)/dble(Jmax)
      end do
      do i = 1, Jmax
         dzeta(i) = zeta(i) - zeta(i - 1)
      end do
   end subroutine
   !lattice for theta_zeta
   subroutine lattice_thete_zeta()
      implicit none
      do k = 0, Kmax
         thete_zeta(k) = pi*dble(k)/dble(Kmax)
      end do
      do i = 1, Kmax
         dthete_zeta(i) = thete_zeta(i) - thete_zeta(i - 1) !問題あり
      end do
      do i = 1, Kmax - 1
         lambda(i) = dthete_zeta(i + 1)/dthete_zeta(i)
      end do
   end subroutine
   !*******************不連続面のr_iおよびtheta_zeta_k
   subroutine r_i()
      r_i(0) = 10d0**99
      !$omp parallel do
      do k = 1, Khat
         r_i(k) = 1d0/dsin(theta_zeta(k))
      end do
      !$omp end parallel do
   end subroutine
   !$omp parallel do
   subroutine theta_zeta_k()
      do i = 0, Imax
         theta_zeta_k(i) = dasin(1d0/r(i))
      end do
   end subroutine
   !theta=piを先に求める
   !Cをあらかじめ求める--------------------------------
   !衝突積分Cを求める---------------------------------
   subroutine compute_K_jklm()!{L1(psi)-L2(psi)のzeta(l),theta_zeta(m)における値}
      implicit none
real(8)::x_gauss(4) = (/-0.861136311594052575224, -0.3399810435848562648027, 0.3399810435848562648027, 0.861136311594052575224/) !,zeta_gauss(4),psi_gauss(4)
    real(8)::w_gauss(4) = (/0.3478548451374538573731, 0.6521451548625461426269, 0.6521451548625461426269, 0.3478548451374538573731/)
      real(8)::zeta_bar_real, theta__zeta_bar_real
      do i = 1, 4 !L1の計算
      do j = 1, 4
      do k = 1, 4

         L1 = L1 + w_gauss(i)*w_gauss(j)*w_gauss(k)*
      end do
      end do
      end do
!関数の定義
!zeta_bar(l=奇数),theta_zeta_bar(m=奇数)
      real(8) function l_kisu_m_kisu(l, m, zeta_bar_real, theta_zeta_bar_real)
         integer l, m
    (zeta_bar_real-zeta_bar(l-1))*(zeta_bar_real-zeta_bar(l+1))*(theta_zeta_bar_real-theta_bar(m-1))*(theta_zeta_bar_real-theta_bar(m+1))/
   $( (zeta_bar(l)-zeta_bar(l-1))*(zeta_bar(l)-zeta_bar(l+1))*(theta_zeta_bar(m)-theta_bar(m-1))*(theta_zeta_bar(m)-theta_bar(m+1)))
      end function

      real(8) function l_kisu_m_gusu(l, m, zeta_bar_real, theta_zeta_bar_real)
         integer l, m
    (zeta_bar_real-zeta_bar(l-1))*(zeta_bar_real-zeta_bar(l+1))*(theta_zeta_bar_real-theta_bar(m-1))*(theta_zeta_bar_real-theta_bar(m+1))/
   $( (zeta_bar(l)-zeta_bar(l-1))*(zeta_bar(l)-zeta_bar(l+1))*(theta_zeta_bar(m)-theta_bar(m-1))*(theta_zeta_bar(m)-theta_bar(m+1)))
      end function

      !vdf -----------------------------------------------------
      subroutine COMPUTE_vdf()
         implicit none
         real(8) :: c1, c2
         real(8) :: a1, a2, a3 !?��̂ĕϐ�?��i?��m?��[?��g?��?��A,B,C)
!a3 =  (zeta(j)*cos(thete_zeta(k))/(dr(i-1)+dr(i))*(1/Cap_lambda(i-1)+2)-(zeta(j)*sin(thete_zeta(k))/(dthete_zeta(k-1)+dthete_zeta(k)*r(i)))*(1/lambda(k-1)+2)-zeta(j)*cos(thete_zeta(k))/r(i)-1/Kn)

         !?��_?��~?��[?��ϐ� Cap_lambda?��?��lambda?��̓O?��?��?��[?��o?��?��?��ŗp?��ӂ�?��?��
         !$omp parallel default(shared), private(i,j,c1,c2,c3)
         !$omp do
         if (houkou=1)
         do k = 2, Khat, houkou
            do j = 2, Jmax, houkou
               do i = 2, Imax, houkou
                  a1 = (zeta(j)*cos(thete_zeta(k))/(dr(i - 1) + dr(i)))
                  a2 = (zeta(j)*sin(thete_zeta(k))/((dthete_zeta(k - 1) + dthete_zeta(k))*r(i)))
                  a3 = a1*(1/Cap_lambda(i - 1) + 2) - a2*(1/lambda(k - 1) + 2) - zeta(j)*cos(thete_zeta(k))/r(i) - 1/Kn
                  ! print*,thete_zeta(k)
                  c1 = (Cap_lambda(i - 1) + 1/Cap_lambda(i - 1) + 2)
                  c2 = (lambda(k - 1) + 1/Cap_lambda(k - 1) + 2)
                  g(i, j, k) = g(i - 1, j, k)*a1*c1/a3 + g(i - 2, j, k)*a1*Cap_lambda(i - 1)/a3 &
                               + g(i, j, k - 1)*a2*c2/a3 + g(i, j, k - 2)*a2*lambda(k - 1)/a3 - uP(i)/Kn
               end do
            end do
         end do
         end if
         elseif (houkou=-1)
         do k = 2, Khat, houkou
            do j = 2, Jmax, houkou
               do i = 2, Imax, houkou
                  a1 = (zeta(j)*cos(thete_zeta(k))/(dr(i - 1) + dr(i)))
                  a2 = (zeta(j)*sin(thete_zeta(k))/((dthete_zeta(k - 1) + dthete_zeta(k))*r(i)))
                  a3 = a1*(1/Cap_lambda(i - 1) + 2) - a2*(1/lambda(k - 1) + 2) - zeta(j)*cos(thete_zeta(k))/r(i) - 1/Kn
                  ! print*,thete_zeta(k)
                  c1 = (Cap_lambda(i - 1) + 1/Cap_lambda(i - 1) + 2)
                  c2 = (lambda(k - 1) + 1/Cap_lambda(k - 1) + 2)
                  g(i, j, k) = g(i - 1, j, k)*a1*c1/a3 + g(i - 2, j, k)*a1*Cap_lambda(i - 1)/a3 &
                               + g(i, j, k - 1)*a2*c2/a3 + g(i, j, k - 2)*a2*lambda(k - 1)/a3 - uP(i)/Kn
               end do
            end do
         end do
         end if
         !$omp end do
         !$omp end parallel

         end subroutine

!uPを計算する
         subroutine COMPUTE_mac()
            use integral_mod
            implicit none

            real(8) :: h(1:Imax, 1:Khat)
            real(8) :: s
            ! call integral(配列の大きさ,配列,刻み幅,積分結果,積分手法)
            !  call integral(size(w,1),size(w,2),w,hx,hy,s,"simpson")
            do i = 1, Imax
               do j = 1, Jmax
                  do k = 1, Khat
                     h(j, k) = zeta(j)**4*exp(-zeta(j)**2)*sin(thete_zeta(k))**2*g(i, j, k)
                  end do
               end do
               call dintegral2d(size(h, 1), size(h, 2), h, dzeta(1), dthete_zeta(1), s, "trapezoid")!motoha integral datta
               uP(i) = 1/spi*s
            end do
         end subroutine

         !BC -----------------------------------------------------
         subroutine COMPUTE_bc(i_input)
            implicit none
            integer, intent(in) :: i_input

            !BC at infinity
            if (i_input == Imax) then
               i = Imax
               do j = 1, Jmax
                  do K = 1, Khat
                     g(i, j, k) = 0d0
                  end do
               end do
            end if
            !BC on the sphere
            if (i_input == 1) then
               i = 1
               do j = 1, Jmax
                  do K = 1, Khat
                     g(i, j, k) = 2
                  end do
               end do
            end if
         end subroutine
         !vdf------------------------------------
         subroutine OUTPUT_vdf()
            ! real(8) :: arg, anal

            ! if( mod(n,save1) == 0 ) then

            !If uP = 0, g can be obtained analytically.
            !   file_name = "vdf_check.dat"
            !   open( unit = 22, file = dir_name//"_running"//"/vdf/"//trim(adjustl(file_name)))
            !     j = Nz/2 !choose any j
            !     do i = 0, Nx
            !       arg  = -( x(i) - 0.5d0 )/z(j)/Kn
            !       if( j < 0 ) then
            !         anal = Kn*( dexp(arg) - 1d0 )
            !       else
            !         anal = Kn*( dexp(arg - 1d0/z(j)/Kn ) - 1d0 )
            !       end if
            !       write(22,*) x(i), g(i,j), anal
            !     end do
            !   close(22)

            implicit none
            character :: file_name*22
            file_name = 'vdf.dat'
            open (unit=22, file="C:\workspace\bgk3\"//"_running"//"/vdf/"//trim(adjustl(file_name)))
            do i = 0, Imax
               write (22, '(a,e15.5,a)') 'zone t = "', r(i), '"'
               do j = 0, Jmax
                  if (j == 0) cycle
                  do k = 0, Kmax
                     write (22, *) zeta(j), g(i, j, k)
                  end do
               end do
            end do
            close (22)

            ! end if

         end subroutine

         !mac------------------------------------
         subroutine OUTPUT_mac()
            implicit none
            character :: file_name*22
            file_name = "mac.dat"
            open (unit=22, file="C:/workspace/bgk3/"//"_running"//"/mac/"//trim(adjustl(file_name)))
            do i = 1, Imax
               write (22, *) r(i), uP(i)
            end do
            close (22)

            !end if

            !   file_name = "MP.dat"
            !   open( unit = 22, file = dir_name//"_running"//"/mac/"//trim(adjustl(file_name)), position = "append" )
            !     write(22,*) n, MP
            !   close(22)

         end subroutine

      end program main
