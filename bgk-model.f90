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
   real(8) :: Kn    !Knudsen number*sqrt(pi)/2
   integer :: Nr = 100    !Nr (maximum value of r)
   integer :: Nthetaz = 100   !Nz (maximum value of theta_zeta)
   integer :: Nzeta = 100  !Nz (maximum value of zeta)

   !------------------------------------------------------------
   !Physical quantities----------------------------------------
   real(8), allocatable :: r(:), dr(:)            !(0:Nr)
   real(8), allocatable :: zeta(:), dzeta(:)            !(0:Nr)
   real(8), allocatable :: thetaz(:), dthetaz(:)            !(0:Nthetaz)
   real(8), allocatable :: Cap_lambda(:), lambda(:)            !(0:Nthetaz)

   real(8), allocatable :: g(:, :, :)                 !(0:Nr,-Nz:Nz,)
   real(8), allocatable :: uP(:)                  !(0:Nr)
   !  real(8)              :: MP, MP_old = 0d0       !flow rate
   !  real(8)              :: err, threshold = 1d-10 !if "err" < threshold, computation finishs
   !  real(8)              :: gamma1 = 1d0           !for BGK
   !  real(8)              :: k0     = -1.01619d0    !for BGK

!main processing
   call allocating()
   call initialize()
   call lattice_r()
   call lattice_zeta()
   call lattice_thetaz()

   call COMPUTE_bc(1)
   call COMPUTE_vdf()
   call COMPUTE_mac()
   call COMPUTE_bc(Nr)

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
      allocate (r(0:Nr), dr(1:Nr))
      allocate (zeta(0:Nzeta), dzeta(1:Nzeta))
      allocate (thetaz(0:Nthetaz), dthetaz(1:Nthetaz))
      allocate (g(0:Nr, 0:Nzeta, 0:Nthetaz))
      allocate (uP(0:Nr))
      allocate (Cap_lambda(0:Nr - 1), lambda(0:Nthetaz))
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
      do i = 0, Nr
         r(i) = Cap_r*dble(i)/dble(Nr)
      end do
      do i = 1, Nr
         dr(i) = r(i) - r(i - 1)
      end do
      do i = 0, Nr - 1
         Cap_lambda(i) = dr(i + 1)/dr(i)
      end do
   end subroutine
   !lattice for zeta
   subroutine lattice_zeta()
      implicit none
      do i = 0, Nzeta
         zeta(i) = Cap_zeta*dble(i)/dble(Nzeta)
      end do
      do i = 1, Nr
         dzeta(i) = zeta(i) - zeta(i - 1)
      end do
   end subroutine
   !lattice for theta_zeta
   subroutine lattice_thetaz()
      implicit none
      do k = 0, Nthetaz
         thetaz(i) = pi*dble(i)/dble(Nthetaz)
      end do
      do i = 1, Nthetaz
         dthetaz(i) = thetaz(i) - thetaz(i - 1)
      end do
      do i = 0, Nthetaz - 1
         lambda(i) = dthetaz(i + 1)/dthetaz(i)
      end do
   end subroutine
   !vdf -----------------------------------------------------
   subroutine COMPUTE_vdf()
      implicit none
      real(8) :: c1, c2
      real(8) :: a1, a2, a3 !?��̂ĕϐ�?��i?��m?��[?��g?��?��A,B,C)
!a3 =  (zeta(j)*cos(thetaz(k))/(dr(i-1)+dr(i))*(1/Cap_lambda(i-1)+2)-(zeta(j)*sin(thetaz(k))/(dthetaz(k-1)+dthetaz(k)*r(i)))*(1/lambda(k-1)+2)-zeta(j)*cos(thetaz(k))/r(i)-1/Kn)

      !?��_?��~?��[?��ϐ� Cap_lambda?��?��lambda?��̓O?��?��?��[?��o?��?��?��ŗp?��ӂ�?��?��
      !$omp parallel default(shared), private(i,j,c1,c2,c3)
      !$omp do
      do k = 1, Nthetaz/2
         do j = 1, Nzeta
            do i = 1, Nr
               a1 = (zeta(j)*cos(thetaz(k))/(dr(i - 1) + dr(i)))
               a2 = (zeta(j)*sin(thetaz(k))/(dthetaz(k - 1) + dthetaz(k)*r(i)))
               a3 = a1*(1/Cap_lambda(i - 1) + 2) - a2*(1/lambda(k - 1) + 2) - zeta(j)*cos(thetaz(k))/r(i) - 1/Kn
               c1 = (Cap_lambda(i - 1) + 1/Cap_lambda(i - 1) + 2)
               c2 = (lambda(k - 1) + 1/Cap_lambda(k - 1) + 2)
               g(i, j, k) = g(i - 1, j, k)*a1*c1/a3 + g(i - 2, j, k)*a1*Cap_lambda(i - 1)/a3 &
                            + g(i, j, k - 1)*a2*c2/a3 + g(i, j, k - 2)*a2*lambda(k - 1)/a3 - uP(i)/Kn
            end do
         end do
      end do
      !$omp end do
      !$omp end parallel
   end subroutine

!uPを計算する
   subroutine COMPUTE_mac()
      use integral_mod
      implicit none

      real(8) :: h(1:Nr, 1:Nthetaz)
      real(8) :: s
      ! call integral(配列の大きさ,配列,刻み幅,積分結果,積分手法)
      !  call integral(size(w,1),size(w,2),w,hx,hy,s,"simpson")
      do i = 1, Nr
         do j = 1, Nzeta
            do k = 1, Nthetaz/2
               h(j, k) = zeta(j)**4*exp(-zeta(j)**2)*sin(thetaz(k))**2*g(i, j, k)
            end do
         end do
         call dintegral2d(size(h, 1), size(h, 2), h, dzeta(1), dthetaz(1), s, "simpson")!motoha integral datta
         uP(i) = 1/spi*s
      end do
   end subroutine

   !BC -----------------------------------------------------
   subroutine COMPUTE_bc(i_input)
      implicit none
      integer, intent(in) :: i_input

      !BC at infinity
      if (i_input == Nr) then
         i = Nr
         do j = 1, Nzeta
            do K = 1, Nthetaz/2
               g(i, j, k) = 0d0
            end do
         end do
      end if
      !BC on the sphere
      if (i_input == 1) then
         i = 1
         do j = 1, Nzeta
            do K = 1, Nthetaz/2
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
      do i = 0, Nr
         write (22, '(a,e15.5,a)') 'zone t = "', r(i), '"'
         do j = 0, Nzeta
            if (j == 0) cycle
            do k = 0, Nthetaz
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
      open (unit=22, file="C:\workspace\bgk3\"//"_running"//"/mac/"//trim(adjustl(file_name)))
      do i = 0, Nr
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
