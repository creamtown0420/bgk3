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
