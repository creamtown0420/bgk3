subroutine dde3d(f, xa, xb, ya, yb, za, zb, eps, s, info)
   ! Integrate for 3d
   implicit none
   real(16), intent(in)::xa, xb, ya, yb, za, zb, eps
   real(16), intent(out)::s
   integer, intent(out)::info

   !============全部にこれいる
   integer, parameter::kmin = 3
   integer, parameter::kmax = 14
   ! max order of DE integration
   real(16), parameter::hr = 6.0_16    ! [-1,1]
   real(16), parameter::c0 = 0.01_16  ! safe value for eps
   !real(16), parameter::pi = acos(-1.0_16)
   real(16), parameter::pi = 3.141592653589793238462643383279502884_16
   real(16), parameter::pi2 = pi/2.0_16
   !============
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
   seps = c0*sqrt(eps)
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
      !$omp parallel do reduction(+:s) private(t,shk,xt,wt,ft)
      do j = 1, nc
         t = dble(2*j - nc - 1)*h
         shk = pi2*sinh(t)

         xt = tanh(shk)
         wt = pi2*cosh(t)/(cosh(shk)*cosh(shk))
         call dde3dsyz(f, ya, yb, za, zb, mba*xt + pba, eps, ft, info)
         s = s + ft*wt
      end do
      !$omp end parallel do
      s = s0*0.5_16 + s*h*mba
      print *, "s=", s

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
   !============全部にこれいる
   integer, parameter::kmin = 3
   integer, parameter::kmax = 14
   ! max order of DE integration
   real(16), parameter::hr = 6.0_16    ! [-1,1]
   real(16), parameter::c0 = 0.01_16  ! safe value for eps
   !real(16), parameter::pi = acos(-1.0_16)
   real(16), parameter::pi = 3.141592653589793238462643383279502884_16
   real(16), parameter::pi2 = pi/2.0_16
   !============
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
      !$omp parallel do reduction(+:s) private(t,shk,xt,wt,ft)
      do j = 1, nc
         t = dble(2*j - nc - 1)*h
         shk = pi2*sinh(t)

         xt = tanh(shk)
         wt = pi2*cosh(t)/(cosh(shk)*cosh(shk))
         call dde3dsz(f, za, zb, xc, mba*xt + pba, eps, ft, info)
         s = s + ft*wt
         !print *, "y no sekibun l=", l, "j=", j
      end do
      !$omp end parallel do
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
   !============全部にこれいる
   integer, parameter::kmin = 3
   integer, parameter::kmax = 14
   ! max order of DE integration
   real(16), parameter::hr = 6.0_16    ! [-1,1]
   real(16), parameter::c0 = 0.01_16  ! safe value for eps
   !real(16), parameter::pi = acos(-1.0_16)
   real(16), parameter::pi = 3.141592653589793238462643383279502884_16
   real(16), parameter::pi2 = pi/2.0_16
   !============
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
      !$omp parallel do reduction(+:s) private(t,shk,xt,wt)
      do j = 1, nc
         t = dble(2*j - nc - 1)*h
         shk = pi2*sinh(t)

         xt = tanh(shk)
         wt = pi2*cosh(t)/(cosh(shk)*cosh(shk))
         s = s + f(xc, yc, mba*xt + pba)*wt
         !print *, "z no sekibun l=", l, "j=", j
      end do
      !$omp end parallel do

      s = s0*0.5_16 + s*h*mba
      as = abs(s); err = abs(s - s0)
      if (as .ge. 1.0_16) err = err/as
      if (err .le. seps .and. l .ge. kmin) exit
      s0 = s
   end do

   if (l .eq. kmax + 1) info = 1

   return
end subroutine dde3dsz
