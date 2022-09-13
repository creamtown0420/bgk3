    implicit none

    !***************************************************************************
    !*** Definition ************************************************************
    !***************************************************************************

    !------------------------------------------------------------
    !Dummy variables --------------------------------------------
    integer :: i, j, k
    !Auxiliary variables for general part -----------------------
    integer :: cputime(1:99, 1:2) = 0       !#1 -> proc, #2 -> 1 is start time, 2 is end time
    real(8) :: cputime_sum(1:99) = 0d0     !total cpu time used by proc
    integer :: CountPerSec, CountMax       !For cpu time measurement
    real(8) :: pi = 4d0*datan(1d0)        !pi
    real(8) :: spi = dsqrt(4d0*datan(1d0)) !sqrt(pi)
    real(8) :: large = 1d99                !some very large value
    real(8) :: small = 1d-16               !some very small value
    integer :: display = 10                !if mod(n,display) = 0, then show the situation on the console
    integer :: proc                        !For cputime measurement

    !------------------------------------------------------------
    !@@@ Physical and numerical parameters ----------------------
    real(8) :: Kn    !Knudsen number*sqrt(pi)/2
    integer :: Nr    !Nr (maximum value of r)
    integer :: Nthetaz    !Nz (maximum value of theta_zeta)
    integer :: Nzeta    !Nz (maximum value of zeta)
    integer :: LS_r  !lattice system for x
    integer :: LS_thetaz  !lattice system for theta_zeta
    integer :: LS_zeta    !lattice system for zeta
    !  integer :: FD_r  !1-> 1st order finite difference, 2-> second order finite difference (Switch for accuracy)
    !  integer :: FD_thetaz  !1-> 1st order finite difference, 2-> second order finite difference (Switch for accuracy)
    integer :: INT_zeta !1-> trapezoidal rule for integration, 2-> simpson rule (Switch for accuracy)
    integer :: INT_thetaz !1-> trapezoidal rule for integration, 2-> simpson rule (Switch for accuracy
    integer :: cap_K !capital Z for the maximum value of thetaz(k) ??
    integer :: cap_zeta !capital Z for the maximum value of zeta ??
    ! integer :: save1 !data save if mod(n,save1) = 0 ??

    !------------------------------------------------------------
    !Physical quantities----------------------------------------
    real(8), allocatable :: r(:), dr(:)            !(0:Nr)
    real(8), allocatable :: zeta(:), dzeta(:)            !(0:Nr)
    real(8), allocatable :: thetaz(:), dthetaz(:)            !(0:Nthetaz)
    real(8), allocatable :: g(:, :, :)                 !(0:Nr,-Nz:Nz,)
    real(8), allocatable :: uP(:)                  !(0:Nr)
    !  real(8)              :: MP, MP_old = 0d0       !flow rate
    !  real(8)              :: err, threshold = 1d-10 !if "err" < threshold, computation finishs
    !  real(8)              :: gamma1 = 1d0           !for BGK
    !  real(8)              :: k0     = -1.01619d0    !for BGK
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
       allocate (g(0:Nr, 0:Nzeta, 0:Nthetaz))
       allocate (uP(0:Nr))
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
          r(i) = 0.5d0*dble(i)/dble(Nr)
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
       do j = 0, Nzeta
          zeta(i) = 0.5d0*dble(i)/dble(Nzeta)
       end do
       do i = 1, Nr
          dzeta(i) = zeta(i) - zeta(i - 1)
       end do
    end subroutine
    !lattice for theta_zeta
    subroutine lattice_thetaz()
       implicit none
       do k = 0, Nthetaz
          zeta(i) = 0.5d0*dble(i)/dble(Nthetaz)
       end do
       do i = 1, Nr
          dthetaz(i) = thetaz(i) - thetaz(i - 1)
       end do
       do i = 0, Nr - 1
          Cap_lambda(i) = dthetaz(i + 1)/dthetaz(i)
       end do
    end subroutine
    !vdf -----------------------------------------------------
    subroutine COMPUTE_vdf()
       implicit none
       real(8) :: c1, c2
       real(8) :: a1, a2, a3 !?��̂ĕϐ�?��i?��m?��[?��g?��?��A,B,C)
       real(8) :: a, b, S0, S1, S2
       integer :: st, en, ip
!a3 =  (zeta(j)*cos(thetaz(k))/(dr(i-1)+dr(i))*(1/Cap_lambda(i-1)+2)-(zeta(j)*sin(thetaz(k))/(dthetaz(k-1)+dthetaz(k)*r(i)))*(1/lambda(k-1)+2)-zeta(j)*cos(thetaz(k))/r(i)-1/Kn)

       !?��_?��~?��[?��ϐ� Cap_lambda?��?��lambda?��̓O?��?��?��[?��o?��?��?��ŗp?��ӂ�?��?��
       !$omp parallel default(shared), private(i,j,c1,c2,c3)
       !$omp do
       do k = 1, Cap_K/2
          do j = dir*1, dir*Nz, dir
             do i = st, en, dir
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
       implicit none
       real(8) :: g0, g1, g2

       do i = 0, Nx
          uP(i) = 0d0
          do j = 0, Nzeta
             do k = 0, Nthetaz
                if (j == -1 .or. j == 0) cycle
                g0 = 1d0/spi*g(i, j + 0)*dexp(-z(j + 0)*z(j + 0))
                g1 = 1d0/spi*g(i, j + 1)*dexp(-z(j + 1)*z(j + 1))
                uP(i) = uP(i) + 0.25d0*dz(j + 1)*(g0 + g1)
             end do
          end do
       end do
    end subroutine
