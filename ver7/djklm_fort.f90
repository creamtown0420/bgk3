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
         !$omp parallel do reduction(+:s) private(t,shk,xt,wt,ft)
         do j = 1, nc
            t = dble(2*j - nc - 1)*h
            shk = pi2*sinh(t)

            xt = tanh(shk)
            wt = pi2*cosh(t)/(cosh(shk)*cosh(shk))
          call dde3dsyz(f, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, ya, yb, za, zb, mba*xt + pba, eps, ft, info)
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
         !$omp parallel do reduction(+:s) private(t,shk,xt,wt,ft)
         do j = 1, nc
            t = dble(2*j - nc - 1)*h
            shk = pi2*sinh(t)

            xt = tanh(shk)
            wt = pi2*cosh(t)/(cosh(shk)*cosh(shk))
            call dde3dsz(f, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, za, zb, xc, mba*xt + pba, eps, ft, info)
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
         !$omp parallel do reduction(+:s) private(t,shk,xt,wt)
         do j = 1, nc
            t = dble(2*j - nc - 1)*h
            shk = pi2*sinh(t)

            xt = tanh(shk)
            wt = pi2*cosh(t)/(cosh(shk)*cosh(shk))
            s = s + f(xc, yc, mba*xt + pba, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)*wt
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

end module DE_integration

module gauss_int
   implicit none
   real(16), private, parameter::pi = 3.141592653589793238462643383279502884_16
   real(16), private, parameter::pi2 = pi/2.0_16
   real(16), private::xt(100)
   data xt/-0.999713726773441233678228469342q0, &
      -0.998491950639595818400163359186q0, &
      -0.996295134733125149186131732241q0, &
      -0.993124937037443459652009892849q0, &
      -0.98898439524299174800441874581q0, &
      -0.983877540706057015496100155511q0, &
      -0.977809358486918288553781088429q0, &
      -0.970785775763706331930897857898q0, &
      -0.96281365425581552729365932603q0, &
      -0.953900782925491742849336930894q0, &
      -0.944055870136255977962774706415q0, &
      -0.933288535043079545924333668131q0, &
      -0.921609298145333952666951328482q0, &
      -0.909029570982529690467126337789q0, &
      -0.89556164497072698669852102243q0, &
      -0.881218679385018415573316825428q0, &
      -0.866014688497164623410739969676q0, &
      -0.84996452787959128429336259142q0, &
      -0.833083879888400823542915833845q0, &
      -0.815389238339176254393988758649q0, &
      -0.79689789239031447638957288218q0, &
      -0.777627909649495475627551386835q0, &
      -0.757598118519707176035667964438q0, &
      -0.73682808980202070551242771482q0, &
      -0.7153381175730564464599671227q0, &
      -0.69314919935580196594864794168q0, &
      -0.670283015603141015802587014323q0, &
      -0.646761908514129279832630304459q0, &
      -0.622608860203707771604190845172q0, &
      -0.597847470247178721264806545149q0, &
      -0.572501932621381191316870443526q0, &
      -0.546597012065094167467994257182q0, &
      -0.520158019881763056646815749455q0, &
      -0.493210789208190933569308793449q0, &
      -0.465781649773358042249216623396q0, &
      -0.437897402172031513108978043622q0, &
      -0.40958529167830154252886840006q0, &
      -0.38087298162462995676336254887q0, &
      -0.35178852637242172097234382955q0, &
      -0.322360343900529151722476582398q0, &
      -0.292617188038471964737555888235q0, &
      -0.262588120371503479168929336255q0, &
      -0.232302481844973969649509963208q0, &
      -0.20178986409573599723604885953q0, &
      -0.171080080538603274887532374707q0, &
      -0.14020313723611397320751460468q0, &
      -0.10918920358006111500342600658q0, &
      -0.078068582813436636694817371202q0, &
      -0.046871682421591631614923912934q0, &
      -0.015628984421543082872216699997q0, &
      0.015628984421543082872216699997q0, &
      0.046871682421591631614923912934q0, &
      0.078068582813436636694817371202q0, &
      0.109189203580061115003426006579q0, &
      0.140203137236113973207514604682q0, &
      0.171080080538603274887532374707q0, &
      0.20178986409573599723604885953q0, &
      0.232302481844973969649509963208q0, &
      0.262588120371503479168929336255q0, &
      0.292617188038471964737555888235q0, &
      0.3223603439005291517224765824q0, &
      0.351788526372421720972343829549q0, &
      0.38087298162462995676336254887q0, &
      0.40958529167830154252886840006q0, &
      0.43789740217203151310897804362q0, &
      0.465781649773358042249216623396q0, &
      0.49321078920819093356930879345q0, &
      0.520158019881763056646815749455q0, &
      0.54659701206509416746799425718q0, &
      0.572501932621381191316870443526q0, &
      0.597847470247178721264806545149q0, &
      0.622608860203707771604190845172q0, &
      0.646761908514129279832630304459q0, &
      0.670283015603141015802587014323q0, &
      0.693149199355801965948647941675q0, &
      0.715338117573056446459967122704q0, &
      0.73682808980202070551242771482q0, &
      0.757598118519707176035667964438q0, &
      0.777627909649495475627551386835q0, &
      0.796897892390314476389572882183q0, &
      0.815389238339176254393988758649q0, &
      0.833083879888400823542915833845q0, &
      0.84996452787959128429336259142q0, &
      0.866014688497164623410739969676q0, &
      0.881218679385018415573316825428q0, &
      0.89556164497072698669852102243q0, &
      0.909029570982529690467126337789q0, &
      0.921609298145333952666951328482q0, &
      0.93328853504307954592433366813q0, &
      0.944055870136255977962774706415q0, &
      0.953900782925491742849336930894q0, &
      0.96281365425581552729365932603q0, &
      0.970785775763706331930897857898q0, &
      0.977809358486918288553781088429q0, &
      0.98387754070605701549610015551q0, &
      0.988984395242991748004418745808q0, &
      0.99312493703744345965200989285q0, &
      0.996295134733125149186131732241q0, &
      0.99849195063959581840016335919q0, &
      0.99971372677344123367822846934q0/
   real(16), private::wt(100)
   data wt/0.000734634490505671730406320658q0, &
      0.001709392653518105239529358371q0, &
      0.002683925371553482419439590429q0, &
      0.003655961201326375182342458728q0, &
      0.004624450063422119351095789083q0, &
      0.005588428003865515157211946348q0, &
      0.00654694845084532276415210333q0, &
      0.007499073255464711578828744016q0, &
      0.0084438714696689714026208349q0, &
      0.009380419653694457951418237661q0, &
      0.010307802574868969585782101728q0, &
      0.011225114023185977117221573366q0, &
      0.012131457662979497407744792449q0, &
      0.01302594789297154228555858376q0, &
      0.013907710703718772687954149108q0, &
      0.01477588452744130176887998752q0, &
      0.015629621077546002723936865954q0, &
      0.016468086176145212643104980088q0, &
      0.017290460568323582439344198367q0, &
      0.01809594072212811666439075142q0, &
      0.018883739613374904552941165882q0, &
      0.019653087494435305865381470245q0, &
      0.020403232646209432766838851658q0, &
      0.021133442112527641542672300441q0, &
      0.0218430024162473863139537413q0, &
      0.022531220256336272701796970932q0, &
      0.023197423185254121622488854183q0, &
      0.023840960265968205962560411902q0, &
      0.02446120270795705271997502335q0, &
      0.025057544481579589703764225621q0, &
      0.0256294029102081160756420098622q0, &
      0.026176219239545676342308741757q0, &
      0.026697459183570962660384664186q0, &
      0.027192613446576880136491567802q0, &
      0.02766119822079238829420415587q0, &
      0.028102755659101173317648330187q0, &
      0.028516854322395097990936762864q0, &
      0.028903089601125203134876228135q0, &
      0.029261084110638276620119023496q0, &
      0.029590488059912642511754510679q0, &
      0.029890979593332830916836806669q0, &
      0.03016226510516914491906868161q0, &
      0.030404079526454820016507859819q0, &
      0.03061618658398044849645944326q0, &
      0.03079837903115259042771390303q0, &
      0.03095047885049098823406346347q0, &
      0.03107233742756651658781017024q0, &
      0.0311638356962099067838183212172q0, &
      0.031224884254849357732376498648q0, &
      0.03125542345386335694764247439q0, &
      0.03125542345386335694764247439q0, &
      0.031224884254849357732376498648q0, &
      0.031163835696209906783818321217q0, &
      0.031072337427566516587810170243q0, &
      0.030950478850490988234063463471q0, &
      0.030798379031152590427713903031q0, &
      0.030616186583980448496459443262q0, &
      0.03040407952645482001650785982q0, &
      0.03016226510516914491906868161q0, &
      0.029890979593332830916836806669q0, &
      0.029590488059912642511754510679q0, &
      0.029261084110638276620119023496q0, &
      0.02890308960112520313487622813q0, &
      0.028516854322395097990936762864q0, &
      0.028102755659101173317648330187q0, &
      0.02766119822079238829420415587q0, &
      0.0271926134465768801364915678q0, &
      0.02669745918357096266038466419q0, &
      0.026176219239545676342308741757q0, &
      0.02562940291020811607564200986q0, &
      0.02505754448157958970376422562q0, &
      0.02446120270795705271997502335q0, &
      0.0238409602659682059625604119023q0, &
      0.02319742318525412162248885418q0, &
      0.022531220256336272701796970932q0, &
      0.0218430024162473863139537413044q0, &
      0.021133442112527641542672300441q0, &
      0.02040323264620943276683885166q0, &
      0.01965308749443530586538147025q0, &
      0.01888373961337490455294116588q0, &
      0.01809594072212811666439075142q0, &
      0.017290460568323582439344198367q0, &
      0.01646808617614521264310498009q0, &
      0.01562962107754600272393686595q0, &
      0.01477588452744130176887998752q0, &
      0.01390771070371877268795414911q0, &
      0.013025947892971542285558583759q0, &
      0.01213145766297949740774479245q0, &
      0.011225114023185977117221573366q0, &
      0.010307802574868969585782101728q0, &
      0.00938041965369445795141823766081q0, &
      0.008443871469668971402620834902q0, &
      0.007499073255464711578828744016q0, &
      0.00654694845084532276415210333q0, &
      0.0055884280038655151572119463484q0, &
      0.00462445006342211935109578908q0, &
      0.00365596120132637518234245873q0, &
      0.00268392537155348241943959043q0, &
      0.001709392653518105239529358371q0, &
      0.00073463449050567173040632066q0/
contains
   subroutine gauss(f, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, xa, xb, ya, yb, za, zb, eps, s, info)
      ! Integrate for 3d
      implicit none
      real(16), intent(in)::xa, xb, ya, yb, za, zb, eps, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta
      real(16), intent(out)::s
      integer, intent(in):: info
      integer::i, j, k
      interface
         function f(x, y, z, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)
            implicit none
            real(16), intent(in)::x, y, z, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta
            real(16)::f
         end function f
      end interface

      real(16)::xmba, xpba, ymba, ypba, zmba, zpba

      xmba = (xb - xa)*0.5_16
      xpba = (xb + xa)*0.5_16
      ymba = (yb - ya)*0.5_16
      ypba = (yb + ya)*0.5_16
      zmba = (zb - za)*0.5_16
      zpba = (zb + za)*0.5_16
      do i = 0, 99
         do j = 0, 99
            do k = 0, 99
               s = s + xmba*ymba*zmba*wt(i)*wt(j)*wt(k) &
                   *f(xmba*xt(i) + xpba, ymba*xt(j) + ypba, zmba*xt(k) + zpba, &
                      zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)
            end do
         end do
      end do
      print *, "s=", s
   end subroutine gauss
end module
program main
   use DE_integration
   use gauss_int
   implicit none
   real(16), parameter::pi2 = 1.5707963267948966_16
   !real(16), parameter::pi = acos(-1.0q0)
   real(16), parameter::pi = 3.141592653589793238462643383279502884_16
   !real(16), parameter::pi2 = pi/2.0q0

   real(16):: eps = 1d-24, s = 0q0!1d-24
   real(16):: alpha, beta, gamma, delta
   real(16):: zeta_l, theta_m, zeta_j, theta_k
   real(16):: val_kiki_1, val_kiki_2
   real(16)::val_gugu_pp_1, val_gugu_pp_2, val_gugu_pm_1, val_gugu_pm_2, val_gugu_mp_1, val_gugu_mp_2, val_gugu_mm_1, val_gugu_mm_2
   real(16)::val_kigu_0p_1, val_kigu_0p_2, val_kigu_0m_1, val_kigu_0m_2, val_guki_p0_1, val_guki_p0_2, val_guki_m0_1, val_guki_m0_2
   real(16):: val1, val2, val3, val4!特異点を避ける積分範囲の値

   !eps = 1d-24; s = 0.0_16
   !CPU time
   integer :: cputime(1:99, 1:2) = 0       !#1 -> proc, #2 -> 1 is start time, 2 is end time
   real(8) :: cputime_sum(1:99) = 0d0     !total cpu time used by proc
   integer :: CountPerSec, CountMax       !For cpu time measurement
   integer :: proc                        !For cputime measurement
   !--------
   integer::info
   real(16) :: Cap_theta = 100.0_16
   real(16) :: Cap_zeta = 100.0_16
   !  integer :: display = 10                !if mod(n,display) = 0, then show the situation on the console
   !  integer :: proc                        !For cputime measurement

   !------------------------------------------------------------
   !@@@ Physical and numerical parameters ----------------------
   integer, parameter :: Kmax = 64  !Nz (maximum value of theta)
   integer, parameter :: Jmax = 64!Nz (maximum value of zeta)
   integer, parameter :: Lmax = 64 !Nz (maximum value of zeta)
   integer, parameter :: Mmax = 64 !Nz (maximum value of zeta)
   !------------------------------------------------------------
   !Physical quantities----------------------------------------
   real(16), allocatable :: zeta(:), theta(:), K_jklm(:, :, :, :)
   !for integral index
   integer::i, j, k, l, m, ios
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
   !並列化がうまくいってない
   !$omp parallel do private(j,k,l,m,s,zeta_l, theta_m, zeta_j, theta_k,val_kiki_1, val_kiki_2&
   , val_gugu_pp_1, val_gugu_pp_2, val_gugu_pm_1, val_gugu_pm_2, val_gugu_mp_1, val_gugu_mp_2, val_gugu_mm_1, val_gugu_mm_2 &
      , val_kigu_0p_1, val_kigu_0p_2, val_kigu_0m_1, val_kigu_0m_2, val_guki_p0_1, val_guki_p0_2, val_guki_m0_1, val_guki_m0_2 &
      , val1, val2, val3, val4, a_l, b_l, c_m, d_m, a_l_p, a_l_m, b_l_p, b_l_m, c_m_p, c_m_m, d_m_p, d_m_m)
   do j = 1, 1!j=0では常にL1,L2は0 1.163269376341041095875714758599207E-0008
      do k = 1, 1!35 3 27 13
         do l = 1, 1! j50 k3 l8 m32 !j55 k3 l4 m22
            do m = 1, 1
               print *, "j", j, "k", k, "l", l, "m", m
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

                           !val1 = 0.0q0
                           val2 = 0.0q0
                           val3 = 0.0q0
                           val4 = 0.0q0
                           ! s = 0.0q0
                           ! call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                           !            zeta_j, zeta(b_l), theta(c_m), theta_k, 0.0q0, 2.0q0*pi, eps, s, info)!右下
                           ! val2 = s
                           ! s = 0.0q0
                           ! call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                           !            zeta(a_l), zeta_j, theta_k, theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)!左上
                           ! val3 = s
                           ! s = 0.0q0
                           ! call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                           !            zeta_j, zeta(b_l), theta_k, theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)!右上
                           ! val4 = s
                           ! s = 0.0q0
                           val_kiki_1 = val1 + val2 + val3 + val4

                        else
                           print *, "no singurality"
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l), zeta(b_l), theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_kiki_1 = s

                        end if
                        ! call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, zeta(a_l), zeta(b_l), &
                        !            theta(c_m), theta(d_m), 0.0q0, 2.0q0*pi, eps, s, info)
                        ! val_kiki_2 = s
                        ! s = 0.0q0
                     else
                        val_kiki_1 = 0.0q0
                        val_kiki_2 = 0.0q0
                     end if
                     K_jklm(j, k, l, m) = val_kiki_1 + val_kiki_2

                     !----------------偶数偶数---------------Nanがある(j=45 k=1 l=0 m=8)ほかにもいろいろある 原因不明　おそらく初期化不良
                  else if (mod(l, 2) .eq. 0 .and. mod(m, 2) .eq. 0) then!偶数偶数
                     a_l_p = l; b_l_p = l + 2; c_m_p = m; d_m_p = m + 2; ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です
                     a_l_m = l - 2; b_l_m = l; c_m_m = m - 2; d_m_m = m  ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です

                     !-----------------------pp---------------------l_p,m_pになる要チェック
                     if (0 .le. a_l_p .and. a_l_p .le. Lmax .and. 0 .le. b_l_p .and. b_l_p .le. Lmax .and. &!pp
                         0 .le. c_m_p .and. c_m_p .le. Mmax .and. 0 .le. d_m_p .and. d_m_p .le. Mmax) then
                        alpha = zeta(l + 1); beta = zeta(l + 2); gamma = theta(m + 1); delta = theta(m + 2)

                        if (zeta(a_l_p) .le. zeta_j .and. zeta_j .le. zeta(b_l_p) &
                            .and. theta(c_m_p) .le. theta_k .and. theta_k .le. theta(d_m_p)) then!特異点があるとき
                           print *, "singurality"
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
                           print *, "no singurality"
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_p), zeta(b_l_p), theta(c_m_p), theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_gugu_pp_1 = s
                           s = 0.0q0
                        end if
                        print *, "L2"
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l_p), zeta(b_l_p), theta(c_m_p), theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_gugu_pp_2 = s
                        s = 0.0q0
                     else
                        print *, "hanigai_pp"
                        val_gugu_pp_1 = 0.0q0
                        val_gugu_pp_2 = 0.0q0
                     end if
                     !-----------------------pm---------------------l_p,m_mになる要チェック
                     if (0 .le. a_l_p .and. a_l_p .le. Lmax .and. 0 .le. b_l_p .and. b_l_p .le. Lmax .and. &!pm
                         0 .le. c_m_m .and. c_m_m .le. Mmax .and. 0 .le. d_m_m .and. d_m_m .le. Mmax) then
                        alpha = zeta(l + 1); beta = zeta(l + 2); gamma = theta(m - 1); delta = theta(m - 2)

                        if (zeta(a_l_p) .le. zeta_j .and. zeta_j .le. zeta(b_l_p) &
                            .and. theta(c_m_m) .le. theta_k .and. theta_k .le. theta(d_m_m)) then!特異点があるとき
                           print *, "singurality"
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
                           print *, "no singurality"
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_p), zeta(b_l_p), theta(c_m_m), theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_gugu_pm_1 = s
                           s = 0.0q0
                        end if
                        print *, "L2"
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l_p), zeta(b_l_p), theta(c_m_m), theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_gugu_pm_2 = s
                        s = 0.0q0
                     else
                        print *, "hanigai_pm"
                        val_gugu_pm_1 = 0.0q0
                        val_gugu_pm_2 = 0.0q0
                     end if
                     !-----------------------mp---------------------l_m,m_pになる要チェック
                     if (0 .le. a_l_m .and. a_l_m .le. Lmax .and. 0 .le. b_l_m .and. b_l_m .le. Lmax .and. &!mp
                         0 .le. c_m_p .and. c_m_p .le. Mmax .and. 0 .le. d_m_p .and. d_m_p .le. Mmax) then
                        alpha = zeta(l - 1); beta = zeta(l - 2); gamma = theta(m + 1); delta = theta(m + 2)

                        if (zeta(a_l_m) .le. zeta_j .and. zeta_j .le. zeta(b_l_m) &
                            .and. theta(c_m_p) .le. theta_k .and. theta_k .le. theta(d_m_p)) then!特異点があるとき
                           print *, "singurality"
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
                           print *, "no singurality"
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_m), zeta(b_l_m), theta(c_m_p), theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_gugu_mp_1 = s
                           s = 0.0q0
                        end if
                        print *, "L2"
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l_m), zeta(b_l_m), theta(c_m_p), theta(d_m_p), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_gugu_mp_2 = s
                        s = 0.0q0
                     else
                        print *, "hanigai_mp"
                        val_gugu_mp_1 = 0.0q0
                        val_gugu_mp_2 = 0.0q0
                     end if
                     !-----------------------mm---------------------l_m,m_mになる要チェック
                     if (0 .le. a_l_m .and. a_l_m .le. Lmax .and. 0 .le. b_l_m .and. b_l_m .le. Lmax .and. &
                         0 .le. c_m_m .and. c_m_m .le. Mmax .and. 0 .le. d_m_m .and. d_m_m .le. Mmax) then
                        alpha = zeta(l - 1); beta = zeta(l - 2); gamma = theta(m - 1); delta = theta(m - 2)

                        if (zeta(a_l_m) .le. zeta_j .and. zeta_j .le. zeta(b_l_m) &
                            .and. theta(c_m_m) .le. theta_k .and. theta_k .le. theta(d_m_m)) then!特異点があるとき
                           print *, "singurality"
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
                           print *, "no singurality"
                           call dde3d(L1, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                      zeta(a_l_m), zeta(b_l_m), theta(c_m_m), theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)
                           val_gugu_mm_1 = s
                           s = 0.0q0
                        end if
                        print *, "L2"
                        call dde3d(L2, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, &
                                   zeta(a_l_m), zeta(b_l_m), theta(c_m_m), theta(d_m_m), 0.0q0, 2.0q0*pi, eps, s, info)
                        val_gugu_mm_2 = s
                        s = 0.0q0
                     else
                        print *, "haingai mm"
                        val_gugu_mm_1 = 0.0q0
                        val_gugu_mm_2 = 0.0q0
                     end if
                     K_jklm(j, k, l, m) = val_gugu_pp_1 + val_gugu_pp_2 + val_gugu_pm_1 + val_gugu_pm_2 &
                                          + val_gugu_mp_1 + val_gugu_mp_2 + val_gugu_mm_1 + val_gugu_mm_2
!-----------------奇数偶数------------------!-----------------奇数偶数------------------!-----------------奇数偶数------------------
                     !-----------------奇数偶数------------------J=77 k=1 l=1 m=16で止まってた あんまりうまくいってなさそう
                  else if (mod(l, 2) .eq. 1 .and. mod(m, 2) .eq. 0) then
                     a_l = l - 1; b_l = l + 1; c_m_p = m; d_m_p = m + 2; ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です
                     c_m_m = m - 2; d_m_m = m  ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です

                     !-----------------------0p---------------------m_pになる要チェック
                     if (0 .le. a_l .and. a_l .le. Lmax .and. 0 .le. b_l .and. b_l .le. Lmax .and. &!pp
                         0 .le. c_m_p .and. c_m_p .le. Mmax .and. 0 .le. d_m_p .and. d_m_p .le. Mmax) then
                        alpha = zeta(l + 1); beta = zeta(l - 1); gamma = theta(m + 1); delta = theta(m + 2)

                        if (zeta(a_l) .le. zeta_j .and. zeta_j .le. zeta(b_l) &
                            .and. theta(c_m_p) .le. theta_k .and. theta_k .le. theta(d_m_p)) then!特異点があるとき
                           print *, "singurality"
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
                           print *, "no singurality"
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
                        print *, "hanigai 0p"
                        val_kigu_0p_1 = 0.0q0
                        val_kigu_0p_2 = 0.0q0
                     end if
                     !-----------------------0m---------------------m_mになる要チェック
                     if (0 .le. a_l .and. a_l .le. Lmax .and. 0 .le. b_l .and. b_l .le. Lmax .and. &!pm
                         0 .le. c_m_m .and. c_m_m .le. Mmax .and. 0 .le. d_m_m .and. d_m_m .le. Mmax) then
                        alpha = zeta(l + 1); beta = zeta(l - 1); gamma = theta(m - 1); delta = theta(m - 2)

                        if (zeta(a_l) .le. zeta_j .and. zeta_j .le. zeta(b_l) &
                            .and. theta(c_m_m) .le. theta_k .and. theta_k .le. theta(d_m_m)) then!特異点があるとき
                           print *, "singurality"
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
                           print *, "no singurality"
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
                        print *, "hanigai 0m"
                        val_kigu_0m_1 = 0.0q0
                        val_kigu_0m_2 = 0.0q0
                     end if
                     K_jklm(j, k, l, m) = val_kigu_0p_1 + val_kigu_0p_2 + val_kigu_0m_1 + val_kigu_0m_2

                     !------------------偶数奇数-------------------!j=9 k=1 l=0 m=11 で0　墓にもいろいろ0
                  else if (mod(l, 2) .eq. 0 .and. mod(m, 2) .eq. 1) then
                     a_l_p = l; b_l_p = l + 2; a_l_m = l - 2; b_l_p = l  ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です
                     c_m = m - 1; d_m = m + 1  ! こいつらは積分区間の中身　マイナスになるときは範囲外なので積分しない　そのためだけの変数です

                     !-----------------------p0---------------------l_pになる要チェック
                     if (0 .le. a_l_p .and. a_l_p .le. Lmax .and. 0 .le. b_l_p .and. b_l_p .le. Lmax .and. &!pp
                         0 .le. c_m .and. c_m .le. Mmax .and. 0 .le. d_m .and. d_m .le. Mmax) then
                        alpha = zeta(l + 1); beta = zeta(l + 2); gamma = theta(m + 1); delta = theta(m - 1)

                        if (zeta(a_l_p) .le. zeta_j .and. zeta_j .le. zeta(b_l_p) &
                            .and. theta(c_m) .le. theta_k .and. theta_k .le. theta(d_m)) then!特異点があるとき
                           print *, "singurality"
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
                           print *, "no singurality"
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
                        print *, "hanigai p0"
                        val_guki_p0_1 = 0.0q0
                        val_guki_p0_2 = 0.0q0
                     end if
                     !-----------------------m0---------------------l_mになる要チェック
                     if (0 .le. a_l_m .and. a_l_m .le. Lmax .and. 0 .le. b_l_m .and. b_l_m .le. Lmax .and. &!pm
                         0 .le. c_m .and. c_m .le. Mmax .and. 0 .le. d_m .and. d_m .le. Mmax) then
                        alpha = zeta(l - 1); beta = zeta(l - 2); gamma = theta(m - 1); delta = theta(m + 1)

                        if (zeta(a_l_m) .le. zeta_j .and. zeta_j .le. zeta(b_l_m) &
                            .and. theta(c_m) .le. theta_k .and. theta_k .le. theta(d_m)) then!特異点があるとき
                           print *, "singurality"

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
                           print *, "no singurality"

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

   open (unit=10, iostat=ios, file='Kjklm64_binary_ver3.dat', action='write', &
       & form='unformatted', status='new')
   ! ファイルが正常に開けたかどうかをチェックする
   if (ios /= 0) then
      write (*, *) 'Failed to open file for output'
      stop
   end if
   write (10) k_jklm
   close (10)
contains
   function F1(zeta_bar, theta_bar, psi, zeta_j, theta_k)!計算できる
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_j, theta_k
      real(16)::F1

      F1 = (zeta_bar)*(zeta_bar) + (zeta_j)*(zeta_j) + (theta_bar)*(theta_bar) + (theta_k)*(theta_k) &
           - 2*((zeta_bar)*(zeta_bar) + (theta_bar)*(theta_k)*cos(psi))
   end function
   function F2(zeta_bar, theta_bar, psi, zeta_j, theta_k)
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_j, theta_k
      real(16)::F2
      F2 = (zeta_j*(zeta_j) + (theta_k)*(theta_k))*((zeta_j)*(zeta_j) + (theta_k)*(theta_k)) &
           - ((zeta_bar)*(zeta_bar) + (theta_bar)*(theta_k)*cos(psi))**2
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

      L1 = 1/sqrt(2.0q0)/pi*cos(psi) &
           /sqrt(F1(zeta_bar, theta_bar, psi, zeta_j, theta_k)) &
           *exp(-zeta_bar**2 + F2(zeta_bar, theta_bar, psi, zeta_j, theta_k) &
                /F1(zeta_bar, theta_bar, psi, zeta_j, theta_k)) &
           *f_for_psi(zeta_bar, theta_bar, zeta_l, theta_m, alpha, beta, gamma, delta)
   end function
   function L1_3variable(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, switch)
      implicit none
      real(16), intent(in)::zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta, switch
      real(16)::L1_3variable
      ! if (mod(l, 2) .eq. 1 .and. mod(m, 2) .eq. 1) then!奇数奇数
      if (switch .eq. 0) then
         L1_3variable = L1(zeta_bar, theta_bar, psi, zeta_j, theta_k, zeta_l, theta_m, alpha, beta, gamma, delta)
      end if
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

      L2 = 1/2/sqrt(2.0q0)/pi*cos(psi) &
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
