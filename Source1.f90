
  !$ use omp_lib  
  implicit none

  !***************************************************************************   
  !*** Definition ************************************************************
  !***************************************************************************
  
  !------------------------------------------------------------
  !Dummy variables --------------------------------------------
  integer :: i, j, n
  integer :: nmax = 10000000    !maximum number of iteration
  logical :: converge = .false. !if this becomes "true", computation finishs
  integer :: mythread           !openmp: my thread id
  integer :: threads            !openmp: total number of threads
    
  !------------------------------------------------------------
  !Auxiliary variables for general part -----------------------
  integer :: cputime(1:99,1:2) = 0       !#1 -> proc, #2 -> 1 is start time, 2 is end time
  real(8) :: cputime_sum(1:99) = 0d0     !total cpu time used by proc
  integer :: CountPerSec, CountMax       !For cpu time measurement
  real(8) :: pi  = 4d0*datan(1d0)        !pi
  real(8) :: spi = dsqrt(4d0*datan(1d0)) !sqrt(pi)
  real(8) :: large = 1d99                !some very large value
  real(8) :: small = 1d-16               !some very small value 
  integer :: display = 10                !if mod(n,display) = 0, then show the situation on the console
  integer :: proc                        !For cputime measurement
  
  !------------------------------------------------------------
  !@@@ Physical and numerical parameters ----------------------
  real(8) :: Kn    !Knudsen number*sqrt(pi)/2
  integer :: Nx    !Nx
  integer :: Nz    !Nz
  integer :: LS_x  !lattice system for x
  integer :: LS_z  !lattice system for z
  integer :: FD_x  !1-> 1st order finite difference, 2-> second order finite difference
  integer :: INT_x !1-> trapezoidal rule for integration, 2-> simpson rule
  integer :: INT_z !1-> trapezoidal rule for integration, 2-> simpson rule
  integer :: cap_Z !capital Z for the maximum value of zeta
  integer :: save1 !data save if mod(n,save1) = 0
    
  !------------------------------------------------------------
  !Physical quantities----------------------------------------
  real(8), allocatable :: x(:), dx(:)            !(0:Nx) 
  real(8), allocatable :: z(:), dz(:)            !(-Nz:Nz) 
  real(8), allocatable :: g(:,:)                 !(0:Nx,-Nz:Nz) 
  real(8), allocatable :: uP(:)                  !(0:Nx) 
  real(8)              :: MP, MP_old = 0d0       !flow rate
  real(8)              :: err, threshold = 1d-10 !if "err" < threshold, computation finishs
  real(8)              :: gamma1 = 1d0           !for BGK
  real(8)              :: k0     = -1.01619d0    !for BGK
    
  !------------------------------------------------------------
  !Variables for I/O ------------------------------------------
  integer       :: ierr
  integer       :: hostnm              !Name of machine
  integer       :: OS                  !1-> windows, 2-> unix
  character(6)  :: id_char, OS_char
  character(10) :: date, time, zone
  integer       :: date_time(8)
  character(99) :: file_name           !Name of output file
  character(38) :: dir_name            !Name of data directory
  character(15) :: comment_char        !comment shown in the finename
  character(99) :: input_name          !Name of input file
  character(99) :: program_name        !Name of f90 file
  character(99) :: host_name           !Name of machine
  character(999):: parameter_set       !Name of parameters contained in one character
  integer       :: process_id          !unix process id

  !*************************************************************************** 
  !*** File read section *****************************************************
  !***************************************************************************

  !Obtain input file name (6 digits number)
  call getarg( 1, input_name )
  id_char = trim(adjustl(input_name))
  if( len(trim(adjustl(input_name))) .ne. 6 ) call terminate(6424)
  input_name = "00input-"//trim(adjustl(input_name))//".txt"

  !Open input file and read parameters
  open( unit = 21, file = trim(adjustl(input_name)), action = "read" )
  read(21,*) OS
  write(OS_char,'(i1.1)') OS
  read(21,*) threads
  !---@@@ Modify as you want ---
  read(21,*) Kn
  read(21,*) Nx
  read(21,*) Nz
  read(21,*) LS_x  
  read(21,*) LS_z  
  read(21,*) Fd_x 
  read(21,*) INT_x
  read(21,*) INT_z
  read(21,*) cap_Z 
  read(21,*) save1
  
  if( INT_z == 2 ) then
    if( mod(Nz,2) == 0 ) then
      write(*,*) "warning! Nz should be odd if INT_z = 2"
      write(*,*) "Nz -> Nz + 1"
      Nz = Nz + 1
    end if
  end if
  
  !Check mythread for openmp
  !$ call omp_set_num_threads(threads)
  !$omp parallel private(mythread)
  !$ mythread = omp_get_thread_num()
  !$ write(*,*) "mythread = ", mythread
  !$omp end parallel

  !Make "data directory" with unique name using date & time
  call date_and_time(date, time, zone, date_time)
  time = time(1:6)
  call getarg( 2, comment_char ) 
  comment_char = trim(adjustl(comment_char))//"----------------"
  dir_name = trim(date)//"_"//trim(time)//"-"//id_char//"-"//comment_char
  call system('mkdir'//" "//dir_name//"_running") !make data directry
  call system('mkdir'//' "'//dir_name//'_running/lattice"')
  call system('mkdir'//' "'//dir_name//'_running/vdf"')
  call system('mkdir'//' "'//dir_name//'_running/mac"')
  call system('mkdir -p'//' MP')

  !Copy input file to "data directory"
  if( OS == 1 ) call system('copy'//" "//trim(adjustl(input_name))//" "//dir_name//"_running")
  if( OS == 2 ) call system('cp'//" "//trim(adjustl(input_name))//" "//dir_name//"_running")
    
  !Copy the source f90 code to "data directory" 
  call get_command_argument( 0, program_name )
  if( OS == 1 ) call system('copy'//" "//trim(adjustl(program_name))//".f90 "//dir_name//"_running")
  if( OS == 1 ) call system('copy'//" "//trim(adjustl(program_name))//".exe "//dir_name//"_running")
  if( OS == 1 ) call system('copy'//" "//"*.f90 "//" "//dir_name//"_running") !VS2019
  !** for windows, exe file cannot be deleted when running
  if( OS == 2 ) call system('cp'//" "//trim(adjustl(program_name))//".f90 "//dir_name//"_running")
  if( OS == 2 ) call system('cp'//" "//trim(adjustl(program_name))//" "//dir_name//"_running")
  !if( OS == 2 ) call system('rm'//" "//trim(adjustl(program_name)))
  
  !Copy lay file
  if( OS == 1 ) call system('copy'//" "//"00matome.lay"//" "//dir_name//"_running")
  if( OS == 2 ) call system('cp'//" "//"00matome.lay"//" "//dir_name//"_running")

  !Write parameters to "parameter.csv"
  open( unit = 22, file = "00parameter.csv", position = "append" )
  parameter_set = ""
  call names_char(dir_name,parameter_set)
  call names_char(id_char ,parameter_set)
  call names_char(OS_char ,parameter_set)
  call names_inte(threads ,parameter_set)
  !--- @@@ Keep the same structure with "Modify as you want" ---
  !--- names_XXXX adds #1 to the end of #2 ---------------------
  !--- This is for csv output ----------------------------------
  call names_real(Kn     ,parameter_set)
  call names_inte(Nx     ,parameter_set)
  call names_inte(Nz     ,parameter_set)
  call names_inte(LS_x   ,parameter_set)
  call names_inte(LS_z   ,parameter_set)
  call names_inte(FD_x   ,parameter_set)
  call names_inte(INT_x  ,parameter_set)
  call names_inte(INT_z  ,parameter_set)
  call names_inte(cap_Z  ,parameter_set)
  call names_inte(save1  ,parameter_set)
  !-----------------------------------------------------------
  ierr = hostnm(host_name) 
  call names_char(host_name,parameter_set)
  write(22,'(a,a)') trim(adjustl(parameter_set))//"running, "
  close(22)

  !*************************************************************************** 
  !*** Main processing *******************************************************
  !***************************************************************************

  call allocating()      !allocate 
  call initialize()      !input initial values
  call lattice_x()       !make x lattice
  call lattice_z()       !make z lattice
  call OUTPUT_lattice( ) !visualize lattice system
    
  !--- @ Solve your problem ---

  do n = 0, nmax
      
         if( converge == .true. ) exit
        
         if( mod(n,display) == 0 ) call show_display()
      
         !proc = 1 is to measure whole cpu time      
         proc = 1
         call system_clock( cputime(proc,1), CountPerSec, CountMax)

         !------------------------------------------------------------
         proc = 2
         call system_clock( cputime(proc,1), CountPerSec, CountMax)
      !boundary condition at x2 = 1/2 and z2 < 0   
      call COMPUTE_bc( Nx )
         call system_clock( cputime(proc,2) )
         cputime_sum(proc) = cputime_sum(proc) + dble(cputime(proc,2) - cputime(proc,1))/dble(CountPerSec)
         !------------------------------------------------------------
         
         !------------------------------------------------------------
         proc = 4
         call system_clock( cputime(proc,1), CountPerSec, CountMax)
      !finite difference for z2 < 0
      call COMPUTE_vdf( -1 )
         call system_clock( cputime(proc,2) )
         cputime_sum(proc) = cputime_sum(proc) + dble(cputime(proc,2) - cputime(proc,1))/dble(CountPerSec)
         !------------------------------------------------------------
         
         !------------------------------------------------------------
         proc = 3
         call system_clock( cputime(proc,1), CountPerSec, CountMax)
      !boundary condition at x2 = 0 and z2 > 0
      call COMPUTE_bc( 0 )
         call system_clock( cputime(proc,2) )
         cputime_sum(proc) = cputime_sum(proc) + dble(cputime(proc,2) - cputime(proc,1))/dble(CountPerSec)
         !------------------------------------------------------------
         
         !------------------------------------------------------------
         proc = 5
         call system_clock( cputime(proc,1), CountPerSec, CountMax)
      !finite difference for z2 > 0   
      call COMPUTE_vdf( 1 )
         call system_clock( cputime(proc,2) )
         cputime_sum(proc) = cputime_sum(proc) + dble(cputime(proc,2) - cputime(proc,1))/dble(CountPerSec)
         !------------------------------------------------------------

         !------------------------------------------------------------
         proc = 7
         call system_clock( cputime(proc,1), CountPerSec, CountMax)
      !integration of g with respect to z2 to obtain uP
      call COMPUTE_mac( )
      !integration of uP with respect to x2 to obtain MP
      call COMPUTE_flowrate( )
         call system_clock( cputime(proc,2) )
         cputime_sum(proc) = cputime_sum(proc) + dble(cputime(proc,2) - cputime(proc,1))/dble(CountPerSec)
         !------------------------------------------------------------
         
         !------------------------------------------------------------
         proc = 6
         call system_clock( cputime(proc,1), CountPerSec, CountMax)
      !visualize g
      call OUTPUT_vdf( )
      !visualize uP and MP
      call OUTPUT_mac( )
         call system_clock( cputime(proc,2) )
         cputime_sum(proc) = cputime_sum(proc) + dble(cputime(proc,2) - cputime(proc,1))/dble(CountPerSec)
         !------------------------------------------------------------
         
         !proc = 1 is to measure whole cpu time
         proc = 1
         call system_clock( cputime(proc,2) )
         cputime_sum(proc) = cputime_sum(proc) + dble(cputime(proc,2) - cputime(proc,1))/dble(CountPerSec)
      
  end do
    
  call OUTPUT_final()

  !-----------------------------

  !*************************************************************************** 
  !*** Finalize section ******************************************************
  !***************************************************************************
  
  !Rename data directory
  write(*,*), trim(adjustl(dir_name))
  if( OS == 1 ) call system('xcopy /D /S /R /Y /I /K'//" "//dir_name//"_running"//" "//dir_name)
  if( OS == 1 ) call system('rd /s/ q'//" "//dir_name//"_running")
  if( OS == 2 ) call system('cp -r'//" "//dir_name//"_running"//" "//dir_name)
  
  !Remove "running" directory
  if( OS == 2 ) call system('rm -r'//" "//dir_name//"_running")
  
  !this is the end of main program
  
  !===========================================================
  !===========================================================
  !===========================================================
  !=== Subroutine ============================================
  !===========================================================
  !===========================================================
  !===========================================================

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
    allocate( x(0:Nx), dx(1:Nx) )
    allocate( z(-Nz:Nz), dz(-Nz+1:Nz) )
    allocate( g(0:Nx,-Nz:Nz) )
    allocate( uP(0:Nx) )
  end subroutine 
  
  !Set initial setting
  subroutine initialize()
  implicit none
    g  = large
    uP = 0d0
  end subroutine 

  !===========================================================
  !===========================================================
  !===========================================================
  !=== Lattice system=========================================
  !===========================================================
  !===========================================================
  !===========================================================

  !lattice for x
  subroutine lattice_x()
  implicit none
  
    if( LS_x == 1 ) then
      do i = 0, Nx
        x(i) = 0.5d0*dble(i)/dble(Nx)
      end do
    end if
    
    if( LS_x == 2 ) then
      do i = 0, Nx
        x(i) = 0d0 !20191218
      end do
    end if
    
    !simpson rule
    if( INT_x == 2 ) then
      do i = 1, Nx - 1, 2
        x(i) = ( x(i-1) + x(i+1) ) * 0.5d0
      end do
    end if    
    
    !delta _x
    do i = 1, Nx
      dx(i) = x(i) - x(i-1)
    end do
  end subroutine
  
  !lattice for z
  subroutine lattice_z()
  implicit none
  integer :: pm
  
    z(0) = large !j=0 must not be used
    
    if( LS_z == 1 ) then
      do pm = -1, 1, 2
        do j = 1, Nz
          z(pm*j) = cap_Z * dble(pm) * dble(j-1)/dble(Nz-1)
          if( abs(j) == 1 ) z(pm*j) = dble(pm) * small 
        end do
      end do
    end if
    
    if( LS_z == 2 ) then
      do pm = -1, 1, 2
        do j = 1, Nz
          z(pm*j) = 0d0 !20191218
          if( abs(j) == 1 ) z(pm*j) = dble(pm) * small 
        end do
      end do
    end if
    
    !simpson rule
    if( INT_z == 2 ) then
      do j = -Nz+1, Nz - 1, 2
        z(j) = ( z(j-1) + z(j+1) ) * 0.5d0
      end do
    end if    
    
    !delta z
    do j = -Nz + 1, Nz
      if( j == 0 .or. j == 1 ) dz(j) = large
      dz(j) = z(j) - z(j-1)
    end do  
  end subroutine 
  
  !===========================================================
  !===========================================================
  !===========================================================
  !=== Computation ===========================================
  !===========================================================
  !===========================================================
  !===========================================================  
  
  !macro -----------------------------------------------------
  subroutine COMPUTE_mac( )
  implicit none
  real(8) :: g0, g1, g2
  
    if( INT_z == 1 ) then
    !$omp parallel default(shared), private(i,j,g0,g1)
    !$omp do  
    do i = 0, Nx
      uP(i) = 0d0
      do j = -Nz, Nz - 1
        if( j == -1 .or. j == 0 ) cycle
        g0    = 1d0/spi * g(i,j+0) * dexp( - z(j+0)*z(j+0) )
        g1    = 1d0/spi * g(i,j+1) * dexp( - z(j+1)*z(j+1) )
        uP(i) = uP(i) + 0.25d0 * dz(j+1) * ( g0 + g1 )
      end do
    end do
    !$omp end do
    !$omp end parallel
    end if
    
    if( INT_z == 2 ) then
    !$omp parallel default(shared), private(i,j,g0,g1,g2)
    !$omp do  
    do i = 0, Nx
      uP(i) = 0d0
      do j = -Nz, Nz - 2, 2
        if( j == -1 .or. j == 0 ) cycle
        uP(i) = 0d0 !20191218
      end do
    end do
    !$omp end do
    !$omp end parallel
    end if
  
  end subroutine
  
  !flow rate---------------------------------------------------
  subroutine COMPUTE_flowrate( )
  implicit none
  real(8) :: g0, g1
  
    MP_old = MP
    MP = 0d0
    if( INT_x == 1 ) then
      do i = 0, Nx - 1
        MP = MP + 0.5d0 * dx(i+1) * ( uP(i) + uP(i+1) )
      end do
    end if
    if( INT_x == 2 ) then
      do i = 0, Nx - 2, 2
        MP = 0d0 !20191218
      end do
    end if
    MP = 2d0*MP


    
    err = dabs( MP_old - MP )/dabs(MP)
    if( err < threshold ) then
      converge = .true.
    end if
  
  end subroutine
 
  !vdf -----------------------------------------------------
  subroutine COMPUTE_vdf( dir )
  implicit none
  integer, intent(in) :: dir
  real(8) :: c1, c2, c3
  real(8) :: a, b, S0, S1, S2
  integer :: st, en, ip

  !determin the direction of finite difference scheme
  if( dir == -1 ) then
    st = Nx - 1
    en = 0
    ip = 1
  else if( dir == 1 ) then
    st = 1
    en = Nx
    ip = 0
  else 
    call terminate(6764)
  end if
  
  if( FD_x == 1 ) then
  !$omp parallel default(shared), private(i,j,c1,c2,c3)
  !$omp do  
  do j = dir*1, dir*Nz, dir
    do i = st, en, dir
      c1     = 1d0 + dble(dir)*dx(i+ip)/z(j)/Kn
      c2     = dble(dir)*dx(i+ip)/z(j)/Kn * 2d0*uP(i)
      c3     = - dble(dir)*dx(i+ip)/z(j)
      g(i,j) = ( g(i-dir,j) + c2 + c3 )/c1
    end do
  end do
  !$omp end do
  !$omp end parallel
  end if

  if( FD_x == 2 ) then
  !$omp parallel default(shared), private(i,j,c1,c2,c3,a,b,S0,S1,S2)
  !$omp do  
  do j = dir*1, dir*Nz, dir
    i = st
    g(i,j) = 0d0 !20191218
    do i = st + dir, en, dir
      g(i,j) = 0d0 !20191218
    end do
  end do
  !$omp end do
  !$omp end parallel
  end if
  
  end subroutine
      
  !BC -----------------------------------------------------
  subroutine COMPUTE_bc( i_input )
  implicit none
  integer, intent(in) :: i_input
  
    !BC on the wall
    if( i_input == Nx ) then
      i = Nx
      do j = -Nz, -1
        g(i,j) = 0d0 !20191218
      end do
    end if

    !BC at x2 = 0  
    if( i_input == 0 ) then
      i = 0
      do j = 1, Nz
        g(i,j) = g(i,-j)
      end do
    end if
    
    if( i_input < Nx .and. i_input > 0 ) then
      call terminate(4847)
    end if
    
  end subroutine
  
  
  !===========================================================
  !===========================================================
  !===========================================================
  !=== Output ================================================
  !===========================================================
  !===========================================================
  !===========================================================  

  !Show computational time-----------------------------------------------------
  subroutine show_display( )
  implicit none
  character(2) :: threads_char
  
  write(threads_char, '(i2.2)') threads
  file_name = "output_omp="//threads_char//".dat"
  open( unit = 99, file = dir_name//"_running"//"/"//trim(adjustl(file_name)), position = "append" )
  
  if( n == 0 ) then
    write(99,'(a)') "/// cpu time ///"
    write(99,'(19a)')    "       n", &
                      "      Total", & 
                      "    BC wall", &
                      "    BC refl", &
                      " vdf dir -1", &
                      " vdf dir  1", &
                      "     output", &
                      "compute mac"
  end if
  write(99,'(i8,18f11.2)') n, cputime_sum(1), cputime_sum(2), &
                              cputime_sum(3), cputime_sum(4), &
                              cputime_sum(5), cputime_sum(6), &
                              cputime_sum(7), cputime_sum(8)
  close(99)

  end subroutine
  
  !lattice system------------------------------------
  subroutine OUTPUT_lattice( )
  implicit none
  
  file_name = "lattice_x.dat"
  open( unit = 22, file = dir_name//"_running"//"/lattice/"//trim(adjustl(file_name)) )
    do i = 0, Nx
      write(22,*) i, x(i)
    end do
  close(22)
  
  file_name = "lattice_z.dat"
  open( unit = 22, file = dir_name//"_running"//"/lattice/"//trim(adjustl(file_name)) )
    do j = -Nz, Nz
      if( j == 0 ) cycle
      write(22,*) j, z(j)
    end do
  close(22)

  end subroutine
  
  !vdf------------------------------------
  subroutine OUTPUT_vdf( )
  implicit none
  real(8) :: arg, anal
  
  if( mod(n,save1) == 0 ) then
  
    !If uP = 0, g can be obtained analytically. 
    file_name = "vdf_check.dat" 
    open( unit = 22, file = dir_name//"_running"//"/vdf/"//trim(adjustl(file_name)) )
      j = Nz/2 !choose any j 
      do i = 0, Nx
        arg  = -( x(i) - 0.5d0 )/z(j)/Kn
        if( j < 0 ) then
          anal = Kn*( dexp(arg) - 1d0 )
        else
          anal = Kn*( dexp(arg - 1d0/z(j)/Kn ) - 1d0 )  
        end if
        write(22,*) x(i), g(i,j), anal
      end do
    close(22)
    
    !
    file_name = "vdf.dat" 
    open( unit = 22, file = dir_name//"_running"//"/vdf/"//trim(adjustl(file_name)) )
      do i = 0, Nx
        write(22,'(a,e15.5,a)') 'zone t = "',x(i),'"'
        do j = -Nz, Nz
          if( j == 0 ) cycle
          write(22,*) z(j), g(i,j)
        end do
      end do
    close(22)

  end if  
    
  end subroutine
  
  !mac------------------------------------
  subroutine OUTPUT_mac( )
  implicit none
  
  if( mod(n,save1) == 0 ) then
      print *, n
      
    file_name = "mac.dat" 
    open( unit = 22, file = dir_name//"_running"//"/mac/"//trim(adjustl(file_name)) )
      do i = 0, Nx
        write(22,*) x(i), uP(i)
      end do
    close(22)
    
  end if
    
    file_name = "MP.dat" 
    open( unit = 22, file = dir_name//"_running"//"/mac/"//trim(adjustl(file_name)), position = "append" )
      write(22,*) n, MP
    close(22)

  end subroutine
  
  !final data------------------------------------
  subroutine OUTPUT_final( )
  implicit none
  real(8) :: MP_theoretical
  
    file_name = "final.dat" 
    open( unit = 22, file = dir_name//"_running"//"/"//trim(adjustl(file_name)) )
      MP_theoretical = 1d0/(2d0*gamma1*Kn) * ( 1d0/6d0 - k0*Kn ) !for Kn << 1
      write(22,*) -MP, MP_theoretical, dabs( -MP - MP_theoretical )/MP_theoretical
      MP_theoretical = 1d0/2d0/dsqrt(pi)*dlog(Kn) !for Kn >> 1
      write(22,*) -MP, MP_theoretical, dabs( -MP - MP_theoretical )/MP_theoretical
    close(22)
    
    file_name = "final"//dir_name//".dat" 
    open( unit = 22, file = "MP"//"/"//trim(adjustl(file_name)) )
    write(22,*) 'variables = "log Kn" "-MP" "relative error" "Nx" "Nz" '
    write(22,*) 'zone t = "',dir_name,'"'
    
    file_name = "MP.dat" 
    open( unit = 23, file = "MP"//"/"//trim(adjustl(file_name)), position = "append", status = "old", iostat = ierr )
    if( ierr /= 0 ) then
      open( unit = 23, file = "MP"//"/"//trim(adjustl(file_name)), position = "append" )
      write(23,*) 'variables = "log Kn" "-MP" "relative error" "Nx" "Nz" '
    end if
    write(23,*) 'zone t = "',dir_name,'"'
    
      MP_theoretical = 1d0/(2d0*gamma1*Kn) * ( 1d0/6d0 - k0*Kn ) !for Kn << 1
      write(22,'(3e25.15,2i5)') dlog10(Kn), -MP, dabs( -MP - MP_theoretical )/MP_theoretical, Nx, Nz
      write(23,'(3e25.15,2i5)') dlog10(Kn), -MP, dabs( -MP - MP_theoretical )/MP_theoretical, Nx, Nz
      
    close(22)

  end subroutine
  
  !===========================================================
  !===========================================================
  !===========================================================
  !=== File treatment=========================================
  !===========================================================
  !===========================================================
  !===========================================================  

  !Convert arg1(integer) to CHARACTER and combine with arg2
  subroutine names_inte(arg1,arg2)
    implicit none
    integer        :: arg1
    character(100) :: arg1_char
    character(999) :: arg2
    write(arg1_char,*) arg1
    arg2 = trim(adjustl(arg2))//trim(adjustl(arg1_char))//","
  end subroutine names_inte

  !Convert arg1(real) to CHARACTER and combine with arg2
  subroutine names_real(arg1,arg2)
    implicit none
    real(8)        :: arg1
    character(100) :: arg1_char
    character(999) :: arg2
    write(arg1_char,*) arg1
    arg2 = trim(adjustl(arg2))//trim(adjustl(arg1_char))//","
  end subroutine names_real

  !Convert arg1(character) to CHARACTER and combine with arg2
  subroutine names_char(arg1,arg2)
    implicit none
    character(*)   :: arg1
    character(999) :: arg2
    arg2 = trim(adjustl(arg2))//trim(adjustl(arg1))//","
  end subroutine names_char

  !Terminate with error code X
  subroutine terminate(arg1)
    implicit none
    integer :: arg1
    write(*,*) "Error! ", arg1
    stop
  end subroutine terminate
  end