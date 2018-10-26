MODULE XSbuilder
implicit none

logical,parameter :: TEST_A = .false. !true for test a, flase for test b
real(8),parameter :: PI = (4*datan(1.0_8))
integer,parameter :: g=7, s=144
integer :: n
real(8) :: width, dw=0.01d0, wsum, wnorm, pin_width
real(8) :: R = 0.54d0, dr=0.01d0, r1, r2, r3, r4
integer :: max_lines
real(8),allocatable :: sigt(:,:), sigs(:,:,:), sigf(:,:), nuf(:,:), chi(:,:), siga(:,:)
real(8),allocatable :: Ox(:), Oy(:), w(:), alpha(:,:,:,:), tau(:,:,:,:), e_tau(:,:,:,:), fp(:), Ap(:), mu(:)
integer,allocatable :: char_lines(:,:,:)
real(8),allocatable :: char_lengths(:,:,:), char_x0(:,:), char_y0(:,:), char_xf(:,:), char_yf(:,:)
real(8),allocatable :: LR(:,:), ApStar(:)

contains
  subroutine buildXS()
  implicit none
  integer :: MOD=1 , UO2=3 , GT=2, MOX43=4, MOX70=5, MOX87=6, i, j, ierror, jj
  integer :: st=2, sa=3, sf=5, nf=6, bigX=7
  real(8) :: XSLib(6,7,7), XSlib_scat(6,7,7), XSsmall(2,7,4), gamma
  integer,allocatable :: geometery(:)

  if (TEST_A) then
    n = 40
    width=1.26d0
  else
    n = 360
    width = 3.78d0
  endif

  ! !formats
  ! 100 format (A , ', ')
  ! 101 format (A, I2, ', ')
  ! 105 format (f5.2, ', ')
  110 format (f11.8,', ')
  !
  !
  ! open(unit=71, file="XSdebugB.out", status='replace', action='write', iostat=ierror)
  ! if (ierror.ne.0) then
  !   write(*,*) "An error occured when opening XSdebug.out"
  ! endif

  !open files
  open(unit=11, file="ModSigs.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening ModSigs.dsv"
  endif
  open(unit=12, file="ModScat.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening ModScat.dsv"
  endif
  open(unit=21, file="UO2Sigs.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening UO2.dsv"
  endif
  open(unit=22, file="UO2Scat.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening UO2Scat.dsv"
  endif
  open(unit=31, file="GuideTubeSigs.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening GuideTubeSigs.dsv"
  endif
  open(unit=32, file="GuideTubeScat.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening GuideTubeScat.dsv"
  endif
  open(unit=41, file="MOX4.3.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening MOX4.3.dsv"
  endif
  open(unit=42, file="MOX4.3Scat.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening MOX4.3Scat.dsv"
  endif
  open(unit=51, file="MOX7.0Sigs.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening MOX7.0Sigs.dsv"
  endif
  open(unit=52, file="MOX7.0Scat.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening MOX7.0Scat.dsv"
  endif
  open(unit=61, file="MOX8.7.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening MOX8.7.dsv"
  endif
  open(unit=62, file="MOX8.7Scat.dsv", status='old', action='read', iostat=ierror)
  if (ierror.ne.0) then
    write(*,*) "An error occured when opening MOX8.7Scat.dsv"
  endif


  !cross section data
  XSlib = 0.0d0
  XSlib_scat = 0.00d0
  do i=1, g
    read(11,*) XSsmall(MOD,i,:)
    read(12,*) XSlib_scat(MOD,i,:)
    read(31,*) XSsmall(GT,i,:)
    read(32,*) XSlib_scat(GT,i,:)
    read(21,*) XSlib(UO2,i,:)
    read(22,*) XSlib_scat(UO2,i,:)
    read(41,*) XSlib(MOX43,i,:)
    read(42,*) XSlib_scat(MOX43,i,:)
    read(51,*) XSlib(MOX70,i,:)
    read(52,*) XSlib_scat(MOX70,i,:)
    read(61,*) XSlib(MOX87,i,:)
    read(62,*) XSlib_scat(MOX87,i,:)
  enddo
  XSlib(1:2,:,1:4) = XSsmall(:,:,:)

  ! do i=1, g
  !   do j=1, 7
  !     write(*,110,advance="no") XSlib_scat(MOD,i,j)
  !   enddo
  !   write(*,*)
  ! enddo

  max_lines = 2*INT((sqrt(2.0d0)*width / dw)) + 1
  !write(*,*) max_lines

  allocate( sigt(g,n), sigs(g,g,n), sigf(g,n), nuf(g,n), chi(g,n), geometery(n), siga(g,n), &
            Ox(s), Oy(s), w(s), mu(s), fp(n), Ap(n), ApStar(n), &
            tau(g,s,max_lines+1,n+1), e_tau(g,s,max_lines+1,n+1), alpha(g,s,max_lines+1,n+1), &
            char_lines(s,max_lines+1,n+1), char_lengths(s,max_lines+1,n+1), &
            char_x0(s,max_lines+1), char_y0(s,max_lines+1), char_xf(s,max_lines+1), char_yf(s,max_lines+1))

  char_lines = -1
  char_lengths = 0.0d0

  !build geometery
  r1 = R * 4.0d0/4
  r2 = R * 3.0d0/4
  r3 = R * 2.0d0/4
  r4 = R * 1.0d0/4

  if (TEST_A) then                           ! Test A
    pin_width = width
    do i=1, 8
      geometery(i) = MOD
      Ap(i) = (pin_width**2 - PI*R**2) / 8.0d0
    enddo
    do i=9, n
      geometery(i) = UO2
    enddo
    do i=9,16
      Ap(i) = (PI*r1**2 - PI*r2**2) / 8.0d0
    enddo
    do i=17,24
      Ap(i) = (PI*r2**2 - PI*r3**2) / 8.0d0
    enddo
    do i=25,32
      Ap(i) = (PI*r3**2 - PI*r4**2) / 8.0d0
    enddo
    do i=33,n
      Ap(i) = (PI*r4**2) / 8.0d0
    enddo

  else                                        ! Test B
    pin_width = 1.26d0
    do i=1, 8                                 ! pin cell 1
      geometery(i) = MOD
    enddo
    do i=9, 40
      geometery(i) = MOX43
      !geometery(i) = UO2
    enddo

    do i=41, 48                               ! pin cell 2
      geometery(i) = MOD
    enddo
    do i=49, 80
      geometery(i) = UO2
    enddo

    do i=81, 88                               ! pin cell 3
      geometery(i) = MOD
    enddo
    do i=89, 120
      geometery(i) = MOX43
      !geometery(i) = UO2
    enddo

    do i=121, 128                            ! pin cell 4
      geometery(i) = MOD
    enddo
    do i=129, 160
      geometery(i) = UO2
    enddo

    do i=161, 168                            ! pin cell 5
      geometery(i) = MOD
    enddo
    do i=169, 200
      geometery(i) = GT
      !geometery(i) = UO2
    enddo

    do i=201, 208                            ! pin cell 6
      geometery(i) = MOD
    enddo
    do i=209, 240
      geometery(i) = UO2
    enddo

    do i=241, 248                            ! pin cell 7
      geometery(i) = MOD
    enddo
    do i=249, 280
      geometery(i) = MOX87
      !geometery(i) = UO2
    enddo

    do i=281, 288                            ! pin cell 8
      geometery(i) = MOD
    enddo
    do i=289, 320
      geometery(i) = UO2
    enddo

    do i=321, 328                            ! pin cell 9
      geometery(i) = MOD
    enddo
    do i=329, 360
      geometery(i) = MOX70
      !geometery(i) = UO2
    enddo


    !Cell area
    do j=0, 320, 40
      !write(*,*) "j = ", j
      do i=j+1, j+8
        Ap(i) = (pin_width**2 - PI*R**2) / 8.0d0
        !write(*,*), "i = ", i
      enddo
      do i=j+9,j+16
        Ap(i) = (PI*r1**2 - PI*r2**2) / 8.0d0
        !write(*,*), "i = ", i
      enddo
      do i=j+17,j+24
        Ap(i) = (PI*r2**2 - PI*r3**2) / 8.0d0
        !write(*,*), "i = ", i
      enddo
      do i=j+25,j+32
        Ap(i) = (PI*r3**2 - PI*r4**2) / 8.0d0
        !write(*,*), "i = ", i
      enddo
      do i=j+33,j+40
        Ap(i) = (PI*r4**2) / 8.0d0
        !write(*,*), "i = ", i
      enddo
    enddo
  endif

  ! write(*,*) 'area'
  ! do i=1, n
  !   write(*, '(A I3, A, F7.5)') "Cell ", i, " = ", Ap(i)
  ! enddo
  ! write(*,*)

  !Build XS mesh
  do i=1, n
    sigt(:,i) = XSlib(geometery(i),:,st)
    siga(:,i) = XSlib(geometery(i),:,sa)
    sigf(:,i) = XSlib(geometery(i),:,sf)
    nuf(:,i) = XSlib(geometery(i),:,nf)
    chi(:,i) = XSlib(geometery(i),:,bigX)
    sigs(:,:,i) = XSlib_scat(geometery(i),:,:)
  end do

  ! do i=1, g
  !   do j=1, 7
  !     write(*,110,advance="no") sigs(i,j,9)
  !   enddo
  !   write(*,*)
  ! enddo

  !Define quadrature set

  !Quadrent 1 (+x, +y)
  Ox(1) = 9.717784813336E-01
  Oy(1) = 1.096881837272E-02
  w(1) = 8.454511187252E-03
  Ox(2) = 7.656319455497E-01
  Oy(2) = 1.160393058611E-02
  w(2) = 8.352354145856E-03
  Ox(3) = 4.445439440056E-01
  Oy(3) = 2.447911451942E-02
  w(3) = 1.460888798152E-02
  Ox(4) = 9.701698603928E-01
  Oy(4) = 5.695764868253E-02
  w(4) = 1.913728513580E-02
  Ox(5) = 7.633693960835E-01
  Oy(5) = 5.995074957044E-02
  w(5) = 1.873220073879E-02
  Ox(6) = 1.483114568272E-01
  Oy(6) = 1.670387759191E-02
  w(6) = 6.404244616724E-03
  Ox(7) = 9.622473153642E-01
  Oy(7) = 1.362124657777E-01
  w(7) = 2.863542971348E-02
  Ox(8) = 7.524467626583E-01
  Oy(8) = 1.419535016004E-01
  w(8) = 2.759429759588E-02
  Ox(9) = 9.410672772109E-01
  Oy(9) = 2.426233944222E-01
  w(9) = 3.648716160597E-02
  Ox(10) = 4.288508824476E-01
  Oy(10) = 1.196054590036E-01
  w(10) = 2.995376809966E-02
  Ox(11) = 7.241384940891E-01
  Oy(11) = 2.488983098158E-01
  w(11) = 3.442681426024E-02
  Ox(12) = 8.997294996538E-01
  Oy(12) = 3.673697853806E-01
  w(12) = 4.244873302980E-02
  Ox(13) = 6.711819639118E-01
  Oy(13) = 3.685670882907E-01
  w(13) = 3.901232700510E-02
  Ox(14) = 1.293388490485E-01
  Oy(14) = 7.447663982495E-02
  w(14) = 1.162080754372E-02
  Ox(15) = 8.335743322378E-01
  Oy(15) = 4.996274255819E-01
  w(15) = 4.642823955812E-02
  Ox(16) = 3.670788892962E-01
  Oy(16) = 2.519357740235E-01
  w(16) = 3.798783310581E-02
  Ox(17) = 5.909368760506E-01
  Oy(17) = 4.869502395267E-01
  w(17) = 4.130171453748E-02
  Ox(18) = 7.417637460141E-01
  Oy(18) = 6.279014865859E-01
  w(18) = 4.841339013884E-02
  Ox(19) = 6.279014865859E-01
  Oy(19) = 7.417637460141E-01
  w(19) = 4.841339013884E-02
  Ox(20) = 4.869502395267E-01
  Oy(20) = 5.909368760506E-01
  w(20) = 4.130171453748E-02
  Ox(21) = 2.519357740235E-01
  Oy(21) = 3.670788892962E-01
  w(21) = 3.798783310581E-02
  Ox(22) = 4.996274255819E-01
  Oy(22) = 8.335743322378E-01
  w(22) = 4.642823955812E-02
  Ox(23) = 7.447663982495E-02
  Oy(23) = 1.293388490485E-01
  w(23) = 1.162080754372E-02
  Ox(24) = 3.685670882907E-01
  Oy(24) = 6.711819639118E-01
  w(24) = 3.901232700510E-02
  Ox(25) = 3.673697853806E-01
  Oy(25) = 8.997294996538E-01
  w(25) = 4.244873302980E-02
  Ox(26) = 2.488983098158E-01
  Oy(26) = 7.241384940891E-01
  w(26) = 3.442681426024E-02
  Ox(27) = 1.196054590036E-01
  Oy(27) = 4.288508824476E-01
  w(27) = 2.995376809966E-02
  Ox(28) = 2.426233944222E-01
  Oy(28) = 9.410672772109E-01
  w(28) = 3.648716160597E-02
  Ox(29) = 1.419535016004E-01
  Oy(29) = 7.524467626583E-01
  w(29) = 2.759429759588E-02
  Ox(30) = 1.362124657777E-01
  Oy(30) = 9.622473153642E-01
  w(30) = 2.863542971348E-02
  Ox(31) = 1.670387759191E-02
  Oy(31) = 1.483114568272E-01
  w(31) = 6.404244616724E-03
  Ox(32) = 5.995074957044E-02
  Oy(32) = 7.633693960835E-01
  w(32) = 1.873220073879E-02
  Ox(33) = 5.695764868253E-02
  Oy(33) = 9.701698603928E-01
  w(33) = 1.913728513580E-02
  Ox(34) = 2.447911451942E-02
  Oy(34) = 4.445439440056E-01
  w(34) = 1.460888798152E-02
  Ox(35) = 1.160393058611E-02
  Oy(35) = 7.656319455497E-01
  w(35) = 8.352354145856E-03
  Ox(36) = 1.096881837272E-02
  Oy(36) = 9.717784813336E-01
  w(36) = 8.454511187252E-03

  do i=1, 36
    gamma = datan(Oy(i)/Ox(i))
    mu(i) = Ox(i) / dcos(gamma)
  enddo

  !Qudarent 2 (gamma + pi/2)
  do i=37, 72
    gamma = datan(Oy(i-36)/Ox(i-36)) + PI/2.0d0
    mu(i) = mu(i-36)
    Ox(i) = dsin(mu(i))*dcos(gamma)
    Oy(i) = dsin(mu(i))*dsin(gamma)
    w(i)  = w(i-36)
  enddo

  !Qudarent 3 (-x, -y)
  do i=73, 108
    gamma = datan(Oy(i-72)/Ox(i-72)) + PI/2.0d0
    mu(i) = mu(i-72)
    Ox(i) = dsin(mu(i))*dcos(gamma)
    Oy(i) = dsin(mu(i))*dsin(gamma)
    w(i)  = w(i-72)
  enddo

  !Qudarent 4 (+x, -y)
  do i=109, s
    gamma = datan(Oy(i-108)/Ox(i-108)) + PI/2.0d0
    mu(i) = mu(i-108)
    Ox(i) = dsin(mu(i))*dcos(gamma)
    Oy(i) = dsin(mu(i))*dsin(gamma)
    w(i)  = w(i-108)
  enddo

  !Normalize w
  wsum = 0.0d0
  do i=1, s
    wsum = wsum + w(i)
  enddo
  wnorm = (4*PI) / (wsum)
  w = wnorm * w

  end subroutine buildXS

END MODULE XSbuilder
