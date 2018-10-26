program FinalProgram
use XSbuilder
use RayTrace
use FluxMap
implicit none


!variables
integer :: maxite = 1000
integer :: i, j, ite, ierror, ii, jj, kk
integer :: numChar, my_cell
real(8) :: Vx, Vy, normV, dx, dy, x0, y0, gamma, gamma2
real(8) :: tol, criteria
real(8) :: phisum, normfactor, top, bot, sigs_total
real(8),allocatable :: phii(:,:), phiprev(:,:)
real(8),allocatable :: qg(:,:), qup(:,:), qf(:,:), qdown(:,:)
real(8),allocatable :: diff(:), rho(:), k(:), diffk(:), rhok(:), avgModPhi(:), avgFuelPhi(:)
logical :: converged, convergedphi, convergedk

call buildXS()

allocate(diff(maxite), rho(maxite), k(0:maxite), diffk(maxite), rhok(maxite))
allocate(phii(g,n), phiprev(g,n), &
        qg(g,n), qup(g,n), qdown(g,n), qf(g,n), LR(g,n) )

open(unit=11, file="final_program.out", status='replace', action='write', iostat=ierror)
if (ierror.ne.0) then
  write(*,*) "An error occured when opening program2.out"
endif

100 format (A , ', ')
101 format (A, I2, ', ')
110 format (f11.8,', ')
111 format (I2,' -> ')

!Define constants for Test 1
k(0) = 1.000
tol = 10.0d0**(-5)
converged = .false.

!Initialize arrays
phii = 1.0 !guess phii = 1.0 everywhere
diff(1)  = 1.0      !first difference is assumed 1
ite = 0

!Ray Tracking

do i=1, s/2
  !determine perpindicular norm vector
  gamma = datan(Oy(i)/Ox(i))
  Vx = dcos(gamma)
  Vy = -( dcos(gamma)**2 / dsin(gamma) )
  normV = SQRT(Vx**2 + Vy**2)
  Vx = Vx / normV
  Vy = Vy / normV

  !determine dx and dy for line spacing
  gamma2 = datan(Vy/Vx)
  dy = dw / dsin(gamma2)
  dx = dw / dcos(gamma2)

  ! write(*,*) "Dir = ", i
  ! write(*,*) "dx = ", dx
  ! write(*,*) "dy = ", dy

  !Draw rays
  numChar = 1
  !bottom boundary
  y0 = 0.0d0
  x0 = 0.0d0
  do
    !ray trace
    ! write(*,*) "Dir = ", i
    ! write(*,*) "Char_line = ", numChar
    ! write(*,*) "x0 = ", x0
    call get_intercept(x0, y0, i, numChar)
    call trace_line(x0, y0, i, numChar)
    numChar = numChar + 1

    x0 = x0 + abs(dx)
    if(x0.GT.width) exit
  enddo

  !side boundary
  if (i .le. s/4) then !positive x
    x0 = 0.0d0
  else              !negative x
    x0 = width
  endif
  y0 = abs(dy) !already did corner line
  do
    !ray trace
    ! write(*,*) "Dir = ", i
    ! write(*,*) "Char_line = ", numChar
    ! write(*,*) "y0 = ", y0
    call get_intercept(x0, y0, i, numChar)
    call trace_line(x0, y0, i, numChar)
    numChar = numChar + 1

    y0 = y0 + abs(dy)
    if(y0.GT.width) exit
  enddo
  char_lines(i,numChar,1) = -1

enddo

!reflcet lines
char_lines(s/2+1:s,:,:) = char_lines(1:s/2,:,:)
char_lengths(s/2+1:s,:,:) = char_lengths(1:s/2,:,:)
char_x0(s/2+1:s,:) = char_xf(1:s/2,:)
char_xf(s/2+1:s,:) = char_x0(1:s/2,:)
char_y0(s/2+1:s,:) = char_yf(1:s/2,:)
char_yf(s/2+1:s,:) = char_y0(1:s/2,:)

! do i=1,s
!   write(*,*) "dir = ", i
!   do j=1, max_lines+1
!     my_cell = char_lines(i, j, 1)
!     if (my_cell .eq. -1) exit
!     do ii=1, n+1
!       my_cell = char_lines(i, j, ii)
!       if (my_cell .eq. -1) exit
!       write(*,111, advance='no') my_cell
!     enddo
!     write(*,*)
!   enddo
!   write(*,*)
! enddo

!correct fp and d_length*
do i=1,n
  fp(i) = ( 4*PI * Ap(i) ) / fp(i)
enddo
do i=1,s
  do j=1, max_lines+1
    do kk=1, n+1
      my_cell = char_lines(i, j, kk)
      if (my_cell .ne. -1) then
        char_lengths(i,j,kk) = char_lengths(i,j,kk) * fp(my_cell)
      endif
    enddo
  enddo
enddo

! Calculate tau and alpha
do kk=1, g !group
  do j=1, s !direction
    do i=1, max_lines+1 !line
      do ii =1, n+1 !track
        my_cell = char_lines(j, i, ii)
        if (my_cell .eq. -1) exit
        tau(kk,j,i,ii) = (sigt(kk,my_cell) * char_lengths(j,i,ii)) / mu(j)
        e_tau(kk,j,i,ii) = dexp(-tau(kk,j,i,ii))
        if (tau(kk,j,i,ii) .lt. 1E-6) then ! taylor expansion for alpha
          alpha(kk,j,i,ii) = 0.50d0 - tau(kk,j,i,ii)/12.0d0
        else
          alpha(kk,j,i,ii) = (1.0d0/tau(kk,j,i,ii)) - ( e_tau(kk,j,i,ii) / ( 1.0d0 - e_tau(kk,j,i,ii) ) )
        endif
        ! if(ISNAN(alpha(kk,j,i,ii)) .or. ISNAN(tau(kk,j,i,ii)) .or. tau(kk,j,i,ii) .eq. 0) then
        !   write(*,*) "Dir, Line, Track = ", j, i, ii
        !   write(*,*) "Cell  = ", my_cell
        !   write(*,*) "Length= ", char_lengths(char_length_index)
        !   write(*,*) "Tau   = ", tau(kk,j,i,ii)
        !   write(*,*) "e_Tau = ", e_tau(kk,j,i,ii)
        !   write(*,*) "Alpha = ", alpha(kk,j,i,ii)
        !   goto 9999
        ! endif
      enddo
    enddo
  enddo
enddo

!Solve loop
do
  ite = ite+1
  !qup
  qup = 0.0d0
  do i=1, g
    if (i.ne.g) then
      do ii = i+1, g
        do j=1, n
          qup(i, j) = qup(i,j) + sigs(ii,i,j)*phii(ii,j)
        enddo
      enddo
    endif
  enddo

  !qf
  qf = 0.0d0
  do i=1, g
    do ii=1,g
      do j=1, n
        qf(i, j) = qf(i,j) + (chi(i,j)/k(ite-1))*(nuf(ii,j)*sigf(ii,j)*phii(ii,j))
      enddo
    enddo
  enddo

  !save phii to phiprev
  phiprev(:,:) = phii(:,:)

  do i=1, g
    !qdown
    qdown = 0.
    if (i.ne.1) then
      do ii = 1, i-1
        do j=1, n
          qdown(i, j) = qdown(i,j) + sigs(ii,i,j)*phii(ii,j)
        enddo
      enddo
    endif

    !qtotal
    do j=1,n
      qg(i,j) = qup(i,j) + qf(i,j) + qdown(i,j)
    enddo

    !transport solve
    CALL transport_solve(i, qg, phii)
  enddo

  !calculate k
  ApStar = 0.0d0
  do j=1, s
    do i=1, max_lines+1
      do ii=1, n+1
        my_cell = char_lines(j,i,ii)
        if (my_cell .eq. -1) exit
        ApStar(my_cell) = ApStar(my_cell) + ((1.0d0/(4*PI)) * w(j) * char_lengths(j,i,ii) * dw )
      enddo
    enddo
  enddo

  top = 0.0d0
  bot = 0.0d0
  do j=1, g
    do i=1, n
      top = top + nuf(j,i)*sigf(j,i)*phii(j,i)*ApStar(i)
      sigs_total = 0.0d0
      do jj=1, g
        sigs_total = sigs_total + sigs(j,jj,i)
      enddo
      bot = bot + (sigt(j,i)-sigs_total)*phii(j,i)*ApStar(i)
      !bot = bot + (sigt(j,i)-sigs_total)*phii(j,i)*ApStar(i) + LR(j,i)
      if (sigt(j,i)-sigs_total .lt. 0) then
        write(*,*) "Negative abs xs: ", sigt(j,i)-sigs_total
        write(*,*) "At group : ", j, " and cell ", i
      endif
    enddo
  enddo
  k(ite) = (top/bot)

  !Normalize scalar flux
  phisum = 0.00d0 !reset
  do j=1,g
    do i=1, n
      phisum = phisum + phii(j,i)*ApStar(i)
    enddo
  enddo

  normfactor = 1.0d0 / (phisum)

  do j=1, g
    do i=1, n
      phii(j,i) = phii(j,i) * normfactor
    enddo
  enddo

  !check convergence
  if (ite.gt.1) then
    !converged k
    diffk(ite) = k(ite) - k(ite-1)
    rhok(ite)  = diffk(ite) / diffk(ite-1)
    criteria   = tol * (1.0 - rhok(ite))
    convergedk = diffk(ite).lt.criteria

    !converged phi
    diff(ite) = maxval(abs(phii-phiprev))
    rho(ite) = diff(ite)/diff(ite-1)
    criteria = tol * (1 - rho(ite))
    convergedphi = diff(ite).lt.criteria

    converged = convergedphi.and.convergedk
  endif

  if( converged ) then
    write(*,*) 'Converged on iteration ', ite
    write(*,*) 'Eigenvalue: ', k(ite)
    exit
  endif
  if( ite.ge.maxite .or. ISNAN(diff(ite))) then
    write(*,*) 'MAX ITERATIONS REACHED, NOT CONVERGED'
    exit
  endif
  if( k(ite) .lt. 0 ) then
    write(*,*) 'Eigenvalue less than zero: ', k(ite)
    exit
  endif
enddo

!write edits

write(11, 100, advance='no') 'ite'
do i=1, ite
  write(11, '(I3, A)', advance='no') i, ', '
enddo
write(11,*)

write(11, 100, advance='no') 'diff'
do i=1, ite
  write(11, 110, advance='no') diff(i)
enddo
write(11,*)

write(11, 100, advance='no') 'rho phi'
do i=1, ite
  write(11, 110, advance='no') rho(i)
enddo
write(11,*)

write(11, 100, advance='no') 'diff k'
do i=1, ite
  write(11, 110, advance='no') diffk(i)
enddo
write(11,*)

write(11, 100, advance='no') 'rho k'
do i=1, ite
  write(11, 110, advance='no') rhok(i)
enddo
write(11,*)

write(11, 100, advance='no') 'k'
do i=1, ite
  write(11, 110, advance='no') k(i)
enddo
write(11,*)

! fuel pin average mod and fuel scalar flux
if (TEST_A) then
  allocate(avgModPhi(1), avgFuelPhi(1))
  do j=1, g
    avgModPhi(1) = 0.0d0
    avgFuelPhi(1) = 0.0d0
    write(11,101,advance='no') "phii ", j
    do i=1, n
      if (mod(i,40) .le. 8) then !mod region
        avgModPhi(1) = avgModPhi(1) + phii(j,i)
      else
        avgFuelPhi(1) = avgFuelPhi(1) + phii(j,i)
      endif
    enddo
    write(11,110,advance='no') avgModPhi(1)/8
    write(11,110,advance='no') avgFuelPhi(1)/32
    write(11,*)
  enddo
else
  allocate(avgModPhi(9), avgFuelPhi(9))
  do j=1, g
    avgModPhi = 0.0d0
    avgFuelPhi = 0.0d0
    write(11,101,advance='no') "phii ", j
    do i=1, n
      ii = (i-1)/40 + 1
      if (mod(i,40) .le. 8) then !mod region
        avgModPhi(ii) = avgModPhi(ii) + phii(j,i)
      else
        avgFuelPhi(ii) = avgFuelPhi(ii) + phii(j,i)
      endif
    enddo
    do i=1,9
      write(11,110,advance='no') avgModPhi(i)/8
      write(11,110,advance='no') avgFuelPhi(i)/32
    enddo
    write(11,*)
  enddo
endif

call drawFluxMap(phii)

9999 end program FinalProgram
