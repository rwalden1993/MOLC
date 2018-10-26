SUBROUTINE transport_solve(my_g, qg, phii)
use XSbuilder
implicit none

integer:: my_g
real(8):: qg(g,n), phii(g,n)

integer :: maxite = 1000
integer :: i, j, ii, ite, ierror, my_cell, end_line, reflect_dir
logical :: convergedphi
real(8) :: my_x0, my_y0, psi_right, psi_left, min_left, min_right, diff_left, diff_right, x_left, x_right
real(8) :: tol, criteria
real(8),allocatable :: qi(:), phiprev(:), phisum_top(:), phisum_bot(:), LR_inner(:)
real(8),allocatable:: diff(:), rho(:)
real(8),allocatable :: psi(:,:,:), psii(:,:,:), psi_edge(:,:), psi_end(:,:)

allocate(qi(n), phiprev(n), phisum_top(n), phisum_bot(n), LR_inner(n))
allocate(diff(maxite), rho(maxite))
allocate(psi(s,max_lines+1,n+1), psii(s,max_lines+1,n+1), psi_edge(s,max_lines+1), psi_end(s,max_lines+1))


tol = 10.0d0**(-6)
ite = 0
diff(1) = 1.0
psi_end = 0.0d0


do
  ite = ite +1
  !update source (qi)
  do i=1, n
    qi(i) = (1.0d0/(4*PI))*((sigs(my_g,my_g,i))*phii(my_g,i) + qg(my_g,i))
  enddo

  !calculate reflective BC
  do j=1, s
    do i=1,max_lines+1
      my_cell = char_lines(j,i,1)
      if (my_cell .eq. -1) exit
      my_x0 = char_x0(j,i)
      my_y0 = char_y0(j,i)
      min_left = width
      min_right = width
      psi_left = 0.0d0
      psi_right = 0.0d0
      x_left = 0.0d0
      x_right = width
      if (my_x0 .eq. 0 .or. my_x0 .eq. width) then !bounce off left or right side
        reflect_dir = mod(abs(72 - j), 144) + 1
        do ii=1, max_lines+1
          my_cell = char_lines(j,ii,1)
          if (my_cell .eq. -1) exit
          diff_right = char_yf(reflect_dir, ii) - my_y0
          diff_left  = my_y0 - char_yf(reflect_dir, ii)
          if (diff_left .lt. min_left .and. diff_left .ge. 0.0d0) then
            min_left = diff_left
            x_left = char_yf(reflect_dir, ii)
            psi_left = psi_end(reflect_dir, ii)
          endif
          if (diff_right .lt. min_right .and. diff_right .gt. 0.0d0) then
            min_right = diff_right
            x_right = char_yf(reflect_dir, ii)
            psi_right = psi_end(reflect_dir, ii)
          endif
        enddo
        psi_edge(j,i) = ((my_y0 - x_left) / (x_right - x_left)) * (psi_right - psi_left) + psi_left
      else                                         !bounce off top or bottom side
        reflect_dir = mod(abs(144 - j), 144) + 1
        do ii=1, max_lines+1
          my_cell = char_lines(j,ii,1)
          if (my_cell .eq. -1) exit
          diff_right = char_xf(reflect_dir, ii) - my_x0
          diff_left  = my_x0 - char_xf(reflect_dir, ii)
          if (diff_left .lt. min_left .and. diff_left .ge. 0) then
            min_left = diff_left
            x_left = char_xf(reflect_dir, ii)
            psi_left = psi_end(reflect_dir, ii)
          endif
          if (diff_right .lt. min_right .and. diff_right .gt. 0) then
            min_right = diff_right
            x_right = char_xf(reflect_dir, ii)
            psi_right = psi_end(reflect_dir, ii)
          endif
        enddo
        psi_edge(j,i) = ((my_x0 - x_left) / (x_right - x_left)) * (psi_right - psi_left) + psi_left
      endif
    enddo
  enddo

  !calculate psi m
  do j=1, s
    do i=1, max_lines+1 !line
      my_cell = char_lines(j, i, 1)
      if( my_cell .eq. -1) exit
      if (j .le. 72) then !top directions
        psi(j, i, 1) = psi_edge(j, i) !reflective
        !psi(j,i,1) = 0.0d0 !vaccum
        do ii=2, n+1 !track
          my_cell = char_lines(j, i, ii-1)
          if (my_cell .eq. -1) exit
          psi(j,i,ii) = ( psi(j,i,ii-1)*e_tau(my_g,j,i,ii-1) ) + &
                        ( (qi(my_cell)/sigt(my_g,my_cell))*(1-e_tau(my_g,j,i,ii-1)) )
        enddo
        psi_end(j,i) =  psi(j,i,ii-1)
      else !bot directions
        do ii=1, n+1 !find track end
          my_cell = char_lines(j, i, ii)
          if ( my_cell .eq. -1 ) then
            end_line = ii
            exit
          endif
        enddo
        psi(j, i, end_line) = psi_edge(j, i) !reflective
        !psi(j,i,end_line) = 0.0d0 !vaccum
        do ii=end_line-1, 1, -1 !track
          my_cell = char_lines(j, i, ii)
          psi(j,i,ii) = ( psi(j,i,ii+1)*e_tau(my_g,j,i,ii) ) + &
                        ( (qi(my_cell)/sigt(my_g,my_cell))*(1-e_tau(my_g,j,i,ii)) )
        enddo
        psi_end(j,i) =  psi(j,i,1)
      endif
    enddo
  enddo


  !calculate cell psi

  do j=1, s
    do i=1, max_lines+1 !line
      do ii=1, n+1 !track
        my_cell = char_lines(j, i, ii)
        if (my_cell .eq. -1) exit
        if (j .le. 72) then
          psii(j,i,ii) = alpha(my_g,j,i,ii)*psi(j,i,ii) + (1 - alpha(my_g,j,i,ii))*psi(j,i,ii+1)
        else
          psii(j,i,ii) = alpha(my_g,j,i,ii)*psi(j,i,ii+1) + (1 - alpha(my_g,j,i,ii))*psi(j,i,ii)
        endif
      enddo
    enddo
  enddo


  !calculate cell average phii and LR
  phiprev(:) = phii(my_g,:) !save to phiprev
  phii(my_g,:) = 0.0d0 !reset array for summation
  phisum_top(:) = 0.0d0
  phisum_bot(:) = 0.0d0
  LR(my_g,:) = 0.0d0
  LR_inner(:) = 0.0d0
  do j=1, s
    do i=1, max_lines+1 !line
      do ii=1, n+1 !track
        my_cell = char_lines(j, i, ii)
        if (my_cell .eq. -1) exit
        phisum_top(my_cell) = phisum_top(my_cell) + &
                              ( w(j) * psii(j,i,ii) * char_lengths(j,i,ii) * dw )
        phisum_bot(my_cell) = phisum_bot(my_cell) + &
                              ((1.0d0/(4*PI)) * w(j) * char_lengths(j,i,ii) * dw )
        if (j .le. 72) then
          LR_inner(my_cell) = LR_inner(my_cell) + dw * (psi(j,i,ii+1)-psi(j,i,ii))
        else
          LR_inner(my_cell) = LR_inner(my_cell) + dw * (psi(j,i,ii)-psi(j,i,ii+1))
        endif
      enddo
    enddo
    LR(my_g,:) = LR(my_g,:) + (w(j) * mu(j) * LR_inner(:))
  enddo
  phii(my_g,:) = phisum_top(:) / phisum_bot(:)


  !converged phii
  if(ite.gt.1) then
    diff(ite) = maxval(abs(phii(my_g,:)-phiprev))
    rho(ite) = diff(ite)/diff(ite-1)
    criteria = tol * (1.0 - rho(ite))
    convergedphi = diff(ite).lt.criteria
  endif

  !write(*,*) "Inner loop diff ", diff(ite)

  if(convergedphi) exit
  if(ite.ge.maxite) then
    write(*,*) "TRANSPORT EQUATION NOT CONVERGED"
    exit
  endif
enddo

end subroutine transport_solve
