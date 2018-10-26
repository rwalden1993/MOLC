MODULE RayTrace
use XSbuilder
implicit none

real(8) :: ray_percision=10E-8

contains
  subroutine trace_line(x0, y0, m, numChar)
  use XSbuilder
  implicit none
  real(8) :: x0, y0, dx, dy, x, y, length, track_edge
  real(8) :: prev_x, prev_y
  integer :: numChar, m, my_cell, prev_cell, num_track
  real(8) :: gamma, cos_gamma

  num_track = 1
  length = 0.00d0
  gamma = datan(dabs(Oy(m)/Ox(m)))
  if (Ox(m) .gt. 0) then
     cos_gamma = dcos(gamma)
  else
    cos_gamma = -dcos(gamma)
  endif
  dx = cos_gamma  * dr
  dy = dsin(gamma) * dr
  x = x0
  y = y0
  my_cell = determine_cell(x,y)

  do
    prev_x = x
    prev_y = y
    x = x + dx
    y = y + dy
    prev_cell = my_cell
    my_cell = determine_cell(x,y)
    if (prev_cell .ne. my_cell) then
      !pusedo binary search for track_edge
      track_edge = bin_search(prev_x, prev_y, Ox(m), Oy(m))
      length = length + track_edge
      x = prev_x + cos_gamma*track_edge
      y = prev_y + dsin(gamma)*track_edge

      !store track
      ! write(*,*) prev_cell, " -> ", my_cell
      ! write(*,*) "prev_x = ", prev_x
      ! write(*,*) "prev_y = ", prev_y
      ! write(*,*) "x = ", x
      ! write(*,*) "y = ", y
      ! write(*,*) "length = ", length
      ! write(*,*) "track_edge = ", track_edge

      fp(prev_cell) = fp(prev_cell) + 2 * w(m) * length * dw !2x -> double for other direction
      char_lines(m, numChar, num_track) = prev_cell
      char_lengths(m, numChar, num_track) = length

      !reset variables
      num_track = num_track + 1
      !length = dr - track_edge
      length = 0.0d0
    else
      length = length + dr
    endif

    !if out of bounds, store and exit
    if(my_cell .eq. -1) then
      char_lines(m, numChar, num_track) = -1
      exit
    endif

  enddo

  end subroutine trace_line


  function determine_cell(x, y) result(cell)
    implicit none
    integer :: cell, pincell, radii, octant, pincell_x, pincell_y
    real(8) :: x, y, my_r, diag_1, diag_2
    real(8) :: pin_bot, pin_top, pin_left, pin_right, pin_mid_x, pin_mid_y

    !find pincel
    if (TEST_A) then
      pincell = 0
      pin_bot = 0.00d0
      pin_top = pin_width
      pin_left = 0.00d0
      pin_right = pin_width
      diag_1 = x
      diag_2 = pin_width - x
    else
      pincell_x = INT(x/pin_width)
      pincell_y = INT(y/pin_width)
      if (x .eq. width) pincell_x = 2
      if (y .eq. width) pincell_y = 2
      pincell = pincell_y*3 + pincell_x
      pin_left = pincell_x * pin_width
      pin_right = (pincell_x + 1) * pin_width
      pin_bot = pincell_y * pin_width
      pin_top = (pincell_y + 1) * pin_width
      diag_1  = x + (pincell_y - pincell_x) * pin_width
      diag_2  = ((pincell_x + pincell_y + 1) * pin_width) - x
    endif

    pin_mid_x = pin_left + pin_width/2.0d0
    pin_mid_y = pin_bot + pin_width/2.0d0

    !find radii
    my_r = sqrt((x-pin_mid_x)**2 + (y-pin_mid_y)**2)
    if (my_r .ge. r1) then
      radii = 0
    else if (my_r .lt. r1 .and. my_r .ge. r2) then
      radii = 1
    else if (my_r .lt. r2 .and. my_r .ge. r3) then
      radii = 2
    else if (my_r .lt. r3 .and. my_r .ge. r4) then
      radii = 3
    else if (my_r .lt. r4) then
      radii = 4
    endif

    !find octant
    if (x .ge. pin_left .and. x .le. pin_mid_x ) then !left side of pin cell
      if (y .ge. pin_bot .and. y .lt. diag_1) then
        octant = 7
      else if (y .ge. diag_1 .and. y .lt. pin_mid_y) then
        octant = 0
      else if (y .ge. pin_mid_y .and. y .lt. diag_2 ) then
        octant = 1
      else if (y .ge. diag_2 .and. y .le. pin_top) then
        octant = 2
      endif
    else !right side of pin cell
      if (y .ge. diag_1 .and. y .le. pin_top) then
        octant = 3
      else if (y .ge. pin_mid_y .and. y .lt. diag_1) then
        octant = 4
      else if (y .ge. diag_2 .and. y .lt. pin_mid_y) then
        octant = 5
      else if (y .ge. pin_bot .and. y .lt. diag_2 ) then
        octant = 6
      endif
    endif

    !find cell
    cell = (pincell * 40) + (radii * 8) + octant + 1

    if (x.lt.0 .or. x.gt.width .or. y.lt.0 .or. y.gt.width) cell = -1 !out of bounds

  end function determine_cell

  function bin_search(x, y, Ox, Oy) result(track_edge)
    implicit none
    real(8) :: x, y, Ox, Oy, track_edge, gamma
    real(8) :: my_dx, my_dy, my_dr, prev_x, prev_y
    integer :: prev_cell, cur_cell

    prev_cell = determine_cell(x,y)
    track_edge = 0.00d0
    my_dr = dr
    gamma = datan(dabs(Oy/Ox))

    do
      prev_x = x
      prev_y = y
      my_dr = my_dr/2.0d0
      if (Ox.gt.0) then
        my_dx = dcos(gamma) * my_dr
      else
        my_dx = -dcos(gamma) * my_dr
      endif
      my_dy = dsin(gamma) * my_dr
      x = x + my_dx
      y = y + my_dy
      cur_cell = determine_cell(x,y)
      if (cur_cell .ne. prev_cell) then
        x = prev_x
        y = prev_y
      else
        track_edge = track_edge + my_dr
      endif
      if (my_dr .lt. ray_percision) exit
    enddo
  end function bin_search

  subroutine get_intercept(x0, y0, dir, numChar)
    real(8) :: x0, y0, m, b, top_i, side_i
    integer :: dir, numChar

    char_x0(dir, numChar) = x0
    char_y0(dir, numChar) = y0
    m = Oy(dir)/Ox(dir)
    b = y0 - m * x0
    top_i = (width - b)/ m
    if (dir .le. s/4) then
      side_i = m * width + b
    else
      side_i = b
    endif
    if(top_i .le. width .and. top_i .gt. 0) then
      char_yf(dir,numChar) = width
      char_xf(dir,numChar) = top_i
    else
      if (dir .le. s/4) then
        char_xf(dir,numChar) = width
      else
        char_xf(dir, numChar) = 0.0d0
      endif
      char_yf(dir,numChar) = side_i
    endif
  end subroutine get_intercept
END MODULE RayTrace
