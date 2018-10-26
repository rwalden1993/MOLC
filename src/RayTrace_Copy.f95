MODULE RayTrace
use XSbuilder
implicit none

contains
  subroutine trace_line(x0, y0, m, numChar)
  use XSbuilder
  implicit none
  real(8) :: x0, y0, dx, dy, x, y, length, track_edge
  real(8) :: prev_x, prev_y
  integer :: numChar, m, my_cell, prev_cell, num_track, char_length_index

  num_track = 1
  length = 0.00d0
  dx = Ox(m) * dr
  dy = Oy(m) * dr
  x = x0
  y = y0
  my_cell = determine_cell(x,y)
  prev_cell = my_cell
  prev_x = x
  prev_y = y

  do
    ! write(*,*) "x = ", x
    ! write(*,*) "y = ", y
    if(my_cell .ne. 1000) prev_cell = my_cell
    my_cell = determine_cell(x,y)
    if (my_cell .ne. 1000 .and. prev_cell .ne. 1000 .and. prev_cell .ne. my_cell) then
      write(*,*) prev_cell, " -> ", my_cell
      !pusedo binary search for track_edge
      track_edge = bin_search(prev_x, prev_y, Ox(m), Oy(m))

      length = length + track_edge
      x = prev_x + Ox(m) * track_edge
      y = prev_y + Oy(m) * track_edge

      !store track
      fp(prev_cell) = fp(prev_cell) + 2 * w(m) * length * dw !2x -> double for other direction
      char_lines(m, numChar, num_track) = prev_cell
      char_length_index = (m-1) * s + (numChar-1) * max_lines + num_track
      char_lengths( char_length_index ) = length

      !reset variables
      num_track = num_track + 1
      length = 0.00d0
    else
      if(my_cell .ne. 1000) then
        prev_x = x
        prev_y = y
      endif
      x = x + dx
      y = y + dy
    endif
    !if out of bounds, store and exit
    if(my_cell .eq. -1) then
      char_lines(m, numChar, num_track) = -1
      exit
    endif
    length = length + dr
  enddo

  end subroutine trace_line

  function determine_cell(x, y) result(cell)
    implicit none
    integer :: cell, pinCell, radii, octant
    real(8) :: x, y, my_r, tol=10E-8
    real(8) :: pin_bot, pin_top, pin_left, pin_right, pin_mid_x, pin_mid_y
    logical :: glance

    glance = .false.

    !find pincel
    pinCell = 0
    pin_bot = 0.00d0
    pin_top = pin_width
    pin_left = 0.00d0
    pin_right = pin_width

    pin_mid_x = pin_left + pin_width/2
    pin_mid_y = pin_bot + pin_width/2

    !find radii
    my_r = sqrt((x-pin_mid_x)**2 + (y-pin_mid_y)**2)
    if (my_r .ge. r1+tol) then
      radii = 0
    else if (my_r .lt. r1-tol .and. my_r .gt. r2+tol) then
      radii = 1
    else if (my_r .lt. r2-tol .and. my_r .gt. r3+tol) then
      radii = 2
    else if (my_r .lt. r3-tol .and. my_r .gt. r4+tol) then
      radii = 3
    else if (my_r .lt. r4-tol) then
      radii = 4
    else
      glance = .true.
    endif

    !find octant
    if (x .ge. pin_left .and. x .lt. pin_mid_x-tol ) then !left side of pin cell
      if (y .ge. pin_bot .and. y .lt. x-tol) then
        octant = 7
      else if (y .gt. x+tol .and. y .lt. pin_mid_y-tol) then
        octant = 0
      else if (y .gt. pin_mid_y+tol .and. y .lt. (pin_width-x)-tol ) then
        octant = 1
      else if (y .gt. (pin_width-x)+tol .and. y .le. pin_top) then
        octant = 2
      else
        glance = .true.
      endif
    else if (x.gt.pin_mid_x+tol) then !right side of pin cell
      if (y .gt. x+tol .and. y .le. pin_top) then
        octant = 3
      else if (y .gt. pin_mid_y+tol .and. y .lt. x-tol) then
        octant = 4
      else if (y .gt. (pin_width-x)+tol .and. y .lt. pin_mid_y-tol) then
        octant = 5
      else if (y .ge. pin_bot .and. y .lt. (pin_width-x)-tol ) then
        octant = 6
      else
        glance = .true.
      endif
    else
      glance = .true.
    endif

    !find cell
    if (glance) then
      cell = 1000
    else
      cell = pincell * (4*8) + radii * 8 + octant + 1
    endif

    if (x.lt.0 .or. x.gt.width .or. y.lt.0 .or. y.gt.width) cell = -1 !out of bounds


  end function determine_cell

  function determine_cell_exact(x, y) result(cell)
    implicit none
    integer :: cell, pinCell, radii, octant
    real(8) :: x, y, my_r
    real(8) :: pin_bot, pin_top, pin_left, pin_right, pin_mid_x, pin_mid_y

    !find pincel
    pinCell = 0
    pin_bot = 0.00d0
    pin_top = pin_width
    pin_left = 0.00d0
    pin_right = pin_width

    pin_mid_x = pin_left + pin_width/2
    pin_mid_y = pin_bot + pin_width/2

    !find radii
    my_r = sqrt((x-pin_mid_x)**2 + (y-pin_mid_y)**2)
    if (my_r .ge. r1) then
      radii = 0
    else if (my_r .lt. r1 .and. my_r .gt. r2) then
      radii = 1
    else if (my_r .lt. r2 .and. my_r .gt. r3) then
      radii = 2
    else if (my_r .lt. r3 .and. my_r .gt. r4) then
      radii = 3
    else if (my_r .lt. r4) then
      radii = 4
    endif

    !find octant
    if (x .ge. pin_left .and. x .lt. pin_mid_x ) then !left side of pin cell
      if (y .ge. pin_bot .and. y .lt. x) then
        octant = 7
      else if (y .gt. x .and. y .lt. pin_mid_y) then
        octant = 0
      else if (y .gt. pin_mid_y .and. y .lt. (pin_width-x) ) then
        octant = 1
      else if (y .gt. (pin_width-x) .and. y .le. pin_top) then
        octant = 2
      endif
    else if (x.gt.pin_mid_x) then !right side of pin cell
      if (y .gt. x .and. y .le. pin_top) then
        octant = 3
      else if (y .gt. pin_mid_y .and. y .lt. x) then
        octant = 4
      else if (y .gt. (pin_width-x) .and. y .lt. pin_mid_y) then
        octant = 5
      else if (y .ge. pin_bot .and. y .lt. (pin_width-x) ) then
        octant = 6
      endif
    endif

    !find cell
    cell = pincell * (4*8) + radii * 8 + octant + 1

    if (x.lt.0 .or. x.gt.width .or. y.lt.0 .or. y.gt.width) cell = -1 !out of bounds


  end function determine_cell_exact

  function bin_search(x, y, Ox, Oy) result(track_edge)
    implicit none
    real(8) :: x, y, Ox, Oy, track_edge
    real(8) :: my_dx, my_dy, my_dr, prev_x, prev_y
    integer :: prev_cell, cur_cell

    prev_cell = determine_cell(x,y)
    track_edge = 0.00d0
    my_dr = dr

    do
      prev_x = x
      prev_y = y
      my_dr = my_dr/2.0d0
      my_dx = Ox * my_dr
      my_dy = Oy * my_dr
      x = x + my_dx
      y = y + my_dx
      cur_cell = determine_cell_exact(x,y)
      if (cur_cell .ne. prev_cell) then
        x = prev_x
        y = prev_y
      else
        track_edge = track_edge + my_dr
      endif
      if (my_dr .lt. 10.0d0**(-8)) exit
    enddo
  end function bin_search
END MODULE RayTrace
