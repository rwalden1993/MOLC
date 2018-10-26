MODULE FluxMap
use XSbuilder
use RayTrace
implicit none

contains
  subroutine drawFluxMap(phii)
    integer :: ierror, my_cell, i
    real(8) :: x, y, dx, dy, map_res, phii(g,n)
    character(14) :: filename

    if (TEST_A) then
      map_res = 0.01d0 !cm
      filename = "flux_map_A.dat"
    else
      map_res = 0.0252d0 !cm
      filename = "flux_map_B.dat"
    endif

    open(unit=11, file=filename, status='replace', action='write', iostat=ierror)
    if (ierror.ne.0) then
      write(*,*) "An error occured when opening program1.out"
    endif

    110 format (f11.8,' ')
    
    x = 0.0d0
    y = 0.0d0
    dx = map_res
    dy = map_res

    do
      x = 0.0d0
      do
        write(11, 110, advance='no') x
        write(11, 110, advance='no') y
        my_cell = determine_cell(x,y)
        do i=1, g
          write(11, 110, advance='no') phii(i,my_cell)
        enddo
        write(11,*)
        x = x + dx
        if(x .gt. width) exit
      enddo
      y = y + dy
      if(y .gt. width) exit
    enddo

  end subroutine drawFluxMap

END MODULE FluxMap
