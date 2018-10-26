program RayTraceTest
use RayTrace
use XSbuilder
implicit none

real(8) :: x, y, expected, act, diff, tol=10E-8
integer :: i, j, kk, ii, end_line, char_length_index
logical :: test_bin_srch, test_ray_trace, test_intercept

call buildXS()

test_bin_srch = .false.
test_ray_trace = .false.
test_intercept = .false.

!bin_search tests
if (test_bin_srch) then
  !test 1
  x = width/2 - 10.0d0**(-5)
  y = 0.045
  expected = 10.0d0**(-5)
  act = bin_search(x,y,1.0d0,0.00d0)
  diff = abs(expected - act)

  if (diff .lt. tol) then
    write(*,*) "TEST 1 PASS"
  else
    write(*,*) "TEST 1 FAIL"
  endif

  write(*,*) "Actual   = ", act
  write(*,*) "Expected = ", expected
  write(*,*) "diff     = ", diff

  !bin_search test 2
  x = width/2 - 2.75E-4
  y = 0.045
  expected = 2.75E-4
  act = bin_search(x,y,1.0d0,0.00d0)
  diff = abs(expected - act)

  if (diff .lt. tol) then
    write(*,*) "TEST 2 PASS"
  else
    write(*,*) "TEST 2 FAIL"
  endif

  write(*,*) "Actual   = ", act
  write(*,*) "Expected = ", expected
  write(*,*) "diff     = ", diff

  !bin_search test 3
  x = width/2 - 7.5E-4
  y = 0.045
  expected = 7.5E-4
  act = bin_search(x,y,1.0d0,0.00d0)
  diff = abs(expected - act)

  if (diff .lt. tol) then
    write(*,*) "TEST 3 PASS"
  else
    write(*,*) "TEST 3 FAIL"
  endif

  write(*,*) "Actual   = ", act
  write(*,*) "Expected = ", expected
  write(*,*) "diff     = ", diff

  !bin_search test 4
  x = width/2 - 5.0E-4
  y = 0.045
  expected = 5.0E-4
  act = bin_search(x,y,1.0d0,0.00d0)
  diff = abs(expected - act)

  if (diff .lt. tol) then
    write(*,*) "TEST 4 PASS"
  else
    write(*,*) "TEST 4 FAIL"
  endif

  write(*,*) "Actual   = ", act
  write(*,*) "Expected = ", expected
  write(*,*) "diff     = ", diff

  !bin_search test 5
  x = width/2 - 10E-4
  y = 0.045
  expected = 10E-4
  act = bin_search(x,y,1.0d0,0.00d0)
  diff = abs(expected - act)

  if (diff .lt. tol) then
    write(*,*) "TEST 5 PASS"
  else
    write(*,*) "TEST 5 FAIL"
  endif

  write(*,*) "Actual   = ", act
  write(*,*) "Expected = ", expected
  write(*,*) "diff     = ", diff

  !bin_search test 6
  x = width/2 - 10E-8
  y = 0.045
  expected = 10E-8
  act = bin_search(x,y,1.0d0,0.00d0)
  diff = abs(expected - act)

  if (diff .lt. tol) then
    write(*,*) "TEST 6 PASS"
  else
    write(*,*) "TEST 6 FAIL"
  endif

  write(*,*) "Actual   = ", act
  write(*,*) "Expected = ", expected
  write(*,*) "diff     = ", diff

  !bin_search test 7
  x = width/2
  y = 0.045
  expected = 0.0d0
  act = bin_search(x,y,1.0d0,0.00d0)
  diff = abs(expected - act)

  if (diff .lt. tol) then
    write(*,*) "TEST 7 PASS"
  else
    write(*,*) "TEST 7 FAIL"
  endif

  write(*,*) "Actual   = ", act
  write(*,*) "Expected = ", expected
  write(*,*) "diff     = ", diff

  !bin_search test 8, boundary cross
  x = width - 2.75E-4
  y = 0.045
  expected = 2.75E-4
  act = bin_search(x,y,1.0d0,0.00d0)
  diff = abs(expected - act)

  if (diff .lt. tol) then
    write(*,*) "TEST 8 PASS"
  else
    write(*,*) "TEST 8 FAIL"
  endif

  write(*,*) "Actual   = ", act
  write(*,*) "Expected = ", expected
  write(*,*) "diff     = ", diff

  !bin_search test 9, right to left
  x = width/2 + 2.75E-4
  y = 0.045
  expected = 2.75E-4
  act = bin_search(x,y,-1.0d0,0.00d0)
  diff = abs(expected - act)

  if (diff .lt. tol) then
    write(*,*) "TEST 9 PASS"
  else
    write(*,*) "TEST 9 FAIL"
  endif

  write(*,*) "Actual   = ", act
  write(*,*) "Expected = ", expected
  write(*,*) "diff     = ", diff

  !bin_search test 10, up
  y = width/2 - 2.75E-4
  x = 0.042901
  expected = 2.75E-4 / sin(atan(abs(Oy(36)/Ox(36))))
  act = bin_search(x,y,Ox(36),Oy(36))
  diff = abs(expected - act)

  if (diff .lt. tol) then
    write(*,*) "TEST 10 PASS"
  else
    write(*,*) "TEST 10 FAIL"
  endif

  write(*,*) "Actual   = ", act
  write(*,*) "Expected = ", expected
  write(*,*) "diff     = ", diff
endif

!Ray trace tests
if(test_ray_trace) then
  write(*,*)
  write(*,*) "Test 1"
  x = 0.0
  y = 0.0
  i = 1
  j = 1
  call trace_line(x, y, i, j)
  write(*,*) "Expected: 1 -> 8 -> 7 -> 6 -> -1"

  write(*,*)
  write(*,*) "Test 2"
  y=width-0.045
  call trace_line(x, y, i, j)
  write(*,*) "Expected: 2 -> 3 -> 4 -> 5 -> -1"

  write(*,*)
  write(*,*) "Test 4"
  y=width/2+10E-4
  call trace_line(x, y, i, j)
  write(*,*) "Expected: 2 -> 10 -> 18 -> 26 -> 34 -> 35 -> 36 -> 37 -> 29 -> 21 -> 13 -> 5 -> -1"

  write(*,*)
  write(*,*) "Test 5"
  y = 0
  i = 18
  call trace_line(x, y, i, j)

  write(*,*)
  write(*,*) "Test 6"
  y = 0
  i = 19
  call trace_line(x, y, i, j)

  !glancing
  write(*,*)
  write(*,*) "Test 7"
  y = 7.1271459988169431E-002
  i = 1
  call trace_line(x, y, i, j)


  write(*,*)
  write(*,*) "Test 8"
  y = 0.35633907594447795
  i = 2
  call trace_line(x, y, i, j)

  write(*,*)
  write(*,*) "Test 9"
  y = 0.64152486270927545
  x = 0
  i = 1
  call trace_line(x, y, i, j)

  write(*,*)
  write(*,*) "Test 10"
  y = 0.64152486270927545
  x = width
  i = 37
  call trace_line(x, y, i, j)

  write(*,*)
  write(*,*) "Test 11"
  y = 0.38752468385781114
  x = width
  i = 37
  call trace_line(x, y, i, j)

  write(*,*)
  write(*,*) "Test 12"
  y = 0
  x = 0.27085788590400617
  i = 54
  call trace_line(x, y, i, j)

  write(*,*)
  write(*,*) "Test 13"
  y = 0.59408100896034999
  x = 0
  i = 12
  j = 1
  call trace_line(x, y, i, j)

  write(*,*)
  write(*,*) "Test 14"
  y = 0.59408100896034999
  x = width
  i = 48
  j = 1
  call trace_line(x, y, i, j)

  write(*,*)
  write(*,*) "Test 15"
  y = 0
  x = 0
  i = 36
  call trace_line(x, y, i, j)

  write(*,*)
  write(*,*) "Test 16 - Backwards"
  do ii=1, n+1 !find track end
    if ( char_lines(i, j, ii) .eq. -1 ) then
      end_line = ii
      exit
    endif
  enddo
  do ii=end_line-1, 1, -1
    write(*,*) char_lines(i, j, ii)
    write(*,*) char_lengths(i,j,ii)
  enddo
endif

if (test_intercept) then

  write(*,*) "Test 1 - side intercept, forwards"
  x = 0.0d0
  y = 1.0d0
  i = 1
  j = 1
  call get_intercept(x, y, i, j)
  write(*,*) "act. xf = ", char_xf(i,j)
  write(*,*) "exp. xf = ", 1.26d0
  write(*,*) "act. yf = ", char_yf(i,j)
  write(*,*) "exp. yf = ", 1.01422d0

  write(*,*)
  write(*,*) "Test 2 - side intercept, backwards"
  x = width
  y = 1.0d0
  i = 72
  j = 1
  call get_intercept(x, y, i, j)
  write(*,*) "act. xf = ", char_xf(i,j)
  write(*,*) "exp. xf = ", 0.0d0
  write(*,*) "act. yf = ", char_yf(i,j)
  write(*,*) "exp. yf = ", 1.01422d0

  write(*,*)
  write(*,*) "Test 3 - top intercept, forewards"
  x = 1.0d0
  y = 0
  i = 36
  j = 1
  call get_intercept(x, y, i, j)
  write(*,*) "act. xf = ", char_xf(i,j)
  write(*,*) "exp. xf = ", 1.01422d0
  write(*,*) "act. yf = ", char_yf(i,j)
  write(*,*) "exp. yf = ", width

  write(*,*)
  write(*,*) "Test 4 - top intercept, backwards"
  x = 1.0d0
  y = 0
  i = 37
  j = 1
  call get_intercept(x, y, i, j)
  write(*,*) "act. xf = ", char_xf(i,j)
  write(*,*) "exp. xf = ", 1 - 0.01422d0
  write(*,*) "act. yf = ", char_yf(i,j)
  write(*,*) "exp. yf = ", width

  write(*,*)
  write(*,*) "Test 5 - side intercept, forwards"
  x = width/2
  y = 0.0d0
  i = 1
  j = 1
  call get_intercept(x, y, i, j)
  write(*,*) "act. xf = ", char_xf(i,j)
  write(*,*) "exp. xf = ", width
  write(*,*) "act. yf = ", char_yf(i,j)
  write(*,*) "exp. yf = ", 7.1110401042564662E-003

  write(*,*)
  write(*,*) "Test 6 - side intercept, backwards"
  x = width/2
  y = 0.0d0
  i = 72
  j = 1
  call get_intercept(x, y, i, j)
  write(*,*) "act. xf = ", char_xf(i,j)
  write(*,*) "exp. xf = ", 0.0d0
  write(*,*) "act. yf = ", char_yf(i,j)
  write(*,*) "exp. yf = ", 7.1110401042564662E-003
endif

end program RayTraceTest
