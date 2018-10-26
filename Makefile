GCC = gfortran
FFLAGS = -Wall -fbacktrace -fbounds-check
LIBDIR_LAPACK = /usr/lib/liblapack
LIBDIR_BLAS = /usr/lib/libblas
LIBS = -L$(LIBDIR_LAPACK) -L$(LIBDIR_BLAS) -llapack -lblas

a: FinalProgramA
b: FinalProgramB
test: RayTraceTest

FinalProgramA: NE723-FinalProgram.o TransportSolver.o XSbuilder.o RayTrace.o FluxMap.o
	$(GCC) $(FFLAGS) NE723-FinalProgram.o TransportSolver.o XSbuilder.o RayTrace.o FluxMap.o -o $@

FinalProgramB: NE723-FinalProgram.o TransportSolver.o XSbuilderB.o
	$(GCC) $(FFLAGS) NE723-FinalProgram.o TransportSolver.o XSbuilderB.o -o $@

RayTraceTest: RayTraceTest.o XSbuilder.o RayTrace.o
	$(GCC) $(FFLAGS) RayTraceTest.o XSbuilder.o RayTrace.o -o $@

NE723-FinalProgram.o: NE723-FinalProgram.f95 TransportSolver.o XSbuilder.o RayTrace.o FluxMap.o
	$(GCC) $(FFLAGS) -c NE723-FinalProgram.f95

TransportSolver.o: TransportSolver.f95 XSbuilder.o
	$(GCC) $(FFLAGS) -c TransportSolver.f95

XSbuilder.o: XSbuilder.f95
	$(GCC) $(FFLAGS) -c XSbuilder.f95

XSbuilderB.o: XSbuilderB.f95
	$(GCC) $(FFLAGS) -c XSbuilderB.f95

RayTrace.o: RayTrace.f95 XSbuilder.o
	$(GCC) $(FFLAGS) -c RayTrace.f95

FluxMap.o: FluxMap.f95 RayTrace.o XSbuilder.o
	$(GCC) $(FFLAGS) -c FluxMap.f95

RayTraceTest.o: RayTraceTest.f95  RayTrace.o XSbuilder.o
	$(GCC) $(FFLAGS) -c RayTraceTest.f95
