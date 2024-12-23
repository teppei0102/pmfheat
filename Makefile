a.out: parameters.f90 get_history.f90
# gfortran -o a.out parameters.f90 get_history.f90
	mpif90 -O0 -g -o a.out parameters.f90 get_history.f90
# mpif90 -check all -warn all -O0 -g -traceback -fpe0 -ftrapuv -check bounds -o a.out parameters.f90 get_history.f90
clean:
	rm parameters.mod *~ *.out
