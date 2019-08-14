all:
	gfortran MCTA.f90 -o MCTA -O3
clean:
	rm result_file*
	rm multiparametric-flow-cytometry
