all:
	gfortran multiparametric-flow-cytometry.f90 -o multiparametric-flow-cytometry -O3
clean:
	rm result_file*
	rm multiparametric-flow-cytometry
