#====================================================
#
#Makefile for symChern by  Zhang Zeying
#
#====================================================
obj =	 module.o readsymhr.o  hamilton.o util.o chern.o main.o
	

f90 =	ifort -fpp
f90 =	mpif90 -fpp -DMPI
f90 =	mpif90 -cpp -DMPI
f90 =	gfortran 

#f90 = gfortran -cpp
#flag=	-O3 -std03 -gen-interfaces
#flag=	-O3 -std03 
#flag=	-O3 -assume realloc_lhs
flag=	-O3 
#flag=	-O3 -std=f2008
libs = -L/usr/lib/ -llapack -lblas

check=	-g -fcheck=all -fbacktrace -Wall -Wextra -Wno-maybe-uninitialized -pedantic
check=	


clean=	rm -f *.o *.mod *__genmod.f90
#clean=	rm -f *.o *.mod 


main :	$(obj)
	$(f90) $(flag) -o chern $(obj) $(libs)
#	$(clean)




$(obj):	%.o:	%.f90
	$(f90) $(check)  $(flag) -c  $< -o $@


clean:	
	$(clean)

#module.o:	../constants.F90
