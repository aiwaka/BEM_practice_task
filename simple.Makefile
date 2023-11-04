FC = gfortran
MATH_LIB = -llapack -lblas
OPTIMIZE=-O3

SOURCE_DIR = .

OBJS=\
constants.o\
subprogram.o\
bem.o

main : $(OBJS) main.o
	$(FC) $(OPTIMIZE) -o main.out $(OBJS) main.o $(MATH_LIB)

main.o : $(SOURCE_DIR)/main.f90
	$(FC) -c $(OPTIMIZE) $(SOURCE_DIR)/main.f90

constants.o : $(SOURCE_DIR)/constants.f90
	$(FC) -c $(OPTIMIZE) $(SOURCE_DIR)/constants.f90

bem.o : $(SOURCE_DIR)/bem.f90
	$(FC) -c $(OPTIMIZE) $(SOURCE_DIR)/bem.f90

subprogram.o : $(SOURCE_DIR)/subprogram.f90
	$(FC) -c $(OPTIMIZE) $(SOURCE_DIR)/subprogram.f90

clean :
	rm  *.o *.mod *.smod main.out
