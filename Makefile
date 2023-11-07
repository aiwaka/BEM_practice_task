FC = $${FORTRAN_FC}
MATH_LIB = $${FORTRAN_MATH_LIB}
OPENMP = $${FORTRAN_OPENMP}
OPTIMIZE=-O3

SOURCE_DIR = .
MOD_DIR = $(SOURCE_DIR)/.mod

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
	rm  *.o *.mod *.smod main.out $(SUB_UTILS_DIR)/*.o $(SUB_UTILS_DIR)/*.out

test: test/test.f90 $(OBJS)
	$(FC) $(OPTIMIZE) -o ./test/test.out $(OBJS) ./test/test.f90 $(MATH_LIB)
	time ./test/test.out

# openmpのスレッド数環境変数がセットされているか確認する.
openmp_validation:
ifdef OPENMP
ifndef OMP_NUM_THREADS
	$(error When using OpenMP, set environmental variable `OMP_NUM_THREADS`.\
	Run `source set_envs.sh`)
endif
endif
	@echo "OMP_NUM_THREADS: $${OMP_NUM_THREADS}"
