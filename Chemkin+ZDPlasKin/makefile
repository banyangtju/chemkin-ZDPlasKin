##########################
### makefile in skle
##########################

FF = gfortran

FFLAGS = -c -o $@ -O3

OBJ = dvodeoriginal.o cklib.o dvode_f90_m.o zdplaskin_m.o main.o
dylib = bolsig_x86_64.so

all : $(OBJ)
	$(FF) -O3 $(OBJ) -o run -L . $(dylib)

####################################################
# Chemkin files
####################################################

dvodeoriginal.o: dvodeoriginal.f
	$(FF) $(FFLAGS) $<

cklib.o: cklib.f ckstrt.h
	$(FF) $(FFLAGS) $<

dvode_f90_m.o: dvode_f90_m.F90
	$(FF) $(FFLAGS) $<
dvode_f90_m.mod : dvode_f90_m.o
	@echo Compiling .\dvode_f90_m.F90
	@$(FF) $(FFLAGS) dvode_f90_m.o 

zdplaskin_m.o: zdplaskin_m.F90 dvode_f90_m.mod
	$(FF) $(FFLAGS) $<
zdplaskin_m.mod : zdplaskin_m.o
	@echo Compiling .\zdplaskin_m.F90
	@$(FF) $(FFLAGS) zdplaskin_m.o 

main.o: main.f90
	$(FF) $(FFLAGS) $<


cleanall:
	rm *.o *.mod run

clean:
	rm $(OBJ) *.mod run
