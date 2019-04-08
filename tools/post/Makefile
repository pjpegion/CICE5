CDF=/apps/netcdf/4.3.0-intel
#####################################################################
# compiler options
# #####################################################################
#FOPT = -O
#FOPT = -C -warn all
FOPT = -C 
#FOPT = -convert big_endian
#FOPT = -p

 F90 = ifort
# #F90 = ifort -warn
######################################################################
# 
#####################################################################
OBJS = fixgriddefs.o gen_fixgrid.o

gengrid: $(OBJS)
	$(F90) $(FOPT) -o gengrid $(OBJS) -L$(CDF)/lib -lnetcdff -lnetcdf 

%.o: %.F90
	$(F90) $(FOPT) $(optall) -c -I$(CDF)/include $<
	cpp $(optall) -I$(CDF)/include $*.F90>$*.i

clean:
	/bin/rm -f makeit *.o *.i *.mod

