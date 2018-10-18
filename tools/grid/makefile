#NETCDF_DIR = /usrx/local/NetCDF/3.6.3
#NETCDF_DIR = /usrx/local/NetCDF/4.2/serial
NETCDF_DIR=/apps/netcdf/4.3.0-intel
NETCDF_INC = -I${NETCDF_DIR}/include
NETCDF_LIB = -L${NETCDF_DIR}/lib -lnetcdff

#FC = ifort
FC = ifort -C
FC_FLAGS = -c -g $(NETCDF_INC)

SRC = icedefs.F90 generate_cice_fix_file.F90 
OBJ = icedefs.o generate_cice_fix_file.o 
EXE = generate_cice_fix_file.x

.SUFFIXES :
.SUFFIXES : .F90 .o

all: $(EXE)

$(EXE): $(OBJ)
	$(FC) -o $(EXE) $(OBJ) $(NETCDF_LIB) $(NETCDF_INC)

.F90.o:
	$(FC) $(FC_FLAGS) $<

clean: 
	rm -f $(OBJ) $(EXE) *.mod
