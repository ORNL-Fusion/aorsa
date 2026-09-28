include compileropts.gnu
# module load  cray-hdf5
# module load cray-netcdf
PGPLOT_PATH=/global/homes/j/jwright/perlmutter-builds/pgplot
PGPLOT_PATH=/global/cfs/cdirs/m77/CompX/pgplot/pgplot.perl.gnu/
FC = ftn
LIBS = $(PGPLOT_PATH)/libpgplot.a

NETCDF_INCLUDE_DIR = ${NETCDF_DIR}/include


LIBS += -L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf
INCLUDE_DIRS +=  -I ${NETCDF_INCLUDE_DIR}

F90FLAGS += -cpp -O2 -std=f2003 -fno-align-commons
