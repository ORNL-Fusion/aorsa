include compileropts.gnu

# module load  cray-hdf5
# module load cray-netcdf
WARNING_FLAGS+=-fallow-argument-mismatch
FC = ftn
#change to your pgplot path if needed
PGPLOT_PATH = /global/cfs/cdirs/m77/pgplot
LIBS = -L $(PGPLOT_PATH) -lpgplot

NETCDF_INCLUDE_DIR = ${NETCDF_DIR}/include

LIBS += -L $(NETCDF_DIR) -lnetcdff -lnetcdf
INCLUDE_DIRS +=  -I ${NETCDF_INCLUDE_DIR}

F90FLAGS += -cpp -O2 -std=f2008 -fno-align-commons

