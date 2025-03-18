include compileropts.gnu

FC = ftn
PGPLOT_PATH = /global/cfs/cdirs/m77/pgplot
LIBS = -L $(PGPLOT_PATH) -lpgplot
NETCDF_INCLUDE_DIR = ${NETCDF_DIR}/include
LIBS += -L $(NETCDF_DIR) -lnetcdff -lnetcdf

INCLUDE_DIRS +=  -I ${NETCDF_INCLUDE_DIR}
#presently only works without -O2 below
F90FLAGS += -Ofast -std=f2003
