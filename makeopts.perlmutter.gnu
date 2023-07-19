include compileropts.gnu

FC = ftn
LIBS = ../pgplot/libpgplot.a
# /global/cfs/projectdirs/m77/pgplot-perlmutter
NETCDF_INCLUDE_DIR =/opt/cray/pe/netcdf/4.9.0.1/gnu/9.1/include
#${CRAY_PARALLEL_NETCDF_DIR}/include

NETCDF_DIR =/opt/cray/pe/netcdf/4.9.0.1/gnu/9.1/lib
#${CRAY_PARALLEL_NETCDF_DIR}/lib

LIBS += $(NETCDF_DIR)/libnetcdff.a -L $(NETCDF_DIR) -lnetcdf
INCLUDE_DIRS +=  -I ${NETCDF_INCLUDE_DIR}
#presently only works without -O2 below
F90FLAGS += -g -std=f2003 -fno-align-commons
