EXEC = xaorsa2d

SRC_DIR = src
OBJ_DIR = obj

FFT_DIR = $(SRC_DIR)/FFTPACK
CQL3D_SETUP_DIR = $(SRC_DIR)/CQL3D_SETUP

MOD_DIR = mod
OBJ_DIR = obj

INCLUDE_DIRS = 
LIBS = 
F90FLAGS = 
WARNING_FLAGS = 

# per file build warning flags
ORBIT_F_WARNING_FLAGS = 

OBJ_FILES = \
 $(OBJ_DIR)/cauchy_mod.o \
 $(OBJ_DIR)/size_mod.o \
 $(OBJ_DIR)/aorsa2din_mod.o \
 $(OBJ_DIR)/swim_global_data_mod.o \
 $(OBJ_DIR)/precision_mod.o \
 $(OBJ_DIR)/profile_mod.o \
 $(OBJ_DIR)/z_erfc_mod.o \
 $(OBJ_DIR)/types_mod.o \
 $(OBJ_DIR)/zfun_hilbert_mod.o \
 $(OBJ_DIR)/qlsum.o \
 $(OBJ_DIR)/ql_myra.o \
 $(OBJ_DIR)/wdot_sum.o \
 $(OBJ_DIR)/wdot_test.o \
 $(OBJ_DIR)/mets2aorsa.o \
 $(OBJ_DIR)/cauchy_ppart.o \
 $(OBJ_DIR)/fftpack5.1d.o \
 $(OBJ_DIR)/vlog.o \
 $(OBJ_DIR)/aorsaSubs.o \
 $(OBJ_DIR)/sigma.o \
 $(OBJ_DIR)/zfunction.o \
 $(OBJ_DIR)/ztable.o \
 $(OBJ_DIR)/current.o \
 $(OBJ_DIR)/mets2aorsa_myra.o \
 $(OBJ_DIR)/slowDown.o \
 $(OBJ_DIR)/fourier.o \
 $(OBJ_DIR)/assert.o \
 $(OBJ_DIR)/setupblacs.o \
 $(OBJ_DIR)/bessel.o \
 $(OBJ_DIR)/check.o \
 $(OBJ_DIR)/rf2x_setup2.o \
 $(OBJ_DIR)/profile_setup.o \
 $(OBJ_DIR)/eqdsk_setup.o \
 $(OBJ_DIR)/orbit.o \
 $(OBJ_DIR)/eqdsk_plot.o \
 $(OBJ_DIR)/fieldws.o \
 $(OBJ_DIR)/scale.o \
 $(OBJ_DIR)/dql_write.o \
 $(OBJ_DIR)/dshell.o \
 $(OBJ_DIR)/aorsa2dMain.o \
 $(OBJ_DIR)/plot.o \
 $(OBJ_DIR)/aorsa2dSum.o 

OBJ_FFT = \
 $(OBJ_DIR)/cfftb1.o \
 $(OBJ_DIR)/cfftf1.o \
 $(OBJ_DIR)/cffti1.o \
 $(OBJ_DIR)/dfftb.o \
 $(OBJ_DIR)/dfftf.o \
 $(OBJ_DIR)/dffti.o \
 $(OBJ_DIR)/passb.o \
 $(OBJ_DIR)/passb2.o \
 $(OBJ_DIR)/passb3.o \
 $(OBJ_DIR)/passb4.o \
 $(OBJ_DIR)/passb5.o \
 $(OBJ_DIR)/passf.o \
 $(OBJ_DIR)/passf2.o \
 $(OBJ_DIR)/passf3.o \
 $(OBJ_DIR)/passf4.o \
 $(OBJ_DIR)/passf5.o \
 $(OBJ_DIR)/zfftb.o \
 $(OBJ_DIR)/zfftf.o \
 $(OBJ_DIR)/zffti.o
 
OBJ_CQL3D_SETUP = \
 $(OBJ_DIR)/global_data_m.o \
 $(OBJ_DIR)/basis_functions_m.o \
 $(OBJ_DIR)/f_expanded_m.o \
 $(OBJ_DIR)/CQL_kinds_m.o \
 $(OBJ_DIR)/vector_write_m.o \
 $(OBJ_DIR)/read_cql3d.o \
 $(OBJ_DIR)/ceez.o \
 $(OBJ_DIR)/cubic_B_splines_v.o \
 $(OBJ_DIR)/cql3d_setup.o


# Determine machine
# -----------------

UNAME_S := $(shell uname -s)
UNAME_R := $(shell uname -r)
LSB_IS := $(shell lsb_release -is)
LSB_RS := $(shell lsb_release -rs)
HOSTNAME := $(shell hostname)

SYSTEM_IDENTIFIED = 0
ifeq ($(LMOD_SYSHOST),perlmutter)
  ifeq ($(PE_ENV),GNU)
    include makeopts.perlmutter.gnu
    $(info System identified as Perlmutter GNU)
    SYSTEM_IDENTIFIED = 1
  endif
  ifeq ($(PE_ENV),AOCC)
    include makeopts.perlmutter.aocc
    $(info System identified as Perlmutter AOCC)
    SYSTEM_IDENTIFIED = 1
  endif
  ifeq ($(PE_ENV),INTEL)
    include makeopts.perlmutter.intel
    $(info System identified as Perlmutter INTEL)
    SYSTEM_IDENTIFIED = 1
  endif
endif
ifeq ($(NERSC_HOST),cori)
  ifeq ($(PE_ENV),GNU)
    include makeopts.cori.gnu
    $(info System identified as Cori GNU)
    SYSTEM_IDENTIFIED = 1
  endif      
  ifeq ($(PE_ENV),INTEL)
    include makeopts.cori.intel
    $(info System identified as Cori Intel)
    SYSTEM_IDENTIFIED = 1
  endif
endif
ifeq ($(UNAME_S),Darwin) # OSX
  #ifeq ($(UNAME_R),18.7.0)
    include makeopts.osx-mojave
    $(info System identified as osx-mojave)
    SYSTEM_IDENTIFIED = 1
  #endif
endif
ifeq ($(LSB_IS),Ubuntu)
  ifeq ($(LSB_RS),20.04)
    include makeopts.ubuntu20.04
    $(info System identified as Ubuntu20.04)
    SYSTEM_IDENTIFIED = 1
  endif
endif

#$(error EXIT)

ifeq ($(SYSTEM_IDENTIFIED),0)
  $(error No build configuration for this system)
endif


F90          = $(FC) -c $(COMMON_OPTION) #$(INCLUDE_DIRS)
F90_NOSAVE   = $(FC) -c $(COMMON_OPTION2) # $(INCLUDE_DIRS)
F90_r4       = $(FC) -c $(COMMON_OPTION3) #$(INCLUDE_DIRS)
F90_4        = $(FC) -c $(COMMON_OPTION4) #$(INCLUDE_DIRS)
F90_LOAD     = $(FC)    $(COMMON_OPTION) #$(INCLUDE_DIRS)

INLINE=
OPTIMIZATION =  
F90FLAGS += $(INLINE) $(OPTIMIZATION) $(MOD_DIR_FLAG)

COMPILE90          = $(F90)          $(F90FLAGS) $(DEFS)
COMPILE90_NOSAVE   = $(F90_NOSAVE)   $(F90FLAGS) $(DEFS)
COMPILE_r4         = $(F90_r4)       $(F90FLAGS) $(DEFS)
COMPILE90_4        = $(F90_4)        $(F90FLAGS) $(DEFS)

LOADFLAGS = $(DFFTPACK)
LOAD = $(F90_LOAD) $(OPTIMIZATION) 

# remove the "SECONDARY" line and life will be very weird
.SECONDARY:

$(EXEC):  make_directories $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP)
	$(info LIBS are $(LIBS))
	$(LOAD) -o $(EXEC) $(OBJ_FILES) $(OBJ_FFT) $(OBJ_CQL3D_SETUP) $(LOADFLAGS) $(LIBS) 

# Dependencies

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f90
	${COMPILE90} -c $< -o $@ ${INCLUDE_DIRS} ${WARNING_FLAGS}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.f
	${COMPILE90} -c $< -o $@ ${INCLUDE_DIRS} ${WARNING_FLAGS}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.F90
	${COMPILE90} -c $< -o $@ ${INCLUDE_DIRS} ${WARNING_FLAGS}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.F
	${COMPILE90} -c $< -o $@ ${INCLUDE_DIRS} ${WARNING_FLAGS}

${OBJ_DIR}/orbit.o: ${SRC_DIR}/orbit.f
	${COMPILE90} -c $< -o $@ ${INCLUDE_DIRS} ${SIGMA_F_WARNING_FLAGS}


$(OBJ_DIR)/rf2x_setup2.o:    $(SRC_DIR)/rf2x_setup2.f 
	                     $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/rf2x_setup2.o \
                             $(SRC_DIR)/rf2x_setup2.f $(INCLUDE_DIRS)

$(OBJ_DIR)/profile_setup.o:  $(SRC_DIR)/profile_setup.f
	                     $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/profile_setup.o \
                             $(SRC_DIR)/profile_setup.f $(INCLUDE_DIRS)

$(OBJ_DIR)/eqdsk_setup.o:    $(SRC_DIR)/eqdsk_setup.f
	                     $(COMPILE90) -o $(OBJ_DIR)/eqdsk_setup.o \
                             $(SRC_DIR)/eqdsk_setup.f $(INCLUDE_DIRS)

$(OBJ_DIR)/orbit.o:          $(SRC_DIR)/orbit.f
	                     $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/orbit.o \
                             $(SRC_DIR)/orbit.f $(INCLUDE_DIRS)			     

$(OBJ_DIR)/eqdsk_plot.o:     $(SRC_DIR)/eqdsk_plot.f90
	                     $(COMPILE_r4) -o $(OBJ_DIR)/eqdsk_plot.o \
                             $(SRC_DIR)/eqdsk_plot.f90 $(INCLUDE_DIRS) $(WARNING_FLAGS)				     

$(OBJ_DIR)/fieldws.o:        $(SRC_DIR)/fieldws.f
	                     $(COMPILE_r4) -o $(OBJ_DIR)/fieldws.o \
                             $(SRC_DIR)/fieldws.f $(INCLUDE_DIRS)

$(OBJ_DIR)/scale.o:          $(SRC_DIR)/scale.f
	                     $(COMPILE_r4) -o $(OBJ_DIR)/scale.o \
                             $(SRC_DIR)/scale.f $(INCLUDE_DIRS)			     

$(OBJ_DIR)/dql_write.o:      $(SRC_DIR)/dql_write.f 
	                     $(COMPILE_r4) -o $(OBJ_DIR)/dql_write.o \
                             $(SRC_DIR)/dql_write.f $(INCLUDE_DIRS) ${WARNING_FLAGS}

$(OBJ_DIR)/plot.o:           $(SRC_DIR)/plot.f
	                     $(COMPILE_r4) -o $(OBJ_DIR)/plot.o \
                             $(SRC_DIR)/plot.f $(INCLUDE_DIRS)				     

### FFTPACK files:

${OBJ_DIR}/%.o: ${FFT_DIR}/%.f
	${COMPILE90} -c $< -o $@ ${INCLUDE_DIRS}


### CQL3D_SETUP files:

$(OBJ_DIR)/basis_functions_m.o: $(CQL3D_SETUP_DIR)/basis_functions_m.f90
	                        $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/basis_functions_m.o \
                                $(CQL3D_SETUP_DIR)/basis_functions_m.f90 $(INCLUDE_DIRS)

$(OBJ_DIR)/f_expanded_m.o:      $(CQL3D_SETUP_DIR)/f_expanded_m.f90
	                        $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/f_expanded_m.o \
                                $(CQL3D_SETUP_DIR)/f_expanded_m.f90 $(INCLUDE_DIRS)

$(OBJ_DIR)/global_data_m.o:     $(CQL3D_SETUP_DIR)/global_data_m.f90
	                        $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/global_data_m.o \
                                $(CQL3D_SETUP_DIR)/global_data_m.f90 $(INCLUDE_DIRS)

$(OBJ_DIR)/CQL_kinds_m.o:       $(CQL3D_SETUP_DIR)/CQL_kinds_m.f90
	                        $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/CQL_kinds_m.o \
                                $(CQL3D_SETUP_DIR)/CQL_kinds_m.f90 $(INCLUDE_DIRS)

$(OBJ_DIR)/vector_write_m.o:    $(CQL3D_SETUP_DIR)/vector_write_m.f90
	                        $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/vector_write_m.o \
                                $(CQL3D_SETUP_DIR)/vector_write_m.f90 $(INCLUDE_DIRS)

$(OBJ_DIR)/read_cql3d.o:        $(CQL3D_SETUP_DIR)/read_cql3d.f90
	                        $(COMPILE90_4) -o $(OBJ_DIR)/read_cql3d.o \
                                $(CQL3D_SETUP_DIR)/read_cql3d.f90 $(INCLUDE_DIRS) ${WARNING_FLAGS}

$(OBJ_DIR)/ceez.o:              $(CQL3D_SETUP_DIR)/ceez.f
	                        $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/ceez.o \
                                $(CQL3D_SETUP_DIR)/ceez.f $(INCLUDE_DIRS)

$(OBJ_DIR)/cubic_B_splines_v.o: $(CQL3D_SETUP_DIR)/cubic_B_splines_v.f90
	                        $(COMPILE90_NOSAVE) -o $(OBJ_DIR)/cubic_B_splines_v.o \
                                $(CQL3D_SETUP_DIR)/cubic_B_splines_v.f90 $(INCLUDE_DIRS)		

$(OBJ_DIR)/cql3d_setup.o:       $(CQL3D_SETUP_DIR)/cql3d_setup.f90
	                        $(COMPILE90_4) -o $(OBJ_DIR)/cql3d_setup.o \
                                $(CQL3D_SETUP_DIR)/cql3d_setup.f90 $(INCLUDE_DIRS)

.phony: make_directories
make_directories: $(OBJ_DIR)/ $(MOD_DIR)/

$(MOD_DIR)/:
	mkdir -p $@

$(OBJ_DIR)/:
	mkdir -p $@

.phony: clean
clean:
	rm -r $(MOD_DIR)/*.mod $(OBJ_DIR)/*.o $(EXEC)

.phony: cleanrepo
cleanrepo:
	cd test/DIIID-helicon
	./cleanrun.sh
	cd ../DIIID-whistler
	./cleanrun.sh
	cd ../../
	cd examples/DIIID-helicon
	./cleanrun.sh
	cd ../DIIID-whister
	./cleanrun.sh
	cd ../../
