# gfortran options
COMMON_OPTION  = -fno-automatic -fdefault-real-8 -fdefault-double-8 
COMMON_OPTION2 = -fdefault-real-8 -fdefault-double-8 
COMMON_OPTION3 = 
COMMON_OPTION4 = -fdefault-real-8 -fdefault-double-8  
WARNING_FLAGS = -Wall -Wuninitialized -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wno-unused-parameter  -fwhole-file  -fcheck=all  -fbacktrace -Wno-implicit-interface -Wno-conversion  -fno-align-commons -fno-automatic -Wno-unused -Wno-unused-dummy-argument #-fallow-argument-mismatch
#-fno-align-commons -Wno-unused -Wno-unused-dummy-argument
# -std=legacy -fdefault-real-8 -fdefault-double-8 
MOD_DIR_FLAG = -J $(MOD_DIR)
F90FLAGS += 
WARNING_FLAGS += -fcheck=all #-ffpe-trap=invalid -fcheck=all

ORBIT_F_WARNING_FLAGS += -ffpe-trap=invalid #-fcheck=all

