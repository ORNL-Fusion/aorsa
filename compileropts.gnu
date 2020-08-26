# gfortran options
COMMON_OPTION  = -g -fno-automatic -fdefault-real-8 -fdefault-double-8 
COMMON_OPTION2 = -fdefault-real-8 -fdefault-double-8 
COMMON_OPTION3 = 
COMMON_OPTION4 = -fdefault-real-8 -fdefault-double-8 
WARNING_FLAGS = -O2 -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -fbacktrace 
MOD_DIR_FLAG = -J $(MOD_DIR)
F90FLAGS += 
WARNING_FLAGS += -ffpe-trap=invalid -fcheck=all

ORBIT_F_WARNING_FLAGS += -ffpe-trap=invalid 

