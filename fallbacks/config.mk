#
# Default configuration parameters obtained from Abinit's Autoconf
#

# ---------------------------------------------------------------------------- #

#
# System information
#

# Architecture
abi_cpu_model   = xeon
abi_cpu_64bits  = yes
abi_cpu_spec    = intel_xeon
abi_sys_spec    = linux-x86_64

# ---------------------------------------------------------------------------- #

#
# Libraries and linking
#

AR             = ar
ARFLAGS        =      rc
ARFLAGS_64BITS = 
ARFLAGS_DEBUG  = 
ARFLAGS_OPTIM  = 
ARFLAGS_EXTRA  = 
ARFLAGS_CMD    = rc
RANLIB         = ranlib

# ---------------------------------------------------------------------------- #

#
# Preprocessing
#

CPP               = mpiicc -E
CPPFLAGS          = 
CPPFLAGS_DEBUG    = 
CPPFLAGS_EXTRA    = 
CPPFLAGS_HINTS    = 
CPPFLAGS_OPTIM    = 

FPP               = 
FPPFLAGS          =     
FPPFLAGS_DEBUG    = 
FPPFLAGS_EXTRA    = 
FPPFLAGS_HINTS    = 
FPPFLAGS_OPTIM    = 

INCLUDES          = 

TRUE_CPP          = cpp -P -std=c99

# ---------------------------------------------------------------------------- #

#
# C compilation
#

cc_info_string    = Intel(R) C Intel(R) 64 Compiler for applications running on Intel(R) 64, Version 18.0.3.222 Build 20180410

abi_cc_vendor     = intel
abi_cc_version    = 18.0

CC                = mpiicc
CFLAGS            =  -g -O2 -vec-report0 
CFLAGS_64BITS     = 
CFLAGS_DEBUG      = -g
CFLAGS_EXTRA      = 
CFLAGS_HINTS      = -vec-report0
CFLAGS_OPTIM      = -O2

CC_LDFLAGS_64BITS = 
CC_LDFLAGS_DEBUG  = 
CC_LDFLAGS_EXTRA  = 
CC_LDFLAGS_HINTS  = -static-libgcc -static-intel
CC_LDFLAGS_OPTIM  = 
CC_LDFLAGS        =    -static-libgcc -static-intel 

CC_LIBS_DEBUG     = 
CC_LIBS_EXTRA     = 
CC_LIBS_HINTS     = 
CC_LIBS_OPTIM     = 
CC_LIBS           =     

# ---------------------------------------------------------------------------- #

#
# C++ compilation
#

cxx_info_string    = icpc (ICC) 18.0.3 20180410

abi_cxx_vendor     = gnu
abi_cxx_version    = 18.0

CXX                = mpiicpc
CXXFLAGS           =  -g -O2 -mtune=native -march=native  
CXXFLAGS_64BITS    = 
CXXFLAGS_DEBUG     = -g
CXXFLAGS_EXTRA     = 
CXXFLAGS_HINTS     = 
CXXFLAGS_OPTIM     = -O2 -mtune=native -march=native

CXX_LDFLAGS_64BITS = 
CXX_LDFLAGS_DEBUG  = 
CXX_LDFLAGS_EXTRA  = 
CXX_LDFLAGS_HINTS  = 
CXX_LDFLAGS_OPTIM  = 
CXX_LDFLAGS        =     

CXX_LIBS_DEBUG     = 
CXX_LIBS_EXTRA     = 
CXX_LIBS_HINTS     = 
CXX_LIBS_OPTIM     = 
CXX_LIBS           =     

# ---------------------------------------------------------------------------- #

#
# Fortran compilation
#

fc_info_string    = Intel(R) Fortran Intel(R) 64 Compiler for applications running on Intel(R) 64, Version 18.0.3.222 Build 20180410

fc_vendor         = intel
fc_version        = 18.0

FC                = mpiifort
FCFLAGS           =  -g -extend-source -noaltparam -nofpscomp  
FCFLAGS_64BITS    = 
FCFLAGS_DEBUG     = -g
FCFLAGS_EXTRA     = 
FCFLAGS_HINTS     = -extend-source -noaltparam -nofpscomp
FCFLAGS_OPTIM     = -O3

FCFLAGS_FIXEDFORM = -fixed
FCFLAGS_FREEFORM  = -free

FC_LDFLAGS_64BITS = 
FC_LDFLAGS_DEBUG  = 
FC_LDFLAGS_EXTRA  = 
FC_LDFLAGS_HINTS  = -static-intel -static-libgcc
FC_LDFLAGS_OPTIM  = 
FC_LDFLAGS        =    -static-intel -static-libgcc 

FC_LIBS_DEBUG     = 
FC_LIBS_EXTRA     = 
FC_LIBS_HINTS     = 
FC_LIBS_OPTIM     = 
FC_LIBS           =       -L/usr/local/boost/1.55.0/lib64 -L/usr/local/fftw/intel-16.0/3.3.4/lib64 -L/opt/intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/lib/release_mt -L/opt/intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/lib -L/usr/local/intel/lib64/intel-mpi -L/opt/intel/compilers_and_libraries_2018.3.222/linux/tbb/lib/intel64/gcc4.7 -L/opt/intel/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64 -L/opt/intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64 -L/usr/local/intel/lib64 -L/opt/intel/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64_lin -L/usr/local/fftw/intel-16.0/3.3.4/lib64/../lib64 -L/usr/local/fftw/intel-16.0/3.3.4/lib64/../lib64/ -L/usr/local/boost/1.55.0/lib64/../lib64 -L/usr/local/boost/1.55.0/lib64/../lib64/ -L/usr/local/intel/lib64/../lib64 -L/usr/local/intel/lib64/../lib64/ -L/tigress/xingz/gcc/lib/gcc/x86_64-pc-linux-gnu/6.3.0/ -L/tigress/xingz/gcc/lib/gcc/x86_64-pc-linux-gnu/6.3.0/../../../../lib64 -L/tigress/xingz/gcc/lib/gcc/x86_64-pc-linux-gnu/6.3.0/../../../../lib64/ -L/lib/../lib64 -L/lib/../lib64/ -L/usr/lib/../lib64 -L/usr/lib/../lib64/ -L/usr/local/fftw/intel-16.0/3.3.4/lib64/ -L/usr/local/boost/1.55.0/lib64/ -L/opt/intel/compilers_and_libraries_2018.3.222/linux/mpi/intel64/lib/ -L/usr/local/intel/lib64/intel-mpi/ -L/opt/intel/compilers_and_libraries_2018.3.222/linux/tbb/lib/intel64/gcc4.7/ -L/opt/intel/compilers_and_libraries_2018.3.222/linux/compiler/lib/intel64/ -L/opt/intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64/ -L/usr/local/intel/lib64/ -L/tigress/xingz/gcc/lib/gcc/x86_64-pc-linux-gnu/6.3.0/../../../ -L/lib64 -L/lib/ -L/usr/lib64 -L/usr/lib -lmpifort -lmpi -lmpigi -ldl -lrt -lpthread -lifport -lifcoremt -limf -lsvml -lm -lipgo -lirc -lirc_s

MODEXT            = mod

# ---------------------------------------------------------------------------- #

#
# Exports
#

# Algorithmic use flags
lib_algo_incs      = 
lib_algo_libs      = 

                    ########################################

# FFT use flags
lib_fft_incs       = 
lib_fft_libs       = -lfftw3

                    ########################################

# Linear algebra use flags
lib_linalg_fcflags = 
lib_linalg_ldflags = 
lib_linalg_incs    = 
lib_linalg_libs    = -lmkl_intel_lp64  -lmkl_intel_lp64 -L/opt/intel/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64 -lmkl_sequential -lmkl_core -lpthread -lm 

                    ########################################

# Mathematical use flags
lib_math_incs      = 
lib_math_libs      = 

                    ########################################

# Timer use flags
lib_timer_incs     = 
lib_timer_libs     =  -lrt

                    ########################################

# ETSF_IO use flags
lib_etsf_io_incs   = 
lib_etsf_io_libs   = 

# NetCDF use flags
lib_netcdf_incs    = 
lib_netcdf_libs    = 

# YAML use flags
lib_yaml_incs    = 
lib_yaml_libs    = 

                    ########################################

# LibXC use flags
lib_libxc_incs     = -I$(fallbacks_instdir)/include
lib_libxc_libs     = -L$(fallbacks_instdir)/lib -lxcf90 -lxc

                    ########################################
# cthyb 
lib_cthyb_incs = @lib_cthyb_incs@
lib_cthyb_libs = @lib_cthyb_libs@

                    ########################################
# triqs
lib_triqs_incs = 
lib_triqs_libs = 

                    ########################################

# Fortran flags
FCFLAGS_ATOMPAW_EXT   =  -g -extend-source -noaltparam -nofpscomp   -O3
FCFLAGS_BIGDFT_EXT    =  -g -extend-source -noaltparam -nofpscomp   -O3
FCFLAGS_ETSF_IO_EXT   =  -g -extend-source -noaltparam -nofpscomp   @fcflags_opt_etsf_io@
FCFLAGS_LIBXC_EXT     =  -g -extend-source -noaltparam -nofpscomp   -O3
FCFLAGS_LINALG_EXT    =  -g -extend-source -noaltparam -nofpscomp   -O3
FCFLAGS_NETCDF_EXT    =  -g -extend-source -noaltparam -nofpscomp   -O3
FCFLAGS_WANNIER90_EXT =  -g -extend-source -noaltparam -nofpscomp   -O3
