/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/*
 * Copyright (C) 2005-2018 ABINIT Group (Yann Pouillon)
 *
 * This file is part of the Abinit software package. For license information,
 * please see the COPYING file in the top-level directory of the Abinit source
 * distribution.
 *
 */

/* Abinit configuration */

#ifndef _ABINIT_CONFIG_H
#define _ABINIT_CONFIG_H

#ifdef __INTEL_COMPILER
#define FC_INTEL 1
#endif



/* Abinit target description. */
#define ABINIT_TARGET "x86_64_linux_intel18.0"

/* Abinit whole version number. */
#define ABINIT_VERSION "8.10.2"

/* Abinit base version number. */
#define ABINIT_VERSION_BASE "8.10"

/* Abinit build date. */
#define ABINIT_VERSION_BUILD "20181112"

/* Abinit major version number. */
#define ABINIT_VERSION_MAJOR "8"

/* Abinit micro version number (patch level). */
#define ABINIT_VERSION_MICRO "2"

/* Abinit minor version number. */
#define ABINIT_VERSION_MINOR "10"

/* Fortran module name mangling macro. */
/* #undef ABI_FC_MOD */

/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* Define to 1 if you are using the COMPAQ C compiler. */
/* #undef CC_COMPAQ */

/* Define to 1 if you are using the GNU C compiler. */
/* #undef CC_GNU */

/* Define to 1 if you are using the IBM XL C compiler. */
/* #undef CC_IBM */

/* Define to 1 if you are using the Intel C compiler. */
#define CC_INTEL 1

/* Define to 1 if you are using the Open64 C compiler. */
/* #undef CC_OPEN64 */

/* Define to 1 if you are using the PathScale C compiler. */
/* #undef CC_PATHSCALE */

/* Define to 1 if you are using the Portland Group C compiler. */
/* #undef CC_PGI */

/* Define to 1 if you are using the Sun C compiler. */
/* #undef CC_SUN */

/* Define to 1 if you are using the COMPAQ C++ compiler. */
/* #undef CXX_COMPAQ */

/* Define to 1 if you are using the GNU C++ compiler. */
#define CXX_GNU 1

/* Define to 1 if you are using the IBM XL C++ compiler. */
/* #undef CXX_IBM */

/* Define to 1 if you are using the Intel C++ compiler. */
/* #undef CXX_INTEL */

/* Define to 1 if you are using the Open64 C++ compiler. */
/* #undef CXX_OPEN64 */

/* Define to 1 if you are using the PathScale C++ compiler. */
/* #undef CXX_PATHSCALE */

/* Define to 1 if you are using the Portland Group C++ compiler. */
/* #undef CXX_PGI */

/* Define to 1 if you are using the Sun C++ compiler. */
/* #undef CXX_SUN */

/* Define to 1 if you want to activate design-by-contract debugging tests. */
/* #undef DEBUG_CONTRACT */

/* Define to 1 to turn on debug mode. */
/* #undef DEBUG_MODE */

/* Define to 1 to turn on verbose debug messages. */
/* #undef DEBUG_VERBOSE */

/* Define to 1 if you are using the ABSOFT Fortran compiler. */
/* #undef FC_ABSOFT */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to 1 if you are using the Fujitsu Fortran compiler. */
/* #undef FC_FUJITSU */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Define to 1 if you are using the G95 Fortran compiler. */
/* #undef FC_G95 */

/* Define to 1 if you are using the GNU Fortran compiler. */
/* #undef FC_GNU */

/* Define to 1 if you are using the Hitachi Fortran compiler. */
/* #undef FC_HITACHI */

/* Define to 1 if you are using the IBM XL Fortran compiler. */
/* #undef FC_IBM */

/* Define to 1 if you are using the Intel Fortran compiler. */
#define FC_INTEL 1

/* Define to 1 if you are using the MIPSpro Fortran compiler. */
/* #undef FC_MIPSPRO */

/* Define to 1 if you are using the NAGWare Fortran 95 compiler. */
/* #undef FC_NAG */

/* Define to 1 if you are using the Open64 Fortran compiler. */
/* #undef FC_OPEN64 */

/* Define to 1 if you are using the PathScale Fortran compiler. */
/* #undef FC_PATHSCALE */

/* Define to 1 if you are using the Portland Group Fortran compiler. */
/* #undef FC_PGI */

/* Define to 1 if you are using the Sun Fortran compiler. */
/* #undef FC_SUN */

/* Define to 1 if you have the `abort' function. */
#define HAVE_ABORT 1

/* Define to 1 if you have the AtomPAW library. */
#define HAVE_ATOMPAW 1

/* Define to 1 if you want to disable vectorization in problematic procedures.
   */
/* #undef HAVE_AVX_SAFE_MODE */

/* Define to 1 if you have the BigDFT library. */
/* #undef HAVE_BIGDFT */

/* Define to 1 if you want to activate Bethe-Salpeter unpacking
   (EXPERIMENTAL). */
/* #undef HAVE_BSE_UNPACKED */

/* Define to 1 if you are working in a Bazaar branch. */
/* #undef HAVE_BZR_BRANCH */

/* Use C clock for timings. */
/* #undef HAVE_CCLOCK */

/* Define to 1 if you have the `clock_gettime' function. */
#define HAVE_CLOCK_GETTIME 1

/* Define to 1 if you have the <cublas.h> header file. */
/* #undef HAVE_CUBLAS_H */

/* Define to 1 if you have the <cuda_runtime_api.h> header file. */
/* #undef HAVE_CUDA_RUNTIME_API_H */

/* Define to 1 if you have the <cufft.h> header file. */
/* #undef HAVE_CUFFT_H */

/* Define to 1 if you have ELPA Fortran 2008 API support */
/* #undef HAVE_ELPA_FORTRAN2008 */

/* Define to 1 if you have the <errno.h> header file. */
#define HAVE_ERRNO_H 1

/* Define to 1 if you have the ETSF_IO library. */
/* #undef HAVE_ETSF_IO */

/* Define to 1 if your Fortran compiler supports allocatable arrays in
   datatypes. */
#define HAVE_FC_ALLOCATABLE_DTARRAYS 1

/* Define to 1 if your Fortran compiler supports the asynchronous attribute.
   */
#define HAVE_FC_ASYNC 1

/* Define to 1 if your Fortran compiler supports BACKTRACE. */
/* #undef HAVE_FC_BACKTRACE */

/* Define to 1 if your Fortran compiler supports GET_COMMAND_ARGUMENT. */
#define HAVE_FC_COMMAND_ARGUMENT 1

/* Define to 1 if your Fortran compiler supports EXECUTE_COMMAND_LINE. */
#define HAVE_FC_COMMAND_LINE 1

/* Define to 1 if your Fortran compiler supports the contiguous attribute. */
#define HAVE_FC_CONTIGUOUS 1

/* Define to 1 if your Fortran compiler supports cpu_time(). */
#define HAVE_FC_CPUTIME 1

/* Define to 1 if your Fortran compiler supports etime(). */
/* #undef HAVE_FC_ETIME */

/* Define to 1 if your Fortran compiler supports exit(). */
#define HAVE_FC_EXIT 1

/* Define to 1 if your Fortran compiler supports flush(). */
#define HAVE_FC_FLUSH 1

/* Define to 1 if your Fortran compiler supports flush_(). */
/* #undef HAVE_FC_FLUSH_ */

/* Define to 1 if your Fortran compiler supports gamma(). */
#define HAVE_FC_GAMMA 1

/* Define to 1 if your Fortran compiler supports getenv(). */
/* #undef HAVE_FC_GETENV */

/* Define to 1 if your Fortran compiler supports getpid(). */
/* #undef HAVE_FC_GETPID */

/* Define to 1 if your Fortran compiler supports IEEE_EXCEPTIONS. */
#define HAVE_FC_IEEE_EXCEPTIONS 1

/* Define to 1 if your Fortran compiler accepts quadruple integers. */
/* #undef HAVE_FC_INT_QUAD */

/* Define to 1 if your Fortran compiler supports IOMSG. */
#define HAVE_FC_IOMSG 1

/* Define to 1 if your Fortran compiler provides the iso_c_binding module. */
#define HAVE_FC_ISO_C_BINDING 1

/* Define to 1 if your Fortran compiler supports 2008 standard in
   ISO_FORTRAN_ENV. */
#define HAVE_FC_ISO_FORTRAN_2008 1

/* Define to 1 if your Fortran compiler supports long lines. */
#define HAVE_FC_LONG_LINES 1

/* Define to 1 if your Fortran compiler supports \newline in a macros. */
/* #undef HAVE_FC_MACRO_NEWLINE */

/* Define to 1 if your Fortran compiler supports MOVE_ALLOC (F2003). */
#define HAVE_FC_MOVE_ALLOC 1

/* Define to 1 if your Fortran compiler supports the private attribute. */
#define HAVE_FC_PRIVATE 1

/* Define to 1 if your Fortran compiler supports the protected attribute. */
#define HAVE_FC_PROTECTED 1

/* Define to 1 if your Fortran compiler supports stream IO. */
#define HAVE_FC_STREAM_IO 1

/* Define to 1 if your Fortran compiler supports SYSTEM. */
/* #undef HAVE_FC_SYSTEM */

/* Define to 1 if you have an optimized FFT library. */
#define HAVE_FFT 1

/* Define to 1 if you want to use the ASL library for FFT. */
/* #undef HAVE_FFT_ASL */

/* Define to 1 if you want to use the DFTI library. */
/* #undef HAVE_FFT_DFTI */

/* Define to 1 if you want to use the threaded DFTI library. */
/* #undef HAVE_FFT_DFTI_THREADS */

/* Define to 1 if you want to use the FFTW2 library. */
/* #undef HAVE_FFT_FFTW2 */

/* Define to 1 if you want to use the threaded FFTW2 library. */
/* #undef HAVE_FFT_FFTW2_THREADS */

/* Define to 1 if you want to use the FFTW3 library. */
#define HAVE_FFT_FFTW3 1

/* Define to 1 if you want to use the threaded FFTW3 library. */
/* #undef HAVE_FFT_FFTW3_MKL */

/* Define to 1 if you want to use the distributed FFTW3 library. */
/* #undef HAVE_FFT_FFTW3_MPI */

/* Define to 1 if you want to use the threaded FFTW3 library. */
/* #undef HAVE_FFT_FFTW3_THREADS */

/* Define to 1 if you want to use the HP MLIB library for FFT. */
/* #undef HAVE_FFT_MLIB */

/* Define to 1 if you have an optimized MPI-parallel FFT library. */
#define HAVE_FFT_MPI 1

/* Define to 1 if you have an optimized serial FFT library. */
#define HAVE_FFT_SERIAL 1

/* Define to 1 if you want to use the SGIMATH library for FFT. */
/* #undef HAVE_FFT_SGIMATH */

/* Define to 1 if your Fortran compiler supports Fortran 2003. */
/* #undef HAVE_FORTRAN2003 */

/* Define to 1 if you have a GPU library. */
/* #undef HAVE_GPU */

/* Define to 1 if you have the Cuda library. */
/* #undef HAVE_GPU_CUDA */

/* Define to 1 if you have a Cuda version < 4. */
/* #undef HAVE_GPU_CUDA3 */

/* Define to 1 if you want to perform double-precision Cuda calculations. */
/* #undef HAVE_GPU_CUDA_DP */

/* Define to 1 if you want to perform single-precision Cuda calculations. */
/* #undef HAVE_GPU_CUDA_SP */

/* Define to 1 if you have a MPI-aware GPU library. */
/* #undef HAVE_GPU_MPI */

/* Define to 1 if you have a serial GPU library. */
/* #undef HAVE_GPU_SERIAL */

/* Define to 1 if you have the GNU Scientific Library. */
/* #undef HAVE_GSL */

/* Define to 1 if you have the <gsl/gsl_sf_gamma.h> header file. */
/* #undef HAVE_GSL_GSL_SF_GAMMA_H */

/* Define to 1 if you want to activate double-precision GW calculations
   (EXPERIMENTAL). */
/* #undef HAVE_GW_DPC */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the Levenberg-Marquardt algorithmic library. */
/* #undef HAVE_LEVMAR */

/* Define to 1 if you have the <levmar.h> header file. */
/* #undef HAVE_LEVMAR_H */

/* Define to 1 if you want to activate internal support for libPAW in ABINIT.
   */
#define HAVE_LIBPAW_ABINIT 1

/* Define to 1 if you want to activate internal support for libtetra(hedron)
   in ABINIT. */
#define HAVE_LIBTETRA_ABINIT 1

/* Define to 1 if you have the `veclib' library (-lveclib). */
/* #undef HAVE_LIBVECLIB */

/* Define to 1 if you have the LibXC library. */
#define HAVE_LIBXC 1

/* Define to 1 if you want to enable support for XML. */
/* #undef HAVE_LIBXML */

/* Define to 1 if you have an optimized linear algebra library. */
#define HAVE_LINALG 1

/* Define to 1 if you have the ASL linear algebra library. */
/* #undef HAVE_LINALG_ASL */

/* Define to 1 if you have an AXPBY BLAS1 extensions. */
/* #undef HAVE_LINALG_AXPBY */

/* Define to 1 if you have an optimized ELPA linear algebra library. */
/* #undef HAVE_LINALG_ELPA */

/* Define to 1 if you have ELPA 2013 API support */
/* #undef HAVE_LINALG_ELPA_2013 */

/* Define to 1 if you have ELPA 2014 API support */
/* #undef HAVE_LINALG_ELPA_2014 */

/* Define to 1 if you have ELPA 2015 API support */
/* #undef HAVE_LINALG_ELPA_2015 */

/* Define to 1 if you have ELPA 2016 API support */
/* #undef HAVE_LINALG_ELPA_2016 */

/* Define to 1 if you have ELPA 2017 API support */
/* #undef HAVE_LINALG_ELPA_2017 */

/* Define to 1 if you have the ESSL linear algebra library. */
/* #undef HAVE_LINALG_ESSL */

/* Define to 1 if you have ?GEMM3M BLAS3 extensions. */
/* #undef HAVE_LINALG_GEMM3M */

/* Define to 1 if you have an optimized GPU-compatible linear algebra library.
   */
/* #undef HAVE_LINALG_GPU */

/* Define to 1 if you have the MAGMA linear algebra library. */
/* #undef HAVE_LINALG_MAGMA */

/* Define to 1 if you have MAGMA >=1.5 API support */
/* #undef HAVE_LINALG_MAGMA_15 */

/* Define to 1 if you have mkl_?imatcopy extensions. */
/* #undef HAVE_LINALG_MKL_IMATCOPY */

/* Define to 1 if you have mkl_?omatadd extensions. */
/* #undef HAVE_LINALG_MKL_OMATADD */

/* Define to 1 if you have mkl_?omatcopy extensions. */
/* #undef HAVE_LINALG_MKL_OMATCOPY */

/* Define to 1 if you have mkl_*threads extensions. */
/* #undef HAVE_LINALG_MKL_THREADS */

/* Define to 1 if you have the HP MLib Library. */
/* #undef HAVE_LINALG_MLIB */

/* Define to 1 if you have an optimized MPI-parallel linear algebra library.
   */
/* #undef HAVE_LINALG_MPI */

/* Define to 1 if you have an optimized PLASMA linear algebra library. */
/* #undef HAVE_LINALG_PLASMA */

/* Define to 1 if you have an optimized ScaLAPACK linear algebra library. */
/* #undef HAVE_LINALG_SCALAPACK */

/* Define to 1 if you have an optimized serial linear algebra library. */
#define HAVE_LINALG_SERIAL 1

/* Define to 1 if you want to activate workaround for bugged ZDOTC and ZDOTU.
   */
/* #undef HAVE_LINALG_ZDOTC_BUG */

/* Define to 1 if you want to activate workaround for bugged ZDOTC and ZDOTU.
   */
/* #undef HAVE_LINALG_ZDOTU_BUG */

/* Define to 1 if you want to activate LOTF functionality (EXPERIMENTAL). */
/* #undef HAVE_LOTF */

/* Define to 1 if you have the `mallinfo' function. */
#define HAVE_MALLINFO 1

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <malloc/malloc.h> header file. */
/* #undef HAVE_MALLOC_MALLOC_H */

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the <mcheck.h> header file. */
#define HAVE_MCHECK_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you want to enable memory profiling. */
/* #undef HAVE_MEM_PROFILING */

/* Define to 1 if you want to enable MPI support. */
#define HAVE_MPI 1

/* Define to 1 if you have a MPI-1 implementation. */
/* #undef HAVE_MPI1 */

/* Define to 1 if you have a MPI-2 implementation. */
#define HAVE_MPI2 1

/* Define to 1 if you want to activate support for MPI_IN_PLACE. */
/* #undef HAVE_MPI2_INPLACE */

/* Define to 1 if your MPI library supports MPI_IALLREDUCE. */
#define HAVE_MPI_IALLREDUCE 1

/* Define to 1 if your MPI library supports MPI_IALLTOALL. */
#define HAVE_MPI_IALLTOALL 1

/* Define to 1 if your MPI library supports MPI_IALLTOALLV. */
#define HAVE_MPI_IALLTOALLV 1

/* Define to 1 if you are using XLF. */
/* #undef HAVE_MPI_INCLUDED_ONCE */

/* Define to 1 if your MPI library supports MPI_INTEGER16. */
#define HAVE_MPI_INTEGER16 1

/* Define to 1 if you want MPI I/O support. */
#define HAVE_MPI_IO 1

/* Define to 1 if you want to use MPI-IO as default IO library (change the
   default value of iomode) (EXPERIMENTAL). */
/* #undef HAVE_MPI_IO_DEFAULT */

/* Define to 1 if your MPI library supports MPI_TYPE_CREATE_STRUCT. */
#define HAVE_MPI_TYPE_CREATE_STRUCT 1

/* Define to 1 if you have the NetCDF library. */
/* #undef HAVE_NETCDF */

/* Define to 1 if you want to use NetCDF as default IO library (change the
   default value of iomode). */
/* #undef HAVE_NETCDF_DEFAULT */

/* Define to 1 if you have the <netcdf.h> header file. */
/* #undef HAVE_NETCDF_H */

/* Define to 1 if you have MPI-IO support in the NetCDF library. */
/* #undef HAVE_NETCDF_MPI */

/* Define to 1 if you have a standard implementation of NumPy. */
#define HAVE_NUMPY 1

/* Set to 1 if OpenMP has a working implementation of COLLAPSE. */
/* #undef HAVE_OMP_COLLAPSE */

/* Define to 1 if you want to activate support for OpenMP (EXPERIMENTAL). */
/* #undef HAVE_OPENMP */

/* Define to 1 if you are using Linux. */
#define HAVE_OS_LINUX 1

/* Define to 1 if you are using MacOS X. */
/* #undef HAVE_OS_MACOSX */

/* Define to 1 if you are using Windows. */
/* #undef HAVE_OS_WINDOWS */

/* Define to 1 if you have the PAPI library. */
/* #undef HAVE_PAPI */

/* Define to 1 if you have the <papi.h> header file. */
/* #undef HAVE_PAPI_H */

/* Define to 1 if you have the PSML library. */
/* #undef HAVE_PSML */

/* Define to 1 if you have the <stdarg.h> header file. */
#define HAVE_STDARG_H 1

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/malloc.h> header file. */
/* #undef HAVE_SYS_MALLOC_H */

/* Define to 1 if you have the <sys/resource.h> header file. */
#define HAVE_SYS_RESOURCE_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <termios.h> header file. */
#define HAVE_TERMIOS_H 1

/* Define to 1 if you want to use the Abinit timer. */
#define HAVE_TIMER_ABINIT 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you want to activate internal support for TRIQS 1.4. */
/* #undef HAVE_TRIQS_v1_4 */

/* Define to 1 if you want to activate internal support for TRIQS 2.0 (This
   option is dominant over the others versions). */
/* #undef HAVE_TRIQS_v2_0 */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the Wannier90 library. */
/* #undef HAVE_WANNIER90 */

/* Define to 1 if you want to activate support for Wannier90 v1.x (default is
   v2.x). */
/* #undef HAVE_WANNIER90_V1 */

/* Define to 1 if you have the <xc_funcs.h> header file. */
/* #undef HAVE_XC_FUNCS_H */

/* Define to 1 if you have the <xc.h> header file. */
/* #undef HAVE_XC_H */

/* Define to 1 if you have the YAML library. */
/* #undef HAVE_YAML */

/* Define to 1 if assertions should be disabled. */
/* #undef NDEBUG */

/* Name of package */
#define PACKAGE "abinit"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://bugs.launchpad.net/abinit/"

/* Define to the full name of this package. */
#define PACKAGE_NAME "ABINIT"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "ABINIT 8.10.2"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "abinit"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "8.10.2"

/* Define to 1 if you want to tell Abinit to read file lists from standard
   input. */
/* #undef READ_FROM_FILE */

/* The size of `char', as computed by sizeof. */
#define SIZEOF_CHAR 1

/* The size of `double', as computed by sizeof. */
#define SIZEOF_DOUBLE 8

/* The size of `float', as computed by sizeof. */
#define SIZEOF_FLOAT 4

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The size of `long double', as computed by sizeof. */
#define SIZEOF_LONG_DOUBLE 16

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `ptrdiff_t', as computed by sizeof. */
#define SIZEOF_PTRDIFF_T 8

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `size_t', as computed by sizeof. */
#define SIZEOF_SIZE_T 8

/* The size of `unsigned int', as computed by sizeof. */
#define SIZEOF_UNSIGNED_INT 4

/* The size of `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG 8

/* The size of `unsigned long long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG_LONG 8

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define to 1 if you want to enable build of macroave (EXPERIMENTAL). */
#define USE_MACROAVE 1

/* Version number of package */
#define VERSION "8.10.2"

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* *** BEGIN sanity checks *** */

/* MPI options */
#if defined HAVE_MPI 

/* Check that one MPI level is actually defined */
#if ! defined HAVE_MPI1 && ! defined HAVE_MPI2
#error "HAVE_MPI1 and HAVE_MPI2 are both undefined"
#endif

/* Check that only one MPI level has been defined */
#if defined HAVE_MPI1 && defined HAVE_MPI2
#error "HAVE_MPI1 and HAVE_MPI2 are both defined"
#endif

#else /* HAVE_MPI */

/* Check that no MPI level is defined */
#if defined HAVE_MPI1 || defined HAVE_MPI2
#error "HAVE_MPI1 and HAVE_MPI2 must be undefined"
#endif

/* Check that MPI-IO is undefined */
#if defined HAVE_MPI_IO
#error "HAVE_MPI_IO must be undefined"
#endif

#endif /* HAVE_MPI */

/* ETSF_IO support */
#if defined HAVE_ETSF_IO

/* Check that NetCDF is defined */
#if ! defined HAVE_NETCDF
#error "HAVE_NETCDF must but defined for ETSF_IO to work"
#endif

#endif /* HAVE_ETSF_IO */

/* *** END sanity checks *** */

#endif /* _ABINIT_CONFIG_H */
