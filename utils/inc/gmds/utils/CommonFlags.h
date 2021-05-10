/*----------------------------------------------------------------------------*/
/*
 * CommonFlags.h
 *
 *  Created on: April 11  2016
 *
 *  Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef GMDS_COMMON_FLAGS_H_
#define GMDS_COMMON_FLAGS_H_
/*----------------------------------------------------------------------------*/
#ifdef WIN32
#ifdef DLLEXPORT
#define EXPORT_GMDS __declspec(dllexport)
#else
#define EXPORT_GMDS __declspec(dllimport)
#endif //DLLEXPORT
#else
#define EXPORT_GMDS
#endif //WIN32
/*----------------------------------------------------------------------------*/
/**
 * \def GMDS_DEBUG
 * \brief This macro is set when compiling in debug mode
 *
 */
#ifdef NDEBUG
#undef GMDS_DEBUG
#else
#define GMDS_DEBUG
#endif
/*----------------------------------------------------------------------------*/
/**
 * \def GMDS_OS_LINUX
 * \brief This macro is set on Linux systems
 *
 * \def GMDS_OS_UNIX
 * \brief This macro is set on Unix systems (Android included).
 *
 * \def GMDS_OS_WINDOWS
 * \brief This macro is set on Windows systems.
 *
 * \def GMDS_OS_APPLE
 * \brief This macro is set on Apple systems.
 */
#if defined(__linux__)
#define GMDS_OS_LINUX
#endif
#elif defined(WIN32) || defined(_WIN64)

#define GMDS_OS_WINDOWS

#if defined(_OPENMP)
#define GMDS_OPENMP
#endif

#if defined(_MSC_VER)
#define GMDS_COMPILER_MSVC
#else
#error "Unsupported compiler"
#endif

#if defined(_WIN64)
#define GMDS_ARCH_64
#else
#define GMDS_ARCH_32
#endif


#if defined(__APPLE__)
#define GMDS_OS_APPLE
#define GMDS_OS_UNIX
#endif


/*----------------------------------------------------------------------------*/
/**
 * \def GMDS_ARCH_32
 * \brief This macro is set if the current system is a 32 bits architecture.
 *
 * \def GMDS_ARCH_64
 * \brief This macro is set if the current system is a 64 bits architecture.
 */
// On linux Works with GCC and ICC
#if defined(__x86_64) || defined(__ppc64__)
#define GMDS_ARCH_64
#else
#define GMDS_ARCH_32
#endif

/*----------------------------------------------------------------------------*/
/**
 * \def GMDS_OPENMP
 * \brief This macro is set if OpenMP is supported on the current system.
 */
#if defined(_OPENMP)
#  define GMDS_OPENMP
#endif
/*----------------------------------------------------------------------------*/
/**
 * \def GMDS_COMPILER_GCC
 * \brief This macro is set if the source code is compiled with GNU's gcc.
 *
 * \def GMDS_COMPILER_INTEL
 * \brief This macro is set if the source code is compiled with Intel's icc.
 *
 * \def GMDS_COMPILER_CLANG
 * \brief This macro is set if the source code is compiled with CLang
 * Visual C++.
 */
#if defined(__INTEL_COMPILER)
#define GMDS_COMPILER_INTEL
#elif defined(__clang__)
#define GMDS_COMPILER_CLANG
#elif defined(__GNUC__)
#define GMDS_COMPILER_GCC
#elif defined(__APPLE__)
#define GMDS_COMPILER_CLANG
#elif defined(__OSX__)
#define GMDS_COMPILER_CLANG
#else
#error "GMDS - Unsupported compiler"
#endif
/*----------------------------------------------------------------------------*/
#endif //GMDS_COMMON_FLAGS_H_
/*----------------------------------------------------------------------------*/

