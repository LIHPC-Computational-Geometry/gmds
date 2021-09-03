/*----------------------------------------------------------------------------*/
/*
 * KFlags.h
 *
 *  Created on: April 11  2016
 *
 *  Author: ledouxf
 */
/*----------------------------------------------------------------------------*/
#ifndef KMDS_FLAGS_H_
#define KMDS_FLAGS_H_
/*----------------------------------------------------------------------------*/
#ifdef WIN32
//#ifdef DLLEXPORT
#define EXPORT_KMDS __declspec(dllexport)
//#else
//#define EXPORT_KMDS __declspec(dllimport)
//#endif //DLLEXPORT
#else
#define EXPORT_KMDS
#endif  // WIN32
/*----------------------------------------------------------------------------*/
/**
 * \def KMDS_DEBUG
 * \brief This macro is set when compiling in debug mode
 *
 */
#ifdef NDEBUG
#undef KMDS_DEBUG
#else
#define KMDS_DEBUG
#endif
/*----------------------------------------------------------------------------*/
/**
 * \def KMDS_OS_LINUX
 * \brief This macro is set on Linux systems
 *
 * \def KMDS_OS_UNIX
 * \brief This macro is set on Unix systems (Android included).
 *
 * \def KMDS_OS_WINDOWS
 * \brief This macro is set on Windows systems.
 *
 * \def KMDS_OS_APPLE
 * \brief This macro is set on Apple systems.
 */
#if defined(__linux__)
#define KMDS_OS_LINUX
#endif
#elif defined(WIN32) || defined(_WIN64)

#define KMDS_OS_WINDOWS

#if defined(_OPENMP)
#define KMDS_OPENMP
#endif

#if defined(_MSC_VER)
#define KMDS_COMPILER_MSVC
#else
#error "Unsupported compiler"
#endif

#if defined(_WIN64)
#define KMDS_ARCH_64
#else
#define KMDS_ARCH_32
#endif

#if defined(__APPLE__)
#define KMDS_OS_APPLE
#define KMDS_OS_UNIX
#endif

/*----------------------------------------------------------------------------*/
/**
 * \def KMDS_ARCH_32
 * \brief This macro is set if the current system is a 32 bits architecture.
 *
 * \def KMDS_ARCH_64
 * \brief This macro is set if the current system is a 64 bits architecture.
 */
// On linux Works with GCC and ICC
#if defined(__x86_64) || defined(__ppc64__)
#define KMDS_ARCH_64
#else
#define KMDS_ARCH_32
#endif

/*----------------------------------------------------------------------------*/
/**
 * \def KMDS_OPENMP
 * \brief This macro is set if OpenMP is supported on the current system.
 */
#if defined(_OPENMP)
#define KMDS_OPENMP
#endif
/*----------------------------------------------------------------------------*/
/**
 * \def KMDS_COMPILER_GCC
 * \brief This macro is set if the source code is compiled with GNU's gcc.
 *
 * \def KMDS_COMPILER_INTEL
 * \brief This macro is set if the source code is compiled with Intel's icc.
 *
 * \def KMDS_COMPILER_CLANG
 * \brief This macro is set if the source code is compiled with CLang
 * Visual C++.
 */
#if defined(__INTEL_COMPILER)
#define KMDS_COMPILER_INTEL
#elif defined(__clang__)
#define KMDS_COMPILER_CLANG
#elif defined(__GNUC__)
#define KMDS_COMPILER_GCC
#elif defined(__APPLE__)
#define KMDS_COMPILER_CLANG
#elif defined(__OSX__)
#define KMDS_COMPILER_CLANG
#else
#error "KMDS - Unsupported compiler"
#endif
/*----------------------------------------------------------------------------*/
#endif  // KMDS_FLAGS_H_
/*----------------------------------------------------------------------------*/
