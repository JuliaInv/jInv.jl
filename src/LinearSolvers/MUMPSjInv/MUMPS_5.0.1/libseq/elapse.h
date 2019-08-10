/*
 *
 *  This file is part of MUMPS 5.0.1, released
 *  on Thu Jul 23 17:08:29 UTC 2015
 *
 *
 *  Copyright 1991-2015 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license:
 *  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
 *
 */

#ifndef MUMPS_CALL
#if defined(_WIN32)
/* Modify/choose between next 2 lines depending
 *  * on your Windows calling conventions */
/* #define MUMPS_CALL __stdcall */
#define MUMPS_CALL
#else
#define MUMPS_CALL
#endif
#endif

#if (defined(_WIN32) && ! defined(__MINGW32__)) || defined(UPPER)
#define mumps_elapse MUMPS_ELAPSE
#elif defined(Add__)
#define mumps_elapse mumps_elapse__
#elif defined(Add_)
#define mumps_elapse mumps_elapse_
#endif

void MUMPS_CALL mumps_elapse(double *val);
