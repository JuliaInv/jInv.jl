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
#ifdef INTSIZE64
#include <inttypes.h>
#define LIBSEQ_INT int64_t
#else
#define LIBSEQ_INT int
#endif

#ifndef MUMPS_MPI_H
#define MUMPS_MPI_H

/* We define all symbols as extern "C" for users who call MUMPS with its
   libseq from a C++ driver. */
#ifdef __cplusplus
extern "C" {
#endif

/* This is the minimum to have the C interface of MUMPS work.
 * Most of the time, users who need this file have no call to MPI functions in
 * their own code. Hence it is not worth declaring all MPI functions here.
 * However if some users come to request some more stub functions of the MPI
 * standards, we may add them. But it is not worth doing it until then. */

typedef LIBSEQ_INT MPI_Comm; /* Simple type for MPI communicator */
static MPI_Comm MPI_COMM_WORLD=(MPI_Comm)0;

LIBSEQ_INT MPI_Init(LIBSEQ_INT *pargc, char ***pargv);
LIBSEQ_INT MPI_Comm_rank(LIBSEQ_INT  comm, LIBSEQ_INT  *rank);
LIBSEQ_INT MPI_Finalize(void);

#ifdef __cplusplus
}
#endif

#endif /* MUMPS_MPI_H */
