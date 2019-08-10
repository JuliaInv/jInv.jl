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
#include <mpi.h>
LIBSEQ_INT MPI_Init(LIBSEQ_INT *pargc, char ***pargv)
{
  return 0;
}

LIBSEQ_INT MPI_Comm_rank( MPI_Comm comm, LIBSEQ_INT *rank)
{
  *rank=0;
  return 0;
}
LIBSEQ_INT MPI_Finalize(void)
{
   return 0;
}

