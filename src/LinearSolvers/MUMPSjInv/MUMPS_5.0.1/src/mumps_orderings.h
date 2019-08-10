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
#ifndef MUMPS_ORDERINGS_H
#define MUMPS_ORDERINGS_H
#include "mumps_common.h"
#include "mumps_c_types.h"
#if defined(pord)
#include <space.h>
MUMPS_INT mumps_pord( MUMPS_INT, MUMPS_INT, MUMPS_INT *, MUMPS_INT *, MUMPS_INT * );
#define MUMPS_PORDF \
    F_SYMBOL(pordf,PORDF)
void MUMPS_CALL
MUMPS_PORDF( MUMPS_INT *nvtx, MUMPS_INT *nedges,
             MUMPS_INT *xadj, MUMPS_INT *adjncy,
             MUMPS_INT *nv, MUMPS_INT *ncmpa );
MUMPS_INT mumps_pord_wnd( MUMPS_INT, MUMPS_INT, MUMPS_INT *, MUMPS_INT *, MUMPS_INT *, MUMPS_INT * );
#define MUMPS_PORDF_WND          \
    F_SYMBOL(pordf_wnd,PORDF_WND)
void MUMPS_CALL
MUMPS_PORDF_WND( MUMPS_INT *nvtx, MUMPS_INT *nedges,
                 MUMPS_INT *xadj, MUMPS_INT *adjncy,
                 MUMPS_INT *nv, MUMPS_INT *ncmpa, MUMPS_INT *totw );
#endif /*PORD*/
#if defined(scotch) || defined(ptscotch)
MUMPS_INT esmumps( const MUMPS_INT n, const MUMPS_INT iwlen, MUMPS_INT * const pe, const MUMPS_INT pfree,
             MUMPS_INT * const len, MUMPS_INT * const iw, MUMPS_INT * const nv, MUMPS_INT * const elen,
             MUMPS_INT * const last);
#define MUMPS_SCOTCH        \
    F_SYMBOL(scotch,SCOTCH)
void MUMPS_CALL
MUMPS_SCOTCH( const MUMPS_INT * const  n,
              const MUMPS_INT * const  iwlen,
              MUMPS_INT * const        petab,
              const MUMPS_INT * const  pfree,
              MUMPS_INT * const        lentab,
              MUMPS_INT * const        iwtab,
              MUMPS_INT * const        nvtab,
              MUMPS_INT * const        elentab,
              MUMPS_INT * const        lasttab,
              MUMPS_INT * const        ncmpa );
#endif /*scotch or ptscotch*/
#if defined(ptscotch)
#include "mpi.h"
#include <stdio.h>
#include "ptscotch.h"
#define MUMPS_DGRAPHINIT \
  F_SYMBOL(dgraphinit,DGRAPHINIT)
void MUMPS_CALL
MUMPS_DGRAPHINIT(SCOTCH_Dgraph *graphptr, MPI_Fint *comm, MPI_Fint *ierr);
#endif /*ptscotch*/
#if defined(parmetis) || defined(parmetis3)
#include "mpi.h"
#include "parmetis.h"
#define MUMPS_PARMETIS \
  F_SYMBOL(parmetis,PARMETIS)
void MUMPS_CALL
MUMPS_PARMETIS(MUMPS_INT *first,      MUMPS_INT *vertloctab, 
               MUMPS_INT *edgeloctab, MUMPS_INT *numflag, 
               MUMPS_INT *options,    MUMPS_INT *order, 
               MUMPS_INT *sizes,      MUMPS_INT *comm,
               MUMPS_INT *ierr);
#endif /*PARMETIS*/
#endif /* MUMPS_ORDERINGS_H */
