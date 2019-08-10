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
/*
 * This file contains interfaces to external ordering packages.
 * At the moment, PORD (J. Schulze) and SCOTCH are interfaced.
 */
#include "mumps_orderings.h"
#include "mumps_c_types.h"
#if defined(pord)
/* Interface to PORD */
/*int mumps_pord( int, int, int *, int *, int * );
#define MUMPS_PORDF   \
F_SYMBOL(pordf,PORDF)*/
void MUMPS_CALL
MUMPS_PORDF( MUMPS_INT *nvtx, MUMPS_INT *nedges,
             MUMPS_INT *xadj, MUMPS_INT *adjncy,
             MUMPS_INT *nv, MUMPS_INT *ncmpa )
{
    *ncmpa = mumps_pord( *nvtx, *nedges, xadj, adjncy, nv );
}
/* Interface to PORD with weighted graph*/
/*int mumps_pord_wnd( int, int, int *, int *, int *, int * );
#define MUMPS_PORDF_WND           \
    F_SYMBOL(pordf_wnd,PORDF_WND)*/
void MUMPS_CALL
MUMPS_PORDF_WND( MUMPS_INT *nvtx, MUMPS_INT *nedges,
                 MUMPS_INT *xadj, MUMPS_INT *adjncy,
                 MUMPS_INT *nv, MUMPS_INT *ncmpa, MUMPS_INT *totw )
{
    *ncmpa = mumps_pord_wnd( *nvtx, *nedges, xadj, adjncy, nv, totw );
}
/************************************************************
 mumps_pord is used in ana_aux.F
        permutation and inverse permutation not set in output,
        but may be printed in default file: "perm_pord" and "iperm_pord",
        if associated part uncommneted.
        But, if uncommetnted a bug occurs in psl_ma41_analysi.F
******************************************************************/
/*********************************************************/
MUMPS_INT mumps_pord
(
   MUMPS_INT nvtx,
   MUMPS_INT nedges,
   MUMPS_INT *xadj_pe,
   MUMPS_INT *adjncy,
   MUMPS_INT *nv
)
{
/**********************************
Argument Comments:
input:
-----
- nvtx          : dimension of the Problem (N)
- nedges        : number of entries (NZ)
- adjncy        : non-zeros entries (IW input)
input/output:
-------------
- xadj_pe       : pointer through beginning of column non-zeros entries (PTRAR)
- on exit, "father array" (PE)
ouput:
------
- nv            : "nfront array" (NV)
*************************************/
  graph_t    *G;
  elimtree_t *T;
  timings_t  cpus[12];
  options_t  options[] = { SPACE_ORDTYPE, SPACE_NODE_SELECTION1,
                    SPACE_NODE_SELECTION2, SPACE_NODE_SELECTION3,
                    SPACE_DOMAIN_SIZE, 0 };
  MUMPS_INT *ncolfactor, *ncolupdate, *parent, *vtx2front;
  MUMPS_INT *first, *link, nfronts, J, K, u, vertex, vertex_root, count;
      /**************************************************
       declaration to uncomment if printing ordering
      ***************************************************
         FILE *fp1, *fp2;
         int  *perm,  *iperm;
      */
/*** decalage des indices couteux dans un premier temps:
****  A modifier dans une version ulterieure de MA41GD  */
  for (u = nvtx; u >= 0; u--)
   {
     xadj_pe[u] = xadj_pe[u] - 1;
   }
   for (K = nedges-1; K >= 0; K--)
   {
      adjncy[K] = adjncy[K] - 1;
   }
 /* initialization of the graph */
   mymalloc(G, 1, graph_t);
   G->xadj   = xadj_pe;
   G->adjncy = adjncy;
   mymalloc(G->vwght, nvtx, MUMPS_INT);
   G->nvtx = nvtx;
   G->nedges = nedges;
   G->type = UNWEIGHTED;
   G->totvwght = nvtx;
   for (u = 0; u < nvtx; u++)
     G->vwght[u] = 1;
  /* main function of the Ordering */
   T = SPACE_ordering(G, options, cpus);
   nfronts = T->nfronts;
   ncolfactor = T->ncolfactor;
   ncolupdate = T->ncolupdate;
   parent = T->parent;
  /*    firstchild = T->firstchild; */
   vtx2front = T->vtx2front;
    /* -----------------------------------------------------------
     store the vertices/columns of a front in a bucket structure
     ----------------------------------------------------------- */
   mymalloc(first, nfronts, MUMPS_INT);
   mymalloc(link, nvtx, MUMPS_INT);
   for (J = 0; J < nfronts; J++)
      first[J] = -1;
   for (u = nvtx-1; u >= 0; u--)
      {
        J = vtx2front[u];
        link[u] = first[J];
        first[J] = u;
      }
  /* -----------------------------------------------------------
     fill the two arrays corresponding to the MUMPS tree structure
     ----------------------------------------------------------- */
  count = 0;
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
     {
       vertex_root = first[K];
       if (vertex_root == -1)
          {
            /* JY: I think this cannot happen */
            printf(" Internal error in mumps_pord (cf JY), %d\n",K);
            exit(-1);
          }
       /* for the principal column of the supervariable */
       if (parent[K] == -1)
          xadj_pe[vertex_root] = 0; /* root of the tree */
       else
          xadj_pe[vertex_root] = - (first[parent[K]]+1);
          nv[vertex_root] = ncolfactor[K] + ncolupdate[K];
          count++;
       for (vertex = link[vertex_root]; vertex != -1; vertex = link[vertex])
        /* for the secondary columns of the supervariable */
       {
         xadj_pe[vertex] = - (vertex_root+1);
         nv[vertex] = 0;
         count++;
        }
  }
  /* ----------------------
     free memory and return
     ---------------------- */
  free(first); free(link);
  free(G->vwght);
  free(G);
  freeElimTree(T);
  return (0);
}
/*********************************************************/
MUMPS_INT mumps_pord_wnd
(
        MUMPS_INT nvtx,
        MUMPS_INT nedges,
        MUMPS_INT *xadj_pe,
        MUMPS_INT *adjncy,
        MUMPS_INT *nv,
        MUMPS_INT *totw
)
{
/**********************************
Argument Comments:
input:
-----
- nvtx   : dimension of the Problem (N)
- nedges : number of entries (NZ)
- adjncy : non-zeros entries (IW input)
- totw   : sum of the weigth of the vertices
input/output:
-------------
- xadj_pe : pointer through beginning of column non-zeros entries (PTRAR)
- on exit, "father array" (PE)
ouput:
------
- nv      : weight of the vertices
- on exit "nfront array" (NV)
*************************************/
        graph_t    *G;
        elimtree_t *T;
        timings_t  cpus[12];
        options_t  options[] = { SPACE_ORDTYPE, SPACE_NODE_SELECTION1,
                    SPACE_NODE_SELECTION2, SPACE_NODE_SELECTION3,
                    SPACE_DOMAIN_SIZE, 0 };
        MUMPS_INT *ncolfactor, *ncolupdate, *parent, *vtx2front;
        MUMPS_INT *first, *link, nfronts, J, K, u, vertex, vertex_root, count;
      /**************************************************
       declaration to uncomment if printing ordering
      ***************************************************
         FILE *fp1, *fp2;
         int  *perm,  *iperm;
      */
/*** decalage des indices couteux dans un premier temps:
****  A modifier dans une version ulterieure de MA41GD  */
        for (u = nvtx; u >= 0; u--)
        {
          xadj_pe[u] = xadj_pe[u] - 1;
        }
        for (K = nedges-1; K >= 0; K--)
        {
          adjncy[K] = adjncy[K] - 1;
        }
 /* initialization of the graph */
        mymalloc(G, 1, graph_t);
        G->xadj  = xadj_pe;
        G->adjncy= adjncy;
        mymalloc(G->vwght, nvtx, MUMPS_INT);
        G->nvtx = nvtx;
        G->nedges = nedges;
        G->type = WEIGHTED;
        G->totvwght = (*totw);
        for (u = 0; u < nvtx; u++)
          G->vwght[u] = nv[u];
  /* main function of the Ordering */
        T = SPACE_ordering(G, options, cpus);
        nfronts = T->nfronts;
        ncolfactor = T->ncolfactor;
        ncolupdate = T->ncolupdate;
        parent = T->parent;
  /*    firstchild = T->firstchild; */
        vtx2front = T->vtx2front;
    /* -----------------------------------------------------------
     store the vertices/columns of a front in a bucket structure
     ----------------------------------------------------------- */
        mymalloc(first, nfronts, MUMPS_INT);
        mymalloc(link, nvtx, MUMPS_INT);
        for (J = 0; J < nfronts; J++)
          first[J] = -1;
        for (u = nvtx-1; u >= 0; u--)
        {
          J = vtx2front[u];
          link[u] = first[J];
          first[J] = u;
        }
  /* -----------------------------------------------------------
     fill the two arrays corresponding to the MUMPS tree structure
     ----------------------------------------------------------- */
  count = 0;
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
     {
        vertex_root = first[K];
        if (vertex_root == -1)
          {
            /* JY: I think this cannot happen */
            printf(" Internal error in mumps_pord (cf JY), %d\n",K);
            exit(-1);
          }
         /* for the principal column of the supervariable */
        if (parent[K] == -1)
          xadj_pe[vertex_root] = 0; /* root of the tree */
        else
          xadj_pe[vertex_root] = - (first[parent[K]]+1);
          nv[vertex_root] = ncolfactor[K] + ncolupdate[K];
          count++;
          for (vertex = link[vertex_root]; vertex != -1; vertex = link[vertex])
          /* for the secondary columns of the supervariable */
            {
              xadj_pe[vertex] = - (vertex_root+1);
              nv[vertex] = 0;
              count++;
        }
  }
  /* ----------------------
     free memory and return
     ---------------------- */
  free(first); free(link);
  free(G->vwght);
  free(G);
  freeElimTree(T);
  return (0);
}
#endif /* pord */
/************************************************************/
#if defined(scotch) || defined(ptscotch)
/*int esmumps( const int n, const int iwlen, int * const pe, const int pfree,
             int * const len, int * const iw, int * const nv, int * const elen,
             int * const last);*/
/* Fortran interface to SCOTCH */
/*#define MUMPS_SCOTCH    \
  F_SYMBOL(scotch,SCOTCH)*/
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
              MUMPS_INT * const        ncmpa )
{
     *ncmpa = esmumps( *n, *iwlen, petab, *pfree,
                       lentab, iwtab, nvtab, elentab, lasttab );
}
#endif /* scotch */
#if defined(ptscotch)
/*#include "mpi.h"
#include <stdio.h>
#include "ptscotch.h"
int mumps_dgraphinit( SCOTCH_Dgraph *, MPI_Fint *, MPI_Fint *);
#define MUMPS_DGRAPHINIT        \
F_SYMBOL(dgraphinit,DGRAPHINIT)*/
void MUMPS_CALL
MUMPS_DGRAPHINIT(SCOTCH_Dgraph *graphptr, MPI_Fint *comm, MPI_Fint *ierr)
{
  MPI_Comm  int_comm;
  int_comm = MPI_Comm_f2c(*comm);
  *ierr = SCOTCH_dgraphInit(graphptr, int_comm);
  return;
}
#endif
#if defined(parmetis) || defined(parmetis3)
/*PARMETIS*/
#include "parmetis.h"
void MUMPS_CALL
MUMPS_PARMETIS(MUMPS_INT *first,      MUMPS_INT *vertloctab, 
               MUMPS_INT *edgeloctab, MUMPS_INT *numflag, 
               MUMPS_INT *options,    MUMPS_INT *order, 
               MUMPS_INT *sizes,      MUMPS_INT *comm,
               MUMPS_INT *ierr)
{
  MPI_Comm  int_comm;
  int iierr;
  int_comm = MPI_Comm_f2c(*comm);
#if defined(parmetis)
  *ierr=0;
  iierr=ParMETIS_V3_NodeND(first, vertloctab, edgeloctab, numflag, options, order, sizes, &int_comm);
  if(iierr != METIS_OK)
    *ierr=1;
#else
  ParMETIS_V3_NodeND(first, vertloctab, edgeloctab, numflag, options, order, sizes, &int_comm);
#endif
  return;
}
#endif
