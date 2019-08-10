C
C  This file is part of MUMPS 5.0.1, released
C  on Thu Jul 23 17:08:29 UTC 2015
C
C
C  Copyright 1991-2015 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license:
C  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
C
      SUBROUTINE ZMUMPS_SOL_R(N, A, LA, IW, LIW, WCB, LWCB,
     &    NRHS,
     &    PTRICB, IWCB, LIWCB, 
     &    RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD, 
     &    NE_STEPS, NA, LNA, STEP,
     &    FRERE, DAD, FILS,
     &    NSTK_S, IPOOL, LPOOL, PTRIST, PTRFAC, MYLEAF, INFO,
     &    KEEP,KEEP8,
     &    PROCNODE_STEPS,
     &    SLAVEF, COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,
     &    RHS_ROOT, LRHS_ROOT, MTYPE, 
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE 
     &    )
      USE ZMUMPS_OOC
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER(8) :: LA
      INTEGER, INTENT(IN) :: N, LIW, LWCB, LPOOL, LIWCB, LNA
      INTEGER SLAVEF, MYLEAF, COMM, MYID
      INTEGER INFO( 40 ), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER NRHS
      COMPLEX(kind=8) A( LA ), WCB( LWCB )
      INTEGER LRHS_ROOT
      COMPLEX(kind=8) RHS_ROOT( LRHS_ROOT )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER NA( LNA ), NE_STEPS( KEEP(28) )
      INTEGER STEP( N ), FRERE( KEEP(28) ), FILS( N ),
     &        DAD( KEEP(28) )
      INTEGER NSTK_S(KEEP(28)), IPOOL( LPOOL )
      INTEGER PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER PTRICB( KEEP(28) ) 
      LOGICAL, intent(in) :: DO_NBSPARSE
      INTEGER, intent(in) :: LRHS_BOUNDS
      INTEGER, intent(in) :: RHS_BOUNDS(LRHS_BOUNDS)
      INTEGER IW( LIW ), IWCB( LIWCB )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER, intent(in) ::  POSINRHSCOMP_FWD(N), LRHSCOMP 
#if defined(RHSCOMP_BYROWS)
      COMPLEX(kind=8), intent(inout) :: RHSCOMP(NRHS,LRHSCOMP)
#else
      COMPLEX(kind=8), intent(inout) :: RHSCOMP(LRHSCOMP,NRHS)
#endif
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MSGTAG, MSGSOU, DUMMY(1)
      LOGICAL FLAG
      INTEGER NBFIN, MYROOT
      INTEGER POSIWCB,POSWCB,PLEFTWCB
      INTEGER INODE
      INTEGER I
      INTEGER III, NBROOT,LEAF
      LOGICAL BLOQ
      EXTERNAL MUMPS_PROCNODE
      INTEGER MUMPS_PROCNODE
      DUMMY(1) = 1
      POSIWCB = LIWCB
      POSWCB  = LWCB
      PLEFTWCB= 1
      DO I = 1, KEEP(28)
        NSTK_S(I)   = NE_STEPS(I)
      ENDDO
      PTRICB = 0
      CALL MUMPS_INIT_POOL_DIST(N, LEAF, MYID,
     &     SLAVEF, NA, LNA, KEEP,KEEP8, STEP,
     &     PROCNODE_STEPS, IPOOL, LPOOL)
      CALL MUMPS_INIT_NROOT_DIST(N, NBROOT, MYROOT, MYID,
     &     SLAVEF, NA, LNA, KEEP, STEP,
     &     PROCNODE_STEPS)
      NBFIN = SLAVEF
      IF ( MYROOT .EQ. 0 ) THEN
        NBFIN = NBFIN - 1
        CALL ZMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &       RACINE_SOLVE, SLAVEF)
        IF (NBFIN.EQ.0) GOTO 260
      END IF
      MYLEAF = LEAF - 1
      III    = 1
   50 CONTINUE
      IF (SLAVEF .EQ. 1) THEN
         CALL ZMUMPS_GET_INODE_FROM_POOL
     &        ( IPOOL(1), LPOOL, III, LEAF, INODE,
     &          KEEP(208) )
        GOTO 60
      ENDIF
      BLOQ = ( ( III .EQ. LEAF )
     &     )
      CALL ZMUMPS_SOLVE_RECV_AND_TREAT( BLOQ, FLAG,
     &     BUFR, LBUFR, LBUFR_BYTES,
     &     MYID, SLAVEF, COMM,
     &     N, NRHS, IPOOL, LPOOL, III, LEAF,
     &     NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &     IWCB, LIWCB,
     &     WCB, LWCB, POSWCB,
     &     PLEFTWCB, POSIWCB,
     &     PTRICB, INFO, KEEP,KEEP8, STEP,
     &     PROCNODE_STEPS,
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &     )
      IF ( INFO( 1 ) .LT. 0 .OR. NBFIN .EQ. 0 ) GOTO 260
      IF (.not. FLAG) THEN
         IF (III .NE. LEAF) THEN
            CALL ZMUMPS_GET_INODE_FROM_POOL
     &           (IPOOL(1), LPOOL, III, LEAF, INODE,
     &           KEEP(208) )
            GOTO 60
         ENDIF                  
      ENDIF                     
      GOTO 50
 60   CONTINUE
      CALL ZMUMPS_SOLVE_NODE( INODE, BUFR, LBUFR, LBUFR_BYTES,
     &        MSGTAG, MSGSOU, MYID, SLAVEF, COMM,  N,
     &        IPOOL, LPOOL, III, LEAF, NBFIN, NSTK_S,
     &        IWCB, LIWCB, WCB, LWCB, A, LA,
     &        IW, LIW, NRHS, 
     &        POSWCB, PLEFTWCB, POSIWCB,
     &        PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     &        FILS, STEP, FRERE, DAD,
     &        MYROOT, INFO, KEEP,KEEP8, RHS_ROOT, MTYPE, 
     &        RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &        , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE 
     &     )
      IF ( INFO( 1 ) .LT. 0 .OR. NBFIN .EQ. 0 ) GOTO 260
      GOTO 50
  260 CONTINUE
      CALL ZMUMPS_FINISH_RECV( MYID,COMM,BUFR,
     &                            LBUFR,LBUFR_BYTES )
      RETURN
      END SUBROUTINE ZMUMPS_SOL_R
