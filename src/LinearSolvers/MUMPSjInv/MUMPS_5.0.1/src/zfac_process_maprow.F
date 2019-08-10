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
      RECURSIVE SUBROUTINE ZMUMPS_MAPLIG( COMM_LOAD, ASS_IRECV,
     &  BUFR, LBUFR, LBUFR_BYTES,
     &
     &  INODE_PERE, ISON, NSLAVES_PERE, LIST_SLAVES_PERE,
     &  NFRONT_PERE, NASS_PERE, NFS4FATHER,LMAP, TROW,
     &  PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &  LRLUS, N, IW,
     &  LIW, A, LA,
     &  PTRIST, PTLUST, PTRFAC,
     &  PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &  IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     &  NBFIN, ICNTL, KEEP,KEEP8,DKEEP,
     &  root, OPASSW, OPELIW,
     &  ITLOC, RHS_MUMPS,
     &  FILS, PTRARW, PTRAIW, INTARR, DBLARR, ND, FRERE,
     &  LPTRAR, NELT, FRTPTR, FRTELT, 
     &
     &  ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &  )
      USE ZMUMPS_COMM_BUFFER
      USE ZMUMPS_LOAD
#if ! defined(NO_FDM_MAPROW)
      USE MUMPS_FAC_MAPROW_DATA_M
#endif
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
#if ! defined(NO_FDM_MAPROW)
#endif
      TYPE (ZMUMPS_ROOT_STRUC ) :: root
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER ICNTL( 40 ), KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION    DKEEP(130)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER BUFR( LBUFR )
      INTEGER SLAVEF, NBFIN
      INTEGER(8) :: LA, IPTRLU, LRLU, LRLUS, POSFAC
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      COMPLEX(kind=8) A( LA )
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER COMP
      INTEGER NSTK( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM, MYID
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER INODE_PERE, ISON
      INTEGER NFS4FATHER
      INTEGER NBROWS_ALREADY_SENT
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE
      INTEGER LIST_SLAVES_PERE( * )
      INTEGER LMAP 
      INTEGER TROW( LMAP )
      DOUBLE PRECISION OPASSW, OPELIW
      COMPLEX(kind=8) DBLARR(max(1,KEEP(13)))
      INTEGER INTARR(max(1,KEEP(14)))
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N+KEEP(253) ), FILS( N )
      COMPLEX(kind=8) :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER NOSLA, I
      INTEGER I_POSMYIDIN_PERE
      INTEGER INDICE_PERE
      INTEGER PDEST, PDEST_MASTER
      LOGICAL :: LOCAL_ASSEMBLY_TO_BE_DONE
      INTEGER NROWS_TO_SEND
      INTEGER PDEST_MASTER_ISON, IPOS_IN_SLAVE
      LOGICAL DESCLU, SLAVE_ISON
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED 
      INTEGER MSGSOU, MSGTAG
      INTEGER LP
      LOGICAL COMPRESSCB
      LOGICAL IS_ERROR_BROADCASTED, IS_ofType5or6
      INTEGER ITYPE, TYPESPLIT
      INTEGER KEEP253_LOC
#if ! defined(NO_FDM_MAPROW)
      INTEGER :: INFO_TMP(2)
#endif
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE, MUMPS_TYPESPLIT
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, MUMPS_TYPESPLIT
      INTEGER LMAP_LOC, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBROW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SLAVES_PERE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAP, PERM
      IS_ERROR_BROADCASTED = .FALSE.
      TYPESPLIT = MUMPS_TYPESPLIT(PROCNODE_STEPS(STEP(INODE_PERE)),
     &                  SLAVEF)
      IS_ofType5or6 = ((TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6))
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
#if ! defined(NO_FDM_MAPROW)
#endif
      ALLOCATE(SLAVES_PERE(0:max(1,NSLAVES_PERE)), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP > 0) write(LP,*) MYID,
     &  ' : PB allocation SLAVES_PERE in ZMUMPS_MAPLIG'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      endif
      IF (NSLAVES_PERE.GT.0) 
     &SLAVES_PERE(1:NSLAVES_PERE) = LIST_SLAVES_PERE(1:NSLAVES_PERE)
      SLAVES_PERE(0) = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(INODE_PERE)),
     &                 SLAVEF )
      ALLOCATE(NBROW(0:NSLAVES_PERE), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP>0) write(LP,*) MYID,
     &  ' : PB allocation NBROW in ZMUMPS_MAPLIG'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 670 
      endif
      LMAP_LOC = LMAP
      ALLOCATE(MAP(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP>0) THEN
        write(LP,*) MYID, ' : PB allocation LMAP in ZMUMPS_MAPLIG'
        ENDIF
        IFLAG  =-13
        IERROR = LMAP
        GOTO 680
      endif
      MAP( 1 : LMAP ) = TROW( 1 : LMAP )
      PDEST_MASTER_ISON = MUMPS_PROCNODE(PROCNODE_STEPS(STEP(ISON)),
     &                    SLAVEF)
      SLAVE_ISON = PDEST_MASTER_ISON .NE. MYID
      IF (SLAVE_ISON) THEN
        IF ( PTRIST(STEP( ISON )) .EQ. 0 ) THEN
          CALL ZMUMPS_TREAT_DESCBAND( ISON, COMM_LOAD,
     &    ASS_IRECV, 
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE, LPTRAR,
     &    NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &    )
          IF ( IFLAG .LT. 0 ) THEN
            IS_ERROR_BROADCASTED = .TRUE.  
            GOTO 670                       
          ENDIF
        END IF
#if ! defined(NO_FDM_MAPROW)
        IF ( ( IW( PTRIST(STEP(ISON)) + 1 + KEEP(IXSZ) ) .NE.
     &       IW( PTRIST(STEP(ISON)) + 3 + KEEP(IXSZ) ) ) .OR.
     &     ( KEEP(50) .NE. 0 .AND.
     &       IW( PTRIST(STEP(ISON)) + 6 + KEEP(IXSZ) ) .NE. 0 ) )
     &  THEN
          INFO_TMP=0 
          CALL MUMPS_FMRD_SAVE_MAPROW(
     &         IW(PTRIST(STEP(ISON))+XXA),
     &         INODE_PERE, ISON, NSLAVES_PERE, NFRONT_PERE,
     &         NASS_PERE, LMAP, NFS4FATHER,
     &         SLAVES_PERE(1:NSLAVES_PERE),
     &         MAP,        
     &         INFO_TMP)
               IF (INFO_TMP(1) < 0) THEN
                 IFLAG = INFO_TMP(1)
                 IERROR = INFO_TMP(2)
               ENDIF
          GOTO 670 
        ELSE
          GOTO 10
        ENDIF
#endif
        DO WHILE (
     &     ( IW( PTRIST(STEP(ISON)) + 1 + KEEP(IXSZ) ) .NE.
     &       IW( PTRIST(STEP(ISON)) + 3 + KEEP(IXSZ) ) ) .OR.
     &     ( KEEP(50) .NE. 0 .AND.
     &       IW( PTRIST(STEP(ISON)) + 6 + KEEP(IXSZ) ) .NE. 0 ) )
          IF ( KEEP(50).eq.0) THEN
#if defined(IBC_TEST)
            MSGSOU = IW( PTRIST(STEP(ISON)) + 7 +  KEEP(IXSZ) )
            MSGTAG = BLOC_FACTO
#else
            MSGSOU = PDEST_MASTER_ISON
            MSGTAG = BLOC_FACTO
#endif
          ELSE
            IF ( IW( PTRIST(STEP(ISON)) + 1 + KEEP(IXSZ) ) .NE.
     &           IW( PTRIST(STEP(ISON)) + 3 + KEEP(IXSZ) ) ) THEN
#if defined(IBC_TEST)
              MSGSOU = IW( PTRIST(STEP(ISON)) + 9 +  KEEP(IXSZ) )
              MSGTAG = BLOC_FACTO_SYM
#else
              MSGSOU = PDEST_MASTER_ISON
              MSGTAG = BLOC_FACTO_SYM
#endif
            ELSE
              MSGSOU = MPI_ANY_SOURCE
              MSGTAG = BLOC_FACTO_SYM_SLAVE
            END IF
          END IF
          BLOCKING = .TRUE.
          SET_IRECV= .FALSE.
          MESSAGE_RECEIVED = .FALSE.
          CALL ZMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MSGSOU, MSGTAG,
     &    STATUS, 
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE, LPTRAR,
     &    NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &    )
          IF ( IFLAG .LT. 0 ) THEN
            IS_ERROR_BROADCASTED = .TRUE.  
            GOTO 670                       
          ENDIF
        END DO
      ENDIF
#if ! defined(NO_FDM_MAPROW)
 10   CONTINUE
#endif
      IF ( NSLAVES_PERE .EQ. 0 ) THEN
        NBROW( 0 ) = LMAP_LOC
      ELSE
        DO I = 0, NSLAVES_PERE
          NBROW( I ) = 0
        END DO
        DO I = 1, LMAP_LOC
          INDICE_PERE = MAP( I )
          CALL MUMPS_BLOC2_GET_ISLAVE(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &           NASS_PERE,
     &           NFRONT_PERE - NASS_PERE,
     &           NSLAVES_PERE,
     &           INDICE_PERE,
     &           NOSLA,
     &           IPOS_IN_SLAVE )
          NBROW( NOSLA ) = NBROW( NOSLA ) + 1
        END DO
        DO I = 1, NSLAVES_PERE
          NBROW(I)=NBROW(I)+NBROW(I-1)
        ENDDO
      ENDIF
      ALLOCATE(PERM(LMAP_LOC), stat=allocok)
      IF (allocok .GT. 0) THEN
          IF (LP.GT.0) THEN
          write(LP,*) MYID,': PB allocation PERM in ZMUMPS_MAPLIG'
          ENDIF
          IFLAG  =-13
          IERROR = LMAP_LOC
          GOTO 670
      ENDIF
      KEEP253_LOC   = 0
      DO I = LMAP_LOC, 1, -1
          INDICE_PERE = MAP( I )
          IF (INDICE_PERE > NFRONT_PERE - KEEP(253)) THEN
             KEEP253_LOC = KEEP253_LOC + 1
          ENDIF
          CALL MUMPS_BLOC2_GET_ISLAVE(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &           NASS_PERE,
     &           NFRONT_PERE - NASS_PERE,
     &           NSLAVES_PERE,
     &           INDICE_PERE,
     &           NOSLA,
     &           IPOS_IN_SLAVE )
          PERM( NBROW( NOSLA ) ) = I
          NBROW( NOSLA ) = NBROW( NOSLA ) - 1
      ENDDO
      DO I = 0, NSLAVES_PERE
          NBROW(I)=NBROW(I)+1
      END DO
      PDEST_MASTER = SLAVES_PERE(0)
      I_POSMYIDIN_PERE = -99999
      LOCAL_ASSEMBLY_TO_BE_DONE = .FALSE.
      DO I = 0, NSLAVES_PERE
        IF (SLAVES_PERE(I) .EQ. MYID) THEN
          I_POSMYIDIN_PERE = I
          LOCAL_ASSEMBLY_TO_BE_DONE = .TRUE.
#if ! defined(NO_FDM_DESCBAND)
          IF (PTRIST(STEP(INODE_PERE)) .EQ. 0
     &      .AND. MYID .NE. PDEST_MASTER) THEN
            CALL ZMUMPS_TREAT_DESCBAND( INODE_PERE, COMM_LOAD,
     &      ASS_IRECV, 
     &      BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &      PTLUST, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &      IFLAG, IERROR, COMM,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF,
     &      NBFIN, MYID, SLAVEF,
     &    
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, PTRARW, PTRAIW,
     &      INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE, LPTRAR,
     &      NELT, FRTPTR, FRTELT, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &      )
            IF ( IFLAG .LT. 0 ) THEN
              IS_ERROR_BROADCASTED = .TRUE. 
              GOTO 600                      
            ENDIF
          ENDIF
#endif
        ENDIF
      END DO
      IF (KEEP(120).NE.0 .AND. LOCAL_ASSEMBLY_TO_BE_DONE) THEN
        CALL ZMUMPS_LOCAL_ASSEMBLY_TYPE2(I_POSMYIDIN_PERE,
     &     SLAVES_PERE(I_POSMYIDIN_PERE),  
     &     MYID, PDEST_MASTER, ISON, INODE_PERE, 
     &     NSLAVES_PERE, NASS_PERE, NFRONT_PERE, NFS4FATHER,
     &     LMAP_LOC, MAP, NBROW, PERM,
     &     IS_ofType5or6, IFLAG, IERROR,
     &     N, SLAVEF, KEEP, NBPROCFILS, IPOOL, LPOOL, STEP,
     &     PROCNODE_STEPS, COMM_LOAD, ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &           KEEP8, IW, LIW, A, LA, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &           PTRIST, PTLUST, PTRAST, PAMASTER, PIMASTER, ND,
     &           NELT, FRTPTR, FRTELT,
     &           OPASSW, OPELIW,
     &           ITLOC, RHS_MUMPS, KEEP253_LOC,
     &           FILS, LPTRAR, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL)
        LOCAL_ASSEMBLY_TO_BE_DONE = .FALSE.
        IF (IFLAG < 0) THEN
          GOTO 600
        ENDIF
      ENDIF
      DO I = NSLAVES_PERE, 0, -1   
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .NE. MYID ) THEN
           DESCLU = .FALSE.
           NBROWS_ALREADY_SENT = 0
           IF (I == NSLAVES_PERE) THEN
             NROWS_TO_SEND=LMAP_LOC-NBROW(I)+1
           ELSE
             NROWS_TO_SEND=NBROW(I+1)-NBROW(I)
           ENDIF
           COMPRESSCB=(IW(PTRIST(STEP(ISON))+XXS).EQ.S_CB1COMP)
           IERR = -1
           DO WHILE (IERR .EQ. -1)
             IF ( IW ( PTRIST(STEP(ISON) )+KEEP(IXSZ) )
     &            .GT. N + KEEP(253) ) THEN
               WRITE(*,*) MYID,': Internal error in Maplig'
               WRITE(*,*) MYID,': PTRIST(STEP(ISON))/N=',
     &                            PTRIST(STEP(ISON)), N
               WRITE(*,*) MYID,': I, NBROW(I)=',I, NBROW(I)
               WRITE(*,*) MYID,': NSLAVES_PERE=',NSLAVES_PERE
               WRITE(*,*) MYID,': ISON, INODE_PERE=',ISON,INODE_PERE
               WRITE(*,*) MYID,': Son header=',
     &         IW(PTRIST(STEP(ISON)): PTRIST(STEP(ISON))+5+KEEP(IXSZ))
               CALL MUMPS_ABORT()
             END IF
             IF (NROWS_TO_SEND .EQ. 0 .AND. PDEST.NE.PDEST_MASTER) THEN
                IERR = 0
                CYCLE
             ENDIF
             CALL ZMUMPS_BUF_SEND_CONTRIB_TYPE2( NBROWS_ALREADY_SENT,
     &       DESCLU, INODE_PERE,
     &       NFRONT_PERE, NASS_PERE, NFS4FATHER,
     &            NSLAVES_PERE, ISON,
     &       NROWS_TO_SEND, LMAP_LOC, MAP,
     &       PERM(min(LMAP_LOC,NBROW(I))),
     &       IW( PTRIST(STEP(ISON))),
     &       A(PTRAST(STEP(ISON))), I, PDEST, PDEST_MASTER, 
     &       COMM, IERR, 
     &
     &       KEEP,KEEP8, STEP, N, SLAVEF,
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE, COMPRESSCB,
     &       KEEP253_LOC )
             IF ( IERR .EQ. -2 ) THEN
               IFLAG  = -17
               IF (LP .GT. 0) THEN
                 WRITE(LP,*)
     &           "FAILURE: SEND BUFFER TOO SMALL IN ZMUMPS_MAPLIG"
               ENDIF
               IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &         NROWS_TO_SEND * IW(PTRIST(STEP(ISON))+KEEP(IXSZ))
     &        * KEEP( 35 )
               GO TO 600
             END IF
             IF ( IERR .EQ. -3 ) THEN
               IF (LP .GT. 0) THEN
                 WRITE(LP,*)
     &           "FAILURE: RECV BUFFER TOO SMALL IN ZMUMPS_MAPLIG"
               ENDIF
               IFLAG  = -20
               IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &         NROWS_TO_SEND * IW(PTRIST(STEP(ISON))+KEEP(IXSZ))
     &         * KEEP( 35 )
               GOTO 600
             ENDIF
             IF (KEEP(219).NE.0) THEN
              IF ( IERR .EQ. -4 ) THEN
                IFLAG  = -13
               IERROR = NFS4FATHER
               IF (LP .GT. 0) THEN
                 WRITE(LP, *)
     & "FAILURE: MAX_ARRAY allocation failed IN ZMUMPS_MAPLIG"
               ENDIF
               GO TO 600
              END IF
             END IF
             IF ( IERR .EQ. -1 ) THEN
               IF (LOCAL_ASSEMBLY_TO_BE_DONE) THEN
                 CALL ZMUMPS_LOCAL_ASSEMBLY_TYPE2(I_POSMYIDIN_PERE,
     &           SLAVES_PERE(I_POSMYIDIN_PERE),  
     &           MYID, PDEST_MASTER, ISON, INODE_PERE, 
     &           NSLAVES_PERE, NASS_PERE, NFRONT_PERE, NFS4FATHER,
     &           LMAP_LOC, MAP, NBROW, PERM,
     &           IS_ofType5or6, IFLAG, IERROR,
     &           N, SLAVEF, KEEP, NBPROCFILS, IPOOL, LPOOL, STEP,
     &           PROCNODE_STEPS, COMM_LOAD, ISTEP_TO_INIV2,
     &           TAB_POS_IN_PERE,
     &           KEEP8, IW, LIW, A, LA, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &           PTRIST, PTLUST, PTRAST, PAMASTER, PIMASTER, ND,
     &           NELT, FRTPTR, FRTELT,
     &           OPASSW, OPELIW,
     &           ITLOC, RHS_MUMPS, KEEP253_LOC,
     &           FILS, LPTRAR, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL)
                 LOCAL_ASSEMBLY_TO_BE_DONE = .FALSE.
                 IF (IFLAG < 0) THEN
                   GOTO 600
                 ENDIF
               ELSE
                 BLOCKING = .FALSE.
                 SET_IRECV = .TRUE.
                 MESSAGE_RECEIVED = .FALSE.
                 CALL ZMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &           ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &           MPI_ANY_SOURCE, MPI_ANY_TAG,
     &           STATUS,
     &           BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &           IWPOS, IWPOSCB, IPTRLU,
     &           LRLU, LRLUS, N, IW, LIW, A, LA,
     &           PTRIST, PTLUST, PTRFAC,
     &           PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &           IFLAG, IERROR, COMM,
     &           NBPROCFILS,
     &           IPOOL, LPOOL, LEAF,
     &           NBFIN, MYID, SLAVEF,
     &
     &           root, OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, 
     &           PTRARW, PTRAIW,
     &           INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE,LPTRAR,
     &           NELT, FRTPTR, FRTELT, 
     &           ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &           )
                 IF ( IFLAG .LT. 0 ) THEN
                   IS_ERROR_BROADCASTED=.TRUE.
                   GOTO 600
                 ENDIF
               END IF
             END IF 
           ENDDO
        ENDIF 
      END DO
      IF (LOCAL_ASSEMBLY_TO_BE_DONE) THEN
        CALL ZMUMPS_LOCAL_ASSEMBLY_TYPE2(I_POSMYIDIN_PERE,
     &     SLAVES_PERE(I_POSMYIDIN_PERE),  
     &     MYID, PDEST_MASTER, ISON, INODE_PERE, 
     &     NSLAVES_PERE, NASS_PERE, NFRONT_PERE, NFS4FATHER,
     &     LMAP_LOC, MAP, NBROW, PERM,
     &     IS_ofType5or6, IFLAG, IERROR,
     &     N, SLAVEF, KEEP, NBPROCFILS, IPOOL, LPOOL, STEP,
     &     PROCNODE_STEPS, COMM_LOAD, ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &           KEEP8, IW, LIW, A, LA, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &           PTRIST, PTLUST, PTRAST, PAMASTER, PIMASTER, ND,
     &           NELT, FRTPTR, FRTELT,
     &           OPASSW, OPELIW,
     &           ITLOC, RHS_MUMPS, KEEP253_LOC,
     &           FILS, LPTRAR, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL)
        LOCAL_ASSEMBLY_TO_BE_DONE = .FALSE.
        IF (IFLAG < 0) THEN
          GOTO 600
        ENDIF
      ENDIF
      ITYPE = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(ISON)), SLAVEF)
      IF (KEEP(214) .EQ. 2) THEN
        CALL ZMUMPS_STACK_BAND( N, ISON,
     &    PTRIST, PTRAST, PTLUST, PTRFAC, IW, LIW, A, LA,
     &    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP,
     &    IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER,
     &    IFLAG, IERROR, SLAVEF, MYID, COMM, KEEP,KEEP8, DKEEP,ITYPE
     &     )
        IF (IFLAG .LT. 0) THEN
          IS_ERROR_BROADCASTED = .TRUE.
          GOTO 600
        ENDIF
      ENDIF
      CALL ZMUMPS_FREE_BAND( N, ISON, PTRIST, PTRAST, IW, LIW,
     &             A, LA, LRLU, LRLUS, IWPOSCB, IPTRLU,
     &             STEP, MYID, KEEP, ITYPE
     &)
 600  CONTINUE
      DEALLOCATE(PERM)
 670  CONTINUE
      DEALLOCATE(MAP)
 680  CONTINUE
      DEALLOCATE(NBROW)
      DEALLOCATE(SLAVES_PERE)
 700  CONTINUE
      IF (IFLAG .LT. 0 .AND. .NOT. IS_ERROR_BROADCASTED) THEN
        CALL ZMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
      ENDIF
      RETURN
      END SUBROUTINE ZMUMPS_MAPLIG
      SUBROUTINE ZMUMPS_MAPLIG_FILS_NIV1( COMM_LOAD, ASS_IRECV, 
     &  BUFR, LBUFR, LBUFR_BYTES,
     &
     &  INODE_PERE, ISON, NSLAVES_PERE, LIST_SLAVES_PERE,
     &  NFRONT_PERE, NASS_PERE, NFS4FATHER, LMAP, TROW,
     &  PROCNODE_STEPS, SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &  LRLUS, N, IW,
     &  LIW, A, LA,
     &  PTRIST, PTLUST, PTRFAC,
     &  PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &  IFLAG, IERROR, MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     &  NBFIN, ICNTL, KEEP,KEEP8,DKEEP, root,
     &  OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &  FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &  ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT, 
     &
     &  ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &  )
      USE ZMUMPS_COMM_BUFFER
      USE ZMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'zmumps_root.h'
      TYPE (ZMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER ICNTL( 40 ), KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION    DKEEP(130)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER SLAVEF, NBFIN
      INTEGER(8) :: LA, IPTRLU, LRLU, LRLUS, POSFAC
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      COMPLEX(kind=8) A( LA )
      INTEGER COMP
      INTEGER IFLAG, IERROR, COMM, MYID
      INTEGER LPOOL, LEAF
      INTEGER INODE_PERE, ISON
      INTEGER NFS4FATHER
      INTEGER NSLAVES_PERE, NFRONT_PERE, NASS_PERE
      INTEGER LIST_SLAVES_PERE(NSLAVES_PERE)
      INTEGER NELIM, LMAP, TROW( LMAP )
      DOUBLE PRECISION OPASSW, OPELIW
      COMPLEX(kind=8) DBLARR(max(1,KEEP(13)))
      INTEGER LPTRAR, NELT
      INTEGER IW( LIW )
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL )
      INTEGER NSTK( KEEP(28) ), ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) )
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28))
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER INTARR(max(1,KEEP(14)))
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER ITLOC( N+KEEP(253) ), FILS( N )
      COMPLEX(kind=8) :: RHS_MUMPS(KEEP(255))
      INTEGER PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER LP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER NOSLA, I, ISTCHK, ISTCHK_LOC
      INTEGER NBROWS_ALREADY_SENT 
      INTEGER INDICE_PERE
      INTEGER INDICE_PERE_ARRAY_ARG(1)
      INTEGER PDEST, PDEST_MASTER, NFRONT
      LOGICAL SAME_PROC, DESCLU
      INTEGER(8) :: APOS, POSROW, ASIZE
      INTEGER NSLSON, NBCOLS, NROW, NROWS_TO_SEND,
     &        NPIV, NROWS_TO_STACK, II, IROW_SON,
     &        IPOS_IN_SLAVE, DECR
      INTEGER NBCOLS_EFF
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      LOGICAL COMPRESSCB
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INTEGER LMAP_LOC, allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NBROW
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SLAVES_PERE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAP, PERM
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      if (NSLAVES_PERE.le.0) then
       write(6,*) ' error 2 in maplig_fils_niv1 ', NSLAVES_PERE
       CALL MUMPS_ABORT()
      endif
      ALLOCATE(NBROW(0:NSLAVES_PERE), stat=allocok)
      IF (allocok .GT. 0) THEN
        IF (LP > 0)
     &  write(LP,*) MYID,
     &  ' : PB allocation NBROW in ZMUMPS_MAPLIG_FILS_NIV1'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      ENDIF
      ALLOCATE(SLAVES_PERE(0:NSLAVES_PERE), stat =allocok)
      IF ( allocok .GT. 0 ) THEN
        IF (LP > 0) write(LP,*) MYID,
     &  ' : PB allocation SLAVES_PERE in ZMUMPS_MAPLIG_FILS_NIV1'
        IFLAG  =-13
        IERROR = NSLAVES_PERE+1
        GOTO 700
      ENDIF
      SLAVES_PERE(1:NSLAVES_PERE) = LIST_SLAVES_PERE(1:NSLAVES_PERE)
      SLAVES_PERE(0) = MUMPS_PROCNODE( 
     &                       PROCNODE_STEPS(STEP(INODE_PERE)),
     &                       SLAVEF )
      LMAP_LOC = LMAP
      ALLOCATE(MAP(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP > 0) write(LP,*) MYID,
     &   ' : PB allocation LMAP in ZMUMPS_MAPLIG_FILS_NIV1'
        IFLAG  =-13
        IERROR = LMAP_LOC
        GOTO 700
      endif
      MAP( 1 : LMAP_LOC ) = TROW( 1 : LMAP_LOC )
      DO I = 0, NSLAVES_PERE
        NBROW( I ) = 0
      END DO
      IF (NSLAVES_PERE == 0) THEN
        NBROW(0) = LMAP_LOC
      ELSE
       DO I = 1, LMAP_LOC
        INDICE_PERE = MAP( I )
        CALL MUMPS_BLOC2_GET_ISLAVE(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &         NASS_PERE,
     &         NFRONT_PERE - NASS_PERE,
     &         NSLAVES_PERE,
     &         INDICE_PERE,
     &         NOSLA,
     &         IPOS_IN_SLAVE )
        NBROW( NOSLA ) = NBROW( NOSLA ) + 1
       END DO
        DO I = 1, NSLAVES_PERE
          NBROW(I)=NBROW(I)+NBROW(I-1)
        ENDDO
      ENDIF
      ALLOCATE(PERM(LMAP_LOC), stat=allocok)
      if (allocok .GT. 0) THEN
        IF (LP > 0) write(LP,*) MYID,
     &  ': PB allocation PERM in ZMUMPS_MAPLIG_FILS_NIV1'
        IFLAG  =-13
        IERROR = LMAP_LOC
        GOTO 700
      endif
        ISTCHK     = PIMASTER(STEP(ISON))
        NBCOLS     = IW(ISTCHK+KEEP(IXSZ))
      DO I = LMAP_LOC, 1, -1
          INDICE_PERE = MAP( I )
          CALL MUMPS_BLOC2_GET_ISLAVE(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &           NASS_PERE,
     &           NFRONT_PERE - NASS_PERE,
     &           NSLAVES_PERE,
     &           INDICE_PERE,
     &           NOSLA,
     &           IPOS_IN_SLAVE )
          PERM( NBROW( NOSLA ) ) = I
          NBROW( NOSLA ) = NBROW( NOSLA ) - 1
      ENDDO
      DO I = 0, NSLAVES_PERE
          NBROW(I)=NBROW(I)+1
      END DO
      PDEST_MASTER = MYID
      IF ( SLAVES_PERE(0) .NE. MYID ) THEN
        WRITE(*,*) 'Error 1 in MAPLIG_FILS_NIV1:',MYID, SLAVES_PERE
        CALL MUMPS_ABORT()
      END IF
      PDEST        = PDEST_MASTER
        I = 0
        ISTCHK     = PIMASTER(STEP(ISON))
        NBCOLS     = IW(ISTCHK+KEEP(IXSZ))
        NELIM      = IW(ISTCHK+1+KEEP(IXSZ))
        NROW       = IW(ISTCHK+2+KEEP(IXSZ))
        NPIV       = IW(ISTCHK+3+KEEP(IXSZ))
        IF (NPIV.LT.0) THEN
         write(6,*) ' Error 2 in ZMUMPS_MAPLIG_FILS_NIV1 ', NPIV
         CALL MUMPS_ABORT()
        ENDIF
        NSLSON     = IW(ISTCHK+5+KEEP(IXSZ))
        NFRONT     = NPIV + NBCOLS
        COMPRESSCB=(IW(PTRIST(STEP(ISON))+XXS) .eq. S_CB1COMP)
        IF (I == NSLAVES_PERE) THEN
          NROWS_TO_STACK=LMAP_LOC-NBROW(I)+1
        ELSE
          NROWS_TO_STACK=NBROW(I+1)-NBROW(I)
        ENDIF
        DECR=1
        NBPROCFILS(STEP(INODE_PERE)) =
     &                           NBPROCFILS(STEP(INODE_PERE)) - DECR
        NBPROCFILS(STEP(ISON))       = NBPROCFILS(STEP(ISON)) - DECR
#if ! defined(NO_XXNBPR)
        IW(PTLUST(STEP(INODE_PERE))+XXNBPR) =
     &  IW(PTLUST(STEP(INODE_PERE))+XXNBPR) - DECR
          CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE_PERE)),
     &                     IW(PTLUST(STEP(INODE_PERE))+XXNBPR))
        IW(PTRIST(STEP(ISON))+XXNBPR) =
     &  IW(PTRIST(STEP(ISON))+XXNBPR) - DECR
          CALL CHECK_EQUAL(NBPROCFILS(STEP(ISON)),
     &                     IW(PTRIST(STEP(ISON))+XXNBPR))
#endif
        DO II = 1,NROWS_TO_STACK
          IROW_SON=PERM(NBROW(I)+II-1)
          INDICE_PERE = MAP(IROW_SON)
          CALL MUMPS_BLOC2_GET_ISLAVE(
     &         KEEP,KEEP8, INODE_PERE, STEP, N, SLAVEF,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &         NASS_PERE,
     &         NFRONT_PERE - NASS_PERE,
     &         NSLAVES_PERE,
     &         INDICE_PERE,
     &         NOSLA,
     &         IPOS_IN_SLAVE )
          INDICE_PERE = IPOS_IN_SLAVE
          IF (COMPRESSCB) THEN
            IF (NELIM.EQ.0) THEN
            POSROW = PAMASTER(STEP(ISON)) +
     &         int(IROW_SON,8)*int(IROW_SON-1,8)/2_8
            ELSE
            POSROW = PAMASTER(STEP(ISON)) +
     &         int(NELIM+IROW_SON,8)*int(NELIM+IROW_SON-1,8)/2_8
            ENDIF
          ELSE
            POSROW = PAMASTER(STEP(ISON)) +
     &             int(NELIM+IROW_SON-1,8)*int(NBCOLS,8)
          ENDIF
          IF (KEEP(50).NE.0) THEN
            NBCOLS_EFF = NELIM + IROW_SON
          ELSE
            NBCOLS_EFF = NBCOLS
          ENDIF
          INDICE_PERE_ARRAY_ARG(1) = INDICE_PERE
          CALL ZMUMPS_ASM_SLAVE_MASTER(N, INODE_PERE, IW, LIW, 
     &    A, LA, ISON, 1, NBCOLS_EFF, INDICE_PERE_ARRAY_ARG,
     &    A(POSROW), PTLUST, PTRAST,
     &    STEP, PIMASTER, OPASSW, IWPOSCB, 
     &    MYID, KEEP,KEEP8,.FALSE.,NBCOLS_EFF)
        ENDDO
        IF (KEEP(219).NE.0) THEN
         IF(NSLAVES_PERE.GT.0 .AND. KEEP(50).EQ.2) THEN
           IF (COMPRESSCB) THEN
             POSROW = PAMASTER(STEP(ISON))
     &          + int(NELIM+NBROW(1),8)*int(NELIM+NBROW(1)-1,8)/2_8
             ASIZE  = int(LMAP_LOC+NELIM,8)*int(NELIM+LMAP_LOC+1,8)/2_8
     &          - int(NELIM+NBROW(1),8)*int(NELIM+NBROW(1)-1,8)/2_8
           ELSE
             POSROW = PAMASTER(STEP(ISON)) +
     &                 int(NELIM+NBROW(1)-1,8)*int(NBCOLS,8)
             ASIZE  = int(LMAP_LOC-NBROW(1)+1,8) * int(NBCOLS,8)
           ENDIF
           CALL ZMUMPS_BUF_MAX_ARRAY_MINSIZE(NFS4FATHER,IERR)
           IF (IERR .NE.0) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     &    ": PB allocation MAX_ARRAY during ZMUMPS_MAPLIG_FILS_NIV1"
              IFLAG=-13
              IERROR=NFS4FATHER
              GOTO 700
           ENDIF
           IF  ( LMAP_LOC-NBROW(1)+1-KEEP(253).GT. 0 ) THEN
           CALL ZMUMPS_COMPUTE_MAXPERCOL(
     &          A(POSROW),ASIZE,NBCOLS,
     &          LMAP_LOC-NBROW(1)+1-KEEP(253),
     &          BUF_MAX_ARRAY,NFS4FATHER,COMPRESSCB,
     &          NELIM+NBROW(1))
           ELSE
                CALL ZMUMPS_SETMAXTOZERO(BUF_MAX_ARRAY,
     &          NFS4FATHER)
           ENDIF
           CALL ZMUMPS_ASM_MAX(N, INODE_PERE, IW, LIW, 
     &          A, LA, ISON, NFS4FATHER,
     &          BUF_MAX_ARRAY, PTLUST, PTRAST,
     &          STEP, PIMASTER, OPASSW,
     &          IWPOSCB,MYID, KEEP,KEEP8)
         ENDIF
        ENDIF 
#if ! defined(NO_XXNBPR)
          CALL CHECK_EQUAL(NBPROCFILS(STEP(ISON)),
     &                     IW(PTRIST(STEP(ISON))+XXNBPR))
          IF (IW(PTRIST(STEP(ISON))+XXNBPR) .EQ. 0
#else
          IF ( NBPROCFILS(STEP(ISON)) .EQ. 0
#endif
     &       ) THEN
               ISTCHK_LOC = PIMASTER(STEP(ISON))
               SAME_PROC= ISTCHK_LOC .LT. IWPOSCB
               IF (SAME_PROC) THEN
                 CALL ZMUMPS_RESTORE_INDICES(N, ISON, INODE_PERE,
     &            IWPOSCB, PIMASTER, PTLUST, IW, LIW, STEP,
     &            KEEP,KEEP8)
               ENDIF
          ENDIF
#if ! defined(NO_XXNBPR)
          CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE_PERE)),
     &                     IW(PTLUST(STEP(INODE_PERE))+XXNBPR))
          IF ( IW(PTLUST(STEP(INODE_PERE))+XXNBPR) .EQ. 0
#else
          IF ( NBPROCFILS(STEP(INODE_PERE)) .EQ. 0
#endif
     &       ) THEN
            CALL ZMUMPS_INSERT_POOL_N( N, IPOOL, LPOOL,
     &        PROCNODE_STEPS,
     &        SLAVEF, KEEP(28), KEEP(76), KEEP(80),
     &        KEEP(47), STEP, INODE_PERE+N )
            IF (KEEP(47) .GE. 3) THEN
              CALL ZMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &       IPOOL, LPOOL, 
     &       PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &       MYID, STEP, N, ND, FILS )
            ENDIF
          END IF
      DO I = 0, NSLAVES_PERE
        PDEST = SLAVES_PERE( I )
        IF ( PDEST .NE. MYID ) THEN
           NBROWS_ALREADY_SENT = 0
 95        CONTINUE
           NFRONT = IW(PIMASTER(STEP(ISON))+KEEP(IXSZ))
           NELIM  = IW(PIMASTER(STEP(ISON))+1+KEEP(IXSZ))
           APOS = PAMASTER(STEP(ISON))
           DESCLU = .TRUE.
           IF (I == NSLAVES_PERE) THEN
             NROWS_TO_SEND=LMAP_LOC-NBROW(I)+1
           ELSE
             NROWS_TO_SEND=NBROW(I+1)-NBROW(I)
           ENDIF
           IF ( NROWS_TO_SEND .EQ. 0) CYCLE
           CALL ZMUMPS_BUF_SEND_CONTRIB_TYPE2(NBROWS_ALREADY_SENT,
     &      DESCLU, INODE_PERE,
     &      NFRONT_PERE, NASS_PERE, NFS4FATHER, 
     &           NSLAVES_PERE,
     &      ISON, NROWS_TO_SEND, LMAP_LOC,
     &      MAP, PERM(min(LMAP_LOC,NBROW(I))),
     &      IW(PIMASTER(STEP(ISON))),
     &      A(APOS), I, PDEST, PDEST_MASTER, COMM, IERR,
     &
     &      KEEP,KEEP8, STEP, N, SLAVEF,
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &      COMPRESSCB, KEEP(253))
            IF ( IERR .EQ. -2 ) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     &": FAILURE, SEND BUFFER TOO SMALL DURING ZMUMPS_MAPLIG_FILS_NIV1"
              IFLAG  = -17
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &        NROWS_TO_SEND *  KEEP( 35 )
              GO TO 700
            END IF
            IF ( IERR .EQ. -3 ) THEN
              IF (LP > 0) WRITE(LP,*) MYID,
     &": FAILURE, RECV BUFFER TOO SMALL DURING ZMUMPS_MAPLIG_FILS_NIV1"
              IFLAG  = -20
              IERROR =  (NROWS_TO_SEND + 3 )* KEEP( 34 ) +
     &        NROWS_TO_SEND *  KEEP( 35 )
              GO TO 700
            ENDIF
            IF (KEEP(219).NE.0) THEN
             IF ( IERR .EQ. -4 ) THEN
               IFLAG  = -13
               IERROR = BUF_LMAX_ARRAY
              IF (LP > 0) WRITE(LP,*) MYID,
     &": FAILURE, MAX_ARRAY ALLOC FAILED DURING ZMUMPS_MAPLIG_FILS_NIV1"
               GO TO 700
             ENDIF
            ENDIF
            IF ( IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV = .TRUE.
              MESSAGE_RECEIVED = .FALSE.
              CALL ZMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &          ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &          MPI_ANY_SOURCE, MPI_ANY_TAG,
     &          STATUS, 
     &          BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &          IWPOS, IWPOSCB, IPTRLU,
     &          LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &          PTLUST, PTRFAC,
     &          PTRAST, STEP, PIMASTER, PAMASTER, NSTK, COMP,
     &          IFLAG, IERROR, COMM,
     &          NBPROCFILS,
     &          IPOOL, LPOOL, LEAF,
     &          NBFIN, MYID, SLAVEF,
     &          root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &          FILS, PTRARW, PTRAIW,
     &          INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE,
     &          LPTRAR, NELT, FRTPTR, FRTELT, 
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &          )
              IF ( IFLAG .LT. 0 ) GOTO 600
              GO TO 95
            END IF
        END IF
      END DO
      ISTCHK = PTRIST(STEP(ISON))
      PTRIST(STEP( ISON )) = -77777777
            IF ( IW(ISTCHK+KEEP(IXSZ)) .GE. 0 ) THEN
              WRITE(*,*) 'error 3 in ZMUMPS_MAPLIG_FILS_NIV1'
              CALL MUMPS_ABORT()
            ENDIF
      CALL ZMUMPS_FREE_BLOCK_CB(.FALSE., MYID, N, ISTCHK,
     &     PAMASTER(STEP(ISON)),
     &     IW, LIW, LRLU, LRLUS, IPTRLU,
     &     IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &     )
 600  CONTINUE
      DEALLOCATE(NBROW)
      DEALLOCATE(MAP)
      DEALLOCATE(PERM)
      DEALLOCATE(SLAVES_PERE)
      RETURN
 700  CONTINUE
      CALL ZMUMPS_BDC_ERROR(MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE ZMUMPS_MAPLIG_FILS_NIV1
      SUBROUTINE ZMUMPS_LOCAL_ASSEMBLY_TYPE2(I, PDEST, MYID,
     &           PDEST_MASTER, ISON, IFATH, NSLAVES_PERE, NASS_PERE,
     &           NFRONT_PERE, NFS4FATHER, LMAP_LOC, MAP,
     &           NBROW, PERM, IS_ofType5or6, IFLAG, IERROR,
     &           N, SLAVEF, KEEP, NBPROCFILS,
     &           IPOOL, LPOOL, STEP,
     &           PROCNODE_STEPS, COMM_LOAD, ISTEP_TO_INIV2,
     &           TAB_POS_IN_PERE,
     &
     &           KEEP8, IW, LIW, A, LA, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &           PTRIST, PTLUST, PTRAST, PAMASTER, PIMASTER, ND,
     &           NELT, FRTPTR, FRTELT,
     &           OPASSW, OPELIW,
     &           ITLOC, RHS_MUMPS, KEEP253_LOC,
     &           FILS, LPTRAR, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL
     &           )
      USE ZMUMPS_COMM_BUFFER, ONLY: ZMUMPS_BUF_MAX_ARRAY_MINSIZE,
     &                              BUF_MAX_ARRAY
      USE ZMUMPS_LOAD, ONLY : ZMUMPS_LOAD_POOL_UPD_NEW_POOL
      INTEGER ICNTL(40)
      INTEGER, intent(in) :: I, PDEST, MYID, PDEST_MASTER, IFATH, ISON
      INTEGER, intent(in) :: N, SLAVEF
      INTEGER, intent(in) :: NSLAVES_PERE, NASS_PERE, NFRONT_PERE
      INTEGER, intent(in) :: NFS4FATHER
      INTEGER, intent(in) :: KEEP(500), STEP(N)
      INTEGER, intent(inout) :: NBPROCFILS( KEEP(28) )
      INTEGER, intent(in) :: LMAP_LOC
      INTEGER, intent(in) :: NBROW(0:NSLAVES_PERE)
      INTEGER, intent(in) :: MAP(LMAP_LOC), PERM(LMAP_LOC)
      INTEGER, intent(inout) :: IFLAG, IERROR
      INTEGER(8), intent(in) :: KEEP8(150)
      INTEGER, intent(in) :: LIW, NELT, LPTRAR
      INTEGER(8), intent(in) :: LA
      INTEGER(8), intent(inout) :: IPTRLU, LRLU, LRLUS
      INTEGER, intent(inout) :: IWPOSCB
      INTEGER, intent(inout) :: IW(LIW)
      COMPLEX(kind=8), intent(inout) :: A( LA )
      INTEGER(8) :: PTRAST(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER    :: PTRIST(KEEP(28)), PIMASTER(KEEP(28)), ND(KEEP(28))
      INTEGER    :: PTLUST(KEEP(28))
      INTEGER, intent(inout) :: ITLOC(N)
      INTEGER, intent(in) :: FRTPTR( N+1 ), FRTELT( NELT )
      DOUBLE PRECISION, intent(inout) :: OPASSW, OPELIW
      COMPLEX(kind=8) :: RHS_MUMPS(KEEP(255))
      INTEGER, intent(in) :: KEEP253_LOC
      INTEGER, intent(in) :: FILS(N)
      INTEGER, intent(in) :: PTRARW( LPTRAR ), PTRAIW( LPTRAR )
      INTEGER, intent(in) :: PROCNODE_STEPS( KEEP(28) ), COMM_LOAD
      INTEGER ISTEP_TO_INIV2(KEEP(71)),
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      COMPLEX(kind=8) DBLARR(max(1,KEEP(13)))
      INTEGER INTARR(max(1,KEEP(14)))
      INTEGER LPOOL
      INTEGER IPOOL( LPOOL )
      LOGICAL, intent(in) :: IS_ofType5or6
      INCLUDE 'mumps_headers.h'
      INCLUDE 'mpif.h'
      INTEGER    :: ISTCHK, ISTCHK_LOC, NBCOLS, NROW, NPIV, NSLSON, 
     &              NFRONT, LDA_SON, NROWS_TO_STACK, II, INDICE_PERE,
     &              NOSLA, COLLIST, IPOS_IN_SLAVE, IROW_SON, ITMP,
     &              NBCOLS_EFF, DECR
      LOGICAL    :: COMPRESSCB, SAME_PROC
      INTEGER(8) :: SIZFR, POSROW, SHIFTCB_SON
      INTEGER    :: IERR, LP
      INTEGER INDICE_PERE_ARRAY_ARG(1)
#if ! defined(NO_XXNBPR)
      INTEGER :: INBPROCFILS_SON
#endif
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
            IF (I == NSLAVES_PERE) THEN
              NROWS_TO_STACK = LMAP_LOC - NBROW(I) + 1
            ELSE
              NROWS_TO_STACK = NBROW(I+1) - NBROW(I)
            ENDIF
            DECR = 1
            IF ( MYID .EQ. PDEST_MASTER ) THEN
              NBPROCFILS(STEP(IFATH)) =
     &            NBPROCFILS(STEP(IFATH)) - DECR
#if ! defined(NO_XXNBPR)
              IW(PTLUST(STEP(IFATH))+XXNBPR) =
     &            IW(PTLUST(STEP(IFATH))+XXNBPR) - DECR
#endif
              IF ( PDEST .EQ. PDEST_MASTER .AND. DECR .NE. 0) THEN
                NBPROCFILS(STEP(ISON)) = NBPROCFILS(STEP(ISON)) - DECR
#if ! defined(NO_XXNBPR)
                IW(PIMASTER(STEP(ISON))+XXNBPR) =
     &             IW(PIMASTER(STEP(ISON))+XXNBPR) - DECR
#endif
              ENDIF
            ENDIF
            ISTCHK     = PTRIST(STEP(ISON))
            NBCOLS     = IW(ISTCHK+KEEP(IXSZ))
            NROW       = IW(ISTCHK+2+KEEP(IXSZ))
            NPIV       = IW(ISTCHK+3+KEEP(IXSZ))
            NSLSON     = IW(ISTCHK+5+KEEP(IXSZ))
            NFRONT     = NPIV + NBCOLS
            COMPRESSCB = (IW(ISTCHK+XXS).EQ.S_CB1COMP)
            CALL MUMPS_GETI8(SIZFR, IW(ISTCHK+XXR))
            IF (IW(ISTCHK+XXS).EQ.S_NOLCBCONTIG) THEN
               LDA_SON     = NBCOLS
               SHIFTCB_SON = int(NPIV,8)*int(NROW,8)
            ELSE IF (IW(ISTCHK+XXS).EQ.S_NOLCLEANED) THEN
               LDA_SON     = NBCOLS
               SHIFTCB_SON = 0_8
            ELSE
               LDA_SON     = NFRONT
               SHIFTCB_SON = int(NPIV,8)
            ENDIF
            IF (PDEST .NE. PDEST_MASTER) THEN
                IF ( KEEP(55) .eq. 0 ) THEN
                  CALL ZMUMPS_ASM_SLAVE_TO_SLAVE_INIT
     &            (N, IFATH, IW, LIW,
     &            A, LA, NROWS_TO_STACK, NBCOLS,
     &            OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &            ITLOC, RHS_MUMPS,
     &            FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &            KEEP,KEEP8, MYID )
                ELSE
                  CALL ZMUMPS_ELT_ASM_S_2_S_INIT(NELT, FRTPTR, FRTELT,
     &            N, IFATH, IW, LIW,
     &            A, LA, NROWS_TO_STACK, NBCOLS, 
     &            OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &            ITLOC, RHS_MUMPS,
     &            FILS, PTRARW, PTRAIW, INTARR, DBLARR, ICNTL,
     &            KEEP, KEEP8, MYID )
                ENDIF
            ENDIF
            DO II = 1,NROWS_TO_STACK
              IROW_SON = PERM(NBROW(I)+II-1)
              INDICE_PERE=MAP(IROW_SON)
              CALL MUMPS_BLOC2_GET_ISLAVE(
     &        KEEP,KEEP8, IFATH, STEP, N, SLAVEF,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &        NASS_PERE,
     &        NFRONT_PERE - NASS_PERE,
     &        NSLAVES_PERE,
     &        INDICE_PERE,
     &        NOSLA,
     &        IPOS_IN_SLAVE )
              INDICE_PERE = IPOS_IN_SLAVE
              IF ( COMPRESSCB ) THEN
                IF (NBCOLS - NROW .EQ. 0 ) THEN
                  ITMP = IROW_SON 
                  POSROW = PTRAST(STEP(ISON))+
     &                     int(ITMP,8) * int(ITMP-1,8) / 2_8
                ELSE
                  ITMP = IROW_SON + NBCOLS - NROW
                  POSROW = PTRAST(STEP(ISON))
     &               + int(ITMP,8) * int(ITMP-1,8) / 2_8
     &               - int(NBCOLS-NROW,8) * int(NBCOLS-NROW+1,8)/2_8
                ENDIF
              ELSE 
                POSROW = PTRAST(STEP(ISON)) + SHIFTCB_SON
     &               +int(IROW_SON-1,8)*int(LDA_SON,8)
              ENDIF
              IF (PDEST == PDEST_MASTER) THEN
                 IF (KEEP(50).NE.0) THEN
                   NBCOLS_EFF = IROW_SON + NBCOLS - NROW
                 ELSE
                   NBCOLS_EFF = NBCOLS
                 ENDIF
                 INDICE_PERE_ARRAY_ARG(1) = INDICE_PERE
                 IF ((IS_ofType5or6).AND.(KEEP(50).EQ.0)) THEN
                   CALL ZMUMPS_ASM_SLAVE_MASTER(N, IFATH, IW, LIW, 
     &             A, LA, ISON, NROWS_TO_STACK, NBCOLS_EFF, 
     &             INDICE_PERE_ARRAY_ARG,
     &             A(POSROW), PTLUST, PTRAST,
     &             STEP, PIMASTER, OPASSW,
     &             IWPOSCB, MYID, KEEP,KEEP8,
     &             IS_ofType5or6, LDA_SON
     &             )
                   EXIT  
                 ELSE IF ( (KEEP(50).NE.0) .AND. 
     &              (.NOT.COMPRESSCB).AND.(IS_ofType5or6) ) THEN
                   CALL ZMUMPS_ASM_SLAVE_MASTER(N, IFATH, IW, LIW, 
     &             A, LA, ISON, NROWS_TO_STACK,
     &             NBCOLS_EFF, INDICE_PERE_ARRAY_ARG,
     &             A(POSROW), PTLUST, PTRAST,
     &             STEP, PIMASTER, OPASSW,
     &             IWPOSCB, MYID, KEEP,KEEP8,
     &             IS_ofType5or6, LDA_SON
     &)
                   EXIT
                 ELSE
                   CALL ZMUMPS_ASM_SLAVE_MASTER(N, IFATH, IW, LIW, 
     &             A, LA, ISON, 1, NBCOLS_EFF, INDICE_PERE_ARRAY_ARG,
     &             A(POSROW), PTLUST, PTRAST,
     &             STEP, PIMASTER, OPASSW,
     &             IWPOSCB, MYID, KEEP,KEEP8,
     &             IS_ofType5or6, LDA_SON
     &)
                 ENDIF
              ELSE
                 ISTCHK  = PTRIST(STEP(ISON))
                 COLLIST = ISTCHK + 6 + KEEP(IXSZ) 
     &                   + IW( ISTCHK + 5 +KEEP(IXSZ)) + NROW + NPIV
                 IF (KEEP(50).NE.0) THEN
                   NBCOLS_EFF = IROW_SON + NBCOLS - NROW
                 ELSE
                   NBCOLS_EFF = NBCOLS
                 ENDIF
                 INDICE_PERE_ARRAY_ARG(1) = INDICE_PERE
                 IF ( (IS_ofType5or6) .AND.
     &                 (
     &                  ( KEEP(50).EQ.0) 
     &                    .OR. 
     &                  ( (KEEP(50).NE.0).and. (.NOT.COMPRESSCB) )
     &                 )
     &               ) THEN
                   CALL ZMUMPS_ASM_SLAVE_TO_SLAVE(N, IFATH,
     &             IW, LIW,
     &             A, LA, NROWS_TO_STACK, NBCOLS, 
     &             INDICE_PERE_ARRAY_ARG,
     &             IW( COLLIST ), A(POSROW),
     &             OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &             ITLOC, RHS_MUMPS,
     &             FILS, ICNTL, KEEP,KEEP8,
     &             MYID, IS_ofType5or6, LDA_SON)
                   NBPROCFILS(STEP(IFATH)) =
     &                       NBPROCFILS(STEP(IFATH)) -
     &                       NROWS_TO_STACK                
#if ! defined(NO_XXNBPR)
                   IW( PTRIST(STEP(IFATH))+XXNBPR) =
     &               IW( PTRIST(STEP(IFATH))+XXNBPR) - NROWS_TO_STACK
#endif
                   EXIT
                 ELSE
                   CALL ZMUMPS_ASM_SLAVE_TO_SLAVE(N, IFATH,
     &             IW, LIW,
     &             A, LA, 1, NBCOLS_EFF, INDICE_PERE_ARRAY_ARG,
     &             IW( COLLIST ), A(POSROW),
     &             OPASSW, OPELIW, STEP, PTRIST, PTRAST,
     &             ITLOC, RHS_MUMPS,
     &             FILS, ICNTL, KEEP,KEEP8,
     &             MYID, IS_ofType5or6, LDA_SON)
                   NBPROCFILS(STEP(IFATH)) =
     &                         NBPROCFILS(STEP(IFATH)) - 1
#if ! defined(NO_XXNBPR)
                   IW( PTRIST(STEP(IFATH))+XXNBPR) =
     &               IW( PTRIST(STEP(IFATH))+XXNBPR) - 1
#endif
                 ENDIF
              ENDIF
            ENDDO
            IF (PDEST.EQ.PDEST_MASTER) THEN 
             IF (KEEP(219).NE.0) THEN
               IF(NSLAVES_PERE.GT.0 .AND. KEEP(50).EQ.2) THEN
                  IF (COMPRESSCB) THEN
                    WRITE(*,*) "Error 1 in PARPIV/ZMUMPS_MAPLIG"
                    CALL MUMPS_ABORT()
                  ELSE
                    POSROW = PTRAST(STEP(ISON))+SHIFTCB_SON+
     &                       int(NBROW(1)-1,8)*int(LDA_SON,8)
                  ENDIF
                  CALL ZMUMPS_BUF_MAX_ARRAY_MINSIZE(NFS4FATHER,IERR)
                  IF (IERR .NE.0) THEN
                    IF (LP .GT. 0) THEN
                      WRITE(LP, *) "MAX_ARRAY allocation failed"
                    ENDIF
                    IFLAG=-13
                    IERROR=NFS4FATHER
                    RETURN
                  ENDIF
                  ITMP=-9999
                  IF ( LMAP_LOC-NBROW(1)+1-KEEP253_LOC .NE. 0 ) THEN
                  CALL ZMUMPS_COMPUTE_MAXPERCOL(
     &                 A(POSROW),
     &       SIZFR-SHIFTCB_SON-int(NBROW(1)-1,8)*int(LDA_SON,8),
     &                 LDA_SON, LMAP_LOC-NBROW(1)+1-KEEP253_LOC,
     &                 BUF_MAX_ARRAY,NFS4FATHER,COMPRESSCB,ITMP)
                  ELSE
                       CALL ZMUMPS_SETMAXTOZERO(
     &                 BUF_MAX_ARRAY, NFS4FATHER)
                  ENDIF
                  CALL ZMUMPS_ASM_MAX(N, IFATH, IW, LIW, 
     &                 A, LA, ISON, NFS4FATHER,
     &                 BUF_MAX_ARRAY, PTLUST, PTRAST,
     &                 STEP, PIMASTER,
     &                 OPASSW,IWPOSCB,MYID, KEEP,KEEP8)
               ENDIF
             ENDIF 
             ISTCHK_LOC = PIMASTER(STEP(ISON))
               SAME_PROC= ISTCHK_LOC .LT. IWPOSCB
#if ! defined(NO_XXNBPR)
               IF ( SAME_PROC ) THEN
                 INBPROCFILS_SON = PTRIST(STEP(ISON))+XXNBPR
                 WRITE(*,*)
     &           "Internal error 0 in ZMUMPS_LOCAL_ASSEMBLY_TYPE2",
     &           INBPROCFILS_SON, PIMASTER(STEP(ISON))
                 CALL MUMPS_ABORT() 
               ELSE
                 INBPROCFILS_SON = PIMASTER(STEP(ISON))+XXNBPR
               ENDIF
#endif
#if ! defined(NO_XXNBPR)
               CALL CHECK_EQUAL( NBPROCFILS(STEP(ISON)),
     &                         IW(INBPROCFILS_SON) )
               IF ( IW(INBPROCFILS_SON) .EQ. 0 ) THEN
#else
               IF ( NBPROCFILS(STEP(ISON)) .EQ. 0 ) THEN
#endif
                 IF (SAME_PROC) THEN
                   CALL ZMUMPS_RESTORE_INDICES(N, ISON, IFATH,
     &               IWPOSCB, PIMASTER, PTLUST, IW, LIW, STEP,
     &               KEEP,KEEP8)
                 ENDIF
                 IF (SAME_PROC) THEN
                   ISTCHK_LOC = PTRIST(STEP(ISON))
                   PTRIST(STEP( ISON) ) = -99999999
                 ELSE
                   PIMASTER(STEP( ISON )) = -99999999
                 ENDIF
                 CALL ZMUMPS_FREE_BLOCK_CB(.FALSE., MYID, N,
     &            ISTCHK_LOC,
     &            PAMASTER(STEP(ISON)),
     &            IW, LIW, LRLU, LRLUS, IPTRLU, IWPOSCB,
     &            LA, KEEP,KEEP8, .FALSE.
     &            )
               ENDIF
#if ! defined(NO_XXNBPR)
             CALL CHECK_EQUAL( NBPROCFILS(STEP(IFATH)),
     &                         IW(PTLUST(STEP(IFATH))+XXNBPR) )
             IF ( IW(PTLUST(STEP(IFATH))+XXNBPR) .EQ. 0
#else
             IF ( NBPROCFILS(STEP(IFATH)) .EQ. 0
#endif
     &       ) THEN
               CALL ZMUMPS_INSERT_POOL_N( N, IPOOL, LPOOL,
     &           PROCNODE_STEPS,
     &           SLAVEF, KEEP(28), KEEP(76), KEEP(80),
     &           KEEP(47), STEP, IFATH+N )
               IF (KEEP(47) .GE. 3) THEN
                 CALL ZMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &          IPOOL, LPOOL, 
     &          PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &          MYID, STEP, N, ND, FILS )
               ENDIF
             END IF
            ELSE
               CALL ZMUMPS_ASM_SLAVE_TO_SLAVE_END
     &         (N, IFATH, IW, LIW,
     &         NBROW(I), STEP, PTRIST, ITLOC, RHS_MUMPS,
     &         KEEP,KEEP8)
            END IF
      RETURN
      END SUBROUTINE ZMUMPS_LOCAL_ASSEMBLY_TYPE2
