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
      RECURSIVE SUBROUTINE SMUMPS_PROCESS_SYM_BLOCFACTO( 
     &   COMM_LOAD, ASS_IRECV,
     &   BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, NBPROCFILS,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM, IFLAG, IERROR, NBFIN,
     &
     &    PTLUST_S, PTRFAC, root, OPASSW, OPELIW,
     &    ITLOC, RHS_MUMPS, FILS,  
     &    PTRARW, PTRAIW, INTARR, DBLARR,
     &    ICNTL,KEEP,KEEP8,DKEEP,IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE  
     &    )
      USE SMUMPS_COMM_BUFFER
      USE SMUMPS_LOAD
      USE SMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INCLUDE 'mumps_headers.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 40 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      REAL    DKEEP(130)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER(8) IPTRLU, LRLU, LRLUS, LA, POSFAC
      INTEGER COMP
      INTEGER IFLAG, IERROR, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28))
      INTEGER(8) PTRAST(KEEP(28)), PTRFAC(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER NBPROCFILS( KEEP(28) ), STEP(N), 
     & PIMASTER(KEEP(28))
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28)),
     &        ITLOC(N+KEEP(253)), FILS(N), ND(KEEP(28))
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER FRERE_STEPS(KEEP(28))
      INTEGER INTARR( max(1,KEEP(14)) )
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      REAL DBLARR( max(1,KEEP(13)) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER PIVI
      INTEGER (8) POSPV1,POSPV2,OFFDAG,LPOS1
      INTEGER J2
      REAL MULT1,MULT2, A11, DETPIV, A22, A12
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER LP
      INTEGER INODE, POSITION, NPIV, IERR
      INTEGER NCOL, LD_BLOCFACTO
      INTEGER(8) LAELL, POSBLOCFACTO
      INTEGER(8) POSELT
      INTEGER IOLDPS, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW, DEST
      INTEGER ICT11
      INTEGER(8) LPOS, LPOS2, DPOS, UPOS
      INTEGER (8) IPOS, KPOS
      INTEGER I, IPIV, FPERE, NSLAVES_TOT,
     &        NSLAVES_FOLLOW, NB_BLOC_FAC
      INTEGER IPOSK, JPOSK, NPIVSENT, Block, IROW, BLSIZE
      INTEGER allocok, TO_UPDATE_CPT_END
      REAL, DIMENSION(:), ALLOCATABLE :: UIP21K
      INTEGER, DIMENSION(:), ALLOCATABLE :: LIST_SLAVES_FOLLOW
      LOGICAL LASTBL
      INTEGER SRC_DESCBAND
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      REAL ONE,ALPHA
      PARAMETER (ONE = 1.0E0, ALPHA=-1.0E0)
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, NextPivDummy
      LOGICAL LAST_CALL
      TYPE(IO_BLOCK) :: MonBloc
      INTEGER LRELAY_INFO
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      POSITION = 0
      TO_UPDATE_CPT_END = -654321
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LASTBL = (NPIV.LE.0)
      IF (LASTBL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NSLAVES_TOT, 1,
     &                 MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NB_BLOC_FAC, 1,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL, 1,
     &                 MPI_INTEGER, COMM, IERR )
         LAELL = int(NPIV,8) * int(NCOL,8)
      IF ( NPIV.GT.0 ) THEN
       IF ( LRLU .LT. LAELL .OR. IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LAELL ) THEN
          IFLAG = -9
          CALL MUMPS_SET_IERROR(LAELL-LRLUS, IERROR)
          IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE IN SMUMPS_PROCESS_SYM_BLOCFACTO,
     & REAL WORKSPACE TOO SMALL"
          GOTO 700
        END IF
        CALL SMUMPS_COMPRE_NEW(N, KEEP(28), IW, LIW, A, LA,
     &       LRLU, IPTRLU,
     &       IWPOS, IWPOSCB, PTRIST, PTRAST,
     &       STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &       KEEP(IXSZ),COMP,DKEEP(97),MYID)
        IF ( LRLU .NE. LRLUS ) THEN
             WRITE(*,*) 'PB compress SMUMPS_PROCESS_SYM_BLOCFACTO,",
     &       " LRLU,LRLUS='
     &       ,LRLU,LRLUS
             IFLAG = -9
             CALL MUMPS_SET_IERROR(LAELL-LRLUS,IERROR)
             GOTO 700
        END IF
        IF ( IWPOS + NPIV - 1 .GT. IWPOSCB ) THEN
          IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE IN SMUMPS_PROCESS_SYM_BLOCFACTO,
     & INTEGER WORKSPACE TOO SMALL"
          IFLAG = -8
          IERROR = IWPOS + NPIV - 1 - IWPOSCB
          GOTO 700
        END IF
       END IF
       LRLU  = LRLU - LAELL
       LRLUS = LRLUS - LAELL
      ENDIF
      KEEP8(67) = min(LRLUS, KEEP8(67))
      POSBLOCFACTO = POSFAC
      POSFAC = POSFAC + LAELL
      CALL SMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                           LA-LRLUS,0_8,LAELL,KEEP,KEEP8,LRLUS)
      IF ( NPIV.EQ.0 ) THEN
        IPIV = 1 
      ELSE
        IPIV = IWPOS
        IWPOS = IWPOS + NPIV
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 IW( IPIV ), NPIV,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &              A(POSBLOCFACTO), NPIV*NCOL, MPI_REAL,
     &              COMM, IERR )
        LD_BLOCFACTO = NCOL
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 LRELAY_INFO, 1,
     &                 MPI_INTEGER, COMM, IERR )
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
        SRC_DESCBAND =
     &      MUMPS_PROCNODE( PROCNODE_STEPS(STEP(INODE)), SLAVEF )
          CALL SMUMPS_TREAT_DESCBAND( INODE, COMM_LOAD, ASS_IRECV,
     &      BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &      PTLUST_S, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &      IFLAG, IERROR, COMM,
     &      NBPROCFILS,
     &      IPOOL, LPOOL, LEAF,
     &      NBFIN, MYID, SLAVEF,
     &
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, PTRARW, PTRAIW,
     &      INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &      LPTRAR, NELT, FRTPTR, FRTELT, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &        )
          IF ( IFLAG .LT. 0 ) GOTO 600
      ENDIF
      IF ( IW( PTRIST(STEP(INODE)) + 3 + KEEP(IXSZ)) .EQ. 0 ) THEN
#if ! defined(NO_XXNBPR)
       CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE)),
     &                  IW(PTRIST(STEP(INODE))+XXNBPR))
       DO WHILE ( IW(PTRIST(STEP(INODE)) + XXNBPR) .NE. 0)
#else
       DO WHILE ( NBPROCFILS( STEP(INODE)) .NE. 0 )
#endif
        BLOCKING = .TRUE.
        SET_IRECV=.FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL SMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, CONTRIB_TYPE2,
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &     )
        IF ( IFLAG .LT. 0 ) GOTO 600
      END  DO
      ENDIF
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL SMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, MPI_ANY_TAG, 
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC, 
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &       )
      IOLDPS = PTRIST(STEP(INODE))
      POSELT = PTRAST(STEP(INODE))
      LCONT1 = IW( IOLDPS + KEEP(IXSZ))
      NASS1  = IW( IOLDPS + 1 + KEEP(IXSZ))
      IF ( NASS1 < 0 ) THEN
        NASS1 = -NASS1
        IW( IOLDPS + 1 + KEEP(IXSZ)) = NASS1
        IF (KEEP(55) .EQ. 0) THEN 
          CALL SMUMPS_ASM_SLAVE_ARROWHEADS(INODE, N, IW, LIW,
     &       IOLDPS, A, LA, POSELT, KEEP, ITLOC, FILS, PTRAIW,
     &       PTRARW, INTARR, DBLARR, RHS_MUMPS)
        ELSE
          CALL SMUMPS_ASM_SLAVE_ELEMENTS(INODE, N, NELT, IW, LIW,
     &       IOLDPS, A, LA, POSELT, KEEP, ITLOC, FILS, PTRAIW,
     &       PTRARW, INTARR, DBLARR, FRTPTR, FRTELT, RHS_MUMPS)
        ENDIF
      ENDIF
      NROW1  = IW( IOLDPS + 2 + KEEP(IXSZ))
      NPIV1  = IW( IOLDPS + 3 + KEEP(IXSZ))
      NSLAV1 = IW( IOLDPS + 5 + KEEP(IXSZ))
      NSLAVES_FOLLOW = NSLAV1 - XTRA_SLAVES_SYM
      HS     = 6 + NSLAV1 + KEEP(IXSZ)
      NCOL1  = LCONT1 + NPIV1
      IF ( LASTBL ) THEN
        TO_UPDATE_CPT_END = ( NSLAVES_TOT - NSLAVES_FOLLOW - 1 ) * 
     &                       NB_BLOC_FAC
      END IF
      IF (NPIV.GT.0) THEN
        IF ( NPIV1 + NCOL .NE. NASS1 ) THEN
          WRITE(*,*) 'SymBLFC Error: NPIV1 + NCOL .NE. NASS1 :',
     &               NPIV1,NCOL,NASS1
          CALL MUMPS_ABORT()
        END IF
        ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
        DO I = 1, NPIV
          PIVI = abs(IW(IPIV+I-1))
          IF (PIVI.EQ.I) CYCLE
          ISW = IW(ICT11+I)
          IW(ICT11+I) = IW(ICT11+PIVI)
          IW(ICT11+PIVI) = ISW
          IPOS = POSELT + int(NPIV1 + I - 1,8)
          KPOS = POSELT + int(NPIV1 + PIVI - 1,8)
          CALL sswap(NROW1, A(IPOS), NCOL1, A(KPOS), NCOL1)
        ENDDO
        ALLOCATE( UIP21K( NPIV * NROW1 ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
            IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR UIP21K IN SMUMPS_PROCESS_SYM_BLOCFACTO"
          IFLAG = -13
          IERROR = NPIV * NROW1
          GOTO 700
        END IF
        IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
          ALLOCATE( LIST_SLAVES_FOLLOW ( NSLAVES_FOLLOW ),
     &            stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR LIST_SLAVES_FOLLOW
     & IN SMUMPS_PROCESS_SYM_BLOCFACTO"
            IFLAG = -13
            IERROR = NSLAVES_FOLLOW
            GOTO 700
          END IF
          LIST_SLAVES_FOLLOW(1:NSLAVES_FOLLOW)=
     &    IW(IOLDPS+6+XTRA_SLAVES_SYM+KEEP(IXSZ):
     &     IOLDPS+5+XTRA_SLAVES_SYM+KEEP(IXSZ)+NSLAVES_FOLLOW)
        END IF
        CALL strsm( 'L', 'U', 'T', 'U', NPIV, NROW1, ONE,
     &               A( POSBLOCFACTO ), LD_BLOCFACTO,
     &               A(POSELT+int(NPIV1,8)), NCOL1 )
         LPOS = POSELT + int(NPIV1,8)
         UPOS = 1_8
         DO I = 1, NROW1
          UIP21K( UPOS: UPOS + int(NPIV-1,8) ) = 
     &                       A(LPOS: LPOS+int(NPIV-1,8))
          LPOS = LPOS + int(NCOL1,8)
          UPOS = UPOS + int(NPIV,8)
         END DO
        LPOS = POSELT + int(NPIV1,8)
        DPOS = POSBLOCFACTO
        I = 1
        DO
          IF(I .GT. NPIV) EXIT
          IF(IW(IPIV+I-1) .GT. 0) THEN
          A11 = ONE/A(DPOS)
            CALL sscal( NROW1, A11, A(LPOS), NCOL1 )
            LPOS = LPOS + 1_8
            DPOS = DPOS + int(LD_BLOCFACTO + 1,8)
            I = I+1
          ELSE
            POSPV1 = DPOS
            POSPV2 = DPOS+ int(LD_BLOCFACTO + 1,8)
            OFFDAG = POSPV1+1_8
            A11 = A(POSPV1)
            A22 = A(POSPV2)
            A12 = A(OFFDAG)
            DETPIV = A11*A22 - A12**2
            A22 = A11/DETPIV
            A11 = A(POSPV2)/DETPIV
            A12 = -A12/DETPIV
            LPOS1 = LPOS
            DO J2 = 1,NROW1
               MULT1 = A11*A(LPOS1)+A12*A(LPOS1+1_8)
               MULT2 = A12*A(LPOS1)+A22*A(LPOS1+1_8)
               A(LPOS1) = MULT1
               A(LPOS1+1_8) = MULT2
               LPOS1 = LPOS1 + int(NCOL1,8)
            ENDDO
            LPOS = LPOS + 2_8
            DPOS = POSPV2 + int(LD_BLOCFACTO + 1,8)
            I = I+2
          ENDIF
        ENDDO
      ENDIF
      IF (KEEP(201).eq.1) THEN
        MonBloc%INODE = INODE
        MonBloc%MASTER = .FALSE.
        MonBloc%Typenode = 2
        MonBloc%NROW = NROW1  
        MonBloc%NCOL = NCOL1  
        MonBloc%NFS  = NASS1
        MonBloc%LastPiv = NPIV1 + NPIV 
        MonBloc%LastPanelWritten_L = -9999 
        MonBloc%LastPanelWritten_U = -9999 
        NULLIFY(MonBloc%INDICES)
        MonBloc%Last = LASTBL
        STRAT = STRAT_TRY_WRITE 
        NextPivDummy      = -8888 
        LIWFAC = IW(IOLDPS+XXI)
        CALL MUMPS_GETI8(LAFAC, IW(IOLDPS+XXR))
        LAST_CALL=.FALSE.
        CALL SMUMPS_OOC_IO_LU_PANEL( STRAT, TYPEF_L, A(POSELT),
     &       LAFAC, MonBloc, NextPivDummy, NextPivDummy,
     &       IW(IOLDPS), LIWFAC, MYID, KEEP8(31), IFLAG,LAST_CALL)
      ENDIF
      IF (NPIV.GT.0) THEN
        LPOS2 = POSELT + int(NPIV1,8)
        UPOS = POSBLOCFACTO+int(NPIV,8)
        LPOS  = LPOS2 + int(NPIV,8)
        CALL sgemm('N','N', NCOL-NPIV,NROW1,NPIV,ALPHA,A(UPOS),NCOL,
     &           A(LPOS2),NCOL1,ONE, A(LPOS),NCOL1)
        DPOS = POSELT + int(NCOL1 - NROW1,8)
        IF ( NROW1 .GT. KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = NROW1
        ENDIF
        IF ( NROW1 .GT. 0 ) THEN
          DO IROW = 1, NROW1, BLSIZE
            Block = min( BLSIZE, NROW1 - IROW + 1 )
            DPOS  = POSELT + int(NCOL1 - NROW1,8)
     &            + int( IROW - 1, 8 ) * int( NCOL1 + 1, 8 )
            LPOS2 = POSELT + int(NPIV1,8)
     &            + int( IROW - 1, 8 ) * int( NCOL1, 8 )
            UPOS  = int( IROW - 1, 8 ) * int(NPIV, 8) + 1_8
            DO I = 1, Block
              CALL sgemv( 'T', NPIV, Block-I+1, ALPHA,
     &                A( LPOS2 + int(I - 1,8) * int(NCOL1,8) ), NCOL1,
     &                UIP21K( UPOS + int(NPIV,8) * int( I - 1, 8 ) ),
     &                1, ONE, A(DPOS+int(NCOL1+1,8)*int(I-1,8)),NCOL1 )
            END DO
           IF ( NROW1-IROW+1-Block .ne. 0 )
     &     CALL sgemm( 'T', 'N', Block, NROW1-IROW+1-Block, NPIV, ALPHA,
     &             UIP21K( UPOS ), NPIV,
     &             A( LPOS2 + int(Block,8) * int(NCOL1,8) ), NCOL1, ONE,
     &             A( DPOS + int(Block,8) * int(NCOL1,8) ), NCOL1 )
          ENDDO
        ENDIF
        FLOP1 = dble(NROW1) * dble(NPIV) *
     &           dble( 2 * NCOL  - NPIV + NROW1 +1 )
        FLOP1 = -FLOP1
        CALL SMUMPS_LOAD_UPDATE( 1, .FALSE., FLOP1, KEEP,KEEP8 )
      ENDIF 
      IW(IOLDPS+KEEP(IXSZ)) = IW(IOLDPS+KEEP(IXSZ)) - NPIV
      IW(IOLDPS + 3+KEEP(IXSZ)) = IW(IOLDPS+3+KEEP(IXSZ)) + NPIV
      IF (LASTBL) IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS + 3+KEEP(IXSZ))
      LRLU  = LRLU + LAELL
      LRLUS = LRLUS + LAELL
      POSFAC = POSFAC - LAELL
      IWPOS = IWPOS - NPIV
      CALL SMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                           LA-LRLUS,0_8,-LAELL,KEEP,KEEP8,LRLUS)
      IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
         IPOSK = NPIV1 + 1
         JPOSK = NCOL1 - NROW1 + 1
           NPIVSENT = NPIV
          IERR = -1
           DO WHILE ( IERR .eq. -1 )
            CALL SMUMPS_BUF_SEND_BLFAC_SLAVE(
     &                    INODE, NPIVSENT, FPERE,
     &                    IPOSK, JPOSK,
     &                    UIP21K, NROW1,
     &                    NSLAVES_FOLLOW,
     &                    LIST_SLAVES_FOLLOW(1),
     &                    COMM, 
     &                    IERR )
            IF (IERR .EQ. -1 ) THEN
              BLOCKING = .FALSE.
              SET_IRECV= .TRUE.
              MESSAGE_RECEIVED = .FALSE.
              CALL SMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &         MPI_ANY_SOURCE, MPI_ANY_TAG,
     &         STATUS, 
     &         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &         IWPOS, IWPOSCB, IPTRLU,
     &         LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &         PTLUST_S, PTRFAC,
     &         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &         IFLAG, IERROR, COMM,
     &         NBPROCFILS,
     &         IPOOL, LPOOL, LEAF,
     &         NBFIN, MYID, SLAVEF,
     &         root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, PTRARW, PTRAIW,
     &         INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &         LPTRAR, NELT, FRTPTR, FRTELT, 
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &           )
             IF ( IFLAG .LT. 0 ) GOTO 600
            END IF
           END DO
#if defined(IBC_TEST)
           WRITE(*,*) MYID,":Send2slave worked"
#endif
           IF ( IERR .eq. -2 ) THEN
              IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE, SEND BUFFER TOO SMALL DURING
     & SMUMPS_PROCESS_SYM_BLOCFACTO"
             WRITE(LP,*) "NPIV=", NPIV, "NROW1=",NROW1
             IFLAG = -17
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           IF ( IERR .eq. -3 ) THEN
              IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE, RECV BUFFER TOO SMALL DURING
     & SMUMPS_PROCESS_SYM_BLOCFACTO"
             IFLAG = -20
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           DEALLOCATE(LIST_SLAVES_FOLLOW)
      END IF
      IF ( NPIV .NE. 0 )  THEN
        IF (allocated(UIP21K)) DEALLOCATE( UIP21K )
      ENDIF
      IOLDPS = PTRIST(STEP(INODE))
      IF (LASTBL) THEN
         IW(IOLDPS+6+KEEP(IXSZ)) =  IW(IOLDPS+6+KEEP(IXSZ)) -
     &                            TO_UPDATE_CPT_END 
         IF ( IW(IOLDPS+6+KEEP(IXSZ) ) .eq. 0
     &        .and. KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 0
     &        .and. NSLAVES_TOT.NE.1)THEN
         DEST = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(INODE)), SLAVEF )
         CALL SMUMPS_BUF_SEND_1INT( INODE, DEST, END_NIV2_LDLT,
     &                             COMM, IERR )
         IF ( IERR .LT. 0 ) THEN
           write(*,*) ' Internal error in PROCESS_SYM_BLOCFACTO.'
           IFLAG = -99
           GOTO 700
         END IF
         ENDIF
      END IF
      IF (LASTBL) THEN 
        IF (IW(IOLDPS+6+KEEP(IXSZ)) .eq. 0 ) THEN 
          CALL SMUMPS_END_FACTO_SLAVE( COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE, 
     &    root,
     &    MYID, COMM,
     &    
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER,
     &    PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, NBPROCFILS,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &    INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &      )
        ENDIF
      ENDIF
 600  CONTINUE
#if defined(IBC_TEST)
      write(6,*) MYID,' :Exiting SMUMPS_PROCESS_SYM_BLOCFACTO for
     &INODE=', INODE
#endif
      RETURN
 700  CONTINUE
      CALL SMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE SMUMPS_PROCESS_SYM_BLOCFACTO
