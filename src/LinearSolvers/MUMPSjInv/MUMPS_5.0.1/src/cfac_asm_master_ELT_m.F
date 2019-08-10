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
      MODULE CMUMPS_FAC_ASM_MASTER_ELT_M
      CONTAINS
      SUBROUTINE CMUMPS_FAC_ASM_NIV1_ELT( COMM_LOAD, ASS_IRECV,
     &    NELT, FRT_PTR, FRT_ELT,
     &    N, INODE, IW, LIW, A, LA, INFO, ND, 
     &    FILS, FRERE, DAD, MAXFRW, root,
     &    OPASSW, OPELIW, PTRIST, PTLUST, PTRFAC, PTRAST, 
     &    STEP, PIMASTER, PAMASTER,PTRARW, 
     &    PTRAIW, ITLOC, RHS_MUMPS, NSTEPS, SON_LEVEL2,
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS, 
     &    ICNTL, KEEP,KEEP8,DKEEP,INTARR,DBLARR, 
     &
     &    NSTK_S,NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM,MYID,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &    PERM, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    )
      USE MUMPS_BUILD_SORT_INDEX_ELT_M
      USE CMUMPS_COMM_BUFFER
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER IZERO 
      PARAMETER (IZERO=0)
      INTEGER NELT,N,LIW,NSTEPS
      INTEGER(8) LA, LRLU, LRLUS, IPTRLU, POSFAC
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      REAL    DKEEP(130)
      INTEGER, INTENT(INOUT) :: INFO(2)
      INTEGER INODE,MAXFRW,
     &        IWPOSCB, COMP
      INTEGER, TARGET :: IWPOS
      INTEGER IDUMMY(1)
      INTEGER, PARAMETER :: LIDUMMY = 1
      INTEGER, TARGET :: IW(LIW)
      INTEGER ITLOC(N+KEEP(253)),
     &        PTRARW(NELT+1), PTRAIW(NELT+1), ND(KEEP(28)), PERM(N), 
     &        FILS(N), FRERE(KEEP(28)), DAD(KEEP(28)),
     &        PTRIST(KEEP(28)), PTLUST(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28))
      COMPLEX :: RHS_MUMPS(KEEP(255))
      INTEGER(8) :: PTRFAC(KEEP(28)), PTRAST(KEEP(28)),
     &              PAMASTER(KEEP(28))
      INTEGER COMM, NBFIN, SLAVEF, MYID
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      LOGICAL SON_LEVEL2
      COMPLEX, TARGET :: A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER        INTARR(max(1,KEEP(14)))
      COMPLEX DBLARR(max(1,KEEP(13)))
      INTEGER LPOOL, LEAF
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER IPOOL( LPOOL )
      INTEGER NBPROCFILS(KEEP(28)), NSTK_S(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER BUFR( LBUFR )
      INTEGER ETATASS
      LOGICAL COMPRESSCB
      INTEGER(8) :: LCB
      INTEGER  MUMPS_TYPENODE
      EXTERNAL MUMPS_TYPENODE
      INCLUDE 'mumps_headers.h'
      INTEGER LP, HS, HF
      LOGICAL LPOK
      INTEGER NBPANELS_L, NBPANELS_U
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFS4FATHER
      INTEGER(8) NFRONT8, LAELL8
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ
      INTEGER LREQ_OOC
      INTEGER(8) APOS, APOS2, LAPOS2
      INTEGER(8) POSELT, POSEL1, ICT12, ICT21
      INTEGER(8) IACHK
      INTEGER(8) JJ2
      INTEGER(8) LSTK8, SIZFR8
#if defined(ALLOW_NON_INIT)
      INTEGER(8) :: JJ8
#endif
      INTEGER SIZFI, NCB
      INTEGER JJ,J1,J2
      INTEGER NCOLS, NROWS, LDA_SON
      INTEGER NELIM,JJ1,J3,
     &        IORG, IBROT
      INTEGER JPOS,ICT11, IJROW
      INTEGER Pos_First_NUMORG,NUMORG,IOLDPS,
     &        NUMELT, ELBEG
      INTEGER AINPUT, 
     &        AII, J
      INTEGER NSLAVES, NSLSON, NPIVS, NPIV_ANA, NPIV
      INTEGER PTRCOL, ISLAVE, PDEST,LEVEL
      LOGICAL LEVEL1, NIV1
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      INTEGER ELTI, SIZE_ELTI
      INTEGER II, I
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER, POINTER :: SON_IWPOS
      INTEGER, POINTER, DIMENSION(:) :: SON_IW
      COMPLEX, POINTER, DIMENSION(:) :: SON_A
      INTEGER NCBSON
      LOGICAL SAME_PROC
      INTRINSIC real
      COMPLEX ZERO
      PARAMETER( ZERO = (0.0E0,0.0E0) )
      LOGICAL  MUMPS_INSSARBR, SSARBR
      EXTERNAL MUMPS_INSSARBR
      DOUBLE PRECISION FLOP1,FLOP1_EFF
      EXTERNAL MUMPS_IN_OR_ROOT_SSARBR
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR
      LP      = ICNTL(1)
      LPOK    = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      NFS4FATHER = -1
      ETATASS    = 0  
      COMPRESSCB=.FALSE.
      IN = INODE
      NBPROCFILS(STEP(IN)) = 0
      LEVEL = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      IF (LEVEL.NE.1) THEN
       write(6,*) 'Error1 in mpi51f_niv1 '
       CALL MUMPS_ABORT()
      END IF
      NSLAVES = 0
      HF = 6 + NSLAVES + KEEP(IXSZ)
      NUMELT = FRT_PTR(INODE+1) - FRT_PTR(INODE)
      IF ( NUMELT .ne. 0 ) THEN
        ELBEG  = FRT_PTR(INODE)
      ELSE
        ELBEG  = 1
      END IF
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      END DO
      NPIV_ANA=NUMORG
      NSTEPS = NSTEPS + 1
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      IF (ISON .NE. 0) THEN
        DO WHILE (ISON .GT. 0)
         NUMSTK = NUMSTK + 1
         SON_IW => IW
         NASS = NASS + SON_IW(PIMASTER(STEP(ISON))+1+KEEP(IXSZ))
         ISON = FRERE(STEP(ISON))
        END DO
      ENDIF
      NFRONT = ND(STEP(INODE)) + NASS + KEEP(253)
      NASS1 = NASS + NUMORG
      LREQ_OOC = 0
      IF (KEEP(201).EQ.1) THEN 
        CALL CMUMPS_OOC_GET_PP_SIZES(KEEP(50), NFRONT, NFRONT, NASS1,
     &       NBPANELS_L, NBPANELS_U, LREQ_OOC)
      ENDIF
      LREQ = HF + 2 * NFRONT + LREQ_OOC   
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL CMUMPS_COMPRE_NEW(N, KEEP(28),
     &        IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &        KEEP(IXSZ), COMP, DKEEP(97), MYID)
          IF (LRLU .NE. LRLUS) THEN
            WRITE( *, * ) 'PB compress CMUMPS_FAC_ASM_NIV1_ELT'
            WRITE( *, * ) 'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 270
          END IF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      END IF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      NIV1 = .TRUE.
        CALL MUMPS_ELT_BUILD_SORT( 
     &        NUMELT, FRT_ELT(ELBEG),
     &        MYID, INODE, N, IOLDPS, HF,
     &        NFRONT, NFRONT_EFF, PERM,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     &        IFSON, STEP, PIMASTER, PTRIST, PTRAIW, NELT, IW, LIW,
     &        INTARR, KEEP(14), ITLOC, FILS, FRERE,
     &        KEEP,
     &        SON_LEVEL2, NIV1, NBPROCFILS, INFO(1),
     &        DAD,PROCNODE_STEPS, SLAVEF,
     &        FRT_PTR, FRT_ELT, Pos_First_NUMORG,
     &        IDUMMY, LIDUMMY
     &        )
      IF (INFO(1).LT.0) GOTO 300
      IF (NFRONT_EFF.NE.NFRONT) THEN
        IF (NFRONT.GT.NFRONT_EFF) THEN
           IF(MUMPS_IN_OR_ROOT_SSARBR(PROCNODE_STEPS(STEP(INODE)),
     &          SLAVEF))THEN
              NPIV=NASS1-(NFRONT_EFF-ND(STEP(INODE)))
              CALL MUMPS_GET_FLOPS_COST(ND(STEP(INODE))+KEEP(253),
     &                                 NPIV,NPIV,
     &                                 KEEP(50),1,FLOP1)             
              NPIV=NPIV_ANA
              CALL MUMPS_GET_FLOPS_COST(ND(STEP(INODE))+KEEP(253),
     &                                 NPIV,NPIV,
     &                                 KEEP(50),1,FLOP1_EFF)
              CALL CMUMPS_LOAD_UPDATE(0,.FALSE.,FLOP1-FLOP1_EFF,
     &             KEEP,KEEP8)
           ENDIF
        IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF))
        NFRONT = NFRONT_EFF
        LREQ = HF + 2 * NFRONT + LREQ_OOC   
        ELSE
         IF (LPOK) THEN
          WRITE(LP,*) 
     &     ' ERROR 1 during ass_niv1_ELT', NFRONT, NFRONT_EFF
         ENDIF
         GOTO 270
        ENDIF
      ENDIF
      IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
        CALL CMUMPS_OOC_PP_SET_PTR(KEEP(50),
     &       NBPANELS_L, NBPANELS_U, NASS1, 
     &       IOLDPS + HF + 2 * NFRONT, IW, LIW)
      ENDIF
      NCB   = NFRONT - NASS1
      MAXFRW = max0(MAXFRW, NFRONT)
      ICT11 = IOLDPS + HF - 1 + NFRONT
      NFRONT8=int(NFRONT,8)
      LAELL8 = NFRONT8*NFRONT8
      IF (LRLU .LT. LAELL8) THEN
        IF (LRLUS .LT. LAELL8) THEN
          GOTO 280
        ELSE
          CALL CMUMPS_COMPRE_NEW(N, KEEP(28), IW, LIW, A, LA,
     &         LRLU, IPTRLU,
     &         IWPOS, IWPOSCB, PTRIST, PTRAST,
     &         STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &         KEEP(IXSZ), COMP, DKEEP(97),MYID)
          IF (LRLU .NE. LRLUS) THEN
           IF (LPOK) THEN
            WRITE(LP, * ) 'PB compress CMUMPS_FAC_ASM_NIV1_ELT'
            WRITE(LP, * ) 'LRLU,LRLUS=',LRLU,LRLUS
           ENDIF
           GOTO 280
          END IF
        END IF
      END IF
      LRLU = LRLU - LAELL8
      LRLUS = LRLUS - LAELL8
      KEEP8(67) = min(LRLUS, KEEP8(67))
      POSELT = POSFAC
      POSFAC = POSFAC + LAELL8
      SSARBR=MUMPS_INSSARBR(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      CALL CMUMPS_LOAD_MEM_UPDATE(SSARBR,.FALSE.,LA-LRLUS,0_8,LAELL8,
     &    KEEP,KEEP8,
     &     LRLUS)
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT + LAELL8 - 1_8
      A(POSELT:LAPOS2) = ZERO
#else
      IF ( KEEP(50) .eq. 0 .OR. NFRONT .LT. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL8 - 1_8
        A(POSELT:LAPOS2) = ZERO
      ELSE
        APOS = POSELT
        DO JJ8 = 0_8, int(NFRONT -1,8)
          A(APOS:APOS+JJ8) = ZERO
          APOS = APOS + NFRONT8
        END DO
      END IF
#endif
      NASS = NASS1
      PTRAST(STEP(INODE)) = POSELT
      PTRFAC(STEP(INODE)) = POSELT
      PTLUST(STEP(INODE)) = IOLDPS
      IW(IOLDPS+XXI) = LREQ  
      CALL MUMPS_STOREI8(LAELL8,IW(IOLDPS+XXR))
      IW(IOLDPS+XXS) = -9999
      IW(IOLDPS+XXN) = -99999
      IW(IOLDPS+XXP) = -99999
      IW(IOLDPS+XXA) = -99999
      IW(IOLDPS+XXF) = -99999
#if defined(NO_XXNBPR)
      IW(IOLDPS+XXNBPR)=-99999  
#else
      CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE)),IW(IOLDPS+XXNBPR))
#endif
      IW(IOLDPS + KEEP(IXSZ))   = NFRONT
      IW(IOLDPS + KEEP(IXSZ)+ 1) = 0
      IW(IOLDPS + KEEP(IXSZ) + 2) = -NASS1
      IW(IOLDPS + KEEP(IXSZ) + 3) = -NASS1
      IW(IOLDPS + KEEP(IXSZ) + 4) = STEP(INODE)
      IW(IOLDPS + KEEP(IXSZ) + 5)   = NSLAVES
      IF (NUMSTK.NE.0) THEN
        ISON = IFSON
        DO 220 IELL = 1, NUMSTK
          ISTCHK    = PIMASTER(STEP(ISON))
          SON_IW    => IW
          SON_IWPOS => IWPOS
          SON_A     => A
          LSTK      = SON_IW(ISTCHK + KEEP(IXSZ))
          LSTK8     = int(LSTK,8)
          NELIM     = SON_IW(ISTCHK + 1+KEEP(IXSZ))
          NPIVS     = SON_IW(ISTCHK + 3+KEEP(IXSZ))
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON    = SON_IW(ISTCHK + KEEP(IXSZ) + 5)
          HS        = 6 + KEEP(IXSZ) + NSLSON 
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LT.SON_IWPOS)
          IF ( SAME_PROC ) THEN
             COMPRESSCB =
     &           ( SON_IW(PTRIST(STEP(ISON))+XXS) .EQ. S_CB1COMP )
          ELSE
             COMPRESSCB =  ( SON_IW(ISTCHK + XXS) .EQ. S_CB1COMP )
          ENDIF
          LEVEL1    = NSLSON.EQ.0
          IF (.NOT.SAME_PROC) THEN
           NROWS = SON_IW( ISTCHK + KEEP(IXSZ) + 2)
          ELSE
           NROWS = NCOLS
          ENDIF
          SIZFI   = HS + NROWS + NCOLS 
          J1 = ISTCHK + HS + NROWS + NPIVS
          IF ( .NOT. LEVEL1 .AND. NELIM.EQ.0 ) GOTO 205
          IF (LEVEL1) THEN
           J2 = J1 + LSTK - 1
           IF (COMPRESSCB) THEN
             SIZFR8 = (LSTK8*(LSTK8+1_8)/2_8)
           ELSE
             SIZFR8 = LSTK8*LSTK8
           ENDIF
          ELSE
           IF ( KEEP(50).eq.0 ) THEN
             SIZFR8 = int(NELIM,8) * LSTK8
           ELSE
             SIZFR8 = int(NELIM,8) * int(NELIM,8)
           END IF
           J2 = J1 + NELIM - 1
          ENDIF
          OPASSW = OPASSW + dble(SIZFR8)
          IACHK = PAMASTER(STEP(ISON))
          IF ( KEEP(50) .eq. 0 ) THEN
            POSEL1 = PTRAST(STEP(INODE)) - NFRONT8
            IF (J2.GE.J1) THEN
              DO 170 JJ = J1, J2
                APOS = POSEL1 + int(SON_IW(JJ),8) * NFRONT8
                DO 160 JJ1 = 1, LSTK
                  JJ2 = APOS + int(SON_IW(J1 + JJ1 - 1) - 1,8)
                  A(JJ2) = A(JJ2) + SON_A(IACHK + int(JJ1 - 1,8))
  160           CONTINUE
                IACHK = IACHK + LSTK8
  170         CONTINUE
            END IF
          ELSE
            IF (LEVEL1) THEN
             LDA_SON = LSTK
            ELSE
             LDA_SON = NELIM
            ENDIF
            IF (COMPRESSCB) THEN
              LCB = SIZFR8
            ELSE
              LCB = int(LDA_SON,8)*int(J2 - J1 + 1,8)
            ENDIF
            IF (LCB .GT. 0) THEN
                CALL CMUMPS_LDLT_ASM_NIV12(A, LA, SON_A(IACHK),
     &           PTRAST(STEP( INODE )), NFRONT, NASS1,
     &           LDA_SON, LCB,
     &           SON_IW( J1 ), J2 - J1 + 1, NELIM, ETATASS, 
     &           COMPRESSCB
     &          )
            ENDIF
          ENDIF
  205     IF (LEVEL1) THEN 
           IF (SAME_PROC) ISTCHK = PTRIST(STEP(ISON))
           IF (SAME_PROC) THEN
             IF (KEEP(50).NE.0) THEN
              J2 = J1 + LSTK - 1
              DO JJ = J1, J2
               SON_IW(JJ) = SON_IW(JJ - NROWS)
              ENDDO
             ELSE
              J2 = J1 + LSTK - 1
              J3 = J1 + NELIM
              DO JJ = J3, J2
               SON_IW(JJ) = SON_IW(JJ - NROWS)
              ENDDO
              IF (NELIM .NE. 0) THEN
                J3 = J3 - 1
                DO JJ = J1, J3
                 JPOS = SON_IW(JJ) + ICT11
                 SON_IW(JJ) = IW(JPOS)
                ENDDO
              ENDIF
             ENDIF
           ENDIF
           IF ( SAME_PROC ) THEN
               PTRIST(STEP( ISON )) = -99999999
           ELSE
               PIMASTER(STEP( ISON )) = -99999999
           ENDIF
             CALL CMUMPS_FREE_BLOCK_CB(SSARBR, MYID, N, ISTCHK,
     &       IACHK,
     &       IW, LIW, LRLU, LRLUS, IPTRLU,
     &       IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &       )
          ELSE
           PDEST = ISTCHK + 6 + KEEP(IXSZ)
           NCBSON  = LSTK - NELIM
           PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
           DO ISLAVE = 0, NSLSON-1
             IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
              CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP,KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON, 
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX = PTRCOL + SHIFT_INDEX
               CALL CMUMPS_MAPLIG( COMM_LOAD, ASS_IRECV, 
     &                           BUFR, LBUFR, LBUFR_BYTES,
     &                           INODE, ISON, NSLAVES, IDUMMY,
     &                           NFRONT, NASS1,NFS4FATHER,
     &         TROW_SIZE, IW( INDX ),
     &         PROCNODE_STEPS,
     &         SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &         LRLUS, N, IW,
     &         LIW, A, LA, 
     &         PTRIST, PTLUST, PTRFAC, PTRAST, STEP,
     &         PIMASTER, PAMASTER, NSTK_S, COMP,
     &         INFO(1), INFO(2), MYID, COMM, NBPROCFILS, IPOOL, LPOOL,
     &         LEAF, NBFIN, ICNTL, KEEP,KEEP8,DKEEP, root,
     &         OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &         INTARR, DBLARR, ND, FRERE,
     &         NELT+1, NELT, FRT_PTR, FRT_ELT,
     &   
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &         )
               IF ( INFO(1) .LT. 0 ) GOTO 500
               EXIT
             ENDIF
           END DO
           IF (PIMASTER(STEP(ISON)).GT.0) THEN
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6 + KEEP(IXSZ)
            CALL  CMUMPS_BUF_SEND_MAPLIG( INODE, NFRONT, 
     &       NASS1, NFS4FATHER,ISON, MYID,
     &       IZERO, IDUMMY, IW(PTRCOL), NCBSON,
     &       COMM, IERR, IW(PDEST), NSLSON, 
     &       SLAVEF, 
     &       KEEP,KEEP8, STEP, N, 
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL CMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &         MPI_ANY_SOURCE, MPI_ANY_TAG,
     &         STATUS, 
     &         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &         IWPOS, IWPOSCB, IPTRLU,
     &         LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &         PTLUST, PTRFAC,
     &         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &         INFO(1), INFO(2), COMM,
     &         NBPROCFILS,
     &         IPOOL, LPOOL, LEAF,
     &         NBFIN, MYID, SLAVEF,
     &         root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, PTRARW, PTRAIW,
     &         INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &         NELT+1, NELT, FRT_PTR, FRT_ELT, 
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &          )
               IF ( INFO(1) .LT. 0 ) GOTO 500
            ENDIF
           END DO
           IF (IERR .EQ. -2) GOTO 290
           IF (IERR .EQ. -3) GOTO 295
           ENDIF
          ENDIF
      ISON = FRERE(STEP(ISON))
  220 CONTINUE
      END IF
      DO IELL=ELBEG,ELBEG+NUMELT-1
        ELTI = FRT_ELT(IELL)
        J1= PTRAIW(ELTI)
        J2= PTRAIW(ELTI+1)-1
        AII = PTRARW(ELTI)
        SIZE_ELTI = J2 - J1 + 1
        DO II=J1,J2
         I = INTARR(II)
         IF (KEEP(50).EQ.0) THEN
          AINPUT    = AII + II - J1
          ICT12 = POSELT + int(I-1,8) * NFRONT8
          DO JJ=J1,J2
           APOS2 = ICT12 + int(INTARR(JJ) - 1,8)
           A(APOS2) = A(APOS2) + DBLARR(AINPUT)
           AINPUT = AINPUT + SIZE_ELTI
          END DO
         ELSE
          ICT12 = POSELT + int(- NFRONT + I - 1,8)
          ICT21 = POSELT + int(I-1,8)*NFRONT8 - 1_8
          DO JJ=II,J2
           J =  INTARR(JJ)
           IF (I.LT.J) THEN
              APOS2 = ICT12 + int(J,8)*NFRONT8
           ELSE
              APOS2 = ICT21 + int(J,8)
           ENDIF
           A(APOS2) = A(APOS2) + DBLARR(AII)
           AII = AII + 1
          END DO
         END IF
        END DO
      END DO
      IF (KEEP(253).GT.0) THEN
       POSELT = PTRAST(STEP(INODE))
       IBROT = INODE
       IJROW = Pos_First_NUMORG  
       DO IORG = 1, NUMORG
        IF (KEEP(50).EQ.0) THEN
          DO JJ=1, KEEP(253)
            APOS = POSELT+
     &           int(IJROW-1,8)               * NFRONT8 +
     &           int(NFRONT-KEEP(253)+JJ-1,8)
            A(APOS) = A(APOS) + RHS_MUMPS( (JJ-1) * KEEP(254) + IBROT )
          ENDDO
        ELSE
          DO JJ=1, KEEP(253)
            APOS = POSELT+
     &           int(NFRONT-KEEP(253)+JJ-1,8)  * NFRONT8 +
     &           int(IJROW-1,8)
            A(APOS) = A(APOS) + RHS_MUMPS( (JJ-1) * KEEP(254) + IBROT )
          ENDDO
        ENDIF
       IBROT = FILS(IBROT)
       IJROW = IJROW+1
       ENDDO
      ENDIF
      GOTO 500
  270 CONTINUE
      INFO(1) = -8
      INFO(2) = LREQ
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE IN INTEGER ALLOCATION DURING CMUMPS_ASM_NIV1_ELT'
      ENDIF
      GOTO 490
  280 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE, WORKSPACE TOO SMALL DURING CMUMPS_ASM_NIV1_ELT'
      ENDIF
      INFO(1) = -9
      CALL MUMPS_SET_IERROR(LAELL8-LRLUS, INFO(2))
      GOTO 500
  290 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &  ' FAILURE, SEND BUFFER TOO SMALL DURING CMUMPS_ASM_NIV1_ELT'
      ENDIF
      INFO(1) = -17
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &  ' FAILURE, RECV BUFFER TOO SMALL DURING CMUMPS_ASM_NIV1_ELT'
      ENDIF
      INFO(1) = -20
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * ) ' FAILURE IN INTEGER',
     &                 ' DYNAMIC ALLOCATION DURING CMUMPS_ASM_NIV1_ELT'
      ENDIF
      INFO(1)   = -13
      INFO(2)  = NUMSTK 
  490 CALL  CMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_FAC_ASM_NIV1_ELT
      SUBROUTINE CMUMPS_FAC_ASM_NIV2_ELT( COMM_LOAD, ASS_IRECV,
     &    NELT, FRT_PTR, FRT_ELT,
     &    N, INODE, IW, LIW, A, LA, INFO,
     &    ND, FILS, FRERE, DAD,
     &    CAND,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    MAXFRW, root,
     &    OPASSW, OPELIW, PTRIST, PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_S,
     &    PTRAIW, ITLOC, RHS_MUMPS, NSTEPS, 
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     &    ICNTL, KEEP, KEEP8,DKEEP,INTARR,DBLARR,
     &    NBPROCFILS, PROCNODE_STEPS, SLAVEF, COMM,MYID,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, LEAF, IPOOL, LPOOL,
     &    PERM, MEM_DISTRIB
     &    )
      USE MUMPS_BUILD_SORT_INDEX_ELT_M
      USE CMUMPS_COMM_BUFFER
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER NELT, N,LIW,NSTEPS, NBFIN
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      REAL       DKEEP(130)
      INTEGER(8) :: LRLUS, LRLU, IPTRLU, POSFAC, LA
      INTEGER JJ
      INTEGER, INTENT(INOUT) :: INFO(2)
      INTEGER INODE,MAXFRW,
     &        LPOOL, LEAF, 
     &        IWPOS, 
     &        IWPOSCB, COMP, SLAVEF
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER IPOOL(LPOOL)
      INTEGER IW(LIW), ITLOC(N+KEEP(253)),
     &        PTRARW(NELT+1), PTRAIW(NELT+1), ND(KEEP(28)),
     &        FILS(N), FRERE(KEEP(28)), STEP(N), DAD (KEEP(28)),
     &        PTRIST(KEEP(28)), PTLUST(KEEP(28)),
     & PIMASTER(KEEP(28)),
     &        NSTK_S(KEEP(28)), PERM(N)
      COMPLEX :: RHS_MUMPS(KEEP(255))
      INTEGER(8) :: PTRFAC(KEEP(28)), PAMASTER(KEEP(28)),
     &              PTRAST(KEEP(28))
      INTEGER   CAND(SLAVEF+1,max(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      COMPLEX A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER        INTARR(max(1,KEEP(14)))
      COMPLEX DBLARR(max(1,KEEP(13)))
      INTEGER MYID, COMM
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER NBPROCFILS(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER BUFR( LBUFR )
      INCLUDE 'mumps_headers.h'
      INTEGER LP, HS, HF, HF_OLD, NCBSON, NSLAVES_OLD
      LOGICAL LPOK
      INTEGER NCBSON_MAX
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      LOGICAL COMPRESSCB
      INTEGER(8) :: LCB
      INTEGER :: IBC_SOURCE
      INTEGER NFS4FATHER
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ
      INTEGER(8) :: LAELL8
      INTEGER LREQ_OOC, NBPANELS_L, NBPANELS_U
      INTEGER NCB
      INTEGER J1,J2
      INTEGER(8) NFRONT8, LAPOS2, POSELT, POSEL1, LDAFS8,
     &           JJ2, IACHK, ICT12, ICT21
      INTEGER(8) :: JJ8
      INTEGER(8) APOS, APOS2
      INTEGER NELIM,LDAFS,JJ1,NPIVS,NCOLS,NROWS,
     &        IORG
      INTEGER LDA_SON, IJROW, IBROT
      INTEGER Pos_First_NUMORG,NBCOL,NUMORG,IOLDPS
      INTEGER AINPUT
      INTEGER NSLAVES, NSLSON
      INTEGER NBLIG, PTRCOL, PTRROW, ISLAVE, PDEST
      INTEGER ELTI, SIZE_ELTI
      INTEGER II, ELBEG, NUMELT, I, J, AII
      LOGICAL SAME_PROC, NIV1, SON_LEVEL2
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      logical :: force_cand
      INTEGER(8) APOSMAX
      REAL  MAXARR
      INTEGER INIV2, SIZE_TMP_SLAVES_LIST, allocok
      INTEGER NUMORG_SPLIT, TYPESPLIT,
     &        NCB_SPLIT, SIZE_LIST_SPLIT, NBSPLIT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TMP_SLAVES_LIST, COPY_CAND
      INTEGER IZERO 
      INTEGER IDUMMY(1)
      INTEGER PDEST1(1)
      INTEGER ETATASS
      PARAMETER( IZERO = 0 )
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE, MUMPS_TYPESPLIT
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, MUMPS_TYPESPLIT
      INTRINSIC real
      COMPLEX ZERO
      REAL RZERO
      PARAMETER( RZERO = 0.0E0 )
      PARAMETER( ZERO = (0.0E0,0.0E0) )
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SONROWS_PER_ROW
      INTEGER :: NMB_OF_CAND, NMB_OF_CAND_ORIG
      LOGICAL :: IS_ofType5or6, SPLIT_MAP_RESTART
      LP      = ICNTL(1)
      LPOK    = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      COMPRESSCB=.FALSE.
      ETATASS = 0  
      IN = INODE
      NBPROCFILS(STEP(IN)) = 0
      NSTEPS = NSTEPS + 1
      NUMELT = FRT_PTR(INODE+1) - FRT_PTR(INODE)
      IF ( NUMELT .NE. 0 ) THEN
        ELBEG = FRT_PTR(INODE)
      ELSE
        ELBEG = 1
      END IF
      NUMORG = 0
      DO WHILE (IN.GT.0)
        NUMORG = NUMORG + 1
        IN = FILS(IN)
      ENDDO
      NUMSTK = 0
      NASS = 0
      IFSON = -IN
      ISON = IFSON
      NCBSON_MAX = 0
      DO WHILE (ISON .GT. 0)
        NUMSTK = NUMSTK + 1
        IF ( KEEP(48)==5 .AND. 
     &       MUMPS_TYPENODE(PROCNODE_STEPS(STEP(ISON)),
     &       SLAVEF) .EQ. 1) THEN
          NCBSON_MAX =
     &      max(NCBSON_MAX,IW(PIMASTER(STEP(ISON))+KEEP(IXSZ)))
        ENDIF
        NASS = NASS + IW(PIMASTER(STEP(ISON)) + 1 + KEEP(IXSZ))
        ISON = FRERE(STEP(ISON))
      ENDDO
      NFRONT = ND(STEP(INODE)) + NASS + KEEP(253)
      MAXFRW = max0(MAXFRW, NFRONT)
      NASS1 = NASS + NUMORG
      NCB   = NFRONT - NASS1
      IF((KEEP(24).eq.0).or.(KEEP(24).eq.1)) then
         force_cand=.FALSE.
      ELSE
         force_cand=(mod(KEEP(24),2).eq.0)
      end if
      TYPESPLIT =  MUMPS_TYPESPLIT (PROCNODE_STEPS(STEP(INODE)), 
     &              SLAVEF)
      IS_ofType5or6 =    (TYPESPLIT.EQ.5 .OR. TYPESPLIT.EQ.6)
      ISTCHK            = PIMASTER(STEP(IFSON))
      PDEST             = ISTCHK + 6 + KEEP(IXSZ)
      NSLSON            = IW(ISTCHK + KEEP(IXSZ) + 5)
      SPLIT_MAP_RESTART = .FALSE.
      IF (force_cand) THEN
         INIV2                = ISTEP_TO_INIV2( STEP( INODE ))
         NMB_OF_CAND          = CAND( SLAVEF+1, INIV2 )
         NMB_OF_CAND_ORIG     = NMB_OF_CAND
         SIZE_TMP_SLAVES_LIST = NMB_OF_CAND
         IF  (IS_ofType5or6) THEN
           DO I=NMB_OF_CAND+1,SLAVEF
            IF ( CAND( I, INIV2 ).LT.0) EXIT
            NMB_OF_CAND = NMB_OF_CAND +1
           ENDDO
           SIZE_TMP_SLAVES_LIST = NSLSON-1
          WRITE(6,*) "NMB_OF_CAND, SIZE_TMP_SLAVES_LIST ", 
     & NMB_OF_CAND, SIZE_TMP_SLAVES_LIST
           IF (INODE.EQ.-999999) THEN
              SPLIT_MAP_RESTART = .TRUE.
           ENDIF
         ENDIF
         IF (IS_ofType5or6.AND.SPLIT_MAP_RESTART) THEN
           TYPESPLIT     = 4
           IS_ofType5or6 = .FALSE.
           SIZE_TMP_SLAVES_LIST = NMB_OF_CAND 
           CAND (SLAVEF+1, INIV2) = SIZE_TMP_SLAVES_LIST
         ENDIF
      ELSE
         INIV2 = 1
         SIZE_TMP_SLAVES_LIST = SLAVEF - 1 
         NMB_OF_CAND          =  SLAVEF - 1
         NMB_OF_CAND_ORIG     =  SLAVEF - 1
      ENDIF
      ALLOCATE(TMP_SLAVES_LIST(SIZE_TMP_SLAVES_LIST),stat=allocok)
      IF (allocok > 0 ) THEN
        GOTO 265
      ENDIF
       TYPESPLIT =  MUMPS_TYPESPLIT (PROCNODE_STEPS(STEP(INODE)), 
     &              SLAVEF)
       IF  ( (TYPESPLIT.EQ.4) 
     &               .OR.(TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6) 
     &     )  THEN
        IF (TYPESPLIT.EQ.4) THEN
         ALLOCATE(COPY_CAND(SLAVEF+1),stat=allocok)
         IF (allocok > 0 ) THEN
           GOTO 245
         ENDIF
         CALL CMUMPS_SPLIT_PREP_PARTITION (
     &      INODE, STEP, N, SLAVEF, 
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      CAND(1,INIV2), ICNTL, COPY_CAND,
     &      NBSPLIT, NUMORG_SPLIT, TMP_SLAVES_LIST(1),
     &      SIZE_TMP_SLAVES_LIST 
     &                                    )
         NCB_SPLIT = NCB-NUMORG_SPLIT
         SIZE_LIST_SPLIT = SIZE_TMP_SLAVES_LIST - NBSPLIT
         CALL CMUMPS_LOAD_SET_PARTITION( NCBSON_MAX, SLAVEF, KEEP,KEEP8,
     &     ICNTL, COPY_CAND,
     &     MEM_DISTRIB(0), NCB_SPLIT, NFRONT, NSLAVES,
     &     TAB_POS_IN_PERE(1,INIV2),
     &     TMP_SLAVES_LIST(NBSPLIT+1),
     &     SIZE_LIST_SPLIT,INODE )
         DEALLOCATE (COPY_CAND)
         CALL CMUMPS_SPLIT_POST_PARTITION (
     &      INODE, STEP, N, SLAVEF, NBSPLIT, NCB,
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      ICNTL, 
     &      TAB_POS_IN_PERE(1,INIV2),
     &      NSLAVES
     &                                    )
         IF (SPLIT_MAP_RESTART) THEN
          IS_ofType5or6 = .TRUE.
          TYPESPLIT =  MUMPS_TYPESPLIT (PROCNODE_STEPS(STEP(INODE)), 
     &              SLAVEF)
          CAND( SLAVEF+1, INIV2 ) = NMB_OF_CAND_ORIG
         ENDIF
        ELSE
         ISTCHK    = PIMASTER(STEP(IFSON))
         PDEST     = ISTCHK + 6 + KEEP(IXSZ)
         NSLSON    = IW(ISTCHK + KEEP(IXSZ) + 5)
         CALL CMUMPS_SPLIT_PROPAGATE_PARTI (
     &      INODE, TYPESPLIT, IFSON, 
     &      CAND(1,INIV2), NMB_OF_CAND_ORIG,
     &      IW(PDEST), NSLSON,
     &      STEP, N, SLAVEF, 
     &      PROCNODE_STEPS, KEEP, DAD, FILS,
     &      ICNTL, ISTEP_TO_INIV2, INIV2,
     &      TAB_POS_IN_PERE, NSLAVES, 
     &      TMP_SLAVES_LIST,
     &      SIZE_TMP_SLAVES_LIST
     &                                    )
        ENDIF
       ELSE
        CALL CMUMPS_LOAD_SET_PARTITION( NCBSON_MAX, SLAVEF, KEEP,KEEP8,
     &     ICNTL, CAND(1,INIV2),
     &     MEM_DISTRIB(0), NCB, NFRONT, NSLAVES,
     &     TAB_POS_IN_PERE(1,INIV2),
     &     TMP_SLAVES_LIST,
     &     SIZE_TMP_SLAVES_LIST,INODE )
       ENDIF
      HF   = NSLAVES + 6 + KEEP(IXSZ)
      LREQ_OOC = 0
      IF (KEEP(201).EQ.1) THEN
        CALL CMUMPS_OOC_GET_PP_SIZES(KEEP(50), NASS1, NFRONT, NASS1,
     &                               NBPANELS_L, NBPANELS_U, LREQ_OOC)
      ENDIF
      LREQ = HF + 2 * NFRONT + LREQ_OOC
      IF ((IWPOS + LREQ -1) .GT. IWPOSCB) THEN
          CALL CMUMPS_COMPRE_NEW(N, KEEP(28),
     &        IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER,
     &        KEEP(216),LRLUS,KEEP(IXSZ),
     &        COMP, DKEEP(97), MYID)
          IF (LRLU .NE. LRLUS) THEN
           IF (LPOK) THEN
            WRITE(LP, * ) 'PB compress CMUMPS_FAC_ASM_NIV2_ELT'
            WRITE(LP, * ) 'LRLU,LRLUS=',LRLU,LRLUS
           ENDIF
           GOTO 270
          ENDIF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      ENDIF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      NIV1 = .FALSE.
      ALLOCATE(SONROWS_PER_ROW(NFRONT-NASS1), stat=allocok)
      IF (allocok > 0) GOTO 275
        CALL MUMPS_ELT_BUILD_SORT(
     &        NUMELT, FRT_ELT(ELBEG),
     &        MYID, INODE, N, IOLDPS, HF,
     &        NFRONT, NFRONT_EFF, PERM,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     &        IFSON, STEP, PIMASTER, PTRIST, PTRAIW, NELT, IW, LIW,
     &        INTARR, KEEP(14), ITLOC, FILS, FRERE,
     &        KEEP, SON_LEVEL2, NIV1, NBPROCFILS, INFO(1),
     &        DAD,PROCNODE_STEPS, SLAVEF,
     &        FRT_PTR, FRT_ELT, Pos_First_NUMORG,
     &        SONROWS_PER_ROW, NFRONT - NASS1)
      IF (INFO(1).LT.0) GOTO 250
      IF ( NFRONT .NE. NFRONT_EFF ) THEN
        IF (
     &        (TYPESPLIT.EQ.5) .OR. (TYPESPLIT.EQ.6)) THEN
          WRITE(6,*) ' Internal error 1 in fac_ass due to splitting ',
     &     ' INODE, NFRONT, NFRONT_EFF =', INODE, NFRONT, NFRONT_EFF 
          WRITE(6,*) ' SPLITTING NOT YET READY FOR THAT'
          CALL MUMPS_ABORT()
        ENDIF
        IF (NFRONT.GT.NFRONT_EFF) THEN
            NCB    = NFRONT_EFF - NASS1
            NSLAVES_OLD = NSLAVES
            HF_OLD      = HF
            IF (TYPESPLIT.EQ.4) THEN
             ALLOCATE(COPY_CAND(SLAVEF+1),stat=allocok)
             IF (allocok > 0 ) THEN
               GOTO 245
             ENDIF
             CALL CMUMPS_SPLIT_PREP_PARTITION (
     &          INODE, STEP, N, SLAVEF, 
     &          PROCNODE_STEPS, KEEP, DAD, FILS,
     &          CAND(1,INIV2), ICNTL, COPY_CAND,
     &          NBSPLIT, NUMORG_SPLIT, TMP_SLAVES_LIST(1),
     &          SIZE_TMP_SLAVES_LIST 
     &                                    )
             NCB_SPLIT = NCB-NUMORG_SPLIT
             SIZE_LIST_SPLIT = SIZE_TMP_SLAVES_LIST - NBSPLIT
             CALL CMUMPS_LOAD_SET_PARTITION( NCBSON_MAX, 
     &         SLAVEF, KEEP,KEEP8,
     &         ICNTL, COPY_CAND,
     &         MEM_DISTRIB(0), NCB_SPLIT, NFRONT_EFF, NSLAVES,
     &         TAB_POS_IN_PERE(1,INIV2),
     &         TMP_SLAVES_LIST(NBSPLIT+1),
     &         SIZE_LIST_SPLIT,INODE )
             DEALLOCATE (COPY_CAND)
             CALL CMUMPS_SPLIT_POST_PARTITION (
     &          INODE, STEP, N, SLAVEF, NBSPLIT, NCB,
     &          PROCNODE_STEPS, KEEP, DAD, FILS,
     &          ICNTL, 
     &          TAB_POS_IN_PERE(1,INIV2),
     &          NSLAVES
     &                                    )
            ELSE
             CALL CMUMPS_LOAD_SET_PARTITION( NCBSON_MAX,
     &       SLAVEF, KEEP, KEEP8, ICNTL,
     &       CAND(1,INIV2),
     &       MEM_DISTRIB(0), NCB, NFRONT_EFF, NSLAVES,
     &       TAB_POS_IN_PERE(1,INIV2),
     &       TMP_SLAVES_LIST, SIZE_TMP_SLAVES_LIST,INODE )
            ENDIF
            HF = NSLAVES + 6 + KEEP(IXSZ)
            IWPOS = IWPOS - ((2*NFRONT)-(2*NFRONT_EFF)) -
     &                   (NSLAVES_OLD - NSLAVES)
            IF (NSLAVES_OLD .NE. NSLAVES) THEN
              IF (NSLAVES_OLD > NSLAVES) THEN
               DO JJ=0,2*NFRONT_EFF-1
                 IW(IOLDPS+HF+JJ)=IW(IOLDPS+HF_OLD+JJ)
               ENDDO
              ELSE
               IF (IWPOS - 1 > IWPOSCB ) GOTO 270
               DO JJ=2*NFRONT_EFF-1, 0, -1
                 IW(IOLDPS+HF+JJ) = IW(IOLDPS+HF_OLD+JJ)
               ENDDO
              END IF
            END IF
            NFRONT = NFRONT_EFF
            LREQ = HF + 2 * NFRONT + LREQ_OOC
        ELSE
          IF (LPOK) THEN
           WRITE(LP,*) ' INTERNAL ERROR 2 during ass_niv2'
          ENDIF
          GOTO 270
        ENDIF
      ENDIF
      NFRONT8=int(NFRONT,8)
      IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
        CALL CMUMPS_OOC_PP_SET_PTR(KEEP(50),
     &       NBPANELS_L, NBPANELS_U, NASS1, 
     &       IOLDPS + HF + 2 * NFRONT, IW, LIW)
      ENDIF
      MAXFRW = max0(MAXFRW, NFRONT)
      PTLUST(STEP(INODE)) = IOLDPS
      IW(IOLDPS+KEEP(IXSZ))     = NFRONT
      IW(IOLDPS + 1+KEEP(IXSZ)) = 0
      IW(IOLDPS + 2+KEEP(IXSZ)) = -NASS1
      IW(IOLDPS + 3+KEEP(IXSZ)) = -NASS1
      IW(IOLDPS + 4+KEEP(IXSZ)) = STEP(INODE)
      IW(IOLDPS+5+KEEP(IXSZ)) = NSLAVES
      IW(IOLDPS+6+KEEP(IXSZ):IOLDPS+5+NSLAVES+KEEP(IXSZ))=
     &                     TMP_SLAVES_LIST(1:NSLAVES)
#if defined(OLD_LOAD_MECHANISM)
#if ! defined (CHECK_COHERENCE)
      IF ( KEEP(73) .EQ. 0 ) THEN
#endif
#endif
        CALL CMUMPS_LOAD_MASTER_2_ALL(MYID, SLAVEF, COMM_LOAD,
     &     TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &     NASS1, KEEP, KEEP8, IW(IOLDPS+6+KEEP(IXSZ)), NSLAVES,INODE)
#if defined(OLD_LOAD_MECHANISM)
#if ! defined (CHECK_COHERENCE) 
      ENDIF
#endif
#endif
      IF(KEEP(86).EQ.1)THEN
         IF(mod(KEEP(24),2).eq.0)THEN
            CALL CMUMPS_LOAD_SEND_MD_INFO(SLAVEF,
     &           CAND(SLAVEF+1,INIV2),
     &           CAND(1,INIV2),
     &           TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &           NASS1, KEEP,KEEP8, TMP_SLAVES_LIST, 
     &           NSLAVES,INODE)
         ELSEIF((KEEP(24).EQ.0).OR.(KEEP(24).EQ.1))THEN
            CALL CMUMPS_LOAD_SEND_MD_INFO(SLAVEF,
     &           SLAVEF-1,
     &           TMP_SLAVES_LIST,
     &           TAB_POS_IN_PERE(1,ISTEP_TO_INIV2(STEP(INODE))),
     &           NASS1, KEEP,KEEP8, TMP_SLAVES_LIST, 
     &           NSLAVES,INODE)
         ENDIF
      ENDIF
      DEALLOCATE(TMP_SLAVES_LIST)
      IF (KEEP(50).EQ.0) THEN
        LAELL8 = int(NASS1,8) * NFRONT8
        LDAFS = NFRONT
        LDAFS8 = NFRONT8
      ELSE
        LAELL8 = int(NASS1,8)*int(NASS1,8)
        IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2)
     &     LAELL8 = LAELL8+int(NASS1,8)
        LDAFS = NASS1
        LDAFS8 = int(NASS1,8)
      ENDIF
      IF (LRLU .LT. LAELL8) THEN
        IF (LRLUS .LT. LAELL8) THEN
          GOTO 280
        ELSE
          CALL CMUMPS_COMPRE_NEW(N, KEEP(28), IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, KEEP(216), LRLUS,
     &        KEEP(IXSZ), COMP, DKEEP(97), MYID)
          IF (LRLU .NE. LRLUS) THEN
           IF (LPOK) THEN
            WRITE(LP, * ) 'PB compress CMUMPS_FAC_ASM_NIV2_ELT'
            WRITE(LP, * ) 'LRLU,LRLUS=',LRLU,LRLUS
           ENDIF
           GOTO 280
          ENDIF
        ENDIF
      ENDIF
      LRLU = LRLU - LAELL8
      LRLUS = LRLUS - LAELL8
      KEEP8(67) = min(LRLUS, KEEP8(67))
      POSELT = POSFAC
      PTRAST(STEP(INODE)) = POSELT
      PTRFAC(STEP(INODE)) = POSELT
      POSFAC = POSFAC + LAELL8
      IW(IOLDPS+XXI)   = LREQ  
      CALL MUMPS_STOREI8(LAELL8,IW(IOLDPS+XXR))
      IW(IOLDPS+XXS) =  -9999
      IW(IOLDPS+XXN) = -99999   
      IW(IOLDPS+XXP) = -99999   
      IW(IOLDPS+XXA) = -99999
      IW(IOLDPS+XXF) = -99999
#if defined(NO_XXNBPR)
      IW(IOLDPS+XXNBPR)=-99999  
#else
      CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE)), IW(IOLDPS+XXNBPR))
#endif
      CALL CMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,LA-LRLUS,0_8,LAELL8,
     & KEEP,KEEP8,
     & LRLUS)
      POSEL1 = POSELT - LDAFS8
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT + LAELL8 - 1_8
      A(POSELT:LAPOS2) = ZERO
#else
      IF ( KEEP(50) .eq. 0 .OR. LDAFS .lt. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL8 - 1_8
        A(POSELT:LAPOS2) = ZERO
      ELSE
        APOS = POSELT
        DO JJ8 = 0_8, LDAFS8 - 1_8
          A(APOS:APOS+JJ8) = ZERO
          APOS = APOS + LDAFS8
        END DO
        IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
          A(APOS:APOS+LDAFS8-1_8)=ZERO
        ENDIF
      END IF
#endif
      IF ((NUMSTK.NE.0).AND.(NASS.NE.0)) THEN
        ISON = IFSON
        DO 220 IELL = 1, NUMSTK
          ISTCHK = PIMASTER(STEP(ISON))
          NELIM = IW(ISTCHK + KEEP(IXSZ) + 1)
          IF (NELIM.EQ.0) GOTO 210
          LSTK    = IW(ISTCHK + KEEP(IXSZ))
          NPIVS   = IW(ISTCHK + KEEP(IXSZ) + 3)
          IF (NPIVS.LT.0) NPIVS=0
          NSLSON  = IW(ISTCHK + KEEP(IXSZ) + 5)
          HS      = 6 + KEEP(IXSZ) + NSLSON 
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LT.IWPOS)
          IF ( SAME_PROC ) THEN
            COMPRESSCB=( IW(PTRIST(STEP(ISON))+XXS) .EQ. S_CB1COMP )
          ELSE
            COMPRESSCB=( IW(ISTCHK + XXS) .EQ. S_CB1COMP )
          ENDIF
          IF (.NOT.SAME_PROC) THEN
           NROWS = IW(ISTCHK + KEEP(IXSZ) + 2)
          ELSE
           NROWS = NCOLS
          ENDIF
          OPASSW = OPASSW + dble(NELIM*LSTK)
          J1 = ISTCHK + HS + NROWS + NPIVS
          J2 = J1 + NELIM - 1
          IACHK = PAMASTER(STEP(ISON))
          IF (KEEP(50).eq.0) THEN
           IF (IS_ofType5or6) THEN
            APOS = POSELT  
            DO JJ8 = 1_8, int(NELIM,8)*int(LSTK,8)
             A(APOS+JJ8-1_8) = A(APOS+JJ8-1_8) + A(IACHK+JJ8-1_8)
            ENDDO
           ELSE
            DO 170 JJ = J1, J2
             APOS = POSEL1 + int(IW(JJ),8) * LDAFS8
             DO 160 JJ1 = 1, LSTK
              JJ2 = APOS + int(IW(J1 + JJ1 - 1),8) - 1_8
              A(JJ2) = A(JJ2) + A(IACHK + int(JJ1,8) - 1_8)
  160        CONTINUE
             IACHK = IACHK + int(LSTK,8)
  170       CONTINUE
           ENDIF
          ELSE
            IF (NSLSON.EQ.0) THEN
             LDA_SON = LSTK
            ELSE
             LDA_SON = NELIM
            ENDIF
            IF (COMPRESSCB) THEN
              LCB = (int(NELIM,8)*int(NELIM+1,8))/2_8
            ELSE
              LCB = int(LDA_SON,8)*int(NELIM,8)
            ENDIF
            IF (LCB .GT. 0) THEN
              CALL CMUMPS_LDLT_ASM_NIV12(A, LA, A(IACHK),
     &           POSELT, LDAFS, NASS1,
     &           LDA_SON, LCB,
     &           IW( J1 ), NELIM, NELIM, ETATASS,
     &           COMPRESSCB
     &          )
            ENDIF
          ENDIF
  210     ISON = FRERE(STEP(ISON))
  220   CONTINUE
      ENDIF
      APOSMAX = POSELT + int(NASS1,8)*int(NASS1,8)
      IF (KEEP(219).NE.0) THEN
        IF (KEEP(50).EQ.2) THEN
          A( APOSMAX: APOSMAX+int(NASS1-1,8))=ZERO
        ENDIF
      ENDIF
      DO IELL=ELBEG,ELBEG+NUMELT-1
        ELTI = FRT_ELT(IELL)
        J1= PTRAIW(ELTI)
        J2= PTRAIW(ELTI+1)-1
        AII = PTRARW(ELTI)
        SIZE_ELTI = J2 - J1 + 1
        DO II=J1,J2
         I = INTARR(II)
         IF (KEEP(50).EQ.0) THEN
          IF (I.LE.NASS1) THEN
           AINPUT    = AII + II - J1
           ICT12 = POSELT + int(I-1,8) * LDAFS8
           DO JJ=J1,J2
            APOS2 = ICT12 + int(INTARR(JJ) - 1,8)
            A(APOS2) = A(APOS2) + DBLARR(AINPUT)
            AINPUT = AINPUT + SIZE_ELTI
           END DO
          ENDIF
         ELSE
          ICT12 = POSELT - LDAFS8 + int(I,8) - 1_8
          ICT21 = POSELT + int(I-1,8)*LDAFS8 - 1_8
          IF ( I .GT. NASS1 ) THEN
           IF (KEEP(219).NE.0 .AND. KEEP(50).EQ.2) THEN
              AINPUT=AII
              DO JJ=II,J2
               J=INTARR(JJ)
               IF (J.LE.NASS1) THEN
                A(APOSMAX+int(J-1,8))=cmplx(
     &              max(real(A(APOSMAX+int(J-1,8))),
     &                  abs(DBLARR(AINPUT))),
     &                 kind=kind(A)
     &                 )
               ENDIF
               AINPUT=AINPUT+1
              ENDDO
           ENDIF
           AII = AII + J2 - II + 1
           CYCLE
          ELSE
            IF (KEEP(219).NE.0) THEN
              MAXARR = RZERO
            ENDIF
            DO JJ=II,J2
              J =  INTARR(JJ)
              IF ( J .LE. NASS1) THEN
                IF (I.LT.J) THEN
                  APOS2 = ICT12 + int(J,8)*LDAFS8
                ELSE
                  APOS2 = ICT21 + int(J,8)
                ENDIF
                A(APOS2) = A(APOS2) + DBLARR(AII)
              ELSE IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
                MAXARR = max(MAXARR,abs(DBLARR(AII)))
              ENDIF
              AII = AII + 1
            END DO
            IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2) THEN
                A(APOSMAX+int(I-1,8)) = cmplx(
     &             max( MAXARR, real(A(APOSMAX+int(I-1,8)))),
     &             kind=kind(A)
     &          )
            ENDIF
          ENDIF 
         END IF 
        END DO
      END DO
      IF (KEEP(253).GT.0) THEN
       POSELT = PTRAST(STEP(INODE))
       IBROT = INODE
       IJROW = Pos_First_NUMORG  
       DO IORG = 1, NUMORG
        IF (KEEP(50).EQ.0) THEN
          DO JJ = 1, KEEP(253)
            APOS = POSELT +
     &             int(IJROW-1,8) * int(LDAFS,8) +
     &             int(LDAFS-KEEP(253)+JJ-1,8)
            A(APOS) = A(APOS) + RHS_MUMPS( (JJ-1) * KEEP(254) + IBROT )
          ENDDO
        ENDIF
        IBROT = FILS(IBROT)
        IJROW = IJROW+1
       ENDDO
      ENDIF
      PTRCOL = IOLDPS + HF + NFRONT 
      PTRROW = IOLDPS + HF + NASS1 
      PDEST  = IOLDPS + 6 + KEEP(IXSZ)
      IBC_SOURCE = MYID
      DO ISLAVE = 1, NSLAVES
              CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB,
     &                NSLAVES, 
     &                NBLIG, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
        IERR = -1
        DO WHILE (IERR .EQ.-1)
         IF ( KEEP(50) .eq. 0 ) THEN
           NBCOL =  NFRONT
           CALL CMUMPS_BUF_SEND_DESC_BANDE( INODE,
     &      sum(SONROWS_PER_ROW(FIRST_INDEX:FIRST_INDEX+NBLIG-1)),
     &      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     &      IZERO, IDUMMY,
     &      IW(PDEST), IBC_SOURCE, NFRONT, COMM, IERR
     &      )
         ELSE
           NBCOL = NASS1+SHIFT_INDEX+NBLIG
           CALL CMUMPS_BUF_SEND_DESC_BANDE( INODE,
     &      sum(SONROWS_PER_ROW(FIRST_INDEX:FIRST_INDEX+NBLIG-1)),
     &      NBLIG, IW(PTRROW), NBCOL, IW(PTRCOL), NASS1,
     &      NSLAVES-ISLAVE,
     &      IW( PTLUST(STEP(INODE))+6+KEEP(IXSZ)+ISLAVE),
     &      IW(PDEST), IBC_SOURCE, NFRONT, COMM, IERR
     &      )
         ENDIF
         IF (IERR.EQ.-1) THEN
          BLOCKING  = .FALSE.
          SET_IRECV = .TRUE.
          MESSAGE_RECEIVED = .FALSE.
          CALL CMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &     BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &     MPI_ANY_SOURCE, MPI_ANY_TAG,
     &     STATUS, BUFR, LBUFR,
     &     LBUFR_BYTES,
     &     PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     &     LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &     PTLUST, PTRFAC,
     &     PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, INFO(1),
     &     INFO(2), COMM,
     &     NBPROCFILS,
     &     IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &     root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &     FILS, PTRARW, PTRAIW,
     &     INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &     NELT+1, NELT, FRT_PTR, FRT_ELT,
     &     ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &       )
          IF ( INFO(1) .LT. 0 ) GOTO 500
          IF (MESSAGE_RECEIVED) THEN
           IOLDPS = PTLUST(STEP(INODE))
           PTRCOL = IOLDPS + HF + NFRONT
           PTRROW = IOLDPS + HF + NASS1 + SHIFT_INDEX
          ENDIF
         ENDIF
        ENDDO
        IF (IERR .EQ. -2) GOTO 300
        IF (IERR .EQ. -3) GOTO 305
        PTRROW = PTRROW + NBLIG
        PDEST  = PDEST + 1
      ENDDO
      DEALLOCATE(SONROWS_PER_ROW)
      IF (NUMSTK.EQ.0) GOTO 500
      ISON = IFSON
      DO IELL = 1, NUMSTK
        ISTCHK = PIMASTER(STEP(ISON))
        NELIM = IW(ISTCHK + 1 + KEEP(IXSZ))
        LSTK    = IW(ISTCHK + KEEP(IXSZ))
        NPIVS   = IW(ISTCHK + 3 + KEEP(IXSZ))
        IF ( NPIVS .LT. 0 ) NPIVS = 0
        NSLSON  = IW(ISTCHK + 5 + KEEP(IXSZ))
        HS      = 6 + NSLSON + KEEP(IXSZ)
        NCOLS     = NPIVS + LSTK
        SAME_PROC     = (ISTCHK.LT.IWPOS)
        IF (.NOT.SAME_PROC) THEN
         NROWS = IW(ISTCHK + 2 + KEEP(IXSZ) )
        ELSE
         NROWS = NCOLS
        ENDIF
        PDEST   = ISTCHK + 6 + KEEP(IXSZ)
        NCBSON  = LSTK - NELIM
        PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
        IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
           NFS4FATHER = NCBSON
           DO I=0,NCBSON-1
              IF(IW(PTRCOL+I) .GT. NASS1) THEN
                 NFS4FATHER = I
                 EXIT
              ENDIF
           ENDDO
           NFS4FATHER = NFS4FATHER + NELIM
        ELSE
          NFS4FATHER = 0
        ENDIF
        IF (NSLSON.EQ.0) THEN
          NSLSON = 1
          PDEST1(1)  = MUMPS_PROCNODE(PROCNODE_STEPS(STEP(ISON)),
     &                                SLAVEF)
          IF (PDEST1(1).EQ.MYID) THEN
            CALL CMUMPS_MAPLIG_FILS_NIV1( COMM_LOAD, ASS_IRECV, 
     &      BUFR, LBUFR, LBUFR_BYTES,
     &      INODE, ISON, NSLAVES, 
     &      IW( PTLUST(STEP(INODE)) + 6 +KEEP(IXSZ)),
     &      NFRONT, NASS1, NFS4FATHER, NCBSON, IW( PTRCOL ),
     &      PROCNODE_STEPS,
     &      SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &      LRLUS, N, IW, LIW, A, LA,
     &      PTRIST, PTLUST, PTRFAC, PTRAST, STEP,
     &      PIMASTER, PAMASTER, NSTK_S, COMP,
     &      INFO(1), INFO(2), MYID, COMM, NBPROCFILS,
     &      IPOOL, LPOOL, LEAF,
     &      NBFIN, ICNTL, KEEP, KEEP8, DKEEP, root,
     &      OPASSW, OPELIW,
     &      ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW, INTARR, DBLARR,
     &      ND, FRERE, NELT+1, NELT, 
     &      FRT_PTR, FRT_ELT, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &      )
           IF ( INFO(1) .LT. 0 ) GOTO 500
          ELSE
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            CALL  CMUMPS_BUF_SEND_MAPLIG( 
     &           INODE, NFRONT,NASS1,NFS4FATHER,
     &           ISON, MYID,
     &      NSLAVES, IW( PTLUST(STEP(INODE)) + 6 +KEEP(IXSZ)),
     &      IW(PTRCOL), NCBSON,
     &      COMM, IERR, PDEST1, NSLSON, SLAVEF, 
     &      KEEP,KEEP8, STEP, N, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL CMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &        MPI_ANY_SOURCE, MPI_ANY_TAG,
     &        STATUS, BUFR, LBUFR, LBUFR_BYTES,
     &        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     &        LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &        PTLUST, PTRFAC,
     &        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, INFO(1),
     &        INFO(2), COMM,
     &        NBPROCFILS,
     &        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &        root,OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS,
     &        PTRARW, PTRAIW,
     &        INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &        NELT+1, NELT, FRT_PTR, FRT_ELT, 
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &        )
              IF ( INFO(1) .LT. 0 ) GOTO 500
            ENDIF
           ENDDO
           IF (IERR .EQ. -2) GOTO 290
           IF (IERR .EQ. -3) GOTO 295
          ENDIF
        ELSE
          IF (PIMASTER(STEP(ISON)).GT.0) THEN
          IERR = -1
          DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6 + KEEP(IXSZ)
            CALL  CMUMPS_BUF_SEND_MAPLIG( 
     &           INODE, NFRONT, NASS1, NFS4FATHER,
     &           ISON, MYID,
     &      NSLAVES, IW(PTLUST(STEP(INODE))+6+KEEP(IXSZ)),
     &      IW(PTRCOL), NCBSON,
     &      COMM, IERR, IW(PDEST), NSLSON, SLAVEF, 
     &      KEEP,KEEP8, STEP, N, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL CMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &        BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &        MPI_ANY_SOURCE, MPI_ANY_TAG,
     &        STATUS, BUFR, LBUFR,
     &        LBUFR_BYTES,
     &        PROCNODE_STEPS, POSFAC, IWPOS, IWPOSCB, IPTRLU,
     &        LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &        PTLUST, PTRFAC,
     &        PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP, INFO(1),
     &        INFO(2), COMM,
     &        NBPROCFILS,
     &        IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &        root,OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &        FILS, PTRARW, PTRAIW,
     &        INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &        NELT+1, NELT, FRT_PTR, FRT_ELT,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &        )
             IF ( INFO(1) .LT. 0 ) GOTO 500
            ENDIF
          ENDDO
          IF (IERR .EQ. -2) GOTO 290
          IF (IERR .EQ. -3) GOTO 295
          ENDIF
          DO ISLAVE = 0, NSLSON-1
            IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
               CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP,KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON,
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX        = PTRCOL + SHIFT_INDEX
              CALL CMUMPS_MAPLIG( COMM_LOAD, ASS_IRECV, 
     &        BUFR, LBUFR, LBUFR_BYTES,
     &        INODE, ISON, NSLAVES, 
     &        IW( PTLUST(STEP(INODE))+6+KEEP(IXSZ)),
     &        NFRONT, NASS1,NFS4FATHER,
     &        TROW_SIZE, IW( INDX ),
     &        PROCNODE_STEPS,
     &        SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &        LRLUS, N, IW,
     &        LIW, A, LA,
     &        PTRIST, PTLUST, PTRFAC, PTRAST, STEP,
     &        PIMASTER, PAMASTER, NSTK_S, COMP, INFO(1), INFO(2),
     &        MYID, COMM, NBPROCFILS, IPOOL, LPOOL, LEAF,
     &        NBFIN, ICNTL, KEEP,KEEP8,DKEEP, root,
     &        OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, PTRARW, PTRAIW,
     &        INTARR, DBLARR, ND, FRERE,
     &        NELT+1, NELT, FRT_PTR, FRT_ELT, 
     & 
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &        )
              IF ( INFO(1) .LT. 0 ) GOTO 500
              EXIT
            ENDIF
          ENDDO
        ENDIF
       ISON = FRERE(STEP(ISON))
      ENDDO
      GOTO 500
  250 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE IN INTEGER DYNAMIC ALLOCATION DURING
     & CMUMPS_FAC_ASM_NIV2_ELT'
      ENDIF
      INFO(1)   = -13
      INFO(2)   = NUMSTK + 1
      GOTO 490
  245 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * ) ' FAILURE ALLOCATING COPY_CAND',
     &                 ' DURING CMUMPS_FAC_ASM_NIV2_ELT'
      ENDIF
      INFO(1)  = -13
      INFO(2)  = SLAVEF+1
      GOTO 490
  265 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * ) ' FAILURE ALLOCATING TMP_SLAVES_LIST',
     &                 ' DURING CMUMPS_FAC_ASM_NIV2_ELT'
      ENDIF
      INFO(1)   = -13
      INFO(2)   = SIZE_TMP_SLAVES_LIST
      GOTO 490
  270 CONTINUE
      INFO(1) = -8
      INFO(2) = LREQ
      IF (LPOK) THEN
        WRITE( LP, * )
     &  ' FAILURE IN INTEGER ALLOCATION DURING CMUMPS_ASM_NIV2_ELT'
      ENDIF
      GOTO 490
  275 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * ) ' FAILURE ALLOCATING SONROWS_PER_ROW',
     &                 ' DURING CMUMPS_ASM_NIV2_ELT'
      ENDIF
      INFO(1)  = -13
      INFO(2)  = NFRONT-NASS1
      GOTO 490
  280 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &  ' FAILURE, WORKSPACE TOO SMALL DURING CMUMPS_ASM_NIV2_ELT'
      ENDIF
      INFO(1) = -9
      CALL MUMPS_SET_IERROR(LAELL8-LRLUS, INFO(2))
      GOTO 490
  290 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE, SEND BUFFER TOO SMALL (1) DURING CMUMPS_ASM_NIV2_ELT'
      ENDIF
      INFO(1) = -17
      LREQ = NCBSON + 6 + NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE, RECV BUFFER TOO SMALL (1) DURING CMUMPS_ASM_NIV2_ELT'
      ENDIF
      INFO(1) = -20
      LREQ = NCBSON + 6 + NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE, SEND BUFFER TOO SMALL (2)',
     &' DURING CMUMPS_FAC_ASM_NIV2_ELT'
      ENDIF
      INFO(1) = -17
      LREQ = NBLIG + NBCOL + 4 + KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  305 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE, RECV BUFFER TOO SMALL (2)',
     &' DURING CMUMPS_FAC_ASM_NIV2_ELT'
      ENDIF
      INFO(1) = -17
      LREQ = NBLIG + NBCOL + 4 + KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  490 CALL CMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_FAC_ASM_NIV2_ELT
      END MODULE CMUMPS_FAC_ASM_MASTER_ELT_M
