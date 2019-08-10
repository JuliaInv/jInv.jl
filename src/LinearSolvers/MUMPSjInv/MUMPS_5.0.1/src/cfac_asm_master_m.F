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
      MODULE CMUMPS_FAC_ASM_MASTER_M
      CONTAINS
      SUBROUTINE CMUMPS_FAC_ASM_NIV1( COMM_LOAD, ASS_IRECV,
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
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, JOBASS, ETATASS
     &    )
      USE MUMPS_BUILD_SORT_INDEX_M
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
      INTEGER N,LIW,NSTEPS
      INTEGER(8) LA, LRLU, LRLUS, IPTRLU, POSFAC
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      REAL    DKEEP(130)
      INTEGER, INTENT(INOUT) :: INFO(2)
      INTEGER INODE,MAXFRW,
     &        IWPOSCB, COMP
      INTEGER, TARGET :: IWPOS
      INTEGER JOBASS,ETATASS 
      LOGICAL SON_LEVEL2
      COMPLEX, TARGET :: A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER COMM, NBFIN, SLAVEF, MYID
      INTEGER LPOOL, LEAF
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER IPOOL( LPOOL )
      INTEGER NBPROCFILS(KEEP(28)), NSTK_S(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER BUFR( LBUFR )
      INTEGER IDUMMY(1)
      INTEGER, PARAMETER :: LIDUMMY = 1
      INTEGER, TARGET :: IW(LIW)
      INTEGER ITLOC(N+KEEP(253)),
     &        PTRARW(N), PTRAIW(N), ND(KEEP(28)), PERM(N), 
     &        FILS(N), FRERE(KEEP(28)), DAD(KEEP(28)),
     &        PTRIST(KEEP(28)), PTLUST(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28))
      COMPLEX :: RHS_MUMPS(KEEP(255))
      INTEGER(8) :: PTRFAC(KEEP(28)), PTRAST(KEEP(28)),
     &              PAMASTER(KEEP(28))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER      INTARR(max(1,KEEP(14)))
      COMPLEX DBLARR(max(1,KEEP(13)))
      INTEGER  MUMPS_TYPENODE
      EXTERNAL MUMPS_TYPENODE
      INCLUDE 'mumps_headers.h'
      INTEGER LP, HS, HF
      LOGICAL LPOK
      INTEGER NBPANELS_L, NBPANELS_U
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER NFS4FATHER
      INTEGER(8) NFRONT8, LAELL8, LAELL_REQ8
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ
      INTEGER LREQ_OOC
      INTEGER(8) :: SIZFR
      INTEGER SIZFI, NCB
      INTEGER J1,J2
      INTEGER NCOLS, NROWS, LDA_SON
      INTEGER(8) :: JJ2, ICT13
#if defined(ALLOW_NON_INIT)
      INTEGER(8) :: NUMROWS, JJ3, JJ8, APOS_ini
#endif
      INTEGER NELIM,JJ,JJ1,J3,
     &        IBROT,IORG
      INTEGER JPOS,ICT11
      INTEGER JK,IJROW,NBCOL,NUMORG,IOLDPS,J4
      INTEGER(8) IACHK, POSELT, LAPOS2, IACHK_ini
      INTEGER(8) APOS, APOS2, APOS3, POSEL1, ICT12
      INTEGER AINPUT
      INTEGER NSLAVES, NSLSON, NPIVS,NPIV_ANA,NPIV
      INTEGER PTRCOL, ISLAVE, PDEST,LEVEL
      INTEGER ISON_IN_PLACE 
      INTEGER ISON_TOP 
      INTEGER(8) SIZE_ISON_TOP8
      LOGICAL RESET_TO_ZERO, RISK_OF_SAME_POS,
     &        RISK_OF_SAME_POS_THIS_LINE, OMP_PARALLEL_FLAG
      LOGICAL LEVEL1, NIV1
      INTEGER TROW_SIZE
      INTEGER INDX, FIRST_INDEX, SHIFT_INDEX
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER, POINTER :: SON_IWPOS
      INTEGER, POINTER, DIMENSION(:) :: SON_IW
      COMPLEX, POINTER, DIMENSION(:) :: SON_A
      INTEGER NCBSON
      LOGICAL SAME_PROC
      INTRINSIC real
      COMPLEX ZERO
      PARAMETER( ZERO = (0.0E0,0.0E0) )
      INTEGER NELT, LPTRAR
      EXTERNAL MUMPS_INSSARBR
      LOGICAL MUMPS_INSSARBR
      LOGICAL SSARBR
      LOGICAL COMPRESSCB
      INTEGER(8) :: LCB
      DOUBLE PRECISION FLOP1,FLOP1_EFF
      EXTERNAL MUMPS_IN_OR_ROOT_SSARBR
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR
      LP      = ICNTL(1)
      LPOK    = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      COMPRESSCB =.FALSE.
      NELT       = 1
      LPTRAR     = N
      NFS4FATHER = -1
      IN         = INODE
      NBPROCFILS(STEP(IN)) = 0
      LEVEL = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      IF (LEVEL.NE.1) THEN
       WRITE(*,*) 'INTERNAL ERROR 1 in CMUMPS_FAC_ASM_NIV1 '
       CALL MUMPS_ABORT()
      END IF
      NSLAVES = 0
      HF = 6 + NSLAVES + KEEP(IXSZ)
      IF (JOBASS.EQ.0) THEN
        ETATASS= 0 
      ELSE
        ETATASS= 2 
        IOLDPS = PTLUST(STEP(INODE)) 
        NFRONT = IW(IOLDPS + KEEP(IXSZ)) 
        NASS1  = iabs(IW(IOLDPS + 2 + KEEP(IXSZ)))
        ICT11 = IOLDPS + HF - 1 + NFRONT
        SSARBR=MUMPS_INSSARBR(PROCNODE_STEPS(STEP(INODE)),
     &                        SLAVEF)
        NUMORG = 0
        DO WHILE (IN.GT.0)
          NUMORG = NUMORG + 1
          IN = FILS(IN)
        ENDDO
        NUMSTK = 0
        IFSON = -IN
        ISON = IFSON
        IF (ISON .NE. 0) THEN
         DO WHILE (ISON .GT. 0)
           NUMSTK = NUMSTK + 1
           ISON = FRERE(STEP(ISON))
         ENDDO
        ENDIF
        GOTO 123
      ENDIF
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
           IF (LPOK) THEN
            WRITE(LP, * ) 'INTERNAL ERROR 2 after compress '
            WRITE(LP, * ) 'IN CMUMPS_FAC_ASM_NIV1 '
            WRITE(LP, * ) 'LRLU,LRLUS=', LRLU,LRLUS
           ENDIF
           GOTO 270
          END IF
          IF ((IWPOS + LREQ -1) .GT. IWPOSCB) GOTO 270
      END IF
      IOLDPS = IWPOS
      IWPOS = IWPOS + LREQ
      ISON_TOP      = -9999
      ISON_IN_PLACE = -9999
      SIZE_ISON_TOP8 = 0_8
      IF (KEEP(234).NE.0) THEN
        IF ( IWPOSCB .NE. LIW ) THEN 
        IF ( IWPOSCB+IW(IWPOSCB+1+XXI).NE.LIW) THEN
          ISON = IW( IWPOSCB + 1 + XXN )
          IF ( DAD( STEP( ISON ) ) .EQ. INODE .AND.
     &    MUMPS_TYPENODE(PROCNODE_STEPS(STEP(ISON)),SLAVEF)
     &    .EQ. 1 )
     &    THEN
            ISON_TOP = ISON
            CALL MUMPS_GETI8(SIZE_ISON_TOP8,IW(IWPOSCB + 1 + XXR))
            IF (LRLU .LT. int(NFRONT,8) * int(NFRONT,8)) THEN
              ISON_IN_PLACE = ISON
            ENDIF
          END IF
        END IF
        END IF
      END IF
      NIV1 = .TRUE.
        CALL MUMPS_BUILD_SORT_INDEX( MYID, INODE, N, IOLDPS, HF,
     &        NFRONT, NFRONT_EFF, PERM, DAD,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     &        IFSON, STEP, PIMASTER, PTRIST, PTRAIW, IW, LIW,
     &        INTARR, ITLOC, FILS, FRERE,
     &        SON_LEVEL2, NIV1, NBPROCFILS, KEEP, KEEP8, INFO(1),
     &        ISON_IN_PLACE, 
     &        PROCNODE_STEPS, SLAVEF, IDUMMY, LIDUMMY
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
         WRITE(LP,*) ' INTERNAL ERROR 3 ',
     &               ' IN CMUMPS_FAC_ASM_NIV1 ',
     &               ' NFRONT, NFRONT_EFF = ',
     &                 NFRONT, NFRONT_EFF
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
      NCB   = NFRONT - NASS1
      MAXFRW = max0(MAXFRW, NFRONT)
      ICT11 = IOLDPS + HF - 1 + NFRONT
      LAELL8 = NFRONT8 * NFRONT8
      LAELL_REQ8 = LAELL8
      IF ( ISON_IN_PLACE > 0 ) THEN
        LAELL_REQ8 = LAELL8 - SIZE_ISON_TOP8
      ENDIF
      IF (LRLU .LT. LAELL_REQ8) THEN
        IF (LRLUS .LT. LAELL_REQ8) THEN
          IF (LPOK) THEN
           WRITE(LP, * ) ' NOT ENOUGH MEMORY during ASSEMBLY ',
     &     ' MEMORY REQUESTED = ', LAELL_REQ8,
     &     ' AVAILABLE =', LRLUS
          ENDIF
          GOTO 280
        ELSE
          CALL CMUMPS_COMPRE_NEW
     &        (N, KEEP(28), IW, LIW, A, LA, LRLU, IPTRLU,
     &         IWPOS, IWPOSCB, PTRIST, PTRAST, STEP, PIMASTER,
     &         PAMASTER,KEEP(216),LRLUS,KEEP(IXSZ),COMP,DKEEP(97),
     &         MYID)
          IF (LRLU .NE. LRLUS) THEN
           IF (LPOK) THEN
            WRITE(LP, * ) 'INTERNAL ERROR 4 after compress '
            WRITE(LP, * ) 'IN CMUMPS_FAC_ASM_NIV1'
            WRITE(LP, * ) 'LRLU,LRLUS=',LRLU,LRLUS
           ENDIF
           GOTO 280
          END IF
        END IF
      END IF
      LRLU = LRLU - LAELL8 
      LRLUS = LRLUS - LAELL8 + SIZE_ISON_TOP8
      KEEP8(67) = min(LRLUS, KEEP8(67))
      POSELT = POSFAC
      POSFAC = POSFAC + LAELL8
      SSARBR=MUMPS_INSSARBR(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
      CALL CMUMPS_LOAD_MEM_UPDATE(SSARBR,.FALSE.,
     &     LA-LRLUS, 
     &     0_8,
     &     LAELL8-SIZE_ISON_TOP8, 
     &     KEEP,KEEP8,
     &     LRLUS)
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = min(POSELT + LAELL8 - 1_8, IPTRLU)
      A(POSELT:LAPOS2) = ZERO
#else
      IF ( KEEP(50) .eq. 0 .OR. NFRONT .LT. KEEP(63) ) THEN
        LAPOS2 = min(POSELT + LAELL8 - 1_8, IPTRLU)
!$OMP PARALLEL DO SCHEDULE(STATIC, 3000)  
        DO JJ8 = POSELT, LAPOS2
           A( JJ8 ) = ZERO
        ENDDO
!$OMP END PARALLEL DO
      ELSE
        IF (ETATASS.EQ.1) THEN
         APOS_ini = POSELT
!$OMP PARALLEL DO PRIVATE(APOS, JJ3)
!$OMP& IF (NFRONT8 - 1_8 > 300_8)
         DO JJ8 = 0_8, NFRONT8 - 1_8
          JJ3 = min(JJ8,int(NASS1-1,8)) 
          APOS = APOS_ini + JJ8 * NFRONT8
          A(APOS:APOS+JJ3) = ZERO
         END DO
!$OMP END PARALLEL DO
        ELSE
          APOS_ini = POSELT
          NUMROWS = min(NFRONT8, (IPTRLU-APOS_ini) / NFRONT8 )
!$OMP PARALLEL DO PRIVATE(APOS)
!$OMP& IF (NUMROWS - 1_8 > 300_8)
          DO JJ8 = 0_8, NUMROWS - 1_8
             APOS = APOS_ini + JJ8 * NFRONT8
             A(APOS:APOS + JJ8) = ZERO
          ENDDO
!$OMP END PARALLEL DO
          IF( NUMROWS .LT. NFRONT8 ) THEN
            APOS = APOS_ini + NFRONT8*NUMROWS
            A(APOS : min(IPTRLU,APOS+NUMROWS)) = ZERO
          ENDIF
        ENDIF
      END IF
#endif
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
 123  CONTINUE  
      IF (NUMSTK.NE.0) THEN
        IF (ISON_TOP > 0) THEN
          ISON = ISON_TOP
        ELSE
          ISON = IFSON
        ENDIF
        DO 220 IELL = 1, NUMSTK
          ISTCHK    = PIMASTER(STEP(ISON))
          SON_IW    => IW
          SON_IWPOS => IWPOS
          SON_A     => A
          LSTK      = SON_IW(ISTCHK + KEEP(IXSZ))
          NELIM     = SON_IW(ISTCHK + KEEP(IXSZ) + 1)
          NPIVS     = SON_IW(ISTCHK + KEEP(IXSZ) + 3)
          IF ( NPIVS .LT. 0 ) NPIVS = 0
          NSLSON    = SON_IW(ISTCHK + KEEP(IXSZ) + 5)
          HS        = 6 + KEEP(IXSZ) + NSLSON 
          NCOLS     = NPIVS + LSTK
          SAME_PROC     = (ISTCHK.LT.SON_IWPOS)
          IF ( SAME_PROC ) THEN
            COMPRESSCB=( SON_IW(PTRIST(STEP(ISON))+XXS) .EQ. S_CB1COMP )
          ELSE
            COMPRESSCB=( SON_IW(ISTCHK + XXS) .EQ. S_CB1COMP )
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
           SIZFR  = int(LSTK,8)*int(LSTK,8)
           IF (COMPRESSCB) SIZFR = (int(LSTK,8)*int(LSTK+1,8))/2_8
          ELSE
           IF ( KEEP(50).eq.0 ) THEN
             SIZFR = int(NELIM,8) * int(LSTK,8)
           ELSE
             SIZFR = int(NELIM,8) * int(NELIM,8)
           END IF
           J2 = J1 + NELIM - 1
          ENDIF
          IF (JOBASS.EQ.0) OPASSW = OPASSW + dble(SIZFR)
          IACHK = PAMASTER(STEP(ISON))
          IF ( KEEP(50) .eq. 0 ) THEN
            POSEL1 = PTRAST(STEP(INODE)) - NFRONT8
            IF (NFRONT .EQ. LSTK.AND. ISON.EQ.ISON_IN_PLACE
     &          .AND.IACHK + SIZFR - 1_8 .EQ. POSFAC - 1_8 ) THEN
               GOTO 205
            ENDIF
            IF (J2.GE.J1) THEN
              RESET_TO_ZERO = (IACHK .LT. POSFAC .AND.
     &                         ISON.EQ.ISON_IN_PLACE)
              RISK_OF_SAME_POS = IACHK + SIZFR - 1_8 .EQ. POSFAC - 1_8
     &        .AND. ISON.EQ.ISON_IN_PLACE
              RISK_OF_SAME_POS_THIS_LINE = .FALSE.
              IACHK_ini = IACHK
              OMP_PARALLEL_FLAG = (RESET_TO_ZERO.EQV..FALSE.).AND.
     &            ((J2-J1).GT.300)
!$OMP PARALLEL IF(OMP_PARALLEL_FLAG) PRIVATE(APOS, JJ1, JJ2,IACHK)
!$OMP& FIRSTPRIVATE(RISK_OF_SAME_POS_THIS_LINE,RESET_TO_ZERO) 
!$OMP DO
              DO 170 JJ = J1, J2
                APOS = POSEL1 + int(SON_IW(JJ),8) * int(NFRONT,8)
                IACHK = IACHK_ini + int(JJ-J1,8)*int(LSTK,8)
                IF (RESET_TO_ZERO) THEN
                  IF (RISK_OF_SAME_POS) THEN
                    IF (JJ.EQ.J2) THEN
                      RISK_OF_SAME_POS_THIS_LINE =
     &                  (ISON .EQ. ISON_IN_PLACE)
     &                  .AND. ( APOS + int(SON_IW(J1+LSTK-1)-1,8).EQ.
     &                          IACHK+int(LSTK-1,8) )
                    ENDIF
                  ENDIF
                  IF ((IACHK .GE. POSFAC).AND.(JJ>J1))THEN
                   RESET_TO_ZERO =.FALSE.
                  ENDIF
                  IF (RISK_OF_SAME_POS_THIS_LINE) THEN
                    DO JJ1 = 1, LSTK
                      JJ2 = APOS + int(SON_IW(J1 + JJ1 - 1) - 1,8)
                      IF ( IACHK+int(JJ1-1,8) .NE. JJ2 ) THEN
                        A(JJ2) = A(IACHK + int(JJ1 - 1,8))
                        A(IACHK + int(JJ1 -1,8)) = ZERO
                      ENDIF
                    ENDDO
                  ELSE
                    DO JJ1 = 1, LSTK
                      JJ2 = APOS + int(SON_IW(J1+JJ1-1),8) - 1_8
                      A(JJ2) = A(IACHK + int(JJ1 - 1,8))
                      A(IACHK + int(JJ1 -1,8)) = ZERO
                    ENDDO
                  ENDIF
                ELSE 
                  DO JJ1 = 1, LSTK
                    JJ2 = APOS + int(SON_IW(J1+JJ1-1),8) - 1_8
                    A(JJ2) = A(JJ2) + SON_A(IACHK + int(JJ1 - 1,8))
                  ENDDO
                ENDIF
  170         CONTINUE
!$OMP END DO
!$OMP END PARALLEL 
            END IF
          ELSE
            IF (LEVEL1) THEN
             LDA_SON = LSTK
            ELSE
             LDA_SON = NELIM
            ENDIF
            IF (COMPRESSCB) THEN
              LCB = SIZFR
            ELSE
              LCB = int(LDA_SON,8)* int(J2-J1+1,8)
            ENDIF
            IF (ISON .EQ. ISON_IN_PLACE) THEN
              CALL CMUMPS_LDLT_ASM_NIV12_IP(A, LA,
     &           PTRAST(STEP( INODE )), NFRONT, NASS1,
     &           IACHK, LDA_SON, LCB,
     &           SON_IW( J1 ), J2 - J1 + 1, NELIM, ETATASS, 
     &           COMPRESSCB)
            ELSE
              IF (LCB .GT. 0) THEN
                CALL CMUMPS_LDLT_ASM_NIV12(A, LA, SON_A(IACHK),
     &           PTRAST(STEP( INODE )), NFRONT, NASS1,
     &           LDA_SON, LCB,
     &           SON_IW( J1 ), J2 - J1 + 1, NELIM, ETATASS, 
     &           COMPRESSCB
     &          )
              ENDIF
            ENDIF
          ENDIF
  205     IF (LEVEL1) THEN 
           IF (SAME_PROC) ISTCHK = PTRIST(STEP(ISON))
           IF ((SAME_PROC).AND.ETATASS.NE.1) THEN
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
           IF (ETATASS.NE.1) THEN
             IF ( SAME_PROC ) THEN 
               PTRIST(STEP(ISON))   = -99999999
             ELSE
               PIMASTER(STEP( ISON )) = -99999999
             ENDIF
               CALL CMUMPS_FREE_BLOCK_CB(SSARBR, MYID, N, ISTCHK,
     &          PAMASTER(STEP(ISON)),
     &          IW, LIW, LRLU, LRLUS, IPTRLU,
     &          IWPOSCB, LA, KEEP,KEEP8,
     &          (ISON .EQ. ISON_TOP)
     &          )
           ENDIF
          ELSE
           PDEST = ISTCHK + 6 + KEEP(IXSZ)
           NCBSON  = LSTK - NELIM
           PTRCOL   = ISTCHK +  HS + NROWS + NPIVS + NELIM
           DO ISLAVE = 0, NSLSON-1
             IF (IW(PDEST+ISLAVE).EQ.MYID) THEN
              CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP, KEEP8, ISON, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE+1, NCBSON, 
     &                NSLSON, 
     &                TROW_SIZE, FIRST_INDEX  )
              SHIFT_INDEX = FIRST_INDEX - 1
              INDX = PTRCOL + SHIFT_INDEX
              CALL CMUMPS_MAPLIG( COMM_LOAD, ASS_IRECV, 
     &             BUFR, LBUFR, LBUFR_BYTES,
     &             INODE, ISON, NSLAVES, IDUMMY,
     &             NFRONT, NASS1,NFS4FATHER,
     &             TROW_SIZE, IW( INDX ),
     &         PROCNODE_STEPS,
     &         SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &         LRLUS, N, IW,
     &         LIW, A, LA,
     &         PTRIST, PTLUST, PTRFAC, PTRAST, STEP,
     &         PIMASTER, PAMASTER, NSTK_S, COMP,
     &         INFO(1), INFO(2), MYID, COMM, NBPROCFILS, IPOOL, LPOOL,
     &         LEAF, NBFIN, ICNTL, KEEP, KEEP8,DKEEP,  root,
     &         OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, PTRARW, PTRAIW,
     &         INTARR, DBLARR, ND, FRERE,
     &         LPTRAR, NELT, IW, IW, 
     &
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE 
     &         )
              IF ( INFO(1) .LT. 0 ) GOTO 500
              EXIT
             ENDIF
           ENDDO
           IF (PIMASTER(STEP(ISON)).GT.0) THEN
           IERR = -1
           DO WHILE (IERR.EQ.-1)
            PTRCOL = PIMASTER(STEP(ISON)) + HS + NROWS + NPIVS + NELIM
            PDEST  = PIMASTER(STEP(ISON)) + 6 + KEEP(IXSZ)
            CALL  CMUMPS_BUF_SEND_MAPLIG( 
     &           INODE, NFRONT, NASS1, NFS4FATHER, 
     &           ISON, MYID,
     &       IZERO, IDUMMY, IW(PTRCOL), NCBSON,
     &       COMM, IERR, IW(PDEST), NSLSON, SLAVEF, 
     &       KEEP, KEEP8, STEP, N, 
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &        )
            IF (IERR.EQ.-1) THEN
             BLOCKING  = .FALSE.
             SET_IRECV = .TRUE.
             MESSAGE_RECEIVED = .FALSE.
             CALL CMUMPS_TRY_RECVTREAT( 
     &         COMM_LOAD, ASS_IRECV,
     &         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &         MPI_ANY_SOURCE, MPI_ANY_TAG,
     &         STATUS,
     &         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &         IWPOS, IWPOSCB, IPTRLU,
     &         LRLU, LRLUS, N, IW, LIW, A, LA,
     &         PTRIST, PTLUST, PTRFAC,
     &         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &         INFO(1), INFO(2), COMM,
     &         NBPROCFILS,
     &         IPOOL, LPOOL, LEAF,
     &         NBFIN, MYID, SLAVEF,
     &         root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, PTRARW, PTRAIW,
     &         INTARR, DBLARR, ICNTL, KEEP, KEEP8,DKEEP, ND, FRERE,
     &         LPTRAR, NELT, IW, IW,
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &          )
               IF ( INFO(1) .LT. 0 ) GOTO 500
            ENDIF
           ENDDO
           IF (IERR .EQ. -2) GOTO 290
           IF (IERR .EQ. -3) GOTO 295
           ENDIF
          ENDIF
        ISON = FRERE(STEP(ISON))
        IF (ISON .LE. 0) THEN
          ISON = IFSON
        ENDIF
  220 CONTINUE
      END IF
      IF (ETATASS.EQ.2) GOTO 500
      POSELT = PTRAST(STEP(INODE))
      IBROT = INODE
      DO 260 IORG = 1, NUMORG
        JK = PTRAIW(IBROT)
        AINPUT = PTRARW(IBROT)
        JJ = JK + 1
        J1 = JJ + 1
        J2 = J1 + INTARR(JK)
        J3 = J2 + 1
        J4 = J2 - INTARR(JJ)
        IJROW = INTARR(J1)
        ICT12 = POSELT + int(IJROW - NFRONT - 1,8)
        DO 240 JJ = J1, J2
           APOS2 = ICT12 + int(INTARR(JJ),8) * NFRONT8
           A(APOS2) = A(APOS2) + DBLARR(AINPUT)
          AINPUT = AINPUT + 1
  240   CONTINUE
        IF (J3 .LE. J4) THEN
          ICT13 = POSELT + int(IJROW - 1,8) * NFRONT8
          NBCOL = J4 - J3 + 1
          DO 250 JJ = 1, NBCOL
            APOS3 = ICT13 + int(INTARR(J3 + JJ - 1) - 1,8)
            A(APOS3) = A(APOS3) + DBLARR(AINPUT + JJ - 1)
  250     CONTINUE
        ENDIF
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
  260 CONTINUE
      GOTO 500
  270 CONTINUE
      INFO(1) = -8
      INFO(2) = LREQ
      IF (LPOK)  THEN
        WRITE( LP, * )
     &' FAILURE IN INTEGER ALLOCATION DURING CMUMPS_FAC_ASM'
      ENDIF
      GOTO 490
  280 CONTINUE
      INFO(1) = -9
      CALL MUMPS_SET_IERROR(LAELL_REQ8 - LRLUS, INFO(2))
      IF (LPOK)  THEN
        WRITE( LP, * )
     &' FAILURE, WORKSPACE TOO SMALL DURING CMUMPS_FAC_ASM'
      ENDIF
      GOTO 490
  290 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &  ' FAILURE, SEND BUFFER TOO SMALL DURING CMUMPS_FAC_ASM'
      ENDIF
      INFO(1) = -17
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &  ' FAILURE, SEND BUFFER TOO SMALL DURING CMUMPS_FAC_ASM'
      ENDIF
      INFO(1) = -17
      LREQ = NCBSON + 6+NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     & ' FAILURE IN INTEGER DYNAMIC ALLOCATION DURING CMUMPS_FAC_ASM'
      ENDIF
      INFO(1)   = -13
      INFO(2)  = NUMSTK + 1
  490 CALL  CMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_FAC_ASM_NIV1
      SUBROUTINE CMUMPS_FAC_ASM_NIV2(COMM_LOAD, ASS_IRECV,
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
      USE MUMPS_BUILD_SORT_INDEX_M
      USE CMUMPS_COMM_BUFFER
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      INCLUDE 'mpif.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N,LIW,NSTEPS, NBFIN
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      REAL       DKEEP(130)
      INTEGER(8) :: LRLUS, LRLU, IPTRLU, POSFAC, LA
      INTEGER INFO(2),INODE,MAXFRW,
     &        LPOOL, LEAF, IWPOS, IWPOSCB, COMP
      COMPLEX A(LA)
      DOUBLE PRECISION  OPASSW, OPELIW
      INTEGER COMM, SLAVEF, MYID,  LBUFR, LBUFR_BYTES
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER IPOOL(LPOOL)
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER IW(LIW), ITLOC(N+KEEP(253)),
     &        PTRARW(N), PTRAIW(N), ND(KEEP(28)),
     &        FILS(N), FRERE(KEEP(28)), DAD (KEEP(28)),
     &        PTRIST(KEEP(28)), PTLUST(KEEP(28)),
     &        STEP(N), 
     & PIMASTER(KEEP(28)),
     &        NSTK_S(KEEP(28)), PERM(N)
      COMPLEX :: RHS_MUMPS(KEEP(255))
      INTEGER   CAND(SLAVEF+1, max(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER NBPROCFILS(KEEP(28)),
     &        PROCNODE_STEPS(KEEP(28)), BUFR(LBUFR)
      INTEGER      INTARR(max(1,KEEP(14)))
      COMPLEX DBLARR(max(1,KEEP(13)))
      INTEGER LP, HS, HF, HF_OLD, NCBSON, NSLAVES_OLD,
     &        NBSPLIT
      LOGICAL LPOK
      INTEGER IN,NUMSTK,NASS,ISON,IFSON,NASS1,IELL
      INTEGER :: IBC_SOURCE
      INTEGER :: MAXWASTEDPROCS
      PARAMETER (MAXWASTEDPROCS=1)
      INTEGER NFS4FATHER,I
      INTEGER NFRONT,NFRONT_EFF,ISTCHK,LSTK,LREQ
      INTEGER(8) :: LAELL8
      INTEGER LREQ_OOC
      LOGICAL COMPRESSCB
      INTEGER(8) :: LCB
      INTEGER NCB
      INTEGER J1,J2,J3,MP
      INTEGER(8) :: JJ8, LAPOS2, JJ2, JJ3
      INTEGER NELIM,JJ,JJ1,NPIVS,NCOLS,NROWS,
     &        IBROT,IORG
      INTEGER LDAFS, LDA_SON
      INTEGER JK,IJROW,NBCOL,NUMORG,IOLDPS,J4, NUMORG_SPLIT
      INTEGER(8) :: ICT13
      INTEGER(8) :: IACHK, APOS, APOS2, POSELT, ICT12, POSEL1
      INTEGER AINPUT
      INTEGER NSLAVES, NSLSON
      INTEGER NBLIG, PTRCOL, PTRROW, ISLAVE, PDEST
      INTEGER PDEST1(1)
      INTEGER TYPESPLIT
      INTEGER ISON_IN_PLACE 
      LOGICAL IS_ofType5or6, SPLIT_MAP_RESTART 
      INTEGER NMB_OF_CAND, NMB_OF_CAND_ORIG
      LOGICAL SAME_PROC, NIV1, SON_LEVEL2
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      INTEGER TROW_SIZE, INDX, FIRST_INDEX, SHIFT_INDEX
      INTEGER IZERO
      INTEGER IDUMMY(1)
      PARAMETER( IZERO = 0 )
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE, MUMPS_TYPESPLIT
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, MUMPS_TYPESPLIT
      COMPLEX ZERO
      REAL RZERO
      PARAMETER( RZERO = 0.0E0 )
      PARAMETER( ZERO = (0.0E0,0.0E0) )
      INTEGER NELT, LPTRAR, NCBSON_MAX
      logical :: force_cand
      INTEGER ETATASS
      INCLUDE 'mumps_headers.h'
      INTEGER (8) :: APOSMAX
      REAL  MAXARR
      INTEGER INIV2, SIZE_TMP_SLAVES_LIST, allocok, 
     &        NCB_SPLIT, SIZE_LIST_SPLIT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TMP_SLAVES_LIST, COPY_CAND
      INTEGER NBPANELS_L, NBPANELS_U
      INTEGER, ALLOCATABLE, DIMENSION(:) :: SONROWS_PER_ROW
      MP      = ICNTL(2)
      LP      = ICNTL(1)
      LPOK    = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      IS_ofType5or6    = .FALSE.
      COMPRESSCB = .FALSE.
      ETATASS    = 0  
      IN         = INODE
      NBPROCFILS(STEP(IN)) = 0
      NSTEPS = NSTEPS + 1
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
      NELT = 1
      LPTRAR = 1
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
     &        KEEP(216),LRLUS,KEEP(IXSZ),COMP,DKEEP(97), MYID)
          IF (LRLU .NE. LRLUS) THEN
           IF (LPOK) THEN
            WRITE(LP, * ) 'PB compress CMUMPS_FAC_ASM_NIV2 ',
     &                    'LRLU,LRLUS=',LRLU,LRLUS
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
        ISON_IN_PLACE = -9999
        CALL MUMPS_BUILD_SORT_INDEX( MYID, INODE, N, IOLDPS, HF,
     &        NFRONT, NFRONT_EFF, PERM, DAD,
     &        NASS1, NASS, NUMSTK, NUMORG, IWPOSCB,
     &        IFSON, STEP, PIMASTER, PTRIST, PTRAIW, IW, LIW,
     &        INTARR, ITLOC, FILS, FRERE,
     &        SON_LEVEL2, NIV1, NBPROCFILS, KEEP,KEEP8, INFO(1),
     &        ISON_IN_PLACE,
     &        PROCNODE_STEPS, SLAVEF, SONROWS_PER_ROW,
     &        NFRONT-NASS1)
      IF (INFO(1).LT.0) GOTO 250
      IF ( NFRONT .NE. NFRONT_EFF ) THEN
        IF (
     &        (TYPESPLIT.EQ.5) .OR. (TYPESPLIT.EQ.6)) THEN
          WRITE(*,*) ' Internal error 1 in fac_ass due to splitting ',
     &     ' INODE, NFRONT, NFRONT_EFF =', INODE, NFRONT, NFRONT_EFF 
          WRITE(*,*) ' SPLITTING NOT YET READY FOR THAT'
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
          WRITE(LP,*) MYID,': INTERNAL ERROR 2 ',
     &     ' IN CMUMPS_FAC_ASM_NIV2 , INODE=', 
     &     INODE, ' NFRONT, NFRONT_EFF=', NFRONT, NFRONT_EFF
         ENDIF
         GOTO 270
        ENDIF
      ENDIF
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
      IW(IOLDPS+6+KEEP(IXSZ):IOLDPS+5+KEEP(IXSZ)+NSLAVES)=
     &             TMP_SLAVES_LIST(1:NSLAVES)
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
        LAELL8 = int(NASS1,8) * int(NFRONT,8)
        LDAFS = NFRONT
      ELSE
        LAELL8 = int(NASS1,8)*int(NASS1,8)
        IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2)
     &     LAELL8 = LAELL8+int(NASS1,8)
        LDAFS = NASS1
      ENDIF
      IF (LRLU .LT. LAELL8) THEN
        IF (LRLUS .LT. LAELL8) THEN
          GOTO 280
        ELSE
         CALL CMUMPS_COMPRE_NEW(N, KEEP(28),
     &      IW, LIW, A, LA,
     &      LRLU, IPTRLU,
     &      IWPOS, IWPOSCB, PTRIST, PTRAST,
     &      STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &      KEEP(IXSZ), COMP, DKEEP(97), MYID)
         IF (LRLU .NE. LRLUS) THEN
          IF (LPOK) THEN
           WRITE(LP, * ) 'INTERNAL ERROR 3 after compress '
           WRITE(LP, * ) 'IN CMUMPS_FAC_ASM_NIV2'
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
      CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE)),IW(IOLDPS+XXNBPR))
#endif
      CALL CMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,LA-LRLUS,0_8,LAELL8,
     &     KEEP,KEEP8,LRLUS)
      POSEL1 = POSELT - int(LDAFS,8)
#if ! defined(ALLOW_NON_INIT)
      LAPOS2 = POSELT + LAELL8 - 1_8
      A(POSELT:LAPOS2) = ZERO
#else
      IF ( KEEP(50) .eq. 0 .OR. LDAFS .lt. KEEP(63) ) THEN
        LAPOS2 = POSELT + LAELL8 - 1_8
        A(POSELT:LAPOS2) = ZERO
      ELSE
        APOS = POSELT
        DO JJ8 = 0_8, int(LDAFS-1,8)
          A(APOS:APOS+JJ8) = ZERO
          APOS = APOS + int(LDAFS,8)
        END DO
        IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
          A(APOS:APOS+int(LDAFS,8)-1_8)=ZERO
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
          NPIVS   = IW(ISTCHK + 3+KEEP(IXSZ))
          IF (NPIVS.LT.0) NPIVS=0
          NSLSON  = IW(ISTCHK + 5+KEEP(IXSZ))
          HS      = 6 + NSLSON  + KEEP(IXSZ)
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
             APOS = POSEL1 + int(IW(JJ),8) * int(LDAFS,8)
             DO 160 JJ1 = 1, LSTK
              JJ2 = APOS + int(IW(J1 + JJ1 - 1),8) - 1_8
              A(JJ2) = A(JJ2) + A(IACHK + int(JJ1 - 1,8))
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
      IBROT = INODE
      APOSMAX = POSELT + int(NASS1,8)*int(NASS1,8)
      DO 260 IORG = 1, NUMORG
        JK = PTRAIW(IBROT)
        AINPUT = PTRARW(IBROT)
        JJ = JK + 1
        J1 = JJ + 1
        J2 = J1 + INTARR(JK)
        J3 = J2 + 1
        J4 = J2 - INTARR(JJ)
        IJROW = INTARR(J1) 
        ICT12 = POSELT + int(IJROW - 1 - LDAFS, 8)
        MAXARR = RZERO
        DO 240 JJ = J1, J2
          IF (KEEP(219).NE.0) THEN
            IF (INTARR(JJ).LE.NASS1) THEN
              APOS2 = ICT12 + int(INTARR(JJ),8) * int(LDAFS,8)
              A(APOS2) = A(APOS2) + DBLARR(AINPUT)
            ELSEIF (KEEP(50).EQ.2) THEN
              MAXARR = max(MAXARR,abs(DBLARR(AINPUT)))
            ENDIF
          ELSE
            IF (INTARR(JJ).LE.NASS1) THEN
              APOS2 = ICT12 + int(INTARR(JJ),8) * int(LDAFS,8)
              A(APOS2) = A(APOS2) + DBLARR(AINPUT)
            ENDIF
          ENDIF
          AINPUT = AINPUT + 1
  240   CONTINUE
        IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2) THEN
           A(APOSMAX+int(IJROW-1,8)) = cmplx(MAXARR,kind=kind(A))
        ENDIF
        IF (J3 .GT. J4) GOTO 255
        ICT13 = POSELT + int(IJROW - 1,8) * int(LDAFS,8)
        NBCOL = J4 - J3 + 1
        DO JJ = 1, NBCOL
          JJ3 = ICT13 + int(INTARR(J3 + JJ - 1),8) - 1_8
          A(JJ3) = A(JJ3) + DBLARR(AINPUT + JJ - 1)
        ENDDO
  255   CONTINUE
        IF (KEEP(50).EQ.0) THEN
          DO JJ = 1, KEEP(253)
            APOS = POSELT +
     &             int(IJROW-1,8) * int(LDAFS,8) +
     &             int(LDAFS-KEEP(253)+JJ-1,8)
            A(APOS) = A(APOS) + RHS_MUMPS( (JJ-1) * KEEP(254) + IBROT )
          ENDDO
        ENDIF
        IBROT = FILS(IBROT)
  260 CONTINUE
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
     &     LPTRAR, NELT, IW, IW,
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
     &      ND, FRERE, LPTRAR, NELT, IW, IW,
     &
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
     &        LPTRAR,
     &        NELT, IW, IW, 
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
     &        LPTRAR, NELT, IW, IW, 
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
     &        INTARR, DBLARR, ND, FRERE, LPTRAR, NELT, IW,
     &        IW, 
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
     & CMUMPS_FAC_ASM_NIV2'
      ENDIF
      INFO(1)   = -13
      INFO(2)   = NUMSTK + 1
      GOTO 490
  245 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * ) ' FAILURE ALLOCATING COPY_CAND',
     &                 ' DURING CMUMPS_FAC_ASM_NIV2'
      ENDIF
      INFO(1)  = -13
      INFO(2)  = SLAVEF+1
      GOTO 490
  265 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * ) ' FAILURE ALLOCATING TMP_SLAVES_LIST',
     &                 ' DURING CMUMPS_FAC_ASM_NIV2'
      ENDIF
      INFO(1)   = -13
      INFO(2)   = SIZE_TMP_SLAVES_LIST
      GOTO 490
  270 CONTINUE
      INFO(1) = -8
      INFO(2) = LREQ
      IF (LPOK) THEN
        WRITE( LP, * )
     &  ' FAILURE IN INTEGER ALLOCATION DURING CMUMPS_FAC_ASM_NIV2'
      ENDIF
      GOTO 490
  275 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * ) ' FAILURE ALLOCATING SONROWS_PER_ROW',
     &                 ' DURING CMUMPS_FAC_ASM_NIV2'
      ENDIF
      INFO(1)  = -13
      INFO(2)  = NFRONT-NASS1
      GOTO 490
  280 CONTINUE
      IF (LPOK) THEN
        WRITE( LP, * )
     &' FAILURE, WORKSPACE TOO SMALL DURING CMUMPS_FAC_ASM_NIV2'
      ENDIF
      INFO(1) = -9
      CALL MUMPS_SET_IERROR(LAELL8-LRLUS, INFO(2))
      GOTO 490
  290 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, SEND BUFFER TOO SMALL (1) DURING CMUMPS_FAC_ASM_NIV2'
      ENDIF
      INFO(1) = -17
      LREQ = NCBSON + 6 + NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  295 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, RECV BUFFER TOO SMALL (1) DURING CMUMPS_FAC_ASM_NIV2'
      ENDIF
      INFO(1) = -20
      LREQ = NCBSON + 6 + NSLSON+KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  300 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, SEND BUFFER TOO SMALL (2)',
     &' DURING CMUMPS_FAC_ASM_NIV2'
      ENDIF
      INFO(1) = -17
      LREQ = NBLIG + NBCOL + 4 + KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  305 CONTINUE
      IF ((ICNTL(1) .GT. 0) .AND. (ICNTL(4) .GE. 1)) THEN
        LP = ICNTL(1)
        WRITE( LP, * )
     &' FAILURE, RECV BUFFER TOO SMALL (2)',
     &' DURING CMUMPS_FAC_ASM_NIV2'
      ENDIF
      INFO(1) = -17
      LREQ = NBLIG + NBCOL + 4 + KEEP(IXSZ)
      INFO(2) =  LREQ  * KEEP( 34 ) 
      GOTO 490
  490 CALL CMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
  500 CONTINUE
      RETURN
      END SUBROUTINE CMUMPS_FAC_ASM_NIV2
      END MODULE CMUMPS_FAC_ASM_MASTER_M
