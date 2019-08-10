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
      MODULE CMUMPS_FAC2_LDLT_M
      CONTAINS
      SUBROUTINE CMUMPS_FAC2_LDLT( COMM_LOAD, ASS_IRECV, 
     &           N, INODE, FPERE, IW, LIW, A, LA,
     &           UU, NOFFW,
     &           NPVW,
     &             COMM, MYID, BUFR, LBUFR,LBUFR_BYTES,NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &             PIMASTER, PAMASTER,
     &             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     &             OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL,
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     &             DKEEP,PIVNUL_LIST,LPN_LIST
     &             )
      USE CMUMPS_FAC_FRONT_AUX_M
      USE CMUMPS_FAC_FRONT_TYPE2_AUX_M
      USE CMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'cmumps_root.h'
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, NOFFW, NPVW
      INTEGER(8) :: LA
      INTEGER, TARGET :: IW( LIW )
      COMPLEX A( LA )
      REAL UU, SEUIL
      TYPE (CMUMPS_ROOT_STRUC) :: root
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER LPTRAR, NELT
      INTEGER ICNTL(40), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER NBFIN, SLAVEF,
     &        IFLAG, IERROR, LEAF, LPOOL
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB, COMP 
      INTEGER NB_BLOC_FAC
      INTEGER FRTPTR( N + 1 ), FRTELT( NELT )
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL),
     &        ITLOC(N+KEEP(253)), FILS(N),
     &        PTRARW(LPTRAR), PTRAIW(LPTRAR),
     &        ND( KEEP(28) ), FRERE( KEEP(28) )
      COMPLEX :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(max(1,KEEP(14)))
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28)),
     &        STEP(N), PIMASTER(KEEP(28)),
     &        NSTK_S(KEEP(28)), NBPROCFILS(KEEP(28)),
     &        PROCNODE_STEPS(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION OPASSW, OPELIW
      COMPLEX DBLARR(max(1,KEEP(13)))
      LOGICAL AVOID_DELAYED
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(130)
      INTEGER(8) :: POSELT
      INTEGER IOLDPS, allocok, K263
      INTEGER INOPV, IFINB, NFRONT, NPIV, IEND_BLOCK
      INTEGER NASS, LDAFS, IBEG_BLOCK
      INTEGER :: IBEG_BLOCK_FOR_IPIV
      LOGICAL LASTBL, LR_ACTIVATED
      INTEGER GROUPLOC
      INTEGER NBLR_ORIG, IBEG_BLR, IEND_BLR, CURRENT_BLR
      LOGICAL RESET_TO_ONE
      INTEGER K109_SAVE
      INTEGER XSIZE, NBKJIB_ORIG
      REAL UUTEMP
      INCLUDE 'mumps_headers.h'
      INTEGER , ALLOCATABLE, DIMENSION ( : ) :: IPIV
      REAL , ALLOCATABLE, DIMENSION ( : )    :: DIAG_ORIG
      INTEGER    :: SIZEDIAG_ORIG
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, TYPEFile, NextPiv2beWritten,
     &        IDUMMY, NELIM
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      INTEGER PP_FIRST2SWAP_L, IFLAG_OOC
      INTEGER PP_LastPIVRPTRFilled 
      EXTERNAL CMUMPS_BDC_ERROR
      LOGICAL STATICMODE
      REAL SEUIL_LOC
      REAL GW_FACTCUMUL
      INTEGER PIVSIZ,IWPOSPIV
      COMPLEX ONE
      PARAMETER (ONE=(1.0E0,0.0E0))
      INOPV = 0
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      IF (AVOID_DELAYED) THEN
        STATICMODE = .TRUE.
        UUTEMP=UU
        SEUIL_LOC = max(SEUIL,epsilon(SEUIL))
      ELSE
        SEUIL_LOC=SEUIL
        UUTEMP=UU
      ENDIF
      RESET_TO_ONE = ((KEEP(110).GT.0).AND.(DKEEP(2).LE.0.0E0))
      IF (RESET_TO_ONE) THEN
        K109_SAVE = KEEP(109)
      ENDIF
      IBEG_BLOCK  = 1
      NB_BLOC_FAC = 0
      XSIZE  = KEEP(IXSZ)
      IOLDPS = PTLUST_S(STEP( INODE ))
      POSELT = PTRAST(STEP( INODE ))
      NFRONT = IW(IOLDPS+XSIZE)
      NASS   = iabs(IW(IOLDPS+2+XSIZE))
      LDAFS  = NASS
      IW(IOLDPS+3+XSIZE) =  -99999
      LR_ACTIVATED= .FALSE. 
      IF (NASS.LT.KEEP(4)) THEN
        NBKJIB_ORIG = NASS
      ELSE IF (NASS .GT. KEEP(3)) THEN
        NBKJIB_ORIG = min( KEEP(6), NASS )
      ELSE
        NBKJIB_ORIG = min( KEEP(5), NASS )
      ENDIF
      IF (.not.LR_ACTIVATED) THEN
          NBLR_ORIG     = KEEP(488)
          IF (KEEP(486) .EQ. -1) THEN
             GROUPLOC = 1
             CURRENT_BLR = 1
             write(*,*) '=============================================='
             write(6,*) ' FR facto with LR grouping not validated yet'
             write(6,*) ' try with KEEP(486) = 0 or 1 '
             write(*,*) '=============================================='
             CALL MUMPS_ABORT()
          ELSE
             GROUPLOC = 0
          ENDIF
      ELSE
          NBLR_ORIG  = -9999 
          GROUPLOC    = KEEP(486)
      ENDIF
      IF (LR_ACTIVATED) THEN
         K263 = 1   
      ELSE
         K263 = KEEP(263)
      ENDIF
      IEND_BLOCK  = 0
      IEND_BLR    = 0
      CURRENT_BLR = 0
      ALLOCATE( IPIV( NASS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) MYID, ' : CMUMPS_FAC2_LDLT failed to allocate ',
     &  NASS, ' integers'
        IFLAG = -13
        IERROR=NASS
        GO TO 490
      END IF
      IF (KEEP(219).GE.3) THEN
       SIZEDIAG_ORIG = NASS
      ELSE
       SIZEDIAG_ORIG = 1
      ENDIF
      ALLOCATE ( DIAG_ORIG(SIZEDIAG_ORIG), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
          WRITE(*,*) MYID,
     &      ' : FACTO_NIV2 failed to allocate ',
     &      NASS, ' REAL/COMPLEX entries'
          IFLAG=-13
          IERROR=NASS
          GO TO 490
      END IF
      IF (KEEP(201).EQ.1) THEN 
        IDUMMY    = -9876
        CALL MUMPS_GETI8(LAFAC,IW(IOLDPS+XXR))
        LIWFAC    = IW(IOLDPS+XXI)
        TYPEFile  = TYPEF_L
        NextPiv2beWritten = 1 
        PP_FIRST2SWAP_L = NextPiv2beWritten 
        MonBloc%LastPanelWritten_L = 0 
        MonBloc%INODE    = INODE
        MonBloc%MASTER   = .TRUE.
        MonBloc%Typenode = 2
        MonBloc%NROW     = NASS
        MonBloc%NCOL     = NASS
        MonBloc%NFS      = NASS
        MonBloc%Last     = .FALSE.   
        MonBloc%LastPiv  = -66666    
        MonBloc%INDICES =>
     &  IW(IOLDPS+6+NFRONT+XSIZE+IW(IOLDPS+5+XSIZE)
     &    :IOLDPS+5+2*NFRONT+XSIZE+IW(IOLDPS+5+XSIZE))
      ENDIF
      LASTBL = .FALSE.
      DO WHILE (IEND_BLR < NASS ) 
        CURRENT_BLR = CURRENT_BLR + 1
        IBEG_BLR = IW(IOLDPS+1+KEEP(IXSZ)) + 1 
        IF (.NOT. LR_ACTIVATED .AND. GROUPLOC .EQ. 0)THEN
          IEND_BLR = min(IEND_BLR + NBLR_ORIG, NASS)
        ENDIF
        DO WHILE (IEND_BLOCK < IEND_BLR ) 
          IBEG_BLOCK = IW(IOLDPS+1+KEEP(IXSZ)) + 1
          IEND_BLOCK = min(IEND_BLOCK + NBKJIB_ORIG, IEND_BLR)
  50      CONTINUE  
            IF (K263.EQ.0) THEN
              IBEG_BLOCK_FOR_IPIV = IBEG_BLOCK
            ELSE
              IBEG_BLOCK_FOR_IPIV = IBEG_BLR
            ENDIF
            CALL CMUMPS_FAC_I_LDLT_NIV2(
     &                DIAG_ORIG, SIZEDIAG_ORIG, GW_FACTCUMUL,
     &                NFRONT,NASS,IBEG_BLOCK_FOR_IPIV,IEND_BLOCK,
     &                NASS, IPIV,
     &                N,INODE,IW,LIW,A,LA,NOFFW,INOPV,
     &                IFLAG,IOLDPS,POSELT,UU, SEUIL_LOC,
     &                KEEP,KEEP8,PIVSIZ,
     &           DKEEP(1),PIVNUL_LIST(1),LPN_LIST,
     &           PP_FIRST2SWAP_L, MonBloc%LastPanelWritten_L,
     &           PP_LastPIVRPTRFilled)
            IF (IFLAG.LT.0) GOTO 490
            IF(KEEP(109).GT. 0) THEN
              IF(PIVNUL_LIST(KEEP(109)).EQ.-1) THEN
                IWPOSPIV = IOLDPS+IW(IOLDPS+1+XSIZE)+6
     &              +IW(IOLDPS+5+XSIZE)
                PIVNUL_LIST(KEEP(109)) = IW(IWPOSPIV+XSIZE)
              ENDIF
            ENDIF
            IF (INOPV.EQ. 1) THEN
              IF (STATICMODE) THEN
                INOPV = -1
                GOTO 50
              ENDIF
             LASTBL = .TRUE.
            ELSE IF (INOPV .LE. 0) THEN 
              NPVW = NPVW + PIVSIZ
              CALL CMUMPS_FAC_MQ_LDLT_NIV2(IEND_BLOCK,
     &             NASS, IW(IOLDPS+1+XSIZE), INODE,A,LA,
     &             LDAFS, POSELT,IFINB,
     &             PIVSIZ,
     &             KEEP(219))
              IF(PIVSIZ .EQ. 2) THEN
                IWPOSPIV = IOLDPS+XSIZE+IW(IOLDPS+1+XSIZE)+6+
     &                     IW(IOLDPS+5+XSIZE)
                IW(IWPOSPIV+NFRONT) = -IW(IWPOSPIV+NFRONT)
              ENDIF
              IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + PIVSIZ
            IF (IFINB.EQ.0) THEN
              GOTO 50 
            ELSE IF (IFINB .EQ. -1) THEN
              LASTBL = .TRUE.
            ENDIF
          ENDIF
          NPIV = IW(IOLDPS+1+XSIZE)
          IF (KEEP(201).EQ.1) THEN
            IF (.NOT.RESET_TO_ONE.OR.K109_SAVE.EQ.KEEP(109)) THEN
              MonBloc%Last   = .FALSE.
              MonBloc%LastPiv= NPIV
              LAST_CALL=.FALSE.
              CALL CMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE,
     &        TYPEFile, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &        LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
              IF (IFLAG_OOC .LT. 0 ) IFLAG = IFLAG_OOC
              IF (IFLAG .LT. 0 ) RETURN
            ENDIF
          ENDIF
          IF (K263.eq.0) THEN
            NELIM = IEND_BLR-NPIV
            CALL CMUMPS_SEND_FACTORED_BLOCK( COMM_LOAD, ASS_IRECV,
     &             N, INODE, FPERE, IW, LIW, 
     &             IOLDPS, POSELT, A, LA, LDAFS,
     &             IBEG_BLOCK, NPIV, IPIV, NASS,LASTBL, NB_BLOC_FAC,
     &             COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &             PIMASTER, PAMASTER,
     &             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     &             OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, 
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE
     & )
            IF ( IFLAG .LT. 0 ) GOTO 500
            IF (RESET_TO_ONE.AND.K109_SAVE.LT.KEEP(109)) THEN
              CALL CMUMPS_RESET_TO_ONE( 
     &        IW(IOLDPS+KEEP(IXSZ)+IW(IOLDPS+5+KEEP(IXSZ))+6),
     &        NPIV, IBEG_BLOCK,
     &        K109_SAVE, KEEP(109), PIVNUL_LIST, LPN_LIST,
     &        A, POSELT, LA, LDAFS)
            ENDIF
            IF (KEEP(201).EQ.1) THEN
              MonBloc%Last  = .FALSE.
              MonBloc%LastPiv= NPIV
              LAST_CALL=.FALSE.
              CALL CMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE,
     &        TYPEFile, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &        LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
              IF (IFLAG_OOC .LT. 0 ) THEN
                IFLAG = IFLAG_OOC
                RETURN
              ENDIF
            ENDIF
          ENDIF
          IF ( IEND_BLR .GT. IEND_BLOCK ) THEN
            CALL CMUMPS_FAC_SQ_LDLT(IBEG_BLOCK,IEND_BLOCK,NPIV,
     &             NASS,NASS,IEND_BLR,INODE,A,LA,
     &             LDAFS, POSELT,
     &             .TRUE.,
     &             KEEP,KEEP8)
          ENDIF
        END DO 
        NPIV   = IW(IOLDPS+1+XSIZE)
        IF (K263.NE.0) THEN
          NELIM = IEND_BLR-NPIV
          CALL CMUMPS_SEND_FACTORED_BLOCK( COMM_LOAD, ASS_IRECV,
     &             N, INODE, FPERE, IW, LIW, 
     &             IOLDPS, POSELT, A, LA, LDAFS,
     &             IBEG_BLR, NPIV, IPIV, NASS,LASTBL, NB_BLOC_FAC,
     &
     &             COMM, MYID, BUFR, LBUFR, LBUFR_BYTES,NBFIN,LEAF,
     &             IFLAG, IERROR, IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP,
     &             PTRIST, PTRAST, PTLUST_S, PTRFAC, STEP,
     &             PIMASTER, PAMASTER,
     &             NSTK_S,NBPROCFILS,PROCNODE_STEPS, root,
     &             OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &             FILS, PTRARW, PTRAIW,
     &             INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, 
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &             )
          IF ( IFLAG .LT. 0 ) GOTO 500
          IF (RESET_TO_ONE.AND.K109_SAVE.LT.KEEP(109)) THEN
              CALL CMUMPS_RESET_TO_ONE( 
     &        IW(IOLDPS+KEEP(IXSZ)+IW(IOLDPS+5+KEEP(IXSZ))+6),
     &        NPIV, IBEG_BLR,
     &        K109_SAVE, KEEP(109), PIVNUL_LIST, LPN_LIST,
     &        A, POSELT, LA, LDAFS)
          ENDIF
          IF (KEEP(201).EQ.1) THEN
              MonBloc%Last  = .FALSE.
              MonBloc%LastPiv= NPIV
              LAST_CALL=.FALSE.
              CALL CMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE,
     &        TYPEFile, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &        LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
              IF (IFLAG_OOC .LT. 0 ) THEN
                IFLAG = IFLAG_OOC
                RETURN
              ENDIF
          ENDIF
        ENDIF
        IF (.NOT. LR_ACTIVATED) THEN
          CALL CMUMPS_FAC_SQ_LDLT(IBEG_BLR,IEND_BLR,NPIV,
     &             NASS,NASS,NASS,INODE,A,LA,
     &             LDAFS, POSELT,
     &             .TRUE.,
     &             KEEP,KEEP8)
        ENDIF
        IF (KEEP(201).EQ.1) THEN
          MonBloc%Last   = .FALSE.
          MonBloc%LastPiv= NPIV
          LAST_CALL=.FALSE.
          CALL CMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE,
     &        TYPEFile, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &        LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
          IF (IFLAG_OOC < 0 ) THEN
              IFLAG = IFLAG_OOC
              GOTO 490
          ENDIF
        ENDIF
      END DO 
      IF (KEEP(201).EQ.1) THEN 
          STRAT        = STRAT_WRITE_MAX   
          MonBloc%Last = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+XSIZE)
          LAST_CALL    = .TRUE.
          CALL CMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEFile, 
     &           A(POSELT), LAFAC, MonBloc,
     &           NextPiv2beWritten, IDUMMY,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC, LAST_CALL )
          IF (IFLAG_OOC .LT. 0 ) THEN
             IFLAG = IFLAG_OOC
             RETURN
          ENDIF
          CALL CMUMPS_OOC_PP_TRYRELEASE_SPACE (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
      GOTO 500
 490  CONTINUE
      CALL CMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
 500  CONTINUE
      IF(allocated(IPIV)) DEALLOCATE( IPIV )
      IF (allocated(DIAG_ORIG)) DEALLOCATE(DIAG_ORIG)
      RETURN
      END SUBROUTINE CMUMPS_FAC2_LDLT
      SUBROUTINE CMUMPS_RESET_TO_ONE(FRONT_INDEX_LIST, NPIV,
     & IBEG_BLOCK, K109_SAVE, K109, PIVNUL_LIST, LPN_LIST,
     & A, POSELT, LA, LDAFS)
      INTEGER, INTENT(IN) :: NPIV, IBEG_BLOCK
      INTEGER, INTENT(IN) :: FRONT_INDEX_LIST(NPIV)
      INTEGER, INTENT(IN) :: K109
      INTEGER, INTENT(INOUT) :: K109_SAVE
      INTEGER, INTENT(IN) :: LPN_LIST
      INTEGER, INTENT(IN) :: PIVNUL_LIST(LPN_LIST)
      INTEGER(8), INTENT(IN) :: POSELT, LA
      INTEGER, INTENT(IN) :: LDAFS
      COMPLEX, INTENT(INOUT) :: A(LA)
      LOGICAL :: TO_UPDATE
      INTEGER :: I, JJ, K
      COMPLEX ONE
      PARAMETER (ONE=(1.0E0,0.0E0))
      DO K = K109_SAVE+1, K109
        TO_UPDATE = .FALSE. 
        I = PIVNUL_LIST(K)  
        DO JJ=IBEG_BLOCK, NPIV
          IF (FRONT_INDEX_LIST(JJ) .EQ.I) THEN
            TO_UPDATE=.TRUE. 
            EXIT
          ENDIF
        ENDDO
        IF (TO_UPDATE) THEN
          A(POSELT+int(JJ,8)+int(LDAFS,8)*int(JJ-1,8))= ONE
          TO_UPDATE=.FALSE. 
        ELSE
          write(*,*) ' Internal error related ', 
     &               'to null pivot row detection'
          CALL MUMPS_ABORT()
        ENDIF
      ENDDO
      K109_SAVE = K109
      RETURN
      END SUBROUTINE CMUMPS_RESET_TO_ONE
      END MODULE CMUMPS_FAC2_LDLT_M
