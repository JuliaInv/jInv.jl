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
      MODULE SMUMPS_FAC1_LDLT_M
      CONTAINS
      SUBROUTINE SMUMPS_FAC1_LDLT( N, INODE, IW, LIW, A, LA,
     &                           IOLDPS, POSELT, IFLAG,
     &                           UU, NNEG, NPVW,
     &                           KEEP,KEEP8,
     &                           MYID, SEUIL, AVOID_DELAYED, ETATASS,
     &     DKEEP,PIVNUL_LIST,LPN_LIST, IWPOS
     &     )
      USE SMUMPS_FAC_FRONT_AUX_M
      USE SMUMPS_OOC
      IMPLICIT NONE
      INTEGER(8) :: LA, POSELT
      INTEGER N, INODE, LIW, IFLAG, NNEG, NPVW
      INTEGER MYID, IOLDPS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      REAL UU, SEUIL
      REAL A( LA )
      INTEGER, TARGET :: IW( LIW )
      LOGICAL AVOID_DELAYED
      INTEGER ETATASS, IWPOS
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(130)
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEG_BLOCK, IEND_BLOCK
      INTEGER NASS, NBKJIB_ORIG, XSIZE
      INTEGER :: LDA
      REAL UUTEMP
      LOGICAL STATICMODE
      REAL SEUIL_LOC
      LOGICAL POSTPONE_COL_UPDATE, IS_MAXFROMM_AVAIL
      REAL MAXFROMM
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, NextPiv2beWritten, IFLAG_OOC,
     &        IDUMMY, PP_FIRST2SWAP_L, PP_LastPIVRPTRFilled
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      LOGICAL LASTBL
      INTEGER GROUPLOC, CURRENT_BLR
      LOGICAL LR_ACTIVATED 
      INCLUDE 'mumps_headers.h'
      INTEGER NBLR_ORIG, IBEG_BLR, IEND_BLR
      INTEGER PIVSIZ,IWPOSP2
      IS_MAXFROMM_AVAIL = .FALSE.
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
        UUTEMP=UU
        SEUIL_LOC = SEUIL
      ENDIF
      POSTPONE_COL_UPDATE = (UUTEMP == 0.0E0 .AND. KEEP(201).NE.1)
      LAFAC  = -9999_8  
      XSIZE  = KEEP(IXSZ)
      NFRONT = IW(IOLDPS+XSIZE)
      LDA    = NFRONT
      NASS   = iabs(IW(IOLDPS+2+XSIZE))
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
             write(6,*) ' FR facto with LR grouping not validated yet'
             write(6,*) ' try with KEEP(486) = 0 or 1 '
             CALL MUMPS_ABORT()
          ELSE
             GROUPLOC = 0
          ENDIF
      ELSE
          NBLR_ORIG  = -9999 
          GROUPLOC    = KEEP(486)
      ENDIF
      IEND_BLOCK  = 0
      IEND_BLR    = 0
      CURRENT_BLR = 0
      LASTBL      = .FALSE.
      IF (KEEP(201).EQ.1) THEN 
          IDUMMY    = -8765
          CALL MUMPS_GETI8(LAFAC,IW(IOLDPS+XXR))
          LIWFAC    = IW(IOLDPS+XXI)
          NextPiv2beWritten = 1 
          PP_FIRST2SWAP_L = NextPiv2beWritten 
          MonBloc%LastPanelWritten_L = 0 
          PP_LastPIVRPTRFilled       = 0
          MonBloc%INODE    = INODE
          MonBloc%MASTER   = .TRUE.
          MonBloc%Typenode = 1
          MonBloc%NROW     = NFRONT
          MonBloc%NCOL     = NFRONT
          MonBloc%NFS      = NASS
          MonBloc%Last     = .FALSE.   
          MonBloc%LastPiv  = -77777    
          MonBloc%INDICES  =>
     &              IW(IOLDPS+6+NFRONT+XSIZE:
     &                 IOLDPS+5+NFRONT+XSIZE+NFRONT)
      ENDIF
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
            CALL SMUMPS_FAC_I_LDLT(NFRONT,NASS,INODE,IEND_BLOCK,
     &                IW,LIW,A,LA,
     &                INOPV, NNEG, IFLAG,IOLDPS,POSELT,UUTEMP,
     &                SEUIL_LOC,KEEP,KEEP8,PIVSIZ,
     &      DKEEP(1),PIVNUL_LIST(1),LPN_LIST, XSIZE,
     &      PP_FIRST2SWAP_L, MonBloc%LastPanelWritten_L,
     &      PP_LastPIVRPTRFilled, MAXFROMM, IS_MAXFROMM_AVAIL)
            IF (IFLAG.LT.0) GOTO 500
            IF(KEEP(109).GT. 0) THEN
              IF(PIVNUL_LIST(KEEP(109)).EQ.-1) THEN
                IWPOSP2 = IOLDPS+IW(IOLDPS+1+XSIZE)+6+XSIZE
     &                  +IW(IOLDPS+5+XSIZE)
                PIVNUL_LIST(KEEP(109)) = IW(IWPOSP2)
              ENDIF
            ENDIF
          IF (INOPV.EQ.1) THEN
            IF(STATICMODE) THEN
              INOPV = -1
              GOTO 50 
            ENDIF
            LASTBL = .TRUE.
          ELSE IF ( INOPV.LE.0 ) THEN 
            NPVW = NPVW + PIVSIZ
            CALL SMUMPS_FAC_MQ_LDLT(IEND_BLOCK,
     &             NFRONT, NASS, IW(IOLDPS+1+XSIZE),
     &             INODE,A,LA,
     &             LDA, POSTPONE_COL_UPDATE, 
     &             POSELT,IFINB,
     &             PIVSIZ, MAXFROMM,
     &             IS_MAXFROMM_AVAIL, (UUTEMP.NE.0.0E0),
     &             KEEP(253) )
            IF(PIVSIZ .EQ. 2) THEN
              IWPOSP2 = IOLDPS+IW(IOLDPS+1+XSIZE)+6
              IW(IWPOSP2+NFRONT+XSIZE) =
     &                              -IW(IWPOSP2+NFRONT+XSIZE)
            ENDIF
            IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + PIVSIZ
            IF (IFINB.EQ.0) THEN
              GOTO 50 
            ELSE IF (IFINB.EQ.-1) THEN
              LASTBL = .TRUE.
            ENDIF
          ENDIF
          IF (KEEP(201).EQ.1.AND. .NOT. POSTPONE_COL_UPDATE) THEN
            MonBloc%Last = LASTBL
            MonBloc%LastPiv= IW(IOLDPS+1+XSIZE)
            LAST_CALL=.FALSE.
            CALL SMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE,
     &        TYPEF_L, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY,
     &        IW(IOLDPS), LIWFAC,
     &        MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
            IF (IFLAG_OOC < 0 ) THEN
              IFLAG=IFLAG_OOC
              GOTO 500
            ENDIF
          ENDIF
          NPIV       =  IW(IOLDPS+1+XSIZE)
          IF ( IEND_BLR .GT. IEND_BLOCK ) THEN
            CALL SMUMPS_FAC_SQ_LDLT(IBEG_BLOCK,IEND_BLOCK,
     &             NPIV, NFRONT,NASS,IEND_BLR,INODE,A,LA,
     &             LDA, POSELT,
     &             POSTPONE_COL_UPDATE,
     &             KEEP,KEEP8)
          ENDIF
        END DO 
        NPIV   =  IW(IOLDPS+1+XSIZE)
        IF (
     &       (.NOT. LR_ACTIVATED) 
     &     )
     &  THEN
          CALL SMUMPS_FAC_SQ_LDLT(IBEG_BLR,IEND_BLR,NPIV,
     &             NFRONT,NASS,NASS,INODE,A,LA,
     &             LDA, POSELT,
     &             POSTPONE_COL_UPDATE,
     &             KEEP,KEEP8)
        ENDIF
        IF (KEEP(201) .EQ. 1 .AND. .NOT. POSTPONE_COL_UPDATE) THEN 
             MonBloc%Last = LASTBL
             MonBloc%LastPiv= NPIV
             LAST_CALL=.FALSE.
             CALL SMUMPS_OOC_IO_LU_PANEL(
     &          STRAT_TRY_WRITE,
     &          TYPEF_L, A(POSELT),
     &          LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &          LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
             IF (IFLAG_OOC < 0 ) THEN
                IFLAG=IFLAG_OOC
                GOTO 500
             ENDIF
        ENDIF
      END DO 
      IF ( LR_ACTIVATED
     &   )  THEN
        IF ( POSTPONE_COL_UPDATE ) THEN
          WRITE(*,*) "Internal error 1 in SMUMPS_FACTO_NIV1"
          CALL MUMPS_ABORT()
        ENDIF
      ELSE
        CALL SMUMPS_FAC_T_LDLT(NFRONT,NASS,IW,LIW,A,LA,
     &         LDA, IOLDPS,POSELT, KEEP,KEEP8,
     &         POSTPONE_COL_UPDATE, ETATASS,
     &         TYPEF_L, LAFAC, MonBloc, NextPiv2beWritten,
     &         LIWFAC, MYID, IFLAG)
      ENDIF
      IF (KEEP(201).EQ.1) THEN 
          STRAT        = STRAT_WRITE_MAX   
          MonBloc%Last = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+XSIZE)
          LAST_CALL    = .TRUE.
          CALL SMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_L, 
     &           A(POSELT), LAFAC, MonBloc,
     &           NextPiv2beWritten, IDUMMY,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
          IF (IFLAG_OOC < 0 ) THEN
            IFLAG=IFLAG_OOC
            GOTO 500
          ENDIF
          CALL SMUMPS_OOC_PP_TRYRELEASE_SPACE (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
 500  CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_FAC1_LDLT
      END MODULE SMUMPS_FAC1_LDLT_M
