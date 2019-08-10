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
      MODULE SMUMPS_FACTO_NIV1_M
      CONTAINS
      SUBROUTINE SMUMPS_FACTO_NIV1(
     &                           N, INODE, IW, LIW, A, 
     &                           LA,
     &                           IOLDPS, POSELT, IFLAG, UU, NOFFW,
     &                           NPVW,
     &                           KEEP,KEEP8, STEP,
     &                           PROCNODE_STEPS, MYID, SLAVEF, SEUIL,
     &                           AVOID_DELAYED, ETATASS,
     &     DKEEP,PIVNUL_LIST,LPN_LIST, 
     &     IWPOS 
     &     )
      USE SMUMPS_FAC_FRONT_AUX_M
      USE SMUMPS_OOC
      IMPLICIT NONE
      INTEGER(8) :: LA, POSELT
      INTEGER N, INODE, LIW, IFLAG, NOFFW, NPVW
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER MYID, SLAVEF, IOLDPS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER PROCNODE_STEPS( KEEP(28) ), STEP(N)
      REAL UU, SEUIL
      LOGICAL AVOID_DELAYED
      INTEGER ETATASS, IWPOS
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(130)
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEG_BLOCK, IEND_BLOCK
      INTEGER NASS, NBKJIB_ORIG, XSIZE
      INTEGER NBLR_ORIG, IBEG_BLR, IEND_BLR
      INTEGER LAST_ROW_TMP
      REAL UUTEMP
      LOGICAL STATICMODE
      REAL SEUIL_LOC
      LOGICAL POSTPONE_COL_UPDATE
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, LNextPiv2beWritten, 
     &        UNextPiv2beWritten, IFLAG_OOC,
     &        PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_LastPIVRPTRFilled_U
      INTEGER TYPEF_LOC
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      INTEGER GROUPLOC, CURRENT_BLR
      LOGICAL LR_ACTIVATED, PIVOTING_IN_LRBLOCK
      INCLUDE 'mumps_headers.h'
      INOPV = 0
      SEUIL_LOC = SEUIL
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
      ENDIF
      POSTPONE_COL_UPDATE = .TRUE.
      LAFAC  = -9999_8  
      XSIZE      = KEEP(IXSZ)
      NFRONT     = IW(IOLDPS+XSIZE)
      NASS       = iabs(IW(IOLDPS+2+XSIZE))
      IW(IOLDPS+3+XSIZE) =  -99999
      LR_ACTIVATED= .FALSE.        
      PIVOTING_IN_LRBLOCK= .FALSE. 
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
      IEND_BLOCK  = 0
      IEND_BLR    = 0
      CURRENT_BLR = 0
      IF (KEEP(201).EQ.1) THEN 
          CALL MUMPS_GETI8(LAFAC,IW(IOLDPS+XXR))
          LIWFAC    = IW(IOLDPS+XXI)
          LNextPiv2beWritten = 1 
          UNextPiv2beWritten = 1 
          PP_FIRST2SWAP_L = LNextPiv2beWritten 
          PP_FIRST2SWAP_U = UNextPiv2beWritten 
          MonBloc%LastPanelWritten_L = 0
          MonBloc%LastPanelWritten_U = 0
          PP_LastPIVRPTRFilled_L = 0 
          PP_LastPIVRPTRFilled_U = 0 
          MonBloc%INODE    = INODE
          MonBloc%MASTER   = .TRUE.
          MonBloc%Typenode = 1
          MonBloc%NROW     = NFRONT
          MonBloc%NCOL     = NFRONT
          MonBloc%NFS      = NASS
          MonBloc%Last     = .FALSE.   
          MonBloc%LastPiv  = -88877    
          NULLIFY(MonBloc%INDICES)   
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
            CALL SMUMPS_FAC_I(NFRONT,NASS,NFRONT,
     &      IBEG_BLOCK,IEND_BLOCK,N,INODE,
     &      IW,LIW,A,LA,INOPV,NOFFW,
     &      IFLAG,IOLDPS,POSELT,UU,SEUIL_LOC,KEEP,KEEP8,
     &      DKEEP(1),PIVNUL_LIST(1),LPN_LIST,
     &      PP_FIRST2SWAP_L,  MonBloc%LastPanelWritten_L,
     &      PP_LastPIVRPTRFilled_L,
     &      PP_FIRST2SWAP_U,  MonBloc%LastPanelWritten_U,
     &      PP_LastPIVRPTRFilled_U)
            IF (IFLAG.LT.0) GOTO 500  
          IF (INOPV.EQ.1) THEN
            IF(STATICMODE) THEN
              INOPV = -1
              GOTO 50 
            ENDIF
          ELSE IF ( INOPV.LE.0 ) THEN 
            CALL SMUMPS_FAC_MQ(IBEG_BLOCK, IEND_BLOCK,
     &              NFRONT, NASS, IW(IOLDPS+1+XSIZE),
     &              A, LA, POSELT, IFINB)
            IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + 1
            IF (IFINB.EQ.0) THEN
              GOTO 50 
            ENDIF
          ENDIF
          IF (KEEP(201).EQ.1) THEN  
            MonBloc%LastPiv= IW(IOLDPS+1+XSIZE)
            STRAT          = STRAT_TRY_WRITE
            LAST_CALL      = .FALSE.
            CALL SMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_U,
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
            IF (IFLAG_OOC < 0 ) THEN
              IFLAG=IFLAG_OOC
              GOTO 500
            ENDIF
          ENDIF
          NPIV       =  IW(IOLDPS+1+XSIZE)
          IF ( IEND_BLR .GT. IEND_BLOCK ) THEN
            CALL SMUMPS_FAC_SQ(IBEG_BLOCK, IEND_BLOCK,
     &            NPIV, NFRONT, IEND_BLR,
     &            A, LA, POSELT, .TRUE.)
          ENDIF
        END DO 
        NPIV   =  IW(IOLDPS+1+XSIZE)
        IF (.NOT. LR_ACTIVATED) THEN
          IF (POSTPONE_COL_UPDATE) THEN
            LAST_ROW_TMP = NASS
          ELSE
            LAST_ROW_TMP = NFRONT
          ENDIF
          IF (IEND_BLR.LT.LAST_ROW_TMP) THEN
            CALL SMUMPS_FAC_SQ(IBEG_BLR, 
     &            IEND_BLR,
     &            NPIV, NFRONT, LAST_ROW_TMP, A, LA,
     &            POSELT, .TRUE.)
          ENDIF
        ENDIF
        IF (KEEP(201).EQ.1) THEN  
             IF (POSTPONE_COL_UPDATE) THEN
               TYPEF_LOC = TYPEF_U
             ELSE
               TYPEF_LOC = TYPEF_BOTH_LU
             ENDIF
             MonBloc%LastPiv= IW(IOLDPS+1+XSIZE)
             STRAT          = STRAT_TRY_WRITE
             LAST_CALL      = .FALSE.
             CALL SMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_LOC,
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
             IF (IFLAG_OOC < 0 ) THEN
                IFLAG=IFLAG_OOC
                GOTO 500
             ENDIF
        ENDIF
      END DO 
      IF ( POSTPONE_COL_UPDATE ) THEN
        IF (LR_ACTIVATED) THEN
          WRITE(*,*) "Internal error 1 in SMUMPS_FACTO_NIV1"
          CALL MUMPS_ABORT()
        ENDIF
        CALL  SMUMPS_FAC_FR_UPDATE_CBROWS( INODE,
     &     NFRONT, NASS, A, LA, LAFAC, POSELT, IW, LIW, IOLDPS,
     &     MonBloc, MYID, NOFFW, LIWFAC, 
     &     PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &     LNextPiv2beWritten, UNextPiv2beWritten, 
     &     PP_LastPIVRPTRFilled_L, PP_LastPIVRPTRFilled_U,
     &     
     &     XSIZE, SEUIL, UU, DKEEP, KEEP8, KEEP, IFLAG)
      ENDIF
      IF (KEEP(201).EQ.1) THEN 
          STRAT            = STRAT_WRITE_MAX   
          MonBloc%Last     = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+XSIZE)
          LAST_CALL    = .TRUE.
          CALL SMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_BOTH_LU,
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC, LAST_CALL )
          IF (IFLAG_OOC < 0 ) THEN
            IFLAG=IFLAG_OOC
            GOTO 500
          ENDIF
          CALL SMUMPS_OOC_PP_TRYRELEASE_SPACE (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
 500  CONTINUE
      NPVW = NPVW + IW(IOLDPS+1+XSIZE)
      RETURN
      END SUBROUTINE SMUMPS_FACTO_NIV1
      END MODULE SMUMPS_FACTO_NIV1_M
