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
      MODULE SMUMPS_FACTO_NIV2_M
      CONTAINS
      SUBROUTINE SMUMPS_FACTO_NIV2( COMM_LOAD, ASS_IRECV, 
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
      USE SMUMPS_FAC_FRONT_AUX_M
      USE SMUMPS_FAC_FRONT_TYPE2_AUX_M
      USE SMUMPS_OOC
      IMPLICIT NONE
      INCLUDE 'smumps_root.h'
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER N, INODE, FPERE, LIW, NOFFW, NPVW
      INTEGER(8) :: LA
      INTEGER IW( LIW )
      REAL A( LA )
      REAL UU, SEUIL
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER COMM, MYID, LBUFR, LBUFR_BYTES
      INTEGER LPTRAR, NELT
      INTEGER ICNTL(40), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER NBFIN, SLAVEF,
     &        IFLAG, IERROR, LEAF, LPOOL
      INTEGER(8) :: POSFAC, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB, COMP 
      INTEGER FRTPTR( N + 1 ), FRTELT( NELT )
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL),
     &        ITLOC(N+KEEP(253)), FILS(N),
     &        PTRARW(LPTRAR), PTRAIW(LPTRAR),
     &        ND( KEEP(28) ), FRERE( KEEP(28) )
      REAL :: RHS_MUMPS(KEEP(255))
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
      REAL DBLARR(max(1,KEEP(13)))
      LOGICAL AVOID_DELAYED
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      REAL DKEEP(130)
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEG_BLOCK, IEND_BLOCK
      INTEGER :: IBEG_BLOCK_FOR_IPIV
      INTEGER NASS, NBKJIB_ORIG, XSIZE
      INTEGER NBLR_ORIG, IBEG_BLR, IEND_BLR
      LOGICAL LASTBL 
      INTEGER(8) :: POSELT
      INTEGER IOLDPS, allocok, K263
      INTEGER idummy 
      REAL    UUTEMP
      LOGICAL STATICMODE
      REAL SEUIL_LOC
      INTEGER , ALLOCATABLE, DIMENSION ( : ) :: IPIV
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, LNextPiv2beWritten, 
     &        UNextPiv2beWritten, IFLAG_OOC,
     &        PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_LastPIVRPTRFilled_U
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      INTEGER GROUPLOC, CURRENT_BLR, NELIM
      LOGICAL LR_ACTIVATED, PIVOTING_IN_LRBLOCK
      LOGICAL COMPRESS_BEFORE_SEND
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
      idummy  = 0
      IOLDPS = PTLUST_S(STEP( INODE ))
      POSELT = PTRAST(STEP( INODE ))
      XSIZE  = KEEP(IXSZ)
      NFRONT = IW(IOLDPS+XSIZE)
      NASS   = iabs(IW(IOLDPS+2+XSIZE))
      IW(IOLDPS+3+XSIZE) =  -99999
      LR_ACTIVATED= .FALSE.        
      PIVOTING_IN_LRBLOCK= .FALSE. 
      COMPRESS_BEFORE_SEND = .TRUE. 
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
        WRITE(*,*) MYID,' : SMUMPS_FACTO_NIV2 :failed to allocate ',
     &  NASS, ' integers'
        IFLAG  = -13
        IERROR =NASS
        GO TO 490
      END IF
      IF (KEEP(201).EQ.1) THEN 
          CALL MUMPS_GETI8(LAFAC,IW(IOLDPS+XXR))
          LIWFAC    = IW(IOLDPS+XXI)
          LNextPiv2beWritten = 1 
          UNextPiv2beWritten = 1 
          PP_FIRST2SWAP_L = LNextPiv2beWritten 
          PP_FIRST2SWAP_U = UNextPiv2beWritten 
          MonBloc%LastPanelWritten_L = 0 
          MonBloc%LastPanelWritten_U = 0        
          MonBloc%INODE    = INODE
          MonBloc%MASTER   = .TRUE.
          MonBloc%Typenode = 2
          MonBloc%NROW     = NASS
          MonBloc%NCOL     = NFRONT
          MonBloc%NFS      = NASS
          MonBloc%Last     = .FALSE.   
          MonBloc%LastPiv  = -68877    
          NULLIFY(MonBloc%INDICES)
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
            CALL SMUMPS_FAC_I(NFRONT,NASS,NASS,
     &      IBEG_BLOCK_FOR_IPIV,IEND_BLOCK,N,INODE,
     &      IW,LIW,A,LA,INOPV,NOFFW,
     &      IFLAG,IOLDPS,POSELT,UU,SEUIL_LOC,KEEP,KEEP8,
     &      DKEEP(1),PIVNUL_LIST(1),LPN_LIST,
     &      PP_FIRST2SWAP_L,  MonBloc%LastPanelWritten_L,
     &      PP_LastPIVRPTRFilled_L,
     &      PP_FIRST2SWAP_U,  MonBloc%LastPanelWritten_U,
     &      PP_LastPIVRPTRFilled_U,
     &      IPIV) 
            IF (IFLAG.LT.0) GOTO 490   
          IF (INOPV.EQ.1) THEN
              IF (STATICMODE) THEN
                INOPV = -1
                GOTO 50
              ENDIF
              LASTBL = .TRUE.
          ELSE IF (INOPV .LE. 0) THEN 
            CALL SMUMPS_FAC_MQ(IBEG_BLOCK, IEND_BLOCK,
     &              NFRONT, NASS, IW(IOLDPS+1+XSIZE),
     &              A, LA, POSELT, IFINB)
            IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + 1
            NPVW = NPVW + 1
            IF (IFINB.EQ.0) THEN
              GOTO 50 
            ELSE IF (IFINB .EQ. -1) THEN
              LASTBL = .TRUE.
            ENDIF
          ENDIF
          NPIV = IW(IOLDPS+1+XSIZE)
          IF (K263.EQ.0) THEN
            NELIM = IEND_BLR - NPIV
            CALL SMUMPS_SEND_FACTORED_BLOCK( COMM_LOAD, ASS_IRECV, 
     &             N, INODE, FPERE, IW, LIW, 
     &             IOLDPS, POSELT, A, LA, NFRONT,
     &             IBEG_BLOCK, NPIV, IPIV, NASS,LASTBL, idummy, 
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
          END IF
          IF ( IFLAG .LT. 0 ) GOTO 500
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
              GOTO 490
            ENDIF
          ENDIF
          NPIV       =  IW(IOLDPS+1+XSIZE)
          IF ( IEND_BLR .GT. IEND_BLOCK ) THEN
            CALL SMUMPS_FAC_SQ(IBEG_BLOCK, IEND_BLOCK,
     &            NPIV, NFRONT, IEND_BLR,
     &            A, LA, POSELT, .TRUE.)
          ENDIF
        END DO 
        NPIV   = IW(IOLDPS+1+XSIZE)
        IF (K263.NE.0) THEN
          NELIM = IEND_BLR - NPIV
          CALL SMUMPS_SEND_FACTORED_BLOCK( COMM_LOAD, ASS_IRECV, 
     &             N, INODE, FPERE, IW, LIW, 
     &             IOLDPS, POSELT, A, LA, NFRONT,
     &             IBEG_BLR, NPIV, IPIV, NASS,LASTBL, idummy, 
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
        ENDIF
        IF (.NOT. LR_ACTIVATED) THEN
          IF (IEND_BLR .LT. NASS) THEN
              CALL SMUMPS_FAC_SQ(IBEG_BLR, 
     &            IEND_BLR,
     &            NPIV, NFRONT, NASS, A, LA,
     &            POSELT, .TRUE.
     &            )
          ENDIF
        ENDIF
        IF (KEEP(201).EQ.1) THEN 
          STRAT            = STRAT_TRY_WRITE
          MonBloc%LastPiv  = NPIV
          LAST_CALL= .FALSE.
          CALL SMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_BOTH_LU, 
     &           A(POSELT), LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
             IF (IFLAG_OOC < 0 ) THEN
                IFLAG=IFLAG_OOC
                GOTO 490
             ENDIF
        ENDIF
      END DO 
      IF (KEEP(201).EQ.1) THEN 
          STRAT        = STRAT_WRITE_MAX   
          MonBloc%Last = .TRUE.
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
            GOTO 490
          ENDIF
          CALL SMUMPS_OOC_PP_TRYRELEASE_SPACE (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
      GOTO 500
 490  CONTINUE
      CALL SMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
 500  CONTINUE
      DEALLOCATE( IPIV )
      RETURN
      END SUBROUTINE SMUMPS_FACTO_NIV2
      END MODULE SMUMPS_FACTO_NIV2_M
