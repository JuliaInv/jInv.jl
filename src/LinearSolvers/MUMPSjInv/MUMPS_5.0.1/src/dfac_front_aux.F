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
      MODULE DMUMPS_FAC_FRONT_AUX_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC_H(NFRONT,NASS,IW,LIW,A,LA,
     &   INOPV,NOFFW,IOLDPS,POSELT,UU,SEUIL,KEEP, DKEEP,
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &     PP_LastPIVRPTRFilled_U)
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER NFRONT,NASS,LIW,INOPV
      INTEGER(8) :: LA
      INTEGER KEEP(500)
      DOUBLE PRECISION    DKEEP(130)
      DOUBLE PRECISION UU, SEUIL
      DOUBLE PRECISION A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION AMROW
      DOUBLE PRECISION RMAX
      DOUBLE PRECISION  SWOP
      INTEGER(8) :: APOS, POSELT
      INTEGER(8) :: J1, J2, J3_8, JJ, IDIAG
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS
      INTEGER NOFFW,NPIV,IPIV
      INTEGER J, J3
      INTEGER NPIVP1,JMAX,ISW,ISWPS1
      INTEGER ISWPS2,KSW,XSIZE
      INTEGER I_PIVRPTR_L, I_PIVR_L, NBPANELS_L
      INTEGER I_PIVRPTR_U, I_PIVR_U, NBPANELS_U
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &        PP_LastPIVRPTRFilled_U
      INTEGER DMUMPS_IXAMAX
      INCLUDE 'mumps_headers.h'
      INTRINSIC max
      DOUBLE PRECISION, PARAMETER :: RZERO = 0.0D0
        NFRONT8 = int(NFRONT,8)
        INOPV   = 0
        XSIZE   = KEEP(IXSZ)
        NPIV    = IW(IOLDPS+1+XSIZE)
        NPIVP1  = NPIV + 1
        IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1) THEN
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR_L, I_PIVR_L, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)
     &              +KEEP(IXSZ),
     &       IW, LIW)
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_U, NBPANELS_U, 
     &       I_PIVRPTR_U, I_PIVR_U, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
        ENDIF
          DO 460 IPIV=NPIVP1,NASS
            APOS = POSELT + NFRONT8*int(NPIV,8) + int(IPIV-1,8)
            JMAX = 1
            AMROW = RZERO
            J1    = APOS
            J3    = NASS -NPIV
            JMAX  = DMUMPS_IXAMAX(J3,A(J1),NFRONT)
            JJ    = J1 + int(JMAX-1,8)*NFRONT8
            AMROW = abs(A(JJ))
            RMAX  = AMROW
            J1    = APOS +  int(NASS-NPIV,8) * NFRONT8
            J3 = NFRONT - NASS - KEEP(253)
            IF (J3.EQ.0) GOTO 370
            DO 360 J=1,J3
              RMAX = max(abs(A(J1)),RMAX)
              J1 = J1 + NFRONT8
  360       CONTINUE
  370       IF (RMAX.EQ.RZERO) GO TO 460
            IDIAG = APOS + int(IPIV - NPIVP1,8)*NFRONT8
            IF (abs(A(IDIAG)).GE.max(UU*RMAX,SEUIL)) THEN
              JMAX = IPIV - NPIV
              GO TO 380
            ENDIF
            IF (AMROW.LT.max(UU*RMAX,SEUIL)) GO TO 460
            NOFFW = NOFFW + 1
  380       CONTINUE
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER(
     &             A(APOS + int(JMAX - 1,8) * NFRONT8 ),
     &             DKEEP(6),
     &             KEEP(259) )
            ENDIF
            IF (IPIV.EQ.NPIVP1) GO TO 400
            KEEP(260)=-KEEP(260)
            J1   = POSELT + int(NPIV,8)
            J3_8 = POSELT + int(IPIV-1,8)
            DO 390 J= 1,NFRONT
              SWOP  = A(J1)
              A(J1) = A(J3_8)
              A(J3_8) = SWOP
              J1 = J1 + NFRONT8
              J3_8 = J3_8 + NFRONT8
  390       CONTINUE
            ISWPS1 = IOLDPS + 5 + NPIVP1 + NFRONT + XSIZE
            ISWPS2 = IOLDPS + 5 + IPIV + NFRONT + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            KEEP(260)=-KEEP(260)
            J1 = POSELT + int(NPIV,8) * NFRONT8
            J2 = POSELT + int(NPIV + JMAX - 1,8) * NFRONT8
            DO 410 KSW=1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + 1_8
              J2 = J2 + 1_8
  410       CONTINUE
            ISWPS1 = IOLDPS + 5 + NPIV + 1 + XSIZE
            ISWPS2 = IOLDPS + 5 + NPIV + JMAX + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460     CONTINUE
       INOPV = 1
       GOTO 430
  420 CONTINUE
              IF (KEEP(201).EQ.1) THEN
                IF (KEEP(251).EQ.0) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, NPIV+JMAX,
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
                ENDIF
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, IPIV,
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
              ENDIF
 430  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_H
      SUBROUTINE DMUMPS_FAC_M(IBEG_BLOCK,
     &     NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &     IOLDPS,POSELT,IFINB,LKJIB,LKJIT,XSIZE)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LIW,INODE,IFINB,LKJIB,IBEG_BLOCK
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    VALPIV
      INTEGER(8) :: APOS, POSELT, UUPOS, LPOS
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS
      INTEGER LKJIT, XSIZE
      DOUBLE PRECISION ONE, ALPHA
      INTEGER NPIV,JROW2
      INTEGER NEL2,NPIVP1,KROW,NEL
      INCLUDE 'mumps_headers.h'
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NFRONT8= int(NFRONT,8)
        NPIV   = IW(IOLDPS+1+XSIZE)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        IFINB  = 0
        IF (IW(IOLDPS+3+XSIZE).LE.0) THEN
          IF (NASS.LT.LKJIT) THEN
           IW(IOLDPS+3+XSIZE) = NASS
          ELSE
           IW(IOLDPS+3+XSIZE) = min0(NASS,LKJIB)
          ENDIF
        ENDIF
        JROW2 = IW(IOLDPS+3+XSIZE)
        NEL2   = JROW2 - NPIVP1
        IF (NEL2.EQ.0) THEN
         IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
         ELSE
          IFINB        = 1
          IW(IOLDPS+3+XSIZE) = min0(JROW2+LKJIB,NASS)
          IBEG_BLOCK = NPIVP1+1
         ENDIF
        ELSE
         APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + NFRONT8
         DO 541 KROW = 1,NEL2
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT8
 541     CONTINUE
         LPOS   = APOS + NFRONT8
         UUPOS  = APOS + 1_8
         CALL dger(NEL,NEL2,ALPHA,A(UUPOS),1,A(LPOS),NFRONT,
     &              A(LPOS+1_8),NFRONT)
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_FAC_M
      SUBROUTINE DMUMPS_FAC_N(NFRONT,NASS,IW,LIW,A,LA,
     &       IOLDPS,POSELT,IFINB,XSIZE)
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      INTEGER NFRONT,NASS,LIW,IFINB
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    ALPHA,VALPIV
      INTEGER(8) :: APOS, POSELT, UUPOS, LPOS, IRWPOS
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS,NPIV,KROW, XSIZE
      INTEGER NEL,ICOL,NEL2
      INTEGER NPIVP1
      DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0
        NFRONT8=int(NFRONT,8)
        NPIV   = IW(IOLDPS+1+XSIZE)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        NEL2   = NASS - NPIVP1
        IFINB  = 0
        IF (NPIVP1.EQ.NASS) IFINB = 1
        APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
        VALPIV = ONE/A(APOS)
        LPOS   = APOS + NFRONT8
        DO 541 KROW = 1,NEL
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT8
 541    CONTINUE
        LPOS   = APOS + NFRONT8
        UUPOS  = APOS + 1_8
        DO 440 ICOL = 1,NEL
             IRWPOS  = LPOS + 1_8
             ALPHA   = -A(LPOS)
             CALL daxpy(NEL2,ALPHA,A(UUPOS),1,A(IRWPOS),1)
             LPOS    = LPOS + NFRONT8
  440   CONTINUE
        RETURN
        END SUBROUTINE DMUMPS_FAC_N
      SUBROUTINE DMUMPS_FAC_P(A,LA,NFRONT,
     &       NPIV,NASS,POSELT)
      IMPLICIT NONE
      INTEGER(8) :: LA,POSELT
      DOUBLE PRECISION    A(LA)
      INTEGER NFRONT, NPIV, NASS
      INTEGER(8) :: LPOS, LPOS1, LPOS2
      INTEGER NEL1,NEL11
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        LPOS2  = POSELT + int(NASS,8)*int(NFRONT,8)
        CALL dtrsm('L','L','N','N',NPIV,NEL1,ONE,A(POSELT),NFRONT,
     &              A(LPOS2),NFRONT)
        LPOS   = LPOS2 + int(NPIV,8)
        LPOS1  = POSELT + int(NPIV,8)
        CALL dgemm('N','N',NEL11,NEL1,NPIV,ALPHA,A(LPOS1),
     &          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
      RETURN
      END SUBROUTINE DMUMPS_FAC_P
      SUBROUTINE DMUMPS_FAC_P_PANEL(A,LAFAC,NFRONT,
     &      NPIV,NASS, IW, LIWFAC,
     &      MonBloc, TYPEFile, MYID, KEEP8,
     &      STRAT, IFLAG_OOC,
     &      LNextPiv2beWritten, UNextPiv2beWritten)
      USE DMUMPS_OOC   
      IMPLICIT NONE
      INTEGER NFRONT, NPIV, NASS
      INTEGER(8) :: LAFAC
      INTEGER  LIWFAC, TYPEFile, MYID, IFLAG_OOC,
     &      LNextPiv2beWritten, UNextPiv2beWritten, STRAT
      DOUBLE PRECISION  A(LAFAC)
      INTEGER  IW(LIWFAC)
      INTEGER(8) KEEP8(150)
      TYPE(IO_BLOCK) :: MonBloc 
      INTEGER(8) :: LPOS2,LPOS1,LPOS
      INTEGER NEL1,NEL11
      DOUBLE PRECISION ALPHA, ONE
      LOGICAL LAST_CALL
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        LPOS2  = 1_8 + int(NASS,8) * int(NFRONT,8)
        CALL dtrsm('L','L','N','N',NPIV,NEL1,ONE,A(1),NFRONT,
     &              A(LPOS2),NFRONT)
        LAST_CALL=.FALSE.
           CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEFile, 
     &           A, LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW, LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
        LPOS   = LPOS2 + int(NPIV,8)
        LPOS1  = int(1 + NPIV,8)
        CALL dgemm('N','N',NEL11,NEL1,NPIV,ALPHA,A(LPOS1),
     &          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
        RETURN
        END SUBROUTINE DMUMPS_FAC_P_PANEL
      SUBROUTINE DMUMPS_FAC_T(A,LA,NPIVB,NFRONT,
     &                             NPIV,NASS,POSELT)
      IMPLICIT NONE
      INTEGER NPIVB,NASS
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER(8) :: APOS, POSELT
      INTEGER NFRONT, NPIV, NASSL
      INTEGER(8) :: LPOS, LPOS1, LPOS2
      INTEGER NEL1, NEL11, NPIVE
      DOUBLE PRECISION    ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        NPIVE  = NPIV - NPIVB
        NASSL  = NASS - NPIVB
        APOS   = POSELT + int(NPIVB,8)*int(NFRONT,8)
     &                  + int(NPIVB,8)
        LPOS2  = APOS + int(NASSL,8)
        CALL dtrsm('R','U','N','U',NEL1,NPIVE,ONE,A(APOS),NFRONT,
     &              A(LPOS2),NFRONT)
        LPOS   = LPOS2 + int(NFRONT,8)*int(NPIVE,8)
        LPOS1  = APOS  + int(NFRONT,8)*int(NPIVE,8)
        CALL dgemm('N','N',NEL1,NEL11,NPIVE,ALPHA,A(LPOS2),
     &          NFRONT,A(LPOS1),NFRONT,ONE,A(LPOS),NFRONT)
        RETURN
        END SUBROUTINE DMUMPS_FAC_T
      SUBROUTINE DMUMPS_FAC_SQ(IBEG_BLOCK, IEND_BLOCK, NPIV,
     &    NFRONT, LAST_ROW, A, LA, POSELT, CALL_GEMM
     &    )
      IMPLICIT NONE
      INTEGER, intent(in)     :: IBEG_BLOCK, IEND_BLOCK
      INTEGER, intent(in)     :: NPIV, NFRONT, LAST_ROW
      INTEGER(8), intent(in)  :: LA
      DOUBLE PRECISION, intent(inout)  :: A(LA)
      INTEGER(8), intent(in)  :: POSELT 
      LOGICAL, intent(in)     :: CALL_GEMM
      INTEGER(8) :: NFRONT8
      INTEGER(8) :: LPOS, LPOS1, LPOS2, POSELT_LOCAL
      INTEGER NELIM, LKJIW, NEL1, NEL11 
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      NFRONT8= int(NFRONT,8)
      NELIM  = IEND_BLOCK - NPIV
      NEL1   = LAST_ROW - IEND_BLOCK
      IF ( NEL1 < 0 ) THEN
        WRITE(*,*)
     &  "Internal error 1 in DMUMPS_FAC_SQ,IEND_BLOCK>LAST_ROW",
     &  IEND_BLOCK, LAST_ROW
        CALL MUMPS_ABORT()
      ENDIF
      LKJIW  = NPIV - IBEG_BLOCK + 1
      NEL11  = NFRONT - NPIV
      IF ((NEL1.NE.0).AND.(LKJIW.NE.0)) THEN
          LPOS2  = POSELT + int(IEND_BLOCK,8)*NFRONT8 +
     &             int(IBEG_BLOCK - 1,8)
          POSELT_LOCAL = POSELT + 
     &           int(IBEG_BLOCK-1,8)*NFRONT8 + int(IBEG_BLOCK - 1,8)
          CALL dtrsm('L','L','N','N',LKJIW,NEL1,ONE,
     &               A(POSELT_LOCAL),NFRONT,
     &               A(LPOS2),NFRONT)
          IF (CALL_GEMM) THEN
           LPOS   = LPOS2 + int(LKJIW,8)
           LPOS1  = POSELT_LOCAL + int(LKJIW,8)
           CALL dgemm('N','N',NEL11,NEL1,LKJIW,ALPHA,A(LPOS1),
     &          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
          END IF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_FAC_SQ
      SUBROUTINE DMUMPS_FAC_MQ(IBEG_BLOCK,IEND_BLOCK,
     &     NFRONT, NASS, NPIV, A, LA, POSELT, IFINB)
      IMPLICIT NONE
      INTEGER, intent(in)    :: IBEG_BLOCK, IEND_BLOCK, NFRONT, 
     &                          NASS, NPIV
      INTEGER, intent(out)   ::  IFINB
      INTEGER(8), intent(in) :: LA, POSELT
      DOUBLE PRECISION, intent(inout) :: A(LA)
      DOUBLE PRECISION    :: VALPIV
      INTEGER(8) :: APOS,  UUPOS, LPOS
      INTEGER(8) :: NFRONT8
      DOUBLE PRECISION    :: ONE, ALPHA
      INTEGER    :: NEL2,NPIVP1,KROW,NEL
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NFRONT8= int(NFRONT,8)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        IFINB  = 0
        NEL2   = IEND_BLOCK - NPIVP1
        IF (NEL2.EQ.0) THEN
         IF (IEND_BLOCK.EQ.NASS) THEN
          IFINB        = -1
         ELSE
          IFINB        = 1
         ENDIF
        ELSE
         APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + NFRONT8
         DO 541 KROW = 1,NEL2
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT8
 541     CONTINUE
         LPOS   = APOS + NFRONT8
         UUPOS  = APOS + 1_8
#if defined(MUMPS_USE_BLAS2)
         CALL dger(NEL,NEL2,ALPHA,A(UUPOS),1,A(LPOS),NFRONT,
     &              A(LPOS+1_8),NFRONT)
#else
         CALL dgemm('N','N',NEL,NEL2,1,ALPHA,A(UUPOS),NEL,
     &               A(LPOS),NFRONT,ONE,A(LPOS+1_8),NFRONT)
#endif
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_FAC_MQ
      SUBROUTINE DMUMPS_FAC_FR_UPDATE_CBROWS( INODE,
     &     NFRONT, NASS, A, LA, LAFAC, POSELT, IW, LIW, IOLDPS,
     &     MonBloc, MYID, NOFFW, LIWFAC,
     &     PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &     LNextPiv2beWritten, UNextPiv2beWritten, 
     &     PP_LastPIVRPTRFilled_L, PP_LastPIVRPTRFilled_U,
     &     
     &     XSIZE, SEUIL, UU, DKEEP, KEEP8, KEEP, IFLAG)
      USE DMUMPS_OOC 
      IMPLICIT NONE
      INTEGER, intent(in)    :: INODE, NFRONT, NASS,
     &                          LIW, MYID, XSIZE, IOLDPS, LIWFAC
      INTEGER(8), intent(in) :: LA, POSELT
      INTEGER, intent(inout) :: NOFFW,
     &                   PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &                   LNextPiv2beWritten, UNextPiv2beWritten,
     &                   PP_LastPIVRPTRFilled_L, PP_LastPIVRPTRFilled_U,
     &                   IFLAG
      INTEGER, intent(inout) :: IW(LIW)
      DOUBLE PRECISION, intent(inout) :: A(LA)
      DOUBLE PRECISION, intent(in)       :: SEUIL, UU, DKEEP(130)
      INTEGER, intent(in)    :: KEEP( 500 ) 
      INTEGER(8), intent(inout) :: LAFAC
      INTEGER(8)             :: KEEP8(150)
      TYPE(IO_BLOCK), intent(inout)  :: MonBloc
      INTEGER  :: NPIV, NEL1, STRAT, TYPEFile, IFLAG_OOC, 
     &            IBEG_BLOCK, IFINB, INOPV
      NPIV   = IW(IOLDPS+1+XSIZE)
      NEL1   = NFRONT - NASS
      IF ((NPIV.GT.0).AND.(NEL1.GT.0)) THEN
        IF (KEEP(201).EQ.1) THEN 
         STRAT          = STRAT_TRY_WRITE
         TYPEFile     = TYPEF_BOTH_LU
         MonBloc%LastPiv= NPIV
         CALL DMUMPS_FAC_P_PANEL(A(POSELT), LAFAC, NFRONT, 
     &      NPIV, NASS, IW(IOLDPS), LIWFAC, 
     &      MonBloc, TYPEFile, MYID, KEEP8,
     &      STRAT, IFLAG_OOC,
     &      LNextPiv2beWritten, UNextPiv2beWritten)
          IF (IFLAG_OOC < 0 ) IFLAG=IFLAG_OOC
        ELSE
          CALL DMUMPS_FAC_P(A,LA,NFRONT, NPIV,NASS,POSELT)
        ENDIF
       ENDIF
        NPIV   = IW(IOLDPS+1+XSIZE)
        IBEG_BLOCK = NPIV
        IF (NASS.EQ.NPIV) GOTO 500
 120    CALL DMUMPS_FAC_H(NFRONT,NASS,IW,LIW,A,LA,
     &     INOPV,NOFFW,IOLDPS,POSELT,UU,SEUIL,
     &     KEEP, DKEEP,
     &     PP_FIRST2SWAP_L,  MonBloc%LastPanelWritten_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U,  MonBloc%LastPanelWritten_U,
     &     PP_LastPIVRPTRFilled_U)
        IF (INOPV.NE.1) THEN
         CALL DMUMPS_FAC_N(NFRONT,NASS,IW,LIW,A,LA,
     &                 IOLDPS,POSELT,IFINB,XSIZE)
         IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + 1
         IF (IFINB.EQ.0) GOTO 120
        ENDIF
        NPIV   = IW(IOLDPS+1+XSIZE)
        NEL1   = NFRONT - NASS
        IF ((NPIV.LE.IBEG_BLOCK).OR.(NEL1.EQ.0)) GO TO 500
        CALL DMUMPS_FAC_T(A,LA,IBEG_BLOCK,
     &                NFRONT,NPIV,NASS,POSELT)
 500  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_FR_UPDATE_CBROWS
      SUBROUTINE DMUMPS_FAC_I(NFRONT,NASS,LAST_ROW,
     &    IBEG_BLOCK_FOR_TIPIV, IEND_BLOCK,
     &    N,INODE,IW,LIW,A,LA,
     &    INOPV,NOFFW,IFLAG,IOLDPS,POSELT,UU,SEUIL,KEEP,KEEP8,
     &     DKEEP,PIVNUL_LIST,LPN_LIST,
     &
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &     PP_LastPIVRPTRFilled_U,
     &     TIPIV)    
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER, intent(in)    :: IBEG_BLOCK_FOR_TIPIV, IEND_BLOCK
      INTEGER, intent(inout), OPTIONAL :: TIPIV(:)
      INTEGER(8), intent(in) :: LA
      DOUBLE PRECISION, intent(inout) :: A(LA)
      INTEGER, intent(in)    :: NFRONT,NASS,N,LIW,INODE,LAST_ROW
      INTEGER, intent(inout) :: IFLAG,INOPV,NOFFW
      DOUBLE PRECISION, intent(in)       :: UU, SEUIL
      INTEGER, intent(inout) :: IW(LIW)
      INTEGER, intent(in)    :: IOLDPS
      INTEGER(8), intent(in) :: POSELT
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in)    :: LPN_LIST
      INTEGER, intent(inout) :: PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION    DKEEP(130)
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &        PP_LastPIVRPTRFilled_U
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION SWOP
      INTEGER XSIZE
      INTEGER(8) :: APOS, IDIAG
      INTEGER(8) :: J1, J2, JJ, J3
      INTEGER(8) :: NFRONT8
      INTEGER ILOC
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      DOUBLE PRECISION RZERO, RMAX, AMROW, ONE
      DOUBLE PRECISION PIVNUL
      DOUBLE PRECISION FIXA, CSEUIL
      INTEGER NPIV,IPIV
      INTEGER NPIVP1,JMAX,J,ISW,ISWPS1
      INTEGER ISWPS2,KSW, HF
      INTEGER DMUMPS_IXAMAX
      INTRINSIC max
      DATA RZERO /0.0D0/
      DATA ONE /1.0D0/
      INTEGER I_PIVRPTR_L, I_PIVR_L, NBPANELS_L
      INTEGER I_PIVRPTR_U, I_PIVR_U, NBPANELS_U
        PIVNUL  = DKEEP(1)
        FIXA    = DKEEP(2)
        CSEUIL  = SEUIL
        NFRONT8=int(NFRONT,8)
        XSIZE   = KEEP(IXSZ)
        NPIV    = IW(IOLDPS+1+XSIZE)
        HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
        NPIVP1  = NPIV + 1
        IF (KEEP(201).EQ.1) THEN
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR_L, I_PIVR_L, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_U, NBPANELS_U, 
     &       I_PIVRPTR_U, I_PIVR_U, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
        ENDIF
        IF ( present(TIPIV) ) THEN
          ILOC    = NPIVP1 - IBEG_BLOCK_FOR_TIPIV + 1
          TIPIV(ILOC) = ILOC
        ENDIF
        IF (INOPV .EQ. -1) THEN
           APOS = POSELT + NFRONT8*int(NPIVP1-1,8) + int(NPIV,8)
           IDIAG = APOS
           IF(abs(A(APOS)).LT.SEUIL) THEN
              IF (dble(A(APOS)) .GE. RZERO) THEN
                 A(APOS) = CSEUIL
              ELSE
                 A(APOS) = -CSEUIL
              ENDIF
              KEEP(98) = KEEP(98)+1
           ELSE IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER(A(APOS), DKEEP(6), KEEP(259))
           ENDIF
           IF (KEEP(201).EQ.1) THEN
             IF (KEEP(251).EQ.0) THEN 
               CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, NPIVP1, 
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
             ENDIF
             CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, NPIVP1, 
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
           ENDIF
           GO TO 420
        ENDIF
        INOPV   = 0
          DO 460 IPIV=NPIVP1,IEND_BLOCK
            APOS = POSELT + NFRONT8*int(IPIV-1,8) + int(NPIV,8)
            JMAX = 1
            IF (UU.GT.RZERO) GO TO 340
            IF (A(APOS).EQ.ZERO) GO TO 630 
            GO TO 380 
  340       AMROW = RZERO
            J1 = APOS
            J2 = APOS +int(- NPIV + NASS - 1,8)
             J     = NASS -NPIV
             JMAX  = DMUMPS_IXAMAX(J,A(J1),1)
             JJ    = J1 + int(JMAX - 1,8)
             AMROW = abs(A(JJ))
            RMAX = AMROW
            J1 = J2 + 1_8
            J2 = APOS +int(- NPIV + NFRONT - 1 - KEEP(253),8)
            IF (J2.LT.J1) GO TO 370
            DO 360 JJ=J1,J2
              RMAX = max(abs(A(JJ)),RMAX)
  360       CONTINUE
  370       IDIAG = APOS + int(IPIV - NPIVP1,8)
            IF ( RMAX .LE. PIVNUL ) THEN
               IF (NFRONT - KEEP(253) .EQ. NASS) THEN 
                 IF (IEND_BLOCK.NE.NASS ) THEN 
                   GOTO 460
                 ENDIF
                 J1 = IDIAG
                 J1=POSELT+int(IPIV-1,8)+NPIV*NFRONT
                 J2=POSELT+int(IPIV-1,8)+(LAST_ROW-1)*NFRONT8
               ELSE
                 J1=POSELT+int(IPIV-1,8)
                 J2=POSELT+int(IPIV-1,8)+int(LAST_ROW-1,8)*NFRONT8
               ENDIF
               DO JJ=J1, J2, NFRONT8
                 IF ( abs(A(JJ)) .GT. PIVNUL ) THEN
                   GOTO 460
                 END IF
               ENDDO
               KEEP(109) = KEEP(109)+1
               ISW = IOLDPS+HF+
     &                      IW(IOLDPS+1+XSIZE)+IPIV-NPIVP1
               PIVNUL_LIST(KEEP(109)) = IW(ISW)
               IF(dble(FIXA).GT.RZERO) THEN
                  IF(dble(A(IDIAG)) .GE. RZERO) THEN
                     A(IDIAG) = FIXA
                  ELSE
                     A(IDIAG) = -FIXA
                  ENDIF
               ELSE
                 J1 = APOS
                 J2 = APOS +int(- NPIV + NFRONT - 1 - KEEP(253),8)
                 DO JJ=J1,J2
                   A(JJ) = ZERO
                 ENDDO
                 A(IDIAG) = -FIXA
               ENDIF
               JMAX = IPIV - NPIV
               GOTO 385   
            ENDIF
            IF (abs(A(IDIAG)).GT.max(UU*RMAX,SEUIL)) THEN
               JMAX = IPIV - NPIV
               GO TO 380
            ENDIF
            IF (AMROW .LE. max(UU*RMAX,SEUIL)) GO TO 460
            NOFFW = NOFFW + 1
  380       CONTINUE
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER( A(APOS+int(JMAX-1,8)),
     &                                 DKEEP(6),
     &                                 KEEP(259))
            ENDIF
  385       CONTINUE
            IF (IPIV.EQ.NPIVP1) GO TO 400
            KEEP(260)=-KEEP(260)
            J1 = POSELT + int(NPIV,8)*NFRONT8
            J2 = J1 + NFRONT8 - 1_8
            J3 = POSELT + int(IPIV-1,8)*NFRONT8
            DO 390 JJ=J1,J2
              SWOP = A(JJ)
              A(JJ) = A(J3)
              A(J3) = SWOP
              J3 = J3 + 1_8
  390       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NPIVP1
            ISWPS2 = IOLDPS + HF - 1 + IPIV
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            KEEP(260)=-KEEP(260)
            IF ( present(TIPIV) ) THEN
              TIPIV(ILOC) = ILOC + JMAX - 1
            ENDIF
            J1 = POSELT + int(NPIV,8)
            J2 = POSELT + int(NPIV + JMAX - 1,8)
            DO 410 KSW=1,LAST_ROW
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + NFRONT8
              J2 = J2 + NFRONT8
  410       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NFRONT + NPIV + 1
            ISWPS2 = IOLDPS + HF - 1 + NFRONT + NPIV + JMAX
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420  
  460     CONTINUE
      IF (IEND_BLOCK.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 430
  630 CONTINUE
      IFLAG = -10
      GOTO 430
  420 CONTINUE
              IF (KEEP(201).EQ.1) THEN
                IF (KEEP(251).EQ.0) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, IPIV, 
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
                ENDIF
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, NPIV+JMAX, 
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
              ENDIF
  430 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_I
      SUBROUTINE DMUMPS_FAC_I_LDLT 
     &   (NFRONT,NASS,INODE,IEND_BLOCK,
     &    IW,LIW, A,LA, INOPV,
     &    NNEG,
     &    IFLAG,IOLDPS,POSELT,UU, SEUIL,KEEP,KEEP8,PIVSIZ,
     &     DKEEP,PIVNUL_LIST,LPN_LIST, XSIZE,
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk,
     &     PP_LastPIVRPTRIndexFilled,MAXFROMM,IS_MAXFROMM_AVAIL)
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER(8) :: POSELT, LA
      INTEGER NFRONT,NASS,LIW,INODE,IFLAG,INOPV,
     &        IOLDPS, NNEG
      INTEGER, intent(in) :: IEND_BLOCK
      INTEGER PIVSIZ,LPIV, XSIZE
      DOUBLE PRECISION A(LA) 
      DOUBLE PRECISION UU, UULOC, SEUIL
      INTEGER IW(LIW)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(130)
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk
      INTEGER PP_LastPIVRPTRIndexFilled
      DOUBLE PRECISION, intent(in) :: MAXFROMM
      LOGICAL, intent(inout) :: IS_MAXFROMM_AVAIL
      include 'mpif.h'
      INTEGER (8) :: POSPV1,POSPV2,OFFDAG,APOSJ
      INTEGER JMAX
      DOUBLE PRECISION RMAX,AMAX,TMAX,TOL
      DOUBLE PRECISION MAXPIV
      DOUBLE PRECISION PIVNUL
      DOUBLE PRECISION FIXA, CSEUIL
      DOUBLE PRECISION PIVOT,DETPIV
      PARAMETER(TOL = 1.0D-20)
      INCLUDE 'mumps_headers.h'
      INTEGER :: J
      INTEGER(8) :: APOS, J1, J2, JJ, NFRONT8, KK, J1_ini, JJ_ini
      INTEGER    :: LDA
      INTEGER(8) :: LDA8
      INTEGER NPIV,IPIV
      INTEGER NPIVP1,K
      INTRINSIC max
      DOUBLE PRECISION ZERO, ONE
      PARAMETER( ZERO = 0.0D0 )
      PARAMETER( ONE = 1.0D0 )
      DOUBLE PRECISION RZERO,RONE
      PARAMETER(RZERO=0.0D0, RONE=1.0D0)
      LOGICAL OMP_FLAG
      INTEGER I_PIVRPTR, I_PIVR, NBPANELS_L
      PIVNUL = DKEEP(1)
      FIXA   = DKEEP(2)
      CSEUIL = SEUIL
      LDA     = NFRONT
      LDA8    = int(LDA,8)
      NFRONT8 = int(NFRONT,8)
      IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1) THEN
             CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR, I_PIVR, IOLDPS+2*NFRONT+6+KEEP(IXSZ),
     &       IW, LIW)
      ENDIF
      UULOC = UU
      PIVSIZ = 1
      NPIV    = IW(IOLDPS+1+XSIZE)
      NPIVP1  = NPIV + 1
      IF(INOPV .EQ. -1) THEN
         APOS = POSELT + (LDA8+1_8) * int(NPIV,8)
         POSPV1 = APOS
         IF(abs(A(APOS)).LT.SEUIL) THEN
            IF(dble(A(APOS)) .GE. RZERO) THEN
               A(APOS) = CSEUIL
            ELSE
               A(APOS) = -CSEUIL
               NNEG = NNEG+1
            ENDIF
            KEEP(98) = KEEP(98)+1
         ELSE IF (KEEP(258) .NE. 0) THEN
            CALL DMUMPS_UPDATEDETER( A(APOS), DKEEP(6), KEEP(259) )
         ENDIF
              IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR), NBPANELS_L,
     &               IW(I_PIVR), NASS, NPIVP1, NPIVP1,
     &               PP_LastPanelonDisk,
     &               PP_LastPIVRPTRIndexFilled)
              ENDIF
         GO TO 420
      ENDIF
      INOPV   = 0
      DO 460 IPIV=NPIVP1,IEND_BLOCK
         APOS = POSELT + LDA8*int(IPIV-1,8) + int(NPIV,8)
         POSPV1 = APOS + int(IPIV - NPIVP1,8)
         PIVOT = A(POSPV1)
         IF (UULOC.EQ.RZERO) THEN 
            IF (abs(A(APOS)).EQ.RZERO) GO TO 630
            IF (A(APOS).LT.RZERO) NNEG = NNEG+1
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER(A(APOS), DKEEP(6), KEEP(259))
            ENDIF
            GO TO 420
         ENDIF
         AMAX = RZERO
         JMAX = 0
         IF ( IS_MAXFROMM_AVAIL ) THEN
            IF ( MAXFROMM > PIVNUL .AND. abs(PIVOT) > TOL) THEN
               IF ( abs(PIVOT) .GT. max(UULOC*MAXFROMM,SEUIL) ) THEN
                 IF (PIVOT .LT. RZERO) NNEG = NNEG+1
                 IF (KEEP(258) .NE. 0) THEN
                   CALL DMUMPS_UPDATEDETER(PIVOT, DKEEP(6), KEEP(259))
                 ENDIF
                 GOTO 415
               ENDIF
            ENDIF
            IS_MAXFROMM_AVAIL = .FALSE.
         ENDIF
         J1 = APOS
         J2 = POSPV1 - 1_8
         DO JJ=J1,J2
            IF(abs(A(JJ)) .GT. AMAX) THEN
               AMAX = abs(A(JJ))
               JMAX = IPIV - int(POSPV1-JJ)
            ENDIF
         ENDDO
         J1 = POSPV1 + LDA8
         DO J=1, IEND_BLOCK - IPIV
            IF(abs(A(J1)) .GT. AMAX) THEN
               AMAX = abs(A(J1))
               JMAX = IPIV + J
            ENDIF
            J1 = J1 + LDA8
         ENDDO
           RMAX = RZERO
           J1_ini = J1
           IF ( (NFRONT - KEEP(253) - IEND_BLOCK).GE.300 ) THEN
             OMP_FLAG = .TRUE.
           ELSE
             OMP_FLAG = .FALSE.
           ENDIF
!$OMP PARALLEL DO PRIVATE(J1) REDUCTION(max:RMAX) IF(OMP_FLAG)
           DO J=1, NFRONT - KEEP(253) - IEND_BLOCK
              J1 = J1_ini + int(J-1,8) * LDA8
              RMAX = max(abs(A(J1)),RMAX)
           ENDDO
!$OMP END PARALLEL DO
         IF (max(AMAX,RMAX,abs(PIVOT)).LE.PIVNUL) THEN
            KEEP(109) = KEEP(109)+1
            PIVNUL_LIST(KEEP(109)) = -1
            IF(dble(FIXA).GT.RZERO) THEN
               IF(dble(PIVOT) .GE. RZERO) THEN
                  A(POSPV1) = FIXA
               ELSE
                  A(POSPV1) = -FIXA
               ENDIF
            ELSE
               J1 = APOS
               J2 = POSPV1 - 1_8
               DO JJ=J1,J2
                  A(JJ) = ZERO
               ENDDO
               J1 = POSPV1 + LDA8
               DO J=1, IEND_BLOCK - IPIV
                  A(J1) = ZERO
                  J1 = J1 + LDA8
               ENDDO
               DO J=1,NFRONT - IEND_BLOCK
                  A(J1) = ZERO
                  J1 = J1 + LDA8
               ENDDO
               A(POSPV1) = ONE
            ENDIF
            PIVOT = A(POSPV1)
            GO TO 415
         ENDIF
         IF ((KEEP(19).EQ.0).AND.(KEEP(110).EQ.0)) THEN
           IF (max(AMAX,RMAX,abs(PIVOT)).LE.TOL) THEN
            IF(SEUIL .GT. epsilon(SEUIL)) THEN
               IF(dble(PIVOT) .GE. RZERO) THEN
                  A(POSPV1) = CSEUIL
               ELSE
                  A(POSPV1) = -CSEUIL
                  NNEG = NNEG+1
               ENDIF
               PIVOT = A(POSPV1)
               KEEP(98) = KEEP(98)+1
               GO TO 415
            ENDIF
           ENDIF
         ENDIF
         IF (max(AMAX,abs(PIVOT)).LE.TOL) then
           GO TO 460
         ENDIF
         IF (abs(PIVOT).GT.max(UULOC*max(RMAX,AMAX),SEUIL)) THEN
               IF (PIVOT .LT. ZERO) NNEG = NNEG+1
               IF (KEEP(258) .NE.0 ) THEN
                 CALL DMUMPS_UPDATEDETER(PIVOT, DKEEP(6), KEEP(259))
               ENDIF
               GO TO 415
         END IF
         IF (  (AMAX.LE.TOL).OR.
     &   (
     &   (KEEP(19).NE.0).AND.(max(AMAX,RMAX,abs(PIVOT)).LE.SEUIL)
     &   ) 
     &      )
     &       THEN
           GO TO 460
         ENDIF
         IF (RMAX.LT.AMAX) THEN
               J1 = APOS
               J2 = POSPV1 - 1_8
               DO JJ=J1,J2
                  IF(int(POSPV1-JJ) .NE. IPIV-JMAX) THEN
                     RMAX = max(RMAX,abs(A(JJ)))
                  ENDIF
               ENDDO
               J1 = POSPV1 + LDA8
               DO J=1,NASS-IPIV
                  IF(IPIV+J .NE. JMAX) THEN
                     RMAX = max(abs(A(J1)),RMAX)
                  ENDIF
                  J1 = J1 + LDA8
               ENDDO
           ENDIF
           APOSJ = POSELT + int(JMAX-1,8)*LDA8 + int(NPIV,8)
           POSPV2 = APOSJ + int(JMAX - NPIVP1,8)
           IF (IPIV.LT.JMAX) THEN
              OFFDAG = APOSJ + int(IPIV - NPIVP1,8)
           ELSE
              OFFDAG = APOS + int(JMAX - NPIVP1,8)
           END IF
           TMAX = RZERO
           IF(JMAX .LT. IPIV) THEN
              JJ_ini = POSPV2
              OMP_FLAG = (NFRONT-JMAX-KEEP(253). GE. 300)
!$OMP PARALLEL DO IF (OMP_FLAG)
!$OMP& PRIVATE(JJ) REDUCTION(max:TMAX)
              DO K = 1, NFRONT - JMAX - KEEP(253)
                 JJ = JJ_ini+ int(K,8)*NFRONT8
                 IF (JMAX+K.NE.IPIV) THEN
                    TMAX=max(TMAX,abs(A(JJ)))
                 ENDIF
              ENDDO
!$OMP END PARALLEL DO
              DO KK =  APOSJ, POSPV2-1_8
                 TMAX = max(TMAX,abs(A(KK)))
              ENDDO
           ELSE
              JJ_ini = POSPV2
              OMP_FLAG = (NFRONT-JMAX-KEEP(253). GE. 300)
!$OMP PARALLEL DO PRIVATE(JJ) REDUCTION(max:TMAX) IF(OMP_FLAG)
              DO K = 1, NFRONT-JMAX-KEEP(253)
                 JJ = JJ_ini + int(K,8)*NFRONT8
                 TMAX=max(TMAX,abs(A(JJ)))
              ENDDO
!$OMP END PARALLEL DO
              DO KK =  APOSJ, POSPV2 - 1_8
                 IF (KK.NE.OFFDAG) THEN
                    TMAX = max(TMAX,abs(A(KK)))
                 ENDIF
              ENDDO
           ENDIF
           DETPIV = A(POSPV1)*A(POSPV2) - A(OFFDAG)**2
           IF (SEUIL.GT.RZERO) THEN
                IF (sqrt(abs(DETPIV)) .LE. SEUIL ) THEN
                  GOTO 460
                ENDIF
           ENDIF
           MAXPIV = max(abs(A(POSPV1)),abs(A(POSPV2)))
           IF (MAXPIV.EQ.RZERO) MAXPIV = RONE
           IF (abs(DETPIV)/MAXPIV.LE.TOL) THEN
             GO TO 460
           ENDIF
           IF ((abs(A(POSPV2))*RMAX+AMAX*TMAX)*UULOC.GT.
     &          abs(DETPIV)) THEN
             GO TO 460
           ENDIF
           IF ((abs(A(POSPV1))*TMAX+AMAX*RMAX)*UULOC.GT.
     &          abs(DETPIV)) THEN
             GO TO 460
           ENDIF
           IF (KEEP(258) .NE.0 ) THEN
             CALL DMUMPS_UPDATEDETER(DETPIV, DKEEP(6), KEEP(259))
           ENDIF
           PIVSIZ = 2
           KEEP(103) = KEEP(103)+1
           IF(DETPIV .LT. RZERO) THEN
             NNEG = NNEG+1
           ELSE IF(A(POSPV2) .LT. RZERO) THEN
             NNEG = NNEG+2
           ENDIF
 415       CONTINUE
           DO K=1,PIVSIZ
              IF (PIVSIZ .EQ. 2) THEN
                IF (K==1) THEN
                  LPIV = min(IPIV,JMAX)
                ELSE
                  LPIV   = max(IPIV,JMAX)
                ENDIF
              ELSE
                LPIV = IPIV
              ENDIF
              IF (LPIV.EQ.NPIVP1) THEN
                 GOTO 416
              ENDIF
              CALL DMUMPS_SWAP( A, LA, IW, LIW,
     &             IOLDPS, NPIVP1, LPIV, POSELT, NASS,
     &             LDA, NFRONT, 1, KEEP(219), KEEP(50),
     &             KEEP(IXSZ)) 
 416          CONTINUE
              IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR), NBPANELS_L,
     &               IW(I_PIVR), NASS, NPIVP1, LPIV, PP_LastPanelonDisk,
     &               PP_LastPIVRPTRIndexFilled)
              ENDIF
              NPIVP1 = NPIVP1 + 1
           ENDDO
           IF(PIVSIZ .EQ. 2) THEN
              A(POSELT+(LDA8+1_8)*int(NPIV,8)+1_8) = DETPIV
           ENDIF
           GOTO 420
  460   CONTINUE
      IF (IEND_BLOCK.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      PIVSIZ = 0
      IFLAG = -10
  420 CONTINUE
      IS_MAXFROMM_AVAIL = .FALSE.
      RETURN
      END SUBROUTINE DMUMPS_FAC_I_LDLT
      SUBROUTINE DMUMPS_FAC_MQ_LDLT(IEND_BLOCK,
     &     NFRONT,NASS,NPIV,INODE,
     &     A,LA,LDA, POSTPONE_COL_UPDATE,
     &     POSELT,IFINB,PIVSIZ,
     &     MAXFROMM, IS_MAXFROMM_AVAIL, IS_MAX_USEFUL,
     &     KEEP253)
      IMPLICIT NONE
      INTEGER, intent(out):: IFINB
      INTEGER, intent(in) :: INODE, NFRONT, NASS, NPIV
      INTEGER, intent(in) :: IEND_BLOCK
      INTEGER, intent(in) :: LDA
      INTEGER(8), intent(in) :: LA
      DOUBLE PRECISION, intent(inout) :: A(LA)
      LOGICAL, intent(in)    :: POSTPONE_COL_UPDATE
      INTEGER(8) :: POSELT
      DOUBLE PRECISION, intent(out) :: MAXFROMM
      LOGICAL, intent(out) :: IS_MAXFROMM_AVAIL
      LOGICAL, intent(in) :: IS_MAX_USEFUL
      INTEGER, INTENT(in) :: KEEP253
      DOUBLE PRECISION    VALPIV
      DOUBLE PRECISION :: MAXFROMMTMP
      INTEGER  NCB1
      INTEGER(8) :: NFRONT8
      INTEGER(8) :: LDA8
      INTEGER(8) :: K1POS
      INTEGER NEL2,NEL
      DOUBLE PRECISION ONE, ZERO
      DOUBLE PRECISION A11,A22,A12
      INTEGER(8) :: APOS, LPOS, LPOS1, LPOS2
      INTEGER(8) :: POSPV1, POSPV2
      INTEGER PIVSIZ,NPIV_NEW,J2,I
      INTEGER(8) :: OFFDAG, OFFDAG_OLD, IBEG, IEND
      INTEGER(8) :: JJ, K1, K2, IROW
      DOUBLE PRECISION SWOP,DETPIV,MULT1,MULT2
      INCLUDE 'mumps_headers.h'
      PARAMETER(ONE  = 1.0D0,
     &          ZERO = 0.0D0)
      LDA8     = int(LDA,8)
      NFRONT8  = int(NFRONT,8)
      NPIV_NEW = NPIV + PIVSIZ
      NEL      = NFRONT - NPIV_NEW
      IFINB    = 0
      IS_MAXFROMM_AVAIL = .FALSE.
      NEL2   = IEND_BLOCK - NPIV_NEW
      IF (NEL2.EQ.0) THEN
        IF (IEND_BLOCK.EQ.NASS) THEN
          IFINB        = -1
        ELSE
          IFINB        = 1
        ENDIF
      ENDIF
      IF(PIVSIZ .EQ. 1) THEN
         APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + LDA8
         MAXFROMM = 0.0D00
         IF (NEL2 > 0) THEN
           IF (.NOT. IS_MAX_USEFUL) THEN
             DO I=1, NEL2
               K1POS = LPOS + int(I-1,8)*LDA8
               A(APOS+int(I,8))=A(K1POS)
               A(K1POS) = A(K1POS) * VALPIV
               DO JJ=1_8, int(I,8)
                 A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
               ENDDO
             ENDDO
           ELSE
             IS_MAXFROMM_AVAIL = .TRUE.
             DO I=1, NEL2
               K1POS = LPOS + int(I-1,8)*LDA8
               A(APOS+int(I,8))=A(K1POS)
               A(K1POS) = A(K1POS) * VALPIV
               A(K1POS+1_8)=A(K1POS+1_8) - A(K1POS) * A(APOS+1_8)
               MAXFROMM=max( MAXFROMM,abs(A(K1POS+1_8)) )
               DO JJ = 2_8, int(I,8)
                 A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
               ENDDO
             ENDDO
           ENDIF
         ENDIF
         IF (POSTPONE_COL_UPDATE) THEN
           NCB1 = NASS   - IEND_BLOCK
         ELSE
           NCB1 = NFRONT - IEND_BLOCK
         ENDIF
         IF (.NOT. IS_MAX_USEFUL) THEN
!$OMP      PARALLEL DO PRIVATE(JJ,K1POS) IF (NCB1 > 300)
           DO I=NEL2+1, NEL2 + NCB1
             K1POS = LPOS+ int(I-1,8)*LDA8
             A(APOS+int(I,8))=A(K1POS)
             A(K1POS) = A(K1POS) * VALPIV
             DO JJ = 1_8, int(NEL2,8)
               A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
             ENDDO
           ENDDO
!$OMP      END PARALLEL DO
         ELSE
           MAXFROMMTMP=0.0D0
!$OMP      PARALLEL DO PRIVATE(JJ,K1POS)
!$OMP&     REDUCTION(max:MAXFROMMTMP) IF (NCB1 - KEEP253 > 300)
           DO I=NEL2+1, NEL2 + NCB1 - KEEP253
             K1POS = LPOS+ int(I-1,8)*LDA8
             A(APOS+int(I,8))=A(K1POS)
             A(K1POS) = A(K1POS) * VALPIV
             IF (NEL2 > 0) THEN
               A(K1POS+1_8) = A(K1POS+1_8) - A(K1POS) * A(APOS+1_8)
               MAXFROMMTMP=max(MAXFROMMTMP, abs(A(K1POS+1_8)))
               DO JJ = 2_8, int(NEL2,8)
                 A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
               ENDDO
             ENDIF
           ENDDO
!$OMP      END PARALLEL DO
           DO I = NEL2 + NCB1 - KEEP253 + 1, NEL2 + NCB1
             K1POS = LPOS+ int(I-1,8)*LDA8
             A(APOS+int(I,8))=A(K1POS)
             A(K1POS) = A(K1POS) * VALPIV
             DO JJ = 1_8, int(NEL2,8)
               A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
             ENDDO
           ENDDO
           MAXFROMM=max(MAXFROMM, MAXFROMMTMP)
         ENDIF
      ELSE
         POSPV1 = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         POSPV2 = POSPV1 + NFRONT8 + 1_8
         OFFDAG_OLD = POSPV2 - 1_8
         OFFDAG = POSPV1 + 1_8
         SWOP = A(POSPV2)
         DETPIV = A(OFFDAG)
          A22 = A(POSPV1)/DETPIV   
          A11 =  SWOP/DETPIV       
          A12 = -A(OFFDAG_OLD)/DETPIV   
          A(OFFDAG)     = A(OFFDAG_OLD)  
          A(OFFDAG_OLD) = ZERO
         LPOS1   = POSPV2 + LDA8 - 1_8
         LPOS2   = LPOS1 + 1_8
         CALL dcopy(NFRONT-NPIV_NEW, A(LPOS1), LDA, A(POSPV1+2_8), 1)
         CALL dcopy(NFRONT-NPIV_NEW, A(LPOS2), LDA, A(POSPV2+1_8), 1)
         JJ = POSPV2 + NFRONT8-1_8  
         IBEG = JJ + 2_8
         IEND = IBEG
         DO J2 = 1,NEL2
           K1 = JJ
           K2 = JJ+1_8
           MULT1 = - (A11*A(K1)+A12*A(K2))
           MULT2 = - (A12*A(K1)+A22*A(K2))
           K1 = POSPV1 + 2_8
           K2 = POSPV2 + 1_8
           DO IROW = IBEG, IEND
              A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
              K1 = K1 + 1_8
              K2 = K2 + 1_8
           ENDDO
           A( JJ       ) = -MULT1
           A( JJ + 1_8 ) = -MULT2
           IBEG = IBEG + NFRONT8
           IEND = IEND + NFRONT8 + 1_8
           JJ = JJ+NFRONT8
         ENDDO
         IEND = IEND-1_8
         DO J2 = IEND_BLOCK+1,NFRONT
           K1 = JJ
           K2 = JJ+1_8
           MULT1 = - (A11*A(K1)+A12*A(K2))
           MULT2 = - (A12*A(K1)+A22*A(K2))
           K1 = POSPV1 + 2_8
           K2 = POSPV2 + 1_8
           DO IROW = IBEG, IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1_8
               K2 = K2 + 1_8
           ENDDO
           A( JJ       ) = -MULT1
           A( JJ + 1_8 ) = -MULT2
           IBEG = IBEG + NFRONT8
           IEND = IEND + NFRONT8
           JJ   = JJ   + NFRONT8
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_FAC_MQ_LDLT
      SUBROUTINE DMUMPS_FAC_SQ_LDLT(IBEG_BLOCK,IEND_BLOCK,NPIV,
     &    NFRONT,NASS,LAST_VAR,INODE,A,LA,
     &    LDA,
     &    POSELT,
     &    POSTPONE_COL_UPDATE,
     &    KEEP,KEEP8 )
      IMPLICIT NONE
      INTEGER, intent(in) :: NPIV
      INTEGER, intent(in) :: NFRONT, NASS, IBEG_BLOCK, IEND_BLOCK
      INTEGER(8), intent(in) :: LA
      DOUBLE PRECISION, intent(inout) :: A(LA)
      INTEGER, intent(in) :: INODE
      INTEGER, intent(in) :: LAST_VAR 
      INTEGER :: KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8), intent(in) :: POSELT
      INTEGER, intent(in) :: LDA
      LOGICAL, intent(in) :: POSTPONE_COL_UPDATE
      INTEGER(8) :: LDA8
      INTEGER NPIV_BLOCK, NEL1
      INTEGER(8) :: LPOS,UPOS,APOS
      INTEGER IROW
      INTEGER Block
      INTEGER BLSIZE
      DOUBLE PRECISION ONE, ALPHA
      INCLUDE 'mumps_headers.h'
      PARAMETER (ONE=1.0D0, ALPHA=-1.0D0)
      LDA8 = int(LDA,8)
      NEL1   = LAST_VAR - IEND_BLOCK
      NPIV_BLOCK  = NPIV - IBEG_BLOCK + 1
      IF (NPIV_BLOCK.EQ.0) GO TO 500
      IF (NEL1.NE.0) THEN
        IF ( LAST_VAR - IEND_BLOCK > KEEP(7) ) THEN
          BLSIZE = KEEP(8)
        ELSE
          BLSIZE = LAST_VAR - IEND_BLOCK
        END IF
        IF ( NASS - IEND_BLOCK .GT. 0 ) THEN
#if defined(SAK_BYROW)
         DO IROW = IEND_BLOCK+1, LAST_VAR, BLSIZE
           Block = min( BLSIZE, NASS - IROW + 1 )
           LPOS = POSELT + int(IROW  - 1,8) * LDA8 +
     &                     int(IBEG_BLOCK - 1,8)
           UPOS = POSELT + int(IBEG_BLOCK - 1,8) * LDA8 +
     &                     int(IROW - 1,8)
           APOS = POSELT + int(IROW  - 1,8) * LDA8 + 
     &            int(IEND_BLOCK,8)
           CALL dgemm( 'N','N', IROW + Block - IEND_BLOCK - 1, 
     &            Block, NPIV_BLOCK,
     &                ALPHA, A( UPOS ), LDA,
     &                A( LPOS ), LDA, ONE, A( APOS ), LDA )
         ENDDO
#else
         DO IROW = IEND_BLOCK+1, LAST_VAR, BLSIZE
          Block = min( BLSIZE, LAST_VAR - IROW + 1 )
           LPOS = POSELT + int( IROW - 1,8) * LDA8 +
     &                     int(IBEG_BLOCK - 1,8)
           UPOS = POSELT + int(IBEG_BLOCK - 1,8) * LDA8 +
     &                     int( IROW - 1,8)
           APOS = POSELT + int( IROW - 1,8) * LDA8 + int( IROW - 1,8)
           CALL dgemm( 'N','N', Block, LAST_VAR - IROW + 1, NPIV_BLOCK,
     &                ALPHA, A( UPOS ), LDA,
     &                A( LPOS ), LDA, ONE, A( APOS ), LDA )
         END DO
#endif
        END IF
       LPOS = POSELT + int(LAST_VAR,8)*LDA8 + int(IBEG_BLOCK - 1,8)
       UPOS = POSELT + int(IBEG_BLOCK-1,8) * LDA8 +
     &                 int(IEND_BLOCK,8)
       APOS = POSELT + int(LAST_VAR,8)*LDA8 + int(IEND_BLOCK,8)
       IF ( POSTPONE_COL_UPDATE ) THEN
         IF (NASS.GT. LAST_VAR) THEN
           CALL dgemm('N', 'N', NEL1, NASS-LAST_VAR, NPIV_BLOCK, ALPHA, 
     &              A(UPOS), LDA, A(LPOS), LDA, ONE, 
     &              A(APOS), LDA)
         ENDIF
       ELSE
         CALL dgemm('N', 'N', NEL1, NFRONT-LAST_VAR, NPIV_BLOCK, ALPHA, 
     &              A(UPOS), LDA, A(LPOS), LDA, ONE, 
     &              A(APOS), LDA)
       END IF
      ENDIF
  500 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_SQ_LDLT
        SUBROUTINE DMUMPS_SWAP( A, LA, IW, LIW,
     &                       IOLDPS, NPIVP1, IPIV, POSELT, NASS,
     &                       LDA, NFRONT, LEVEL, K219, K50, XSIZE )
        IMPLICIT NONE
      INTEGER(8) :: POSELT, LA
      INTEGER LIW, IOLDPS, NPIVP1, IPIV
      INTEGER LDA, NFRONT, NASS, LEVEL, K219, K50, XSIZE
      DOUBLE PRECISION A( LA )
      INTEGER IW( LIW )
      INCLUDE 'mumps_headers.h'
      INTEGER ISW, ISWPS1, ISWPS2, HF
      INTEGER(8) :: IDIAG, APOS
      INTEGER(8) :: LDA8
      DOUBLE PRECISION SWOP
            LDA8 = int(LDA,8)
            APOS = POSELT + LDA8*int(IPIV-1,8) + int(NPIVP1-1,8)
            IDIAG = APOS + int(IPIV - NPIVP1,8)
            HF = 6 + IW( IOLDPS + 5 + XSIZE) + XSIZE
            ISWPS1 = IOLDPS + HF + NPIVP1 - 1
            ISWPS2 = IOLDPS + HF + IPIV - 1
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            ISW = IW(ISWPS1+NFRONT)
            IW(ISWPS1+NFRONT) = IW(ISWPS2+NFRONT)
            IW(ISWPS2+NFRONT) = ISW
            IF ( LEVEL .eq. 2 ) THEN
              CALL dswap( NPIVP1 - 1,
     &            A( POSELT + int(NPIVP1-1,8) ), LDA,
     &            A( POSELT + int(IPIV-1,8)   ), LDA )
            END IF
            CALL dswap( NPIVP1-1,
     &           A( POSELT+int(NPIVP1-1,8) * LDA8 ), 1,
     &           A( POSELT + int(IPIV-1,8) * LDA8 ), 1 )
             CALL dswap( IPIV - NPIVP1 - 1,
     &           A( POSELT+int(NPIVP1,8) * LDA8 + int(NPIVP1-1,8) ),
     &           LDA, A( APOS + 1_8 ), 1 )
            SWOP = A(IDIAG)
            A(IDIAG) = A( POSELT+int(NPIVP1-1,8)*LDA8+int(NPIVP1-1,8) )
            A( POSELT + int(NPIVP1-1,8)*LDA8 + int(NPIVP1-1,8) ) = SWOP
            CALL dswap( NASS - IPIV, A( APOS + LDA8 ), LDA,
     &                  A( IDIAG + LDA8 ), LDA )
            IF ( LEVEL .eq. 1 ) THEN
              CALL dswap( NFRONT - NASS,
     &        A( APOS  + int(NASS-IPIV+1,8) * LDA8 ), LDA,
     &        A( IDIAG + int(NASS-IPIV+1,8) * LDA8 ), LDA )
            END IF
            IF (K219.NE.0 .AND.K50.EQ.2) THEN
             IF ( LEVEL .eq. 2) THEN
              APOS                 = POSELT+LDA8*LDA8-1_8
              SWOP                 = A(APOS+int(NPIVP1,8))
              A(APOS+int(NPIVP1,8))= A(APOS+int(IPIV,8))
              A(APOS+int(IPIV,8))  = SWOP
             ENDIF
            ENDIF
        RETURN
        END SUBROUTINE DMUMPS_SWAP
      SUBROUTINE DMUMPS_FAC_T_LDLT(NFRONT,NASS,
     &    IW,LIW,A,LA,
     &    LDA,
     &    IOLDPS,POSELT,KEEP,KEEP8,
     &    POSTPONE_COL_UPDATE, ETATASS,
     &    TYPEFile, LAFAC, MonBloc, NextPiv2beWritten,
     &    LIWFAC, MYID, IFLAG
     &    )
      USE DMUMPS_OOC
      IMPLICIT NONE
      INTEGER NFRONT, NASS,LIW
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW) 
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: POSELT
      INTEGER LDA
      INTEGER IOLDPS, ETATASS
      LOGICAL POSTPONE_COL_UPDATE
      INTEGER(8) :: LAFAC
      INTEGER TYPEFile, NextPiv2beWritten
      INTEGER LIWFAC, MYID, IFLAG
      TYPE(IO_BLOCK):: MonBloc
      INTEGER IDUMMY
      LOGICAL LAST_CALL
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: UPOS, APOS, LPOS
      INTEGER(8) :: LDA8
      INTEGER BLSIZE, BLSIZE2, Block, IROW, NPIV, I, IROWEND
      INTEGER I2, I2END, Block2
      DOUBLE PRECISION  ONE, ALPHA, BETA, ZERO, VALPIV
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      LDA8 = int(LDA,8)
      IF (ETATASS.EQ.1) THEN
        BETA = ZERO
      ELSE
        BETA = ONE
      ENDIF
      IF ( NFRONT - NASS > KEEP(57) ) THEN
        BLSIZE = KEEP(58)
      ELSE
        BLSIZE = NFRONT - NASS
      END IF
      BLSIZE2 = KEEP(218)
      NPIV = IW( IOLDPS + 1 + KEEP(IXSZ))
      IF ( NFRONT - NASS .GT. 0 ) THEN
       IF ( POSTPONE_COL_UPDATE ) THEN
         CALL dtrsm( 'L', 'U', 'T', 'U',
     &               NPIV, NFRONT-NPIV, ONE,
     &               A( POSELT ), LDA,
     &               A( POSELT + LDA8 * int(NPIV,8) ), LDA )
       ENDIF
       DO IROWEND = NFRONT - NASS, 1, -BLSIZE
        Block = min( BLSIZE, IROWEND )
        IROW  = IROWEND - Block + 1
        LPOS = POSELT + int(NASS,8)*LDA8 + int(IROW-1,8) * LDA8
        APOS = POSELT + int(NASS,8)*LDA8 + int(IROW-1,8) * LDA8 +
     &                  int(NASS + IROW - 1,8)
        UPOS = POSELT + int(NASS,8)
        IF (.NOT. POSTPONE_COL_UPDATE) THEN
          UPOS = POSELT + int(NASS + IROW - 1,8)
        ENDIF
        IF (POSTPONE_COL_UPDATE) THEN
         DO I = 1, NPIV
          VALPIV = ONE/A(POSELT+(LDA8+1_8)*int(I-1,8))
          CALL dcopy( Block, A( LPOS+int(I-1,8) ), LDA,
     &                       A( UPOS+int(I-1,8)*LDA8 ), 1 )
          CALL dscal( Block, VALPIV,
     &                A( LPOS + int(I - 1,8) ), LDA )
         ENDDO
        ENDIF
        DO I2END = Block, 1, -BLSIZE2
          Block2 = min(BLSIZE2, I2END)
          I2 = I2END - Block2+1
          CALL dgemm('N', 'N', Block2, Block-I2+1, NPIV, ALPHA,
     &               A(UPOS+int(I2-1,8)), LDA,
     &               A(LPOS+int(I2-1,8)*LDA8), LDA,
     &               BETA,
     &               A(APOS + int(I2-1,8) + int(I2-1,8)*LDA8), LDA)
          IF (KEEP(201).EQ.1) THEN
            IF (NextPiv2beWritten.LE.NPIV) THEN
              LAST_CALL=.FALSE.
              CALL DMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE, TYPEFile,
     &        A(POSELT), LAFAC, MonBloc,
     &        NextPiv2beWritten, IDUMMY,
     &        IW(IOLDPS), LIWFAC, MYID,
     &        KEEP8(31),
     &        IFLAG,LAST_CALL )
              IF (IFLAG .LT. 0 ) RETURN
            ENDIF
          ENDIF
        ENDDO
        IF ( NFRONT - NASS - IROW + 1 - Block > 0 ) THEN
        CALL dgemm( 'N', 'N', Block, NFRONT-NASS-Block-IROW+1, NPIV,
     &              ALPHA,  A( UPOS ), LDA,
     &              A( LPOS + LDA8 * int(Block,8) ), LDA,
     &              BETA,
     &              A( APOS + LDA8 * int(Block,8) ), LDA )
        ENDIF
       END DO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_FAC_T_LDLT
      SUBROUTINE DMUMPS_STORE_PERMINFO( PIVRPTR, NBPANELS, PIVR, NASS,
     &                                  K, P, LastPanelonDisk,
     &                                  LastPIVRPTRIndexFilled )
      IMPLICIT NONE
      INTEGER, intent(in) :: NBPANELS, NASS, K, P
      INTEGER, intent(inout) :: PIVRPTR(NBPANELS), PIVR(NASS)
      INTEGER LastPanelonDisk, LastPIVRPTRIndexFilled
      INTEGER I
      IF ( LastPanelonDisk+1 > NBPANELS ) THEN
           WRITE(*,*) "INTERNAL ERROR IN DMUMPS_STORE_PERMINFO!"
           WRITE(*,*) "NASS=",NASS,"PIVRPTR=",PIVRPTR(1:NBPANELS)
           WRITE(*,*) "K=",K, "P=",P, "LastPanelonDisk=",LastPanelonDisk
           WRITE(*,*) "LastPIVRPTRIndexFilled=", LastPIVRPTRIndexFilled
           CALL MUMPS_ABORT()
      ENDIF
      PIVRPTR(LastPanelonDisk+1) = K + 1
      IF (LastPanelonDisk.NE.0) THEN
        PIVR(K - PIVRPTR(1) + 1) = P
        DO I = LastPIVRPTRIndexFilled + 1, LastPanelonDisk
          PIVRPTR(I)=PIVRPTR(LastPIVRPTRIndexFilled)
        ENDDO
      ENDIF
      LastPIVRPTRIndexFilled = LastPanelonDisk + 1
      RETURN
      END SUBROUTINE DMUMPS_STORE_PERMINFO
      END MODULE DMUMPS_FAC_FRONT_AUX_M
