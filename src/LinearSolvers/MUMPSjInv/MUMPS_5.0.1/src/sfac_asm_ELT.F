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
      SUBROUTINE SMUMPS_ELT_ASM_S_2_S_INIT(
     &    NELT, FRT_PTR, FRT_ELT,
     &    N, INODE, IW, LIW, A, LA, 
     &    NBROWS, NBCOLS,
     &    OPASSW, OPELIW, STEP, PTRIST, PTRAST, ITLOC,
     &    RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW, INTARR, DBLARR, 
     &    ICNTL, KEEP, KEEP8, MYID)
      IMPLICIT NONE
      INTEGER NELT, N,LIW
      INTEGER(8) :: LA
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      INTEGER INODE, MYID
      INTEGER NBROWS, NBCOLS 
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER IW(LIW), ITLOC(N + KEEP(253)), STEP(N),
     &        PTRIST(KEEP(28)),
     &        FILS(N), PTRARW(NELT+1), 
     &        PTRAIW(NELT+1)
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(max(1,KEEP(14)))
      INTEGER FRT_PTR(N+1), FRT_ELT(NELT)
      REAL A(LA),
     &        DBLARR(max(1,KEEP(13)))
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER(8) :: POSELT
      INTEGER IOLDPS, NBCOLF, NBROWF, NSLAVES, HF,
     &        K1,K2,K,J,JPOS,NASS
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INCLUDE 'mumps_headers.h'
      IOLDPS  = PTRIST(STEP(INODE))
      POSELT  = PTRAST(STEP(INODE))
      NBCOLF  = IW(IOLDPS+KEEP(IXSZ))
      NBROWF  = IW(IOLDPS+2+KEEP(IXSZ))
      NASS    = IW(IOLDPS+1+KEEP(IXSZ))
      NSLAVES = IW(IOLDPS+5+KEEP(IXSZ))
      HF      = 6 + NSLAVES+KEEP(IXSZ)
      IF (NASS.LT.0) THEN
          NASS         = -NASS
          IW(IOLDPS+1+KEEP(IXSZ)) = NASS
          CALL SMUMPS_ASM_SLAVE_ELEMENTS( INODE, N, NELT, IW, LIW,
     &    IOLDPS, A, LA, POSELT, KEEP, ITLOC, FILS, PTRAIW, PTRARW,
     &    INTARR, DBLARR, FRT_PTR, FRT_ELT, RHS_MUMPS)
      END IF
      IF (NBROWS.GT.0) THEN
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = JPOS
           JPOS     = JPOS + 1
          END DO
      END IF
      RETURN
      END SUBROUTINE SMUMPS_ELT_ASM_S_2_S_INIT
      SUBROUTINE SMUMPS_ASM_SLAVE_ELEMENTS( INODE, N, NELT, IW, LIW,
     &IOLDPS, A, LA, POSELT, KEEP, ITLOC, FILS, PTRAIW, PTRARW,
     &INTARR, DBLARR, FRT_PTR, FRT_ELT, RHS_MUMPS)
      IMPLICIT NONE
      INTEGER, intent(in)    :: N, NELT, LIW, IOLDPS, INODE
      INTEGER(8), intent(in) :: LA, POSELT
      INTEGER, intent(in)    :: IW(LIW)
      INTEGER, intent(in)    :: KEEP(500)
      INTEGER, intent(inout) :: ITLOC(N+KEEP(253))
      REAL, intent(inout) :: A(LA)
      REAL, intent(in)    :: RHS_MUMPS(KEEP(255))
      INTEGER, intent(in)    :: INTARR(max(1,KEEP(14)))
      REAL, intent(in)    :: DBLARR(max(1,KEEP(13)))
      INTEGER, intent(in)    :: FRT_PTR(N+1), FRT_ELT(NELT)
      INTEGER, intent(in)    :: FILS(N), PTRAIW(NELT+1), PTRARW(NELT+1)
      INCLUDE 'mumps_headers.h'
      INTEGER    :: HF, NBROWF, NBCOLF, NASS, NSLAVES
      INTEGER    :: ILOC, IELL, ELTI, ELBEG, NUMELT, SIZE_ELTI
      INTEGER    :: I, II, J, JJ, J1, J2, K, K1, K2
      INTEGER    :: IPOS, IPOS1, IPOS2, JPOS, IJROW
      INTEGER    :: IN
      INTEGER    :: AINPUT, AII
      INTEGER(8) :: APOS, APOS2, ICT12
      INTEGER    :: K1RHS, K2RHS, JFirstRHS
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      NBCOLF  = IW(IOLDPS+KEEP(IXSZ))
      NBROWF  = IW(IOLDPS+2+KEEP(IXSZ))
      NASS    = IW(IOLDPS+1+KEEP(IXSZ))
      A(POSELT:POSELT+int(NBROWF,8)*int(NBCOLF,8)-1_8) = ZERO
      NSLAVES= IW(IOLDPS+5 + KEEP(IXSZ))
      HF      = 6 + NSLAVES + KEEP(IXSZ)
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          JPOS = 1
          DO K = K1, K2
           J        = IW(K)
           ITLOC(J) = -JPOS
           JPOS     = JPOS + 1
          END DO
          K1 = IOLDPS + HF 
          K2 = K1 + NBROWF - 1
          JPOS = 1
          IF ((KEEP(253).GT.0).AND.(KEEP(50).NE.0)) THEN
           K1RHS = 0
           K2RHS = -1
           DO K = K1, K2
            J        = IW(K)
            ITLOC(J) = -ITLOC(J)*NBCOLF + JPOS
            IF ((K1RHS.EQ.0).AND.(J.GT.N)) THEN
             K1RHS = K
             JFirstRHS=J-N 
            ENDIF
            JPOS     = JPOS + 1
           ENDDO
           IF (K1RHS.GT.0) K2RHS=K2
           IF ( K2RHS.GE.K1RHS ) THEN
             IN = INODE
             DO WHILE (IN.GT.0) 
               IJROW = -ITLOC(IN)  
               DO K = K1RHS, K2RHS
                J    = IW(K)       
                I    = ITLOC(J)    
                ILOC = mod(I,NBCOLF) 
              APOS = POSELT+int(ILOC-1,8)*int(NBCOLF,8) + 
     &               int(IJROW-1,8) 
              A(APOS) = A(APOS) + RHS_MUMPS(
     &                 (JFirstRHS+(K-K1RHS)-1)*KEEP(254)+ IN)
             ENDDO
             IN = FILS(IN)
            ENDDO
           ENDIF
          ELSE  
           DO K = K1, K2
            J        = IW(K)
            ITLOC(J) = -ITLOC(J)*NBCOLF + JPOS
            JPOS     = JPOS + 1
           END DO
          ENDIF
          ELBEG  = FRT_PTR(INODE)
          NUMELT = FRT_PTR(INODE+1) - ELBEG
          DO IELL=ELBEG,ELBEG+NUMELT-1
           ELTI = FRT_ELT(IELL)
           J1= PTRAIW(ELTI)
           J2= PTRAIW(ELTI+1)-1
           AII = PTRARW(ELTI)
           SIZE_ELTI = J2 - J1 + 1
           DO II=J1,J2
            I = ITLOC(INTARR(II))
            IF (KEEP(50).EQ.0) THEN
             IF (I.LE.0) CYCLE
             AINPUT    = AII + II - J1
             IPOS = mod(I,NBCOLF)
             ICT12 = POSELT + int(IPOS-1,8) * int(NBCOLF,8)
             DO JJ = J1, J2
              JPOS = ITLOC(INTARR(JJ))
              IF (JPOS.LE.0) THEN 
                   JPOS = -JPOS
              ELSE
                   JPOS = JPOS/NBCOLF
              END IF
              APOS2    = ICT12 + int(JPOS - 1,8)
              A(APOS2) = A(APOS2) +  DBLARR(AINPUT)
              AINPUT   = AINPUT + SIZE_ELTI
             END DO
            ELSE
              IF ( I .EQ. 0 ) THEN 
               AII = AII + J2 - II + 1
               CYCLE
              ENDIF
              IF ( I .LE. 0 ) THEN 
               IPOS1 = -I
               IPOS2 = 0
              ELSE 
               IPOS1 = I/NBCOLF
               IPOS2 = mod(I,NBCOLF)
              END IF
              ICT12 =  POSELT + int(IPOS2-1,8)*int(NBCOLF,8)
              DO JJ=II,J2
               AII = AII + 1
               J = ITLOC(INTARR(JJ))
               IF ( J .EQ. 0 ) CYCLE
               IF ( IPOS2.EQ.0 .AND. J.LE.0) CYCLE
               IF ( J .LE. 0 ) THEN
                JPOS = -J
               ELSE
                JPOS = J/NBCOLF
               END IF
               IF ( (IPOS1.GE.JPOS) .AND. (IPOS2.GT.0) ) THEN
                 APOS2 = ICT12  + int(JPOS - 1,8)
                 A(APOS2) = A(APOS2) +  DBLARR(AII-1)
               END IF
               IF ( (IPOS1.LT.JPOS) .AND. (J.GT.0) ) THEN
                 IPOS = mod(J,NBCOLF)
                 JPOS = IPOS1
                 APOS2 = POSELT + int(IPOS-1,8)*int(NBCOLF,8)
     &                          + int(JPOS - 1,8)
                 A(APOS2) = A(APOS2) +  DBLARR(AII-1)
               END IF
              END DO
            END IF
           END DO
          END DO
          K1 = IOLDPS + HF + NBROWF
          K2 = K1 + NBCOLF - 1
          DO K = K1, K2
           J = IW(K)
           ITLOC(J) = 0
          END DO
      END SUBROUTINE SMUMPS_ASM_SLAVE_ELEMENTS
