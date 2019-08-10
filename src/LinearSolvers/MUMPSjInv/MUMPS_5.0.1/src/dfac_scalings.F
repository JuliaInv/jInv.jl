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
       SUBROUTINE DMUMPS_FAC_A(N, NZ, NSCA, 
     &      ASPK, IRN, ICN, COLSCA, ROWSCA, WK, LWK, WK_REAL,
     &      LWK_REAL, ICNTL, INFO)
       IMPLICIT NONE
      INTEGER N, NZ, NSCA
      INTEGER IRN(NZ), ICN(NZ)
      INTEGER ICNTL(40), INFO(40)
      DOUBLE PRECISION    ASPK(NZ)
      DOUBLE PRECISION COLSCA(*), ROWSCA(*)
      INTEGER LWK, LWK_REAL
      DOUBLE PRECISION    WK(LWK)
      DOUBLE PRECISION WK_REAL(LWK_REAL)
      INTEGER MPG,LP
      INTEGER IWNOR
      INTEGER I
      LOGICAL PROK
      DOUBLE PRECISION ONE
      PARAMETER( ONE = 1.0D0 )
      LP      = ICNTL(1)
      MPG     = ICNTL(2)
      MPG    = ICNTL(3)
      PROK   = ((MPG.GT.0).AND.(ICNTL(4).GE.2))
      IF (PROK) WRITE(MPG,101)
 101    FORMAT(/' ****** SCALING OF ORIGINAL MATRIX '/)
        IF (NSCA.EQ.1) THEN
         IF (PROK)
     &    WRITE (MPG,*) ' DIAGONAL SCALING '
        ELSEIF (NSCA.EQ.3) THEN
         IF (PROK)
     &   WRITE (MPG,*) ' COLUMN SCALING'
        ELSEIF (NSCA.EQ.4) THEN
         IF (PROK)
     &   WRITE (MPG,*) ' ROW AND COLUMN SCALING (1 Pass)'
        ENDIF
        DO 10 I=1,N
            COLSCA(I) = ONE
            ROWSCA(I) = ONE
 10     CONTINUE
        IF (5*N.GT.LWK_REAL) GOTO 410
        IWNOR = 1
          IF (NSCA.EQ.1) THEN
            CALL DMUMPS_FAC_V(N,NZ,ASPK,IRN,ICN,
     &        COLSCA,ROWSCA,MPG)
          ELSEIF (NSCA.EQ.3) THEN
            CALL DMUMPS_FAC_Y(N,NZ,ASPK,IRN,ICN,WK_REAL(IWNOR),
     &      COLSCA, MPG)
          ELSEIF (NSCA.EQ.4) THEN
            CALL DMUMPS_ROWCOL(N,NZ,IRN,ICN,ASPK,
     &      WK_REAL(IWNOR),WK_REAL(IWNOR+N),COLSCA,ROWSCA,MPG)
          ENDIF
      GOTO 500
      INFO(2) = NZ-LWK
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1))
     &     WRITE(LP,*) '*** ERROR: Not enough space to scale matrix'
      GOTO 500
 410  INFO(1) = -5
      INFO(2) = 5*N-LWK_REAL
      IF ((LP.GT.0).AND.(ICNTL(4).GE.1))
     & WRITE(LP,*) '*** ERROR: Not enough space to scale matrix'
      GOTO 500
 500  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_A
      SUBROUTINE DMUMPS_ROWCOL(N,NZ,IRN,ICN,VAL,
     &    RNOR,CNOR,COLSCA,ROWSCA,MPRINT)
      INTEGER N, NZ
      DOUBLE PRECISION    VAL(NZ)
      DOUBLE PRECISION    RNOR(N),CNOR(N)
      DOUBLE PRECISION    COLSCA(N),ROWSCA(N)
      DOUBLE PRECISION    CMIN,CMAX,RMIN,ARNOR,ACNOR
      INTEGER IRN(NZ), ICN(NZ)
      DOUBLE PRECISION    VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      DOUBLE PRECISION ZERO, ONE
      PARAMETER(ZERO=0.0D0, ONE=1.0D0)
      DO 50 J=1,N
       CNOR(J)   = ZERO
       RNOR(J)   = ZERO
  50  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I.LE.0).OR.(I.GT.N).OR.
     &        (J.LE.0).OR.(J.GT.N)) GOTO 100
            VDIAG = abs(VAL(K))
            IF (VDIAG.GT.CNOR(J)) THEN
              CNOR(J) =     VDIAG
            ENDIF
            IF (VDIAG.GT.RNOR(I)) THEN
              RNOR(I) =     VDIAG
            ENDIF
 100   CONTINUE
      IF (MPRINT.GT.0) THEN
       CMIN = CNOR(1)
       CMAX = CNOR(1)
       RMIN = RNOR(1)
       DO 111 I=1,N
        ARNOR = RNOR(I)
        ACNOR = CNOR(I)
        IF (ACNOR.GT.CMAX) CMAX=ACNOR
        IF (ACNOR.LT.CMIN) CMIN=ACNOR
        IF (ARNOR.LT.RMIN) RMIN=ARNOR
 111   CONTINUE
       WRITE(MPRINT,*) '**** STAT. OF MATRIX PRIOR ROW&COL SCALING'
       WRITE(MPRINT,*) ' MAXIMUM NORM-MAX OF COLUMNS:',CMAX
       WRITE(MPRINT,*) ' MINIMUM NORM-MAX OF COLUMNS:',CMIN
       WRITE(MPRINT,*) ' MINIMUM NORM-MAX OF ROWS   :',RMIN
      ENDIF
      DO 120 J=1,N
       IF (CNOR(J).LE.ZERO) THEN
         CNOR(J)   = ONE
       ELSE
         CNOR(J)   = ONE / CNOR(J)
       ENDIF
 120  CONTINUE
      DO 130 J=1,N
       IF (RNOR(J).LE.ZERO) THEN
         RNOR(J)   = ONE
       ELSE
         RNOR(J)   = ONE / RNOR(J)
       ENDIF
 130  CONTINUE
       DO 110 I=1,N
        ROWSCA(I) = ROWSCA(I) * RNOR(I)
        COLSCA(I) = COLSCA(I) * CNOR(I)
 110   CONTINUE
      IF (MPRINT.GT.0)
     &  WRITE(MPRINT,*) ' END OF SCALING BY MAX IN ROW AND COL'
      RETURN
      END SUBROUTINE DMUMPS_ROWCOL
      SUBROUTINE DMUMPS_FAC_Y(N,NZ,VAL,IRN,ICN,
     &       CNOR,COLSCA,MPRINT)
      INTEGER N,NZ
      DOUBLE PRECISION VAL(NZ)
      DOUBLE PRECISION CNOR(N)
      DOUBLE PRECISION COLSCA(N)
      INTEGER IRN(NZ), ICN(NZ)
      DOUBLE PRECISION VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      DOUBLE PRECISION ZERO, ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DO 10 J=1,N
       CNOR(J)   = ZERO
  10  CONTINUE
      DO 100 K=1,NZ
        I = IRN(K)
        J = ICN(K)
        IF ((I.LE.0).OR.(I.GT.N).OR.
     &      (J.LE.0).OR.(J.GT.N)) GOTO 100
        VDIAG = abs(VAL(K))
        IF (VDIAG.GT.CNOR(J)) THEN
           CNOR(J) =     VDIAG
        ENDIF
 100  CONTINUE
      DO 110 J=1,N
       IF (CNOR(J).LE.ZERO) THEN
         CNOR(J)   = ONE
       ELSE
         CNOR(J)   = ONE/CNOR(J)
       ENDIF
 110  CONTINUE
       DO 215 I=1,N
        COLSCA(I) = COLSCA(I) * CNOR(I)
 215   CONTINUE
      IF (MPRINT.GT.0) WRITE(MPRINT,*) ' END OF COLUMN SCALING'
      RETURN
      END SUBROUTINE DMUMPS_FAC_Y
      SUBROUTINE DMUMPS_FAC_V(N,NZ,VAL,IRN,ICN,
     &      COLSCA,ROWSCA,MPRINT)
      INTEGER   N, NZ
      DOUBLE PRECISION  VAL(NZ)
      DOUBLE PRECISION ROWSCA(N),COLSCA(N)
      INTEGER   IRN(NZ),ICN(NZ)
      DOUBLE PRECISION      VDIAG
      INTEGER   MPRINT,I,J,K
      INTRINSIC sqrt
      DOUBLE PRECISION ZERO, ONE
      PARAMETER(ZERO=0.0D0, ONE=1.0D0)
      DO 10 I=1,N
       ROWSCA(I)   = ONE
  10  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          IF ((I.GT.N).OR.(I.LE.0)) GOTO 100
          J = ICN(K)
          IF (I.EQ.J) THEN
            VDIAG = abs(VAL(K))
            IF (VDIAG.GT.ZERO) THEN
              ROWSCA(J) = ONE/(sqrt(VDIAG))
            ENDIF
          ENDIF
 100   CONTINUE
       DO 110 I=1,N
        COLSCA(I) = ROWSCA(I)
 110   CONTINUE
      IF (MPRINT.GT.0) WRITE(MPRINT,*) ' END OF DIAGONAL SCALING'
      RETURN
      END SUBROUTINE DMUMPS_FAC_V
      SUBROUTINE DMUMPS_FAC_X(NSCA,N,NZ,IRN,ICN,VAL,
     &    RNOR,ROWSCA,MPRINT)
      INTEGER N, NZ, NSCA
      INTEGER IRN(NZ), ICN(NZ)
      DOUBLE PRECISION VAL(NZ)
      DOUBLE PRECISION RNOR(N)
      DOUBLE PRECISION ROWSCA(N)
      DOUBLE PRECISION VDIAG
      INTEGER MPRINT
      INTEGER I,J,K
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
      DO 50 J=1,N
       RNOR(J)   = ZERO
  50  CONTINUE
      DO 100 K=1,NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I.LE.0).OR.(I.GT.N).OR.
     &        (J.LE.0).OR.(J.GT.N)) GOTO 100
            VDIAG = abs(VAL(K))
            IF (VDIAG.GT.RNOR(I)) THEN
              RNOR(I) =  VDIAG
            ENDIF
 100   CONTINUE
      DO 130 J=1,N
       IF (RNOR(J).LE.ZERO) THEN
         RNOR(J)   = ONE
       ELSE
         RNOR(J)   = ONE/RNOR(J)
       ENDIF
 130  CONTINUE
      DO 110 I=1,N
        ROWSCA(I) = ROWSCA(I)* RNOR(I)
 110  CONTINUE
      IF ( (NSCA.EQ.4) .OR. (NSCA.EQ.6) ) THEN
        DO 150 K=1,NZ
          I   = IRN(K)
          J   = ICN(K)
          IF (min(I,J).LT.1 .OR. I.GT.N .OR. J.GT.N) GOTO 150
          VAL(K) = VAL(K) * RNOR(I)
 150    CONTINUE
      ENDIF
      IF (MPRINT.GT.0)
     &  WRITE(MPRINT,'(A)') '  END OF ROW SCALING'
      RETURN
      END SUBROUTINE DMUMPS_FAC_X
      SUBROUTINE DMUMPS_ANORMINF( id,  ANORMINF, LSCAL )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER MASTER, IERR
      PARAMETER( MASTER = 0 )
      TYPE(DMUMPS_STRUC), TARGET :: id
      DOUBLE PRECISION, INTENT(OUT) :: ANORMINF
      LOGICAL :: LSCAL
      INTEGER, DIMENSION (:), POINTER :: KEEP,INFO
      INTEGER(8), DIMENSION (:), POINTER :: KEEP8
      LOGICAL :: I_AM_SLAVE
      DOUBLE PRECISION DUMMY(1)
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0)
      DOUBLE PRECISION, ALLOCATABLE :: SUMR(:), SUMR_LOC(:)
      INTEGER :: allocok, MTYPE, I
      INFO =>id%INFO
      KEEP =>id%KEEP
      KEEP8 =>id%KEEP8
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &             ( id%MYID .eq. MASTER .AND.
     &               KEEP(46) .eq. 1 ) )
      IF (id%MYID .EQ. MASTER) THEN
       ALLOCATE( SUMR( id%N ), stat =allocok )
       IF (allocok .GT.0 ) THEN
        id%INFO(1)=-13
        id%INFO(2)=id%N
        RETURN
       ENDIF
      ENDIF
      IF ( KEEP(54) .eq. 0 ) THEN
          IF (id%MYID .EQ. MASTER) THEN
            IF (KEEP(55).EQ.0) THEN
             IF (.NOT.LSCAL) THEN
              CALL DMUMPS_SOL_X(id%A(1),
     &          id%NZ, id%N,
     &          id%IRN(1), id%JCN(1),
     &          SUMR, KEEP(1),KEEP8(1) )
             ELSE
              CALL DMUMPS_SCAL_X(id%A(1),
     &          id%NZ, id%N,
     &          id%IRN(1), id%JCN(1), 
     &          SUMR, KEEP(1), KEEP8(1),
     &          id%COLSCA(1))
             ENDIF
            ELSE
             MTYPE = 1
             IF (.NOT.LSCAL) THEN
              CALL DMUMPS_SOL_X_ELT(MTYPE, id%N,
     &           id%NELT, id%ELTPTR(1),
     &           id%LELTVAR, id%ELTVAR(1),
     &           id%NA_ELT, id%A_ELT(1),
     &           SUMR, KEEP(1),KEEP8(1) )
             ELSE
              CALL DMUMPS_SOL_SCALX_ELT(MTYPE, id%N,
     &           id%NELT, id%ELTPTR(1),
     &           id%LELTVAR, id%ELTVAR(1),
     &           id%NA_ELT, id%A_ELT(1),
     &           SUMR, KEEP(1),KEEP8(1), id%COLSCA(1))
             ENDIF
            ENDIF
          ENDIF
      ELSE
          ALLOCATE( SUMR_LOC( id%N ), stat =allocok )
          IF (allocok .GT.0 ) THEN
             id%INFO(1)=-13
             id%INFO(2)=id%N
             RETURN
          ENDIF
          IF ( I_AM_SLAVE .and.
     &           id%NZ_loc .NE. 0 ) THEN
           IF (.NOT.LSCAL) THEN
              CALL DMUMPS_SOL_X(id%A_loc(1),
     &          id%NZ_loc, id%N,
     &          id%IRN_loc(1), id%JCN_loc(1), 
     &          SUMR_LOC, id%KEEP(1),id%KEEP8(1) )
           ELSE
              CALL DMUMPS_SCAL_X(id%A_loc(1),
     &          id%NZ_loc, id%N,
     &          id%IRN_loc(1), id%JCN_loc(1), 
     &          SUMR_LOC, id%KEEP(1),id%KEEP8(1),
     &          id%COLSCA(1))
           ENDIF
          ELSE
           SUMR_LOC = ZERO
          ENDIF
          IF ( id%MYID .eq. MASTER ) THEN
              CALL MPI_REDUCE( SUMR_LOC, SUMR,
     &        id%N, MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MASTER,id%COMM, IERR)
          ELSE
              CALL MPI_REDUCE( SUMR_LOC, DUMMY,
     &        id%N, MPI_DOUBLE_PRECISION,
     &        MPI_SUM,MASTER,id%COMM, IERR)
          END IF
        DEALLOCATE (SUMR_LOC)
      ENDIF
      IF ( id%MYID .eq. MASTER ) THEN
       ANORMINF = dble(ZERO)
        IF (LSCAL) THEN
         DO I = 1, id%N
          ANORMINF = max(abs(id%ROWSCA(I) * SUMR(I)), 
     &                  ANORMINF)
         ENDDO
        ELSE
         DO I = 1, id%N
          ANORMINF = max(abs(SUMR(I)), 
     &                  ANORMINF)
         ENDDO
        ENDIF
      ENDIF
      CALL MPI_BCAST(ANORMINF, 1,
     &              MPI_DOUBLE_PRECISION, MASTER,
     &              id%COMM, IERR )
      IF (id%MYID .eq. MASTER) DEALLOCATE (SUMR)
      RETURN
      END SUBROUTINE DMUMPS_ANORMINF
