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
      SUBROUTINE DMUMPS_MV_ELT( N, NELT, ELTPTR, ELTVAR, A_ELT,
     &                          X, Y, K50, MTYPE )
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C  To perform the matrix vector product
C      A_ELT X = Y    if MTYPE = 1
C      A_ELT^T X = Y  if MTYPE = 0
C
C  If K50 is different from 0, then the elements are
C  supposed to be in symmetric packed storage; the
C  lower part is stored by columns.
C  Otherwise, the element is square, stored by columns.
C
C  Note
C  ====
C
C  A_ELT is processed entry by entry and this code is not
C  optimized. In particular, one could gather/scatter
C  X / Y for each element to improve performance.
C
C  Arguments
C  =========
C
      INTEGER N, NELT, K50, MTYPE
      INTEGER ELTPTR( NELT + 1 ), ELTVAR( * )
      DOUBLE PRECISION A_ELT( * ), X( N ), Y( N )
C
C  Local variables
C  ===============
C
      INTEGER IEL, I , J, K, SIZEI, IELPTR
      DOUBLE PRECISION TEMP
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
C
C
C     Executable statements
C     =====================
C
      Y = ZERO
      K = 1
C     --------------------
C     Process the elements
C     --------------------
      DO IEL = 1, NELT
        SIZEI  = ELTPTR( IEL + 1 ) - ELTPTR( IEL )
        IELPTR = ELTPTR( IEL ) - 1
        IF ( K50 .eq. 0 ) THEN
C         -------------------
C         Unsymmetric element
C         stored by columns
C         -------------------
          IF ( MTYPE .eq. 1 ) THEN
C           -----------------
C           Compute A_ELT x X
C           -----------------
            DO J = 1, SIZEI
              TEMP = X( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                Y( ELTVAR( IELPTR + I ) ) =
     &          Y( ELTVAR( IELPTR + I ) ) +
     &             A_ELT( K ) * TEMP
                K = K + 1
              END DO
            END DO
          ELSE
C           -------------------
C           Compute A_ELT^T x X
C           -------------------
            DO J = 1, SIZEI
              TEMP = Y( ELTVAR( IELPTR + J ) )
              DO I = 1, SIZEI
                TEMP = TEMP + 
     &          A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
                K = K + 1
              END DO
              Y( ELTVAR( IELPTR + J ) ) = TEMP
            END DO
          END IF
        ELSE
C         -----------------
C         Symmetric element
C         L stored by cols
C         -----------------
          DO J = 1, SIZEI
C           Diagonal counted once
            Y( ELTVAR( IELPTR + J ) ) =
     &      Y( ELTVAR( IELPTR + J ) ) +
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
            K = K + 1
            DO I = J+1, SIZEI
C             Off diagonal + transpose
              Y( ELTVAR( IELPTR + I ) ) =
     &        Y( ELTVAR( IELPTR + I ) ) +
     &           A_ELT( K ) * X( ELTVAR( IELPTR + J ) )
              Y( ELTVAR( IELPTR + J ) ) =
     &        Y( ELTVAR( IELPTR + J ) ) +
     &           A_ELT( K ) * X( ELTVAR( IELPTR + I ) )
              K = K + 1
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE DMUMPS_MV_ELT
      SUBROUTINE DMUMPS_LOC_MV
     &( N, NZ_loc, IRN_loc, JCN_loc, A_loc, X, Y_loc,
     &  LDLT, MTYPE)
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Perform a distributed matrix vector product.
C        Y_loc <- A X   if MTYPE = 1
C        Y_loc <- A^T X if MTYPE = 0
C
C     Notes:
C     =====
C
C     1) assembly of all Y_loc still has to be done on exit.
C     2) X should be available on all processors.
C
C     Arguments:
C     =========
C
      INTEGER N, NZ_loc
      INTEGER IRN_loc( NZ_loc ), JCN_loc( NZ_loc )
      DOUBLE PRECISION A_loc( NZ_loc ), X( N ), Y_loc( N )
      INTEGER LDLT, MTYPE
C
C     Locals variables:
C     ================
C
      INTEGER I, J, K
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      Y_loc = ZERO
      IF ( LDLT .eq. 0 ) THEN
C       Unsymmetric
        IF ( MTYPE .eq. 1 ) THEN
C         No transpose
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &          (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(I) = Y_loc(I) + A_loc(K) * X(J)
        ENDDO
        ELSE
C         Transpose
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N)
     &        .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(J) = Y_loc(J) + A_loc(K) * X(I)
        ENDDO
        END IF
      ELSE
C       Lower (or upper) part of symmetric
C       matrix was provided (LDLT facto)
        DO K = 1, NZ_loc
          I = IRN_loc(K)
          J = JCN_loc(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &        (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(I) = Y_loc(I) + A_loc(K) * X(J)
          IF (J.NE.I) THEN
            Y_loc(J) = Y_loc(J) + A_loc(K) * X(I)
          ENDIF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_LOC_MV
      SUBROUTINE DMUMPS_MV( N, NZ, IRN, ICN, ASPK, X, Y,
     &                      LDLT, MTYPE, MAXTRANS, PERM,
     &                      IFLAG, IERROR )
C
C     Purpose:
C     =======
C
C     Perform matrix-vector product
C        Y <- A X if MTYPE = 1
C        Y <- A^T X if MTYPE = 0
C
C
C     Note:
C     ====
C
C     MAXTRANS should be set to 1 if a column permutation
C     was applied on A and we still want the matrix vector
C     product wrt the original matrix.
C
C     Arguments:
C     =========
C
      INTEGER N, NZ, LDLT, MTYPE, MAXTRANS
      INTEGER IRN( NZ ), ICN( NZ ) 
      INTEGER PERM( N )
      DOUBLE PRECISION ASPK( NZ ), X( N ), Y( N )
      INTEGER, intent(out) :: IFLAG, IERROR
C
C     Local variables
C     ===============
C
      INTEGER K, I, J
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PX
      DOUBLE PRECISION ZERO
      INTEGER :: allocok
      PARAMETER( ZERO = 0.0D0 )
      Y = ZERO
      ALLOCATE(PX(N), stat=allocok)
      IF (allocok < 0) THEN
        IFLAG  = -13
        IERROR = N
        RETURN
      ENDIF
C
C     --------------------------------------
C     Permute X if A has been permuted
C     with some max-trans column permutation
C     --------------------------------------
      IF ( MAXTRANS .eq. 1 .and. MTYPE .eq. 1) THEN
        DO I = 1, N
          PX(I) = X( PERM( I ) )
        END DO
      ELSE
        PX = X
      END IF
      IF ( LDLT .eq. 0 ) THEN
C
C     Complete unsymmetric matrix was provided (LU facto)
       IF (MTYPE .EQ. 1) THEN
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y(I) = Y(I) + ASPK(K) * PX(J)
        ENDDO
       ELSE
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y(J) = Y(J) + ASPK(K) * PX(I)
        ENDDO
       ENDIF
C
      ELSE
C
C       Lower (or upper) part of symmetric
C       matrix was provided (LDLT facto)
        DO K = 1, NZ
          I = IRN(K)
          J = ICN(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y(I) = Y(I) + ASPK(K) * PX(J)
          IF (J.NE.I) THEN
            Y(J) = Y(J) + ASPK(K) * PX(I)
          ENDIF
        ENDDO
      END IF
      IF ( MAXTRANS .EQ. 1 .AND. MTYPE .eq. 0 ) THEN
      PX = Y
      DO I = 1, N
        Y( PERM( I ) ) = PX( I )
      END DO
      END IF
      DEALLOCATE(PX)
      RETURN
      END SUBROUTINE DMUMPS_MV
C
C
      SUBROUTINE DMUMPS_LOC_OMEGA1
     &( N, NZ_loc, IRN_loc, JCN_loc, A_loc, X, Y_loc,
     &  LDLT, MTYPE)
      IMPLICIT NONE
C
C     Purpose:
C     =======
C     Compute
C        * If MTYPE = 1
C            Y_loc(i) = Sum | Aij | | Xj |
C                        j
C        * If MTYPE = 0
C            Y_loc(j) = Sum | Aij | | Xi |
C
C
C     Notes:
C     =====
C
C     1) assembly of all Y_loc still has to be done.
C     2) X should be available on all processors.
C
C     Arguments:
C     =========
C
      INTEGER N, NZ_loc
      INTEGER IRN_loc( NZ_loc ), JCN_loc( NZ_loc )
      DOUBLE PRECISION A_loc( NZ_loc ), X( N )
      DOUBLE PRECISION Y_loc( N )
      INTEGER LDLT, MTYPE
C
C     Local variables:
C     ===============
C
      INTEGER I, J, K
      DOUBLE PRECISION RZERO
      PARAMETER( RZERO = 0.0D0 )
C
      Y_loc = RZERO
      IF ( LDLT .eq. 0 ) THEN
C       Unsymmetric
        IF ( MTYPE .eq. 1 ) THEN
C         No transpose
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &          (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
            Y_loc(I) = Y_loc(I) + abs( A_loc(K) * X(J) )
          ENDDO
        ELSE
C         Transpose
          DO K = 1, NZ_loc
            I = IRN_loc(K)
            J = JCN_loc(K)
            IF ((I .LE. 0) .OR. (I .GT. N)
     &        .OR. (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(J) = Y_loc(J) + abs( A_loc(K) * X(I) )
          ENDDO
        END IF
      ELSE
C       Lower (or upper) part of symmetric
C       matrix was provided (LDLT facto)
        DO K = 1, NZ_loc
          I = IRN_loc(K)
          J = JCN_loc(K)
          IF ((I .LE. 0) .OR. (I .GT. N) .OR.
     &        (J .LE. 0) .OR. (J .GT. N)
     &        ) CYCLE
          Y_loc(I) = Y_loc(I) + abs( A_loc(K) * X(J) )
          IF (J.NE.I) THEN
            Y_loc(J) = Y_loc(J) + abs( A_loc(K) * X(I) )
          ENDIF
        ENDDO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_LOC_OMEGA1
