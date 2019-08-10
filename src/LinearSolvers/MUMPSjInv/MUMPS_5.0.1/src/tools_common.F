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
      SUBROUTINE MUMPS_MAKE1ROOT( N, FRERE, FILS, NFSIZ, THEROOT )
      IMPLICIT NONE
      INTEGER, intent( in    )  :: N
      INTEGER, intent( in    )  :: NFSIZ( N )
      INTEGER, intent( inout )  :: FRERE( N ), FILS( N )
      INTEGER, intent( out   )  :: THEROOT
      INTEGER INODE, IROOT, IFILS, IN, IROOTLAST, SIZE
      IROOT = -9999
      SIZE  = 0
      DO INODE = 1, N
        IF ( FRERE( INODE ) .EQ. 0 )  THEN
          IF ( NFSIZ( INODE ) .GT. SIZE ) THEN
            SIZE  = NFSIZ( INODE )
            IROOT = INODE
          END IF
        ENDIF
      END DO
      IN = IROOT
      DO WHILE ( FILS( IN ) .GT. 0 )
        IN = FILS( IN )
      END DO
      IROOTLAST = IN
      IFILS     = - FILS ( IN )
      DO INODE = 1, N
        IF ( FRERE( INODE ) .eq. 0 .and. INODE .ne. IROOT ) THEN
          IF ( IFILS .eq. 0 ) THEN
            FILS( IROOTLAST ) = - INODE
            FRERE( INODE )    = -IROOT
            IFILS             = INODE
          ELSE
            FRERE( INODE ) = -FILS( IROOTLAST )
            FILS( IROOTLAST ) = - INODE
          END IF
        END IF
      END DO
      THEROOT = IROOT
      RETURN
      END SUBROUTINE MUMPS_MAKE1ROOT
      INTEGER FUNCTION MUMPS_TYPENODE(PROCINFO_INODE, SLAVEF)
      IMPLICIT NONE
      INTEGER SLAVEF 
      INTEGER PROCINFO_INODE, TPN
      IF (PROCINFO_INODE <= SLAVEF ) THEN
        MUMPS_TYPENODE = 1
      ELSE
        TPN = (PROCINFO_INODE-1+2*SLAVEF)/SLAVEF - 1
        IF ( TPN .LT. 1 ) TPN = 1
        IF (TPN.EQ.4.OR.TPN.EQ.5.OR.TPN.EQ.6) TPN = 2
        MUMPS_TYPENODE = TPN
      END IF
      RETURN 
      END FUNCTION MUMPS_TYPENODE
      INTEGER FUNCTION MUMPS_PROCNODE(PROCINFO_INODE, SLAVEF)
      IMPLICIT NONE
      INTEGER SLAVEF 
      INTEGER PROCINFO_INODE
      IF (SLAVEF == 1) THEN
        MUMPS_PROCNODE = 0
      ELSE
        MUMPS_PROCNODE=mod(2*SLAVEF+PROCINFO_INODE-1,SLAVEF)
      END IF
      RETURN
      END FUNCTION MUMPS_PROCNODE
      INTEGER FUNCTION MUMPS_TYPESPLIT (PROCINFO_INODE, SLAVEF)
      IMPLICIT NONE
      INTEGER, intent(in) ::  SLAVEF 
      INTEGER PROCINFO_INODE, TPN
      IF (PROCINFO_INODE <= SLAVEF ) THEN
         MUMPS_TYPESPLIT = 1
      ELSE
        TPN = (PROCINFO_INODE-1+2*SLAVEF)/SLAVEF - 1
        IF ( TPN .LT. 1 ) TPN = 1
         MUMPS_TYPESPLIT = TPN
      ENDIF
      RETURN
      END FUNCTION MUMPS_TYPESPLIT
      LOGICAL FUNCTION MUMPS_ROOTSSARBR( PROCINFO_INODE, SLAVEF )
      IMPLICIT NONE
      INTEGER SLAVEF
      INTEGER TPN, PROCINFO_INODE
      TPN = (PROCINFO_INODE-1+2*SLAVEF)/SLAVEF - 1
      MUMPS_ROOTSSARBR = ( TPN .eq. 0 )
      RETURN
      END FUNCTION MUMPS_ROOTSSARBR
      LOGICAL FUNCTION MUMPS_INSSARBR( PROCINFO_INODE, SLAVEF )
      IMPLICIT NONE
      INTEGER SLAVEF
      INTEGER TPN, PROCINFO_INODE
      TPN = (PROCINFO_INODE-1+SLAVEF+SLAVEF)/SLAVEF - 1
      MUMPS_INSSARBR = ( TPN .eq. -1 )
      RETURN 
      END FUNCTION MUMPS_INSSARBR
      LOGICAL FUNCTION MUMPS_IN_OR_ROOT_SSARBR
     &        ( PROCINFO_INODE, SLAVEF )
      IMPLICIT NONE
      INTEGER SLAVEF
      INTEGER TPN, PROCINFO_INODE
      TPN = (PROCINFO_INODE-1+SLAVEF+SLAVEF)/SLAVEF - 1
      MUMPS_IN_OR_ROOT_SSARBR =
     &           ( TPN .eq. -1 .OR. TPN .eq. 0 )
      RETURN
      END FUNCTION MUMPS_IN_OR_ROOT_SSARBR
      LOGICAL FUNCTION MUMPS_I_AM_CANDIDATE( MYID, SLAVEF, INODE,
     &                 NMB_PAR2, ISTEP_TO_INIV2 , K71, STEP, N, 
     &                 CANDIDATES, KEEP24 )
      IMPLICIT NONE
      INTEGER MYID, SLAVEF, INODE, NMB_PAR2, KEEP24, I
      INTEGER K71, N
      INTEGER ISTEP_TO_INIV2 ( K71 ), STEP ( N )
      INTEGER CANDIDATES(SLAVEF+1, max(NMB_PAR2,1))
      INTEGER NCAND, POSINODE
      MUMPS_I_AM_CANDIDATE = .FALSE.
      IF (KEEP24 .eq. 0) RETURN
      POSINODE = ISTEP_TO_INIV2 ( STEP (INODE) )
      NCAND = CANDIDATES( SLAVEF+1, POSINODE )
      DO I = 1, NCAND
        IF (MYID .EQ. CANDIDATES( I, POSINODE ))
     &     MUMPS_I_AM_CANDIDATE = .TRUE.
      END DO
      RETURN
      END FUNCTION MUMPS_I_AM_CANDIDATE
      SUBROUTINE MUMPS_SECDEB(T)
      DOUBLE PRECISION T
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
      T=MPI_WTIME()
      RETURN
      END SUBROUTINE MUMPS_SECDEB
      SUBROUTINE MUMPS_SECFIN(T)
      DOUBLE PRECISION T
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
      T=MPI_WTIME()-T
      RETURN
      END SUBROUTINE MUMPS_SECFIN
      SUBROUTINE MUMPS_SORT_DOUBLES( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      DOUBLE PRECISION VAL( N )
      INTEGER I, ISWAP
      DOUBLE PRECISION SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
          IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
            DONE = .FALSE.
            ISWAP = ID( I )
            ID ( I ) = ID ( I + 1 )
            ID ( I + 1 ) = ISWAP
            SWAP = VAL( I )
            VAL( I ) = VAL( I + 1 )
            VAL( I + 1 ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_SORT_DOUBLES
      SUBROUTINE MUMPS_SORT_DOUBLES_DEC( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      DOUBLE PRECISION VAL( N )
      INTEGER I, ISWAP
      DOUBLE PRECISION SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
          IF ( VAL( I ) .LT. VAL( I + 1 ) ) THEN
            DONE = .FALSE.
            ISWAP = ID( I )
            ID ( I ) = ID ( I + 1 )
            ID ( I + 1 ) = ISWAP
            SWAP = VAL( I )
            VAL( I ) = VAL( I + 1 )
            VAL( I + 1 ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_SORT_DOUBLES_DEC
#if defined (PESSL)
      SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT,
     &                     LLD, INFO )
      INTEGER            ICSRC, ICTXT, INFO, IRSRC, LLD, M, MB, N, NB
      INTEGER            DESC( * )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     &                   LLD_, MB_, M_, NB_, N_, RSRC_
# if defined(DESC8)
      PARAMETER          ( DLEN_ = 8, DTYPE_ = 1,
     &                     CTXT_ = 7, M_ = 1, N_ = 2, MB_ = 3, NB_ = 4,
     &                     RSRC_ = 5, CSRC_ = 6, LLD_ = 8 )
# else
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     &                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     &                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
# endif
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
      EXTERNAL           blacs_gridinfo, PXERBLA
      INTEGER            NUMROC
      EXTERNAL           NUMROC
      INTRINSIC          max, min
      CALL blacs_gridinfo( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( MB.LT.1 ) THEN
         INFO = -4
      ELSE IF( NB.LT.1 ) THEN
         INFO = -5
      ELSE IF( IRSRC.LT.0 .OR. IRSRC.GE.NPROW ) THEN
         INFO = -6
      ELSE IF( ICSRC.LT.0 .OR. ICSRC.GE.NPCOL ) THEN
         INFO = -7
      ELSE IF( NPROW.EQ.-1 ) THEN
         INFO = -8
      ELSE IF( LLD.LT.max( 1, numroc( M, MB, MYROW, IRSRC,
     &                                NPROW ) ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 )
     &   CALL PXERBLA( ICTXT, 'DESCINIT', -INFO )
# ifndef DESC8
      DESC( DTYPE_ ) = BLOCK_CYCLIC_2D
# endif
      DESC( M_ )  = max( 0, M )
      DESC( N_ )  = max( 0, N )
      DESC( MB_ ) = max( 1, MB )
      DESC( NB_ ) = max( 1, NB )
      DESC( RSRC_ ) = max( 0, min( IRSRC, NPROW-1 ) )
      DESC( CSRC_ ) = max( 0, min( ICSRC, NPCOL-1 ) )
      DESC( CTXT_ ) = ICTXT
      DESC( LLD_ )  = max( LLD, max( 1, numroc( DESC( M_ ), DESC( MB_ ),
     &                              MYROW, DESC( RSRC_ ), NPROW ) ) )
      RETURN
      END SUBROUTINE DESCINIT
      SUBROUTINE PXERBLA( ICTXT, SRNAME, INFO )
      INTEGER            ICTXT, INFO
      CHARACTER*(*)      SRNAME
      INTEGER            MYCOL, MYROW, NPCOL, NPROW
      EXTERNAL           blacs_gridinfo
      CALL blacs_gridinfo( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      WRITE( *, FMT = 9999 ) MYROW, MYCOL, SRNAME, INFO
 9999 FORMAT( '{', I5, ',', I5, '}:  On entry to ', A,
     &        ' parameter number', I4, ' had an illegal value' )
      END SUBROUTINE PXERBLA
#endif
      SUBROUTINE MUMPS_MEM_CENTRALIZE(MYID, COMM, INFO, INFOG, IRANK)
      IMPLICIT NONE
      INTEGER MYID, COMM, IRANK, INFO, INFOG(2)
      INCLUDE 'mpif.h'
      INTEGER IERR_MPI, MASTER
      INTEGER TEMP1(2), TEMP2(2)
      PARAMETER( MASTER = 0 )
      CALL MPI_REDUCE( INFO, INFOG(1), 1, MPI_INTEGER,
     &                 MPI_MAX, MASTER, COMM, IERR_MPI )
      CALL MPI_REDUCE( INFO, INFOG(2), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, COMM, IERR_MPI )
      TEMP1(1) = INFO
      TEMP1(2) = MYID
      CALL MPI_REDUCE( TEMP1, TEMP2, 1, MPI_2INTEGER,
     &                 MPI_MAXLOC, MASTER, COMM, IERR_MPI )
      IF ( MYID.eq. MASTER ) THEN
        IF ( INFOG(1) .ne. TEMP2(1) ) THEN
          write(*,*) 'Error in MUMPS_MEM_CENTRALIZE'
          CALL MUMPS_ABORT()
        END IF
        IRANK    = TEMP2(2)
      ELSE
        IRANK    = -1
      END IF
      RETURN
      END SUBROUTINE MUMPS_MEM_CENTRALIZE
      INTEGER FUNCTION MUMPS_GET_POOL_LENGTH
     &        (MAX_ACTIVE_NODES,KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER MAX_ACTIVE_NODES
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      MUMPS_GET_POOL_LENGTH = MAX_ACTIVE_NODES + 1 + 3
      RETURN
      END FUNCTION MUMPS_GET_POOL_LENGTH
      SUBROUTINE MUMPS_INIT_POOL_DIST(N, LEAF,
     &           MYID_NODES,
     &           SLAVEF, NA, LNA, KEEP,KEEP8, STEP,
     &           PROCNODE_STEPS, IPOOL, LPOOL)
      IMPLICIT NONE
      INTEGER N, LEAF, MYID_NODES,
     &        SLAVEF, LPOOL, LNA
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28)), NA(LNA),
     &        IPOOL(LPOOL)
      INTEGER NBLEAF, INODE, I
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      NBLEAF = NA(1)
      LEAF = 1
      DO I = 1, NBLEAF
        INODE = NA(I+2)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),SLAVEF)
     &   .EQ.MYID_NODES) THEN
           IPOOL(LEAF) = INODE
           LEAF        = LEAF + 1
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_INIT_POOL_DIST
      SUBROUTINE MUMPS_INIT_NROOT_DIST(N, NBROOT,
     &           NROOT_LOC, MYID_NODES,
     &           SLAVEF, NA, LNA, KEEP, STEP,
     &           PROCNODE_STEPS)
      INTEGER, INTENT( OUT ) :: NROOT_LOC 
      INTEGER, INTENT( OUT ) :: NBROOT 
      INTEGER, INTENT( IN ) :: KEEP( 500 )
      INTEGER, INTENT( IN ) :: SLAVEF
      INTEGER, INTENT( IN ) :: N
      INTEGER, INTENT( IN ) :: STEP(N)
      INTEGER, INTENT( IN ) :: LNA
      INTEGER, INTENT( IN ) :: NA(LNA)
      INTEGER, INTENT( IN ) :: PROCNODE_STEPS(KEEP(28))
      INTEGER, INTENT( IN ) :: MYID_NODES
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INTEGER :: INODE, I, NBLEAF
      NBLEAF = NA(1)
      NBROOT = NA(2)
      NROOT_LOC = 0
      DO I = 1, NBROOT
        INODE = NA(I+2+NBLEAF)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &    SLAVEF).EQ.MYID_NODES) THEN
            NROOT_LOC = NROOT_LOC + 1
        END IF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_INIT_NROOT_DIST
      LOGICAL FUNCTION MUMPS_COMPARE_TAB(TAB1,TAB2,LEN1,LEN2)
      IMPLICIT NONE
      INTEGER LEN1 , LEN2 ,I
      INTEGER TAB1(LEN1)
      INTEGER TAB2(LEN2)
      MUMPS_COMPARE_TAB=.FALSE.
      IF(LEN1 .NE. LEN2) THEN
         RETURN
      ENDIF
      DO I=1 , LEN1
         IF(TAB1(I) .NE. TAB2(I)) THEN
            RETURN
         ENDIF
      ENDDO
      MUMPS_COMPARE_TAB=.TRUE.
      RETURN
      END FUNCTION MUMPS_COMPARE_TAB
      SUBROUTINE MUMPS_SORT_INT( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      INTEGER VAL( N )
      INTEGER I, ISWAP
      INTEGER SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
           IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
              DONE = .FALSE.
              ISWAP = ID( I )
              ID ( I ) = ID ( I + 1 )
              ID ( I + 1 ) = ISWAP
              SWAP = VAL( I )
              VAL( I ) = VAL( I + 1 )
              VAL( I + 1 ) = SWAP
           END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_SORT_INT
      SUBROUTINE MUMPS_SORT_INT_DEC( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      INTEGER VAL( N )
      INTEGER I, ISWAP
      INTEGER SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
           IF ( VAL( I ) .LT. VAL( I + 1 ) ) THEN
              DONE = .FALSE.
              ISWAP = ID( I )
              ID ( I ) = ID ( I + 1 )
              ID ( I + 1 ) = ISWAP
              SWAP = VAL( I )
              VAL( I ) = VAL( I + 1 )
              VAL( I + 1 ) = SWAP
           END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_SORT_INT_DEC
      SUBROUTINE MUMPS_ABORT()
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR, IERRCODE
      IERRCODE = -99
      CALL MPI_ABORT(MPI_COMM_WORLD, IERRCODE, IERR)
      RETURN
      END SUBROUTINE MUMPS_ABORT
      SUBROUTINE MUMPS_GET_PERLU(KEEP12,ICNTL14,
     &     KEEP50,KEEP54,ICNTL6,ICNTL8)
      IMPLICIT NONE
      INTEGER, intent(out)::KEEP12
      INTEGER, intent(in)::ICNTL14,KEEP50,KEEP54,ICNTL6,ICNTL8
      KEEP12 = ICNTL14 
      IF(ICNTL6.EQ.0 .AND. ICNTL8.EQ.0) RETURN
      IF ( (KEEP54.NE.0).AND. (KEEP50.NE.1)
     &     .AND. (KEEP12 .GT. 0) ) KEEP12= KEEP12+5
      RETURN
      END SUBROUTINE MUMPS_GET_PERLU
      SUBROUTINE MUMPS_BCAST_I8( I8_VALUE, ROOT, MYID, COMM, IERR)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ROOT, MYID, COMM, IERR
      INTEGER(8) :: I8_VALUE
      DOUBLE PRECISION :: DBLE_VALUE
      IF (MYID .EQ. ROOT) THEN
        DBLE_VALUE = dble(I8_VALUE)
      ENDIF
      CALL MPI_BCAST( DBLE_VALUE, 1, MPI_DOUBLE_PRECISION,
     &                ROOT,  COMM, IERR )
      I8_VALUE = int( DBLE_VALUE,8)
      RETURN
      END SUBROUTINE MUMPS_BCAST_I8
      SUBROUTINE MUMPS_REDUCEI8( IN, OUT, MPI_OP, ROOT, COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ROOT, COMM, MPI_OP
      INTEGER(8) IN, OUT
      INTEGER IERR
      DOUBLE PRECISION DIN, DOUT
      DIN =dble(IN)
      DOUT=0.0D0
      CALL MPI_REDUCE(DIN, DOUT, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_OP, ROOT, COMM, IERR)
      OUT=int(DOUT,kind=8)
      RETURN
      END SUBROUTINE MUMPS_REDUCEI8
      SUBROUTINE MUMPS_ALLREDUCEI8( IN, OUT, MPI_OP, COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, MPI_OP
      INTEGER(8) IN, OUT
      INTEGER IERR
      DOUBLE PRECISION DIN, DOUT
      DIN =dble(IN)
      DOUT=0.0D0
      CALL MPI_ALLREDUCE(DIN, DOUT, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_OP, COMM, IERR)
      OUT=int(DOUT,kind=8)
      RETURN
      END SUBROUTINE MUMPS_ALLREDUCEI8
      SUBROUTINE MUMPS_IREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      INTEGER, POINTER             :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE, MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      INTEGER, POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+MINSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_IREALLOC
      SUBROUTINE MUMPS_SREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      REAL(kind(1.E0)), POINTER    :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE, MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      REAL(kind(1.E0)), POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+MINSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_SREALLOC
      SUBROUTINE MUMPS_DREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      REAL(kind(1.D0)), POINTER    :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE, MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      REAL(kind(1.D0)), POINTER    :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+MINSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_DREALLOC
      SUBROUTINE MUMPS_CREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      COMPLEX(kind((1.E0,1.E0))), POINTER             :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE, MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      COMPLEX(kind((1.E0,1.E0))), POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+MINSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_CREALLOC
      SUBROUTINE MUMPS_ZREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      COMPLEX(kind((1.D0,1.D0))), POINTER             :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE, MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      COMPLEX(kind((1.D0,1.D0))), POINTER             :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      CHARACTER(len=60)            :: FMTA, FMTD
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+MINSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-size(ARRAY)
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+MINSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_ZREALLOC
      SUBROUTINE MUMPS_SETI8TOI4(I8, I4)
      IMPLICIT NONE
      INTEGER   , INTENT(OUT) :: I4
      INTEGER(8), INTENT(IN)  :: I8
      IF ( I8 .GT. int(huge(I4),8) ) THEN
        I4 = -int(I8/1000000_8,kind(I4))
      ELSE
        I4 = int(I8,kind(I4))
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_SETI8TOI4
      SUBROUTINE MUMPS_ABORT_ON_OVERFLOW(I8, STRING)
      IMPLICIT NONE
      INTEGER(8), INTENT(IN) :: I8
      CHARACTER(*), INTENT(IN) :: STRING
      INTEGER I4
      IF ( I8 .GT. int(huge(I4),8)) THEN
        WRITE(*,*) STRING
        CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_ABORT_ON_OVERFLOW
      SUBROUTINE MUMPS_SET_IERROR( SIZE8, IERROR  )
      INTEGER(8), INTENT(IN) :: SIZE8
      INTEGER, INTENT(OUT) :: IERROR
      CALL MUMPS_SETI8TOI4(SIZE8, IERROR)
      RETURN
      END SUBROUTINE MUMPS_SET_IERROR
      SUBROUTINE MUMPS_STOREI8(I8, INT_ARRAY)
      IMPLICIT NONE
      INTEGER(8), intent(in)  :: I8
      INTEGER,    intent(out) :: INT_ARRAY(2)
      INTEGER(kind(0_4)) :: I32
      INTEGER(8) :: IDIV, IPAR
      PARAMETER (IPAR=int(huge(I32),8))
      PARAMETER (IDIV=IPAR+1_8)
      IF ( I8 .LT. IDIV ) THEN
        INT_ARRAY(1) = 0
        INT_ARRAY(2) = int(I8)
      ELSE
        INT_ARRAY(1) = int(I8 / IDIV)
        INT_ARRAY(2) = int(mod(I8,IDIV))
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_STOREI8
      SUBROUTINE MUMPS_GETI8(I8, INT_ARRAY)
      IMPLICIT NONE
      INTEGER(8), intent(out)  :: I8
      INTEGER,    intent(in)  :: INT_ARRAY(2)
      INTEGER(kind(0_4)) :: I32
      INTEGER(8) :: IDIV, IPAR
      PARAMETER (IPAR=int(huge(I32),8))
      PARAMETER (IDIV=IPAR+1_8)
      IF ( INT_ARRAY(1) .EQ. 0 ) THEN
        I8=int(INT_ARRAY(2),8)
      ELSE
        I8=int(INT_ARRAY(1),8)*IDIV+int(INT_ARRAY(2),8)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_GETI8
      SUBROUTINE MUMPS_ADDI8TOARRAY( INT_ARRAY, I8 )
      IMPLICIT NONE
      INTEGER(8), intent(in) :: I8
      INTEGER, intent(inout) :: INT_ARRAY(2)
      INTEGER(8) :: I8TMP
      CALL MUMPS_GETI8(I8TMP, INT_ARRAY)
      I8TMP = I8TMP + I8
      CALL MUMPS_STOREI8(I8TMP, INT_ARRAY)
      RETURN
      END SUBROUTINE MUMPS_ADDI8TOARRAY
      SUBROUTINE MUMPS_SUBTRI8TOARRAY( INT_ARRAY, I8 )
      IMPLICIT NONE
      INTEGER(8), intent(in) :: I8
      INTEGER, intent(inout) :: INT_ARRAY(2)
      INTEGER(8) :: I8TMP
      CALL MUMPS_GETI8(I8TMP, INT_ARRAY)
      I8TMP = I8TMP - I8
      CALL MUMPS_STOREI8(I8TMP, INT_ARRAY)
      RETURN
      END SUBROUTINE MUMPS_SUBTRI8TOARRAY
      FUNCTION MUMPS_SEQANA_AVAIL(ICNTL7)
      LOGICAL :: MUMPS_SEQANA_AVAIL
      INTEGER, INTENT(IN) :: ICNTL7
      LOGICAL :: SCOTCH=.FALSE.
      LOGICAL :: METIS =.FALSE.
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
      METIS = .TRUE.
#endif
#if defined(scotch) || defined(ptscotch)
      SCOTCH = .TRUE.
#endif
      IF ( ICNTL7 .LT. 0 .OR. ICNTL7 .GT. 7 ) THEN
        MUMPS_SEQANA_AVAIL = .FALSE.
      ELSE
        MUMPS_SEQANA_AVAIL = .TRUE.
      ENDIF
      IF ( ICNTL7 .EQ. 5 ) MUMPS_SEQANA_AVAIL = METIS
      IF ( ICNTL7 .EQ. 3 ) MUMPS_SEQANA_AVAIL = SCOTCH
      RETURN
      END FUNCTION MUMPS_SEQANA_AVAIL
      FUNCTION MUMPS_PARANA_AVAIL(WHICH)
      LOGICAL :: MUMPS_PARANA_AVAIL
      CHARACTER :: WHICH*(*)
      LOGICAL :: PTSCOTCH=.FALSE., PARMETIS=.FALSE.
#if defined(ptscotch)
      PTSCOTCH = .TRUE.
#endif
#if defined(parmetis) || defined(parmetis3)
      PARMETIS = .TRUE.
#endif
      SELECT CASE(WHICH)
      CASE('ptscotch','PTSCOTCH')
         MUMPS_PARANA_AVAIL = PTSCOTCH
      CASE('parmetis','PARMETIS')
         MUMPS_PARANA_AVAIL = PARMETIS
      CASE('both','BOTH')
         MUMPS_PARANA_AVAIL = PTSCOTCH .AND. PARMETIS
      CASE('any','ANY')
         MUMPS_PARANA_AVAIL = PTSCOTCH .OR. PARMETIS
      CASE default
         write(*,'("Invalid input in MUMPS_PARANA_AVAIL")')
      END SELECT
      RETURN
      END FUNCTION MUMPS_PARANA_AVAIL
      SUBROUTINE MUMPS_SORT_STEP(N,FRERE,STEP,FILS,
     &     NA,LNA,NE,ND,DAD,LDAD,USE_DAD,
     &     NSTEPS,INFO,LP,
     &     PROCNODE,SLAVEF
     &     )
      IMPLICIT NONE
      INTEGER N, NSTEPS, LNA, LP,LDAD
      INTEGER FRERE(NSTEPS), FILS(N), STEP(N)
      INTEGER NA(LNA), NE(NSTEPS), ND(NSTEPS)
      INTEGER DAD(LDAD)
      LOGICAL USE_DAD
      INTEGER INFO(40)
      INTEGER SLAVEF,PROCNODE(NSTEPS)
      INTEGER  POSTORDER,TMP_SWAP
      INTEGER, DIMENSION (:), ALLOCATABLE :: STEP_TO_NODE
      INTEGER, DIMENSION (:), ALLOCATABLE :: IPOOL,TNSTK
      INTEGER I,II,allocok
      INTEGER NBLEAF,NBROOT,LEAF,IN,INODE,IFATH
      EXTERNAL MUMPS_TYPENODE
      INTEGER MUMPS_TYPENODE
      POSTORDER=1
      NBLEAF = NA(1)
      NBROOT = NA(2)
      ALLOCATE( IPOOL(NBLEAF), TNSTK(NSTEPS), stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in CMUMPS_SORT_STEP'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      DO I=1,NSTEPS
         TNSTK(I) = NE(I)
      ENDDO
      ALLOCATE(STEP_TO_NODE(NSTEPS),stat=allocok)
      IF (allocok > 0) THEN
         IF ( LP .GT. 0 )
     &        WRITE(LP,*)'Memory allocation error in
     &CMUMPS_REORDER_TREE'
         INFO(1)=-7
         INFO(2)=NSTEPS
         RETURN
      ENDIF
      DO I=1,N
         IF(STEP(I).GT.0)THEN
            STEP_TO_NODE(STEP(I))=I
         ENDIF
      ENDDO
      IPOOL(1:NBLEAF)=NA(3:2+NBLEAF)
      LEAF = NBLEAF + 1
 91   CONTINUE
      IF (LEAF.NE.1) THEN
         LEAF = LEAF -1
         INODE = IPOOL(LEAF)
      ENDIF
 96   CONTINUE
      IF (USE_DAD) THEN
         IFATH = DAD( STEP(INODE) )
      ELSE
         IN = INODE
 113     IN = FRERE(IN)
         IF (IN.GT.0) GO TO 113
         IFATH = -IN
      ENDIF
      TMP_SWAP=FRERE(STEP(INODE))
      FRERE(STEP(INODE))=FRERE(POSTORDER)
      FRERE(POSTORDER)=TMP_SWAP
      TMP_SWAP=ND(STEP(INODE))
      ND(STEP(INODE))=ND(POSTORDER)
      ND(POSTORDER)=TMP_SWAP
      TMP_SWAP=NE(STEP(INODE))
      NE(STEP(INODE))=NE(POSTORDER)
      NE(POSTORDER)=TMP_SWAP
      TMP_SWAP=PROCNODE(STEP(INODE))
      PROCNODE(STEP(INODE))=PROCNODE(POSTORDER)
      PROCNODE(POSTORDER)=TMP_SWAP
      IF(USE_DAD)THEN
         TMP_SWAP=DAD(STEP(INODE))
         DAD(STEP(INODE))=DAD(POSTORDER)
         DAD(POSTORDER)=TMP_SWAP
      ENDIF
      TMP_SWAP=TNSTK(STEP(INODE))
      TNSTK(STEP(INODE))=TNSTK(POSTORDER)
      TNSTK(POSTORDER)=TMP_SWAP
      II=STEP_TO_NODE(POSTORDER)
      TMP_SWAP=STEP(INODE)
      STEP(STEP_TO_NODE(POSTORDER))=TMP_SWAP
      STEP(INODE)=POSTORDER
      STEP_TO_NODE(POSTORDER)=INODE
      STEP_TO_NODE(TMP_SWAP)=II
      IN=II
 101  IN = FILS(IN)
      IF (IN .GT. 0 ) THEN
         STEP(IN)=-STEP(II)
         GOTO 101
      ENDIF
      IN=INODE
 102  IN = FILS(IN)
      IF (IN .GT. 0 ) THEN
         STEP(IN)=-STEP(INODE)
         GOTO 102
      ENDIF
      POSTORDER = POSTORDER + 1
      IF (IFATH.EQ.0) THEN
         NBROOT = NBROOT - 1
         IF (NBROOT.EQ.0) GOTO 116
         GOTO 91
      ENDIF
      TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
      IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN      
         INODE = IFATH
         GOTO 96
      ELSE
         GOTO 91
      ENDIF
 116  CONTINUE
      DEALLOCATE(STEP_TO_NODE)
      DEALLOCATE(IPOOL,TNSTK)
      RETURN
      END SUBROUTINE MUMPS_SORT_STEP
#if ! defined(NO_XXNBPR)
      SUBROUTINE CHECK_EQUAL(NBPR, IWNBPR)
      IMPLICIT NONE
      INTEGER, intent(in) :: NBPR, IWNBPR
      IF (NBPR .NE. IWNBPR) THEN
        WRITE(*,*) " NBPROCFILS(...), IW(..+XXNBPR_ = ", NBPR, IWNBPR
#if ! defined(IBC_TEST)
        CALL MUMPS_ABORT()
#endif
      ENDIF
      RETURN
      END SUBROUTINE CHECK_EQUAL
#endif
