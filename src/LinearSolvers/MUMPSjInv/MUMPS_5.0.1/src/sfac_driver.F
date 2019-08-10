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
      SUBROUTINE SMUMPS_FAC_DRIVER( id)
      USE SMUMPS_COMM_BUFFER
      USE SMUMPS_LOAD
      USE SMUMPS_OOC
      USE SMUMPS_STRUC_DEF
      USE MUMPS_FRONT_DATA_MGT_M
#if ! defined(NO_FDM_DESCBAND)
      USE MUMPS_FAC_DESCBAND_DATA_M
#endif
#if ! defined(NO_FDM_MAPROW)
      USE MUMPS_FAC_MAPROW_DATA_M
#endif
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C  Performs scaling, sorting in arrowhead, then
C  distributes the matrix, and perform
C  factorization.
C
C
      INTERFACE
C     Explicit interface needed because
C     of "id" derived datatype argument
      SUBROUTINE SMUMPS_ANORMINF(id, ANORMINF, LSCAL)
      USE SMUMPS_STRUC_DEF
      TYPE (SMUMPS_STRUC), TARGET :: id
      REAL, INTENT(OUT) :: ANORMINF
      LOGICAL :: LSCAL
      END SUBROUTINE SMUMPS_ANORMINF
C
      END INTERFACE
C
C  Parameters
C  ==========
C
      TYPE(SMUMPS_STRUC), TARGET :: id
C
C  MPI
C  ===
C
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      INTEGER, PARAMETER :: MASTER = 0
C
C  Local variables
C  ===============
C
      INCLUDE 'mumps_headers.h'
      INTEGER NSEND, NSEND_TOT, LDPTRAR, NELT
      INTEGER NLOCAL, NLOCAL_TOT, KEEP13_SAVE, ITMP
      INTEGER(8) K67
      INTEGER(8) ITMP8
      INTEGER  MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INTEGER MP, LP, MPG, allocok
      LOGICAL PROK, PROKG, LSCAL, LPOK
      INTEGER SMUMPS_LBUF, SMUMPS_LBUFR_BYTES, SMUMPS_LBUF_INT
      INTEGER(8) SMUMPS_LBUFR_BYTES8, SMUMPS_LBUF8
      INTEGER PTRIST, PTRWB, MAXELT_SIZE,
     &     ITLOC, IPOOL, NSTEPS, K28, LPOOL, LIW
      INTEGER IRANK, ID_ROOT
      INTEGER KKKK, NZ_locMAX
      INTEGER(8) MEMORY_MD_ARG
      INTEGER(8) MAXS_BASE8, MAXS_BASE_RELAXED8
      REAL CNTL4
      INTEGER MIN_PERLU, MAXIS_ESTIM
      INTEGER   MAXIS
      INTEGER(8) :: MAXS
      DOUBLE PRECISION TIME, TIMEET
      REAL ZERO, ONE, MONE
      PARAMETER( ZERO = 0.0E0, ONE = 1.0E0, MONE = -1.0E0)
      REAL CZERO
      PARAMETER( CZERO = 0.0E0 )
      INTEGER PERLU, TOTAL_MBYTES, K231, K232, K233
      INTEGER COLOUR, COMM_FOR_SCALING ! For Simultaneous scaling
      INTEGER LIWK, LWK, LWK_REAL
C     SLAVE: used to determine if proc has the role of a slave
      LOGICAL I_AM_SLAVE, PERLU_ON, WK_USER_PROVIDED
C     WK_USER_PROVIDED is set to true when workspace WK_USER is provided by user
      REAL :: ANORMINF, SEUIL, SEUIL_LDLT_NIV2
      REAL :: CNTL1, CNTL3, CNTL5, CNTL6, EPS
      INTEGER N, LPN_LIST,POSBUF
      INTEGER, DIMENSION (:), ALLOCATABLE :: ITMP2
      INTEGER I,K
C
C  Workspace.
C
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWK
      REAL, DIMENSION(:), ALLOCATABLE :: WK
      REAL, DIMENSION(:), ALLOCATABLE :: WK_REAL
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IWK8
      INTEGER, DIMENSION(:), ALLOCATABLE :: BURP
      INTEGER, DIMENSION(:), ALLOCATABLE :: BUCP
      INTEGER, DIMENSION(:), ALLOCATABLE :: BURS
      INTEGER, DIMENSION(:), ALLOCATABLE :: BUCS
      INTEGER BUREGISTRE(12)
      INTEGER BUINTSZ, BURESZ, BUJOB
      INTEGER BUMAXMN, M, SCMYID, SCNPROCS
      REAL    SCONEERR, SCINFERR
C
C  Parameters arising from the structure
C  =====================================
C
      INTEGER, POINTER ::  JOB, NZ
*     Control parameters: see description in SMUMPSID
      REAL,DIMENSION(:),POINTER::RINFO, RINFOG
      REAL,DIMENSION(:),POINTER::    CNTL
      INTEGER,DIMENSION(:),POINTER::INFO, INFOG, KEEP
      INTEGER, DIMENSION(:), POINTER :: MYIRN_loc, MYJCN_loc
      REAL, DIMENSION(:), POINTER :: MYA_loc
      INTEGER, TARGET :: DUMMYIRN_loc(1), DUMMYJCN_loc(1)
      REAL, TARGET :: DUMMYA_loc(1)
      INTEGER(8),DIMENSION(:),POINTER::KEEP8
      INTEGER,DIMENSION(:),POINTER::ICNTL
      EXTERNAL MUMPS_GET_POOL_LENGTH
      INTEGER MUMPS_GET_POOL_LENGTH
      INTEGER(8) TOTAL_BYTES
      INTEGER(8) :: I8TMP
C
C  External references
C  ===================
      INTEGER numroc
      EXTERNAL numroc
C  Fwd in facto:
      REAL, DIMENSION(:), POINTER :: RHS_MUMPS
      LOGICAL :: RHS_MUMPS_ALLOCATED
      INTEGER :: NB_ACTIVE_FRONTS_ESTIM
C
C 
      JOB=>id%JOB
      NZ=>id%NZ
      RINFO=>id%RINFO
      RINFOG=>id%RINFOG
      CNTL=>id%CNTL
      INFO=>id%INFO
      INFOG=>id%INFOG
      KEEP=>id%KEEP
      KEEP8=>id%KEEP8
      ICNTL=>id%ICNTL
      IF (id%NZ_loc .NE. 0) THEN
        MYIRN_loc=>id%IRN_loc
        MYJCN_loc=>id%JCN_loc
        MYA_loc=>id%A_loc
      ELSE
        MYIRN_loc=>DUMMYIRN_loc
        MYJCN_loc=>DUMMYJCN_loc
        MYA_loc=>DUMMYA_loc
      ENDIF
      N = id%N
      EPS = epsilon ( ZERO )
C     TIMINGS: reset to 0
      id%DKEEP(92)=0.0E0
      id%DKEEP(93)=0.0E0
      id%DKEEP(94)=0.0E0
      id%DKEEP(97)=0.0E0
      id%DKEEP(98)=0.0E0
C     Related to forward in facto functionality (referred to as "Fwd in facto")
      NULLIFY(RHS_MUMPS)
      RHS_MUMPS_ALLOCATED = .FALSE.
C     -----------------------------------------------------------------------
C     Set WK_USER_PROVIDED to true when workspace WK_USER is provided by user
C     We can accept WK_USER to be provided on only one proc and 
C     different values of WK_USER per processor
C     
      IF (KEEP8(24).GT.0_8) THEN 
C     IMPORTANT: if on entry KEEP8(24) is non zero then
C                it means that WK_USER has been used in a previous phase.
C                In this case S should have already been nullified.
C                We nullify S so that later when we test
C                if (associated(S) we can free space and reallocate it).
           NULLIFY(id%S)
      ENDIF
C
C     --  KEEP8(24) can now then be reset safely
      WK_USER_PROVIDED = (id%LWK_USER.NE.0)
      IF (WK_USER_PROVIDED) THEN
          IF (id%LWK_USER.GT.0) THEN
            KEEP8(24) = int(id%LWK_USER,8)
          ELSE
            KEEP8(24) = -int(id%LWK_USER,8)* 1000000_8 
          ENDIF
      ELSE
          KEEP8(24) = 0_8
      ENDIF
C
C     KEEP(13) might be modified
C       (elmement entry format) 
C       but need be restore for 
C       furture factorisation
C       with different scaling option
C 
      KEEP13_SAVE = KEEP(13)
C     In case of loop on factorization with
C     different scaling options, initialize
C     DKEEP(4:5) to 0.
      id%DKEEP(4)=-1.0E0
      id%DKEEP(5)=-1.0E0
C  Mapping information used during solve. In case of several facto+solve
C  it has to be recomputed. In case of several solves with the same
C  facto, it is not recomputed.
      IF (associated(id%IPTR_WORKING)) THEN
        DEALLOCATE(id%IPTR_WORKING)
        NULLIFY(id%IPTR_WORKING)
      END IF
      IF (associated(id%WORKING)) THEN 
        DEALLOCATE(id%WORKING)
        NULLIFY(id%WORKING)
      END IF
C
C  Units for printing
C  MP: diagnostics
C  LP: errors
C
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LPOK    = ((LP.GT.0).AND.(id%ICNTL(4).GE.1))
      PROK    = ((MP.GT.0).AND.(id%ICNTL(4).GE.2))
      PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG   = (PROKG.AND.(id%ICNTL(4).GE.2))
      IF ( PROK ) WRITE( MP, 130 )
      IF ( PROKG ) WRITE( MPG, 130 )
      IF ( PROKG .and. KEEP(53).GT.0 ) THEN
        WRITE(MPG,'(/A,I3)') ' Null space option :', KEEP(19)
        IF ( KEEP(21) .ne. N ) THEN 
          WRITE( MPG, '(A,I10)') ' Max deficiency    : ', KEEP(21)
        END IF
        IF ( KEEP(22) .ne. 0 ) THEN 
          WRITE( MPG, '(A,I10)') ' Min deficiency    : ', KEEP(22)
        END IF
      END IF
C     -------------------------------------
C     Depending on the type of parallelism,
C     the master can now (soon) potentially
C     have the role of a slave
C     -------------------------------------
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &             ( id%MYID .eq. MASTER .AND.
     &               KEEP(46) .eq. 1 ) )
C
C  Prepare work for out-of-core
C
        IF (id%MYID .EQ. MASTER .AND. KEEP(201) .NE. -1) THEN
C         Note that if KEEP(201)=-1, then we have decided
C         at analysis phase that factors will not be stored
C         (neither in memory nor on disk). In that case,
C         ICNTL(22) is ignored.
C         -- ICNTL(22) must be set before facto phase 
C            (=1 OOC on; =0 OOC off)
C            and cannot be changed for subsequent solve phases.
          KEEP(201)=id%ICNTL(22)
          IF (KEEP(201) .NE. 0) THEN !Later: .GT. to allow ICNTL(22)=-1
#           if defined(OLD_OOC_NOPANEL)
              KEEP(201)=2
#           else
              KEEP(201)=1
#           endif
          ENDIF
        ENDIF
C       ----------------------
C       Broadcast KEEP options
C       defined for facto:
C       ----------------------
        CALL MPI_BCAST( KEEP(12), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(19), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(21), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(201), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        IF (id%MYID.EQ.MASTER) THEN
          IF (KEEP(217).GT.2.OR.KEEP(217).LT.0) THEN
            KEEP(217)=0
          ENDIF
          KEEP(214)=KEEP(217)
          IF (KEEP(214).EQ.0) THEN
            IF (KEEP(201).NE.0) THEN ! OOC or no factors
              KEEP(214)=1
            ELSE
              KEEP(214)=2
            ENDIF
          ENDIF
        ENDIF
        CALL MPI_BCAST( KEEP(214), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        IF (KEEP(201).NE.0) THEN
C         -- Low Level I/O strategy
          CALL MPI_BCAST( KEEP(99), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(205), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(211), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        ENDIF
C 
C
C       KEEP(50)  case
C       ==============
C
C          KEEP(50)  = 0 : matrix is unsymmetric
C          KEEP(50) /= 0 : matrix is symmetric
C          KEEP(50) = 1 : Ask L L^T on the root. Matrix is PSD.
C          KEEP(50) = 2 : Ask for L U on the root
C          KEEP(50) = 3 ... L D L^T ??
C
C       ---------------------------------------
C       For symmetric (non general) matrices 
C       set (directly) CNTL(1) = 0.0
C       ---------------------------------------
        IF ( KEEP(50) .eq. 1 ) THEN
          IF (id%CNTL(1) .ne. ZERO ) THEN
            IF ( MPG .GT. 0 ) THEN
              WRITE(MPG,'(A)')
     &' ** Warning : SPD solver called, resetting CNTL(1) to 0.0E0'
            END IF
          END IF
          id%CNTL(1) = ZERO
        END IF
      IF (KEEP(219).NE.0) THEN
       CALL SMUMPS_BUF_MAX_ARRAY_MINSIZE(max(KEEP(108),1),IERR)
       IF (IERR .NE. 0) THEN
C      ------------------------
C      Error allocating SMUMPS_BUF
C      ------------------------
          INFO(1) = -13
          INFO(2) = max(KEEP(108),1)
       END IF
      ENDIF
C     Fwd in facto: explicitly forbid
C     sparse RHS and A-1 computation
      IF (id%KEEP(252).EQ.1 .AND. id%MYID.EQ.MASTER) THEN
        IF (id%ICNTL(20).EQ.1) THEN ! out-of-range => 0
C         NB: in doc ICNTL(20) only accessed during solve
C         In practice, will have failed earlier if RHS not allocated.
C         Still it looks safer to keep this test.
          id%INFO(1)=-43
          id%INFO(2)=20
          IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: Sparse RHS is incompatible with forward',
     &       ' performed during factorization (ICNTL(32)=1)'
        ELSE IF (id%ICNTL(30).NE.0) THEN ! out-of-range => 1
          id%INFO(1)=-43
          id%INFO(2)=30
          IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: A-1 functionality incompatible with forward',
     &       ' performed during factorization (ICNTL(32)=1)'
        ELSE IF (id%ICNTL(9) .NE. 1) THEN
          id%INFO(1)=-43
          id%INFO(2)=9
          IF (PROKG) WRITE(MPG,'(A)')
     &       ' ERROR: Transpose system (ICNTL(9).NE.0) not ',
     &       ' compatible with forward performed during',
     &       ' factorization (ICNTL(32)=1)'
        ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
C
      IF (INFO(1).LT.0) GOTO 530
      IF ( PROKG ) THEN
          WRITE( MPG, 172 ) id%NSLAVES, id%ICNTL(22),
     &    KEEP8(111), KEEP(126), KEEP(127), KEEP(28),
     &    KEEP8(4)/1000000_8
          IF (KEEP(252).GT.0) 
     &    WRITE(MPG,173) KEEP(253)
      ENDIF
      IF (KEEP(201).LE.0) THEN
C       In-core version or no factors
        KEEP(IXSZ)=XSIZE_IC
      ELSE IF (KEEP(201).EQ.2) THEN
C       OOC version, no panels
        KEEP(IXSZ)=XSIZE_OOC_NOPANEL
      ELSE IF (KEEP(201).EQ.1) THEN
C     Panel versions:
        IF (KEEP(50).EQ.0) THEN
          KEEP(IXSZ)=XSIZE_OOC_UNSYM
        ELSE
          KEEP(IXSZ)=XSIZE_OOC_SYM
        ENDIF
      ENDIF
C
*     **********************************
*     Begin intializations regarding the
*     computation of the determinant
*     **********************************
      IF (id%MYID.EQ.MASTER) KEEP(258)=ICNTL(33)
      CALL MPI_BCAST(KEEP(258), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      IF (KEEP(258) .NE. 0) THEN
        KEEP(259) = 0      ! Initial exponent of the local determinant
        KEEP(260) = 1      ! Number of permutations
        id%DKEEP(6)  = 1.0E0  ! real part of the local determinant
      ENDIF
*     ********************************
*     End intializations regarding the
*     computation of the determinant
*     ********************************
C
*     **********************
*     Begin of Scaling phase
*     **********************
C
C     SCALING MANAGEMENT
C     * Options 1, 3, 4 centralized only
C  
C     * Options 7, 8  : also works for distributed matrix
C
C     At this point, we have the scaling arrays allocated
C     on the master. They have been allocated on the master
C     inside the main MUMPS driver.
C
      CALL MPI_BCAST(KEEP(52), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      LSCAL = ((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 8))
      IF (LSCAL) THEN
C
        IF ( id%MYID.EQ.MASTER ) THEN
          CALL MUMPS_SECDEB(TIMEET)
        ENDIF
C       -----------------------
C       Retrieve parameters for
C       simultaneous scaling
C       -----------------------
        IF (KEEP(52) .EQ. 7) THEN 
C       -- Cheap setting of SIMSCALING (it is the default in 4.8.4)
           K231= KEEP(231)
           K232= KEEP(232)
           K233= KEEP(233)
        ELSEIF (KEEP(52) .EQ. 8) THEN
C       -- More expensive setting of SIMSCALING (it was the default in 4.8.1,2,3)
           K231= KEEP(239)
           K232= KEEP(240)
           K233= KEEP(241)
        ENDIF
        CALL MPI_BCAST(id%DKEEP(3),1,MPI_REAL,MASTER,
     &       id%COMM,IERR)
C
        IF ( ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)) .AND. 
     &       KEEP(54).NE.0 ) THEN
C         ------------------------------
C         Scaling for distributed matrix
C         We need to allocate scaling
C         arrays on all processors, not
C         only the master.
C         ------------------------------
           IF ( id%MYID .NE. MASTER ) THEN
              IF ( associated(id%COLSCA))
     &             DEALLOCATE( id%COLSCA )
              IF ( associated(id%ROWSCA))
     &             DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=N
            ENDIF
            ALLOCATE( id%ROWSCA(N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=N
            ENDIF
         ENDIF
         M = N
         BUMAXMN=M
         IF(N > BUMAXMN) BUMAXMN = N
         LIWK = 4*BUMAXMN
         ALLOCATE (IWK(LIWK),BURP(M),BUCP(N),
     &            BURS(2* (id%NPROCS)),BUCS(2* (id%NPROCS)),
     &            stat=allocok)
         IF (allocok > 0) THEN
            INFO(1)=-13
            INFO(2)=LIWK+M+N+4* (id%NPROCS)
         ENDIF
C        --- Propagate enventual error
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &        id%COMM, id%MYID )
         IF (INFO(1).LT.0) GOTO 530
C        -- estimation of memory and construction of partvecs
         BUJOB = 1
C        -- LWK not used
         LWK_REAL   = 1
         ALLOCATE(WK_REAL(LWK_REAL))
         CALL SMUMPS_SIMSCALEABS(
     &        MYIRN_loc(1), MYJCN_loc(1), MYA_loc(1),
     &        id%NZ_loc,
     &        M, N,  id%NPROCS, id%MYID, id%COMM,
     &        BURP, BUCP,
     &        BURS, BUCS, BUREGISTRE,
     &        IWK, LIWK,
     &        BUINTSZ, BURESZ, BUJOB,
     &        id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &        id%KEEP(50),
     &        K231, K232, K233, 
     &        id%DKEEP(3),
     &        SCONEERR, SCINFERR)
         IF(LIWK < BUINTSZ) THEN
            DEALLOCATE(IWK)
            LIWK = BUINTSZ
            ALLOCATE(IWK(LIWK), stat=allocok)
            IF (allocok > 0) THEN
               INFO(1)=-13
               INFO(2)=LIWK
            ENDIF
         ENDIF
         LWK_REAL = BURESZ
         DEALLOCATE(WK_REAL)
         ALLOCATE (WK_REAL(LWK_REAL), stat=allocok)
         IF (allocok > 0) THEN
            INFO(1)=-13
            INFO(2)=LWK_REAL
         ENDIF
C        --- Propagate enventual error
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &        id%COMM, id%MYID )
         IF (INFO(1).LT.0) GOTO 530
C        -- estimation of memory and construction of partvecs
         BUJOB = 2   
         CALL SMUMPS_SIMSCALEABS(
     &        MYIRN_loc(1), MYJCN_loc(1), MYA_loc(1),
     &        id%NZ_loc,
     &        M, N,  id%NPROCS, id%MYID, id%COMM,
     &        BURP, BUCP,
     &        BURS, BUCS, BUREGISTRE,
     &        IWK, LIWK,
     &        BUINTSZ, BURESZ, BUJOB,
     &        id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &        id%KEEP(50),
     &        K231, K232, K233, 
     &        id%DKEEP(3),
     &        SCONEERR, SCINFERR)
         id%DKEEP(4) = SCONEERR
         id%DKEEP(5) = SCINFERR
CXXXX 
         DEALLOCATE(IWK, WK_REAL,BURP,BUCP,BURS, BUCS)
        ELSE IF ( KEEP(54) .EQ. 0 ) THEN
C         ------------------
C         Centralized matrix
C         ------------------
          IF ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8))  THEN
C             -------------------------------
C             Create a communicator of size 1
C             -------------------------------
              IF (id%MYID.EQ.MASTER) THEN
                COLOUR = 0
              ELSE
                COLOUR = MPI_UNDEFINED
              ENDIF
              CALL MPI_COMM_SPLIT( id%COMM, COLOUR, 0,
     &             COMM_FOR_SCALING, IERR )
              IF (id%MYID.EQ.MASTER) THEN
                 M = N
                 BUMAXMN=N
CXXXX 
                 IF(N > BUMAXMN) BUMAXMN = N
                 LIWK = 1
                 ALLOCATE (IWK(LIWK),BURP(1),BUCP(1),
     &                BURS(1),BUCS(1),
     &                stat=allocok)
                 LWK_REAL = M + N  
                 ALLOCATE (WK_REAL(LWK_REAL), stat=allocok)
                 IF (allocok > 0) THEN
                    INFO(1)=-13
                    INFO(2)=1
                 ENDIF
                 IF (INFO(1) .LT. 0) GOTO 400
                 CALL MPI_COMM_RANK(COMM_FOR_SCALING, SCMYID, IERR)
                 CALL MPI_COMM_SIZE(COMM_FOR_SCALING, SCNPROCS, IERR)
                 BUJOB = 1
                 CALL SMUMPS_SIMSCALEABS(
     &                id%IRN(1), id%JCN(1), id%A(1),
     &                id%NZ,
     &                M, N,  SCNPROCS, SCMYID, COMM_FOR_SCALING,
     &                BURP, BUCP,
     &                BURS, BUCS, BUREGISTRE,
     &                IWK, LIWK,
     &                BUINTSZ, BURESZ, BUJOB,
     &                id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &                id%KEEP(50),
     &                K231, K232, K233, 
     &                id%DKEEP(3),
     &                SCONEERR, SCINFERR)
                 IF(LWK_REAL < BURESZ) THEN
                    ! internal error since LWK=BURESZ=M+N
                    INFO(1) = -136
                    GOTO 400
                 ENDIF
                 BUJOB = 2
                 CALL SMUMPS_SIMSCALEABS(id%IRN(1),
     &                id%JCN(1), id%A(1),
     &                id%NZ,
     &                M, N,  SCNPROCS, SCMYID, COMM_FOR_SCALING,
     &                BURP, BUCP,
     &                BURS, BUCS, BUREGISTRE,
     &                IWK, LIWK,
     &                BUINTSZ, BURESZ, BUJOB,
     &                id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &                id%KEEP(50),
     &                K231, K232, K233, 
     &                id%DKEEP(3),
     &                SCONEERR, SCINFERR)
                 id%DKEEP(4) = SCONEERR
                 id%DKEEP(5) = SCINFERR
CXXXX 
                 DEALLOCATE(WK_REAL)                 
                 DEALLOCATE (IWK,BURP,BUCP,
     &                BURS,BUCS)
              ENDIF
C             Centralized matrix: make DKEEP(4:5) available to all processors
              CALL MPI_BCAST( id%DKEEP(4),2,MPI_REAL,
     &                        MASTER, id%COMM, IERR )
  400         CONTINUE
              IF (id%MYID.EQ.MASTER) THEN
C               Communicator should only be
C               freed on the master process
                CALL MPI_COMM_FREE(COMM_FOR_SCALING, IERR)
              ENDIF
              CALL MUMPS_PROPINFO(ICNTL(1), INFO(1), id%COMM, id%MYID)
              IF (INFO(1).LT.0) GOTO 530
          ELSE IF (id%MYID.EQ.MASTER) THEN
C           ----------------------------------
C           Centralized scaling, options 1 to 6
C           ----------------------------------
            IF (KEEP(52).GT.0 .AND. KEEP(52).LE.6) THEN
C             ---------------------
C             Allocate temporary
C             workspace for scaling
C             ---------------------
              IF ( KEEP(52) .eq. 5 .or. 
     &          KEEP(52) .eq. 6 ) THEN
C               We have an explicit copy of the original
C               matrix in complex format which should probably
C               be avoided (but do we want to keep all
C               those old scaling options ?)
                LWK = NZ
              ELSE
                LWK = 1
              END IF
              LWK_REAL = 5 * N
              ALLOCATE( WK_REAL( LWK_REAL ), stat = IERR )
              IF ( IERR .GT. 0 ) THEN
                INFO(1) = -13
                INFO(2) = LWK_REAL
                GOTO 137
              END IF
              ALLOCATE( WK( LWK ), stat = IERR )
              IF ( IERR .GT. 0 ) THEN
                INFO(1) = -13
                INFO(2) = LWK
                GOTO 137
              END IF
              CALL SMUMPS_FAC_A(N, NZ, KEEP(52), id%A(1),
     &             id%IRN(1), id%JCN(1),
     &             id%COLSCA(1), id%ROWSCA(1),
     &             WK, LWK, WK_REAL, LWK_REAL, ICNTL(1), INFO(1) )
              DEALLOCATE( WK_REAL )
              DEALLOCATE( WK )
            ENDIF
          ENDIF
        ENDIF ! Scaling distributed matrices or centralized
        IF (id%MYID.EQ.MASTER) THEN
            CALL MUMPS_SECFIN(TIMEET)
            id%DKEEP(92)=real(TIMEET)
C         Print inf-norm after last KEEP(233) iterations of
C         scaling option KEEP(52)=7 or 8 (SimScale)
C
          IF (PROKG.AND.(KEEP(52).EQ.7.OR.KEEP(52).EQ.8) 
     &             .AND. (K233+K231+K232).GT.0) THEN
           IF (K232.GT.0) WRITE(MPG, 166) id%DKEEP(4)
          ENDIF
        ENDIF
      ENDIF ! LSCAL
C
C       scaling might also be provided by the user
        LSCAL = (LSCAL .OR. (KEEP(52) .EQ. -1) .OR. KEEP(52) .EQ. -2)
        IF (LSCAL .AND. KEEP(258).NE.0 .AND. id%MYID .EQ. MASTER) THEN
          DO I = 1, id%N
            CALL SMUMPS_UPDATEDETER_SCALING(id%ROWSCA(I),
     &           id%DKEEP(6),    ! determinant
     &           KEEP(259))   ! exponent of the determinant
          ENDDO
          IF (KEEP(50) .EQ. 0) THEN ! unsymmetric
            DO I = 1, id%N
              CALL SMUMPS_UPDATEDETER_SCALING(id%COLSCA(I),
     &           id%DKEEP(6),    ! determinant
     &           KEEP(259))   ! exponent of the determinant
            ENDDO
          ELSE
C           -----------------------------------------
C           In this case COLSCA = ROWSCA
C           Since determinant was initialized to 1,
C           compute square of the current determinant
C           rather than going through COLSCA.
C           -----------------------------------------
            CALL SMUMPS_DETER_SQUARE(id%DKEEP(6), KEEP(259))
          ENDIF
C         Now we should have taken the
C         inverse of the scaling vectors
          CALL SMUMPS_DETER_SCALING_INVERSE(id%DKEEP(6), KEEP(259))
        ENDIF
C
C       ********************
C       End of Scaling phase
C       At this point: either (matrix is distributed and KEEP(52)=7 or 8)
C       in which case scaling arrays are allocated on all processors,
C       or scaling arrays are only on the host processor.
C       In case of distributed matrix input, we will free the scaling
C       arrays on procs with MYID .NE. 0 after the all-to-all distribution
C       of the original matrix.
C       ********************
C
 137  CONTINUE
C     Fwd in facto: in case of repeated factorizations
C     with different Schur options we prefer to free
C     systematically this array now than waiting for
C     the root node. We rely on the fact that it is
C     allocated or not during the solve phase so if
C     it was allocated in a 1st call to facto and not
C     in a second, we don't want the solve to think
C     it was allocated in the second call.
      IF (associated(id%root%RHS_CNTR_MASTER_ROOT)) THEN
        DEALLOCATE (id%root%RHS_CNTR_MASTER_ROOT)
        NULLIFY (id%root%RHS_CNTR_MASTER_ROOT)
      ENDIF
C     Fwd in facto: check that id%NRHS has not changed
      IF ( id%MYID.EQ.MASTER.AND. KEEP(252).EQ.1 .AND.
     &      id%NRHS .NE. id%KEEP(253) ) THEN
C         Error: NRHS should not have
C         changed since the analysis
          id%INFO(1)=-42
          id%INFO(2)=id%KEEP(253)
      ENDIF
C     Fwd in facto: allocate and broadcast RHS_MUMPS
C     to make it available on all processors.
      IF (id%KEEP(252) .EQ. 1) THEN
          IF ( id%MYID.NE.MASTER ) THEN
            id%KEEP(254) = N              ! Leading dimension
            id%KEEP(255) = N*id%KEEP(253) ! Tot size
            ALLOCATE(RHS_MUMPS(id%KEEP(255)),stat=IERR)
            IF (IERR > 0) THEN
               INFO(1)=-13
               INFO(2)=id%KEEP(255)
               IF (LPOK)
     &         WRITE(LP,*) 'ERREUR while allocating RHS on a slave'
               NULLIFY(RHS_MUMPS)
            ENDIF
            RHS_MUMPS_ALLOCATED = .TRUE.
          ELSE 
C           Case of non working master
            id%KEEP(254)=id%LRHS              ! Leading dimension
            id%KEEP(255)=id%LRHS*(id%KEEP(253)-1)+id%N ! Tot size
            RHS_MUMPS=>id%RHS
            RHS_MUMPS_ALLOCATED = .FALSE.
            IF (LSCAL) THEN
C             Scale before broadcast: apply row
C             scaling (remark that we assume no
C             transpose).
              DO K=1, id%KEEP(253)
                DO I=1, N
                  RHS_MUMPS( id%KEEP(254) * (K-1) + I )
     &          = RHS_MUMPS( id%KEEP(254) * (K-1) + I )
     &          * id%ROWSCA(I)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
      ELSE
          id%KEEP(255)=1
          ALLOCATE(RHS_MUMPS(1))
          RHS_MUMPS_ALLOCATED = .TRUE.
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF ( INFO(1).lt.0 ) GOTO 530
      IF (KEEP(252) .EQ. 1) THEN
C
C         Broadcast the columns of the right-hand side
C         one by one. Leading dimension is keep(254)=N
C         on procs with MYID > 0 but may be larger on
C         the master processor.
          DO I= 1, id%KEEP(253)
            CALL MPI_BCAST(RHS_MUMPS((I-1)*id%KEEP(254)+1), N,
     &           MPI_REAL, MASTER,id%COMM,IERR)
          END DO
      ENDIF
C     Keep a copy of ICNTL(24) and make it
C     available on all working processors.
      KEEP(110)=id%ICNTL(24)
      CALL MPI_BCAST(KEEP(110), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
C     KEEP(110) defaults to 0 for out of range values
      IF (KEEP(110).NE.1) KEEP(110)=0
*
*
C     -----------------------------------------------
C     Depending on the option used for 
C       -detecting null pivots (ICNTL(24)/KEEP(110))
C         CNTL(3) is used to set DKEEP(1)
C               ( A row is considered as null if ||row|| < DKEEP(1) )
C         CNTL(5) is then used to define if a large 
C                 value is set on the diagonal or if a 1 is set
C                 and other values in the row are reset to zeros.
C       -Rank revealing on the Schur (ICNTL(16)/KEEP(19))
C         CNTL(6) is used to set SEUIL and SEUIL_LDLT_NIV2
C         SEUIL* corresponds to the minimum required 
C                absolute value of pivot.
C         SEUIL_LDLT_NIV2 is used only in the 
C                case of SYM=2 within a niv2 node for which 
C                we have only a partial view of the fully summed rows.
C        Note that SEUIL* might be reset later in this routine
C                but only when static pivoting is on 
C                which will be excluded if null pivots or 
C                rank-revealing (RR) is on
C     -----------------------------------------------
      IF (id%MYID .EQ. MASTER) CNTL3 = id%CNTL(3)
      CALL MPI_BCAST(CNTL3, 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL5 = id%CNTL(5)
      CALL MPI_BCAST(CNTL5, 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL6 = id%CNTL(6)
      CALL MPI_BCAST(CNTL6, 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL1 = id%CNTL(1)
      CALL MPI_BCAST(CNTL1, 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
C     -------------------------------------------------------
      ANORMINF = ZERO
C        We compute, when needed, 
C        the infinite norm of Rowsca *A*Colsca
C        and make it available on all working processes.
C     -------------------------------------------------------
C     -------------------------------
C     -- Rank revealing on the root; 
C     -------------------------------
C          -- SEUIL will then be used to postpone
C          -- pivots smaller than SEUIL.
C          -- ||A|| always applied to scale threshold
C          -- in order to be less matrix dependent.
C     -- Set both SEUIL and SEUIL_LDLT_NIV2 
C     -------------------------------
      IF (KEEP(19).EQ.0) THEN 
C        -- RR is off
         SEUIL = ZERO
         id%DKEEP(9) = ZERO
      ELSE
C        -- RR is on
C      July 2012
C      CNTL(3) is the threshold used in the following 
C        to compute the SEUIL used for postponing pivots to root
C      SEUIL*CNTL(6) is then the treshold for null pivot detection 
C      (with 0< CNTL(6) <= 1)
         CALL SMUMPS_ANORMINF(  id , ANORMINF, LSCAL )
         IF (CNTL3 .LT. ZERO) THEN
           SEUIL = abs(CNTL(3))
         ELSE IF  (CNTL3 .GT. ZERO) THEN
           SEUIL = CNTL3*ANORMINF
         ELSE  !  (CNTL(3) .EQ. ZERO) THEN
           SEUIL = N*EPS*ANORMINF  ! standard articles
         ENDIF
         IF (PROKG) WRITE(MPG,*)
     &   ' ABSOLUTE PIVOT THRESHOLD for rank revealing =',SEUIL
      ENDIF
C     After QR with pivoting of root or SVD, diagonal entries 
C     need be analysed to determine null space vectors.
C     Two strategies are provided :
      id%DKEEP(9) = SEUIL 
      IF (id%DKEEP(10).LT.MONE) THEN
         id%DKEEP(10)=MONE
      ELSEIF((id%DKEEP(10).LE.ONE).AND.(id%DKEEP(10).GE.ZERO)) THEN
         id%DKEEP(10)=1000.0E0
      ENDIF
      SEUIL_LDLT_NIV2 = SEUIL
C     -------------------------------
C     -- Null pivot row detection
C     -------------------------------
      IF (KEEP(110).EQ.0) THEN 
C        Initialize DKEEP(1) to a negative value
C        in order to avoid detection of null pivots
C        (test max(AMAX,RMAX,abs(PIVOT)).LE.PIVNUL
C        in SMUMPS_FAC_I, where PIVNUL=DKEEP(1))
         id%DKEEP(1) = -1.0E0
         id%DKEEP(2) = ZERO
      ELSE
         IF (ANORMINF.EQ.ZERO) 
     &       CALL SMUMPS_ANORMINF(  id , ANORMINF, LSCAL )
         IF (KEEP(19).NE.0) THEN
C     RR postponing considers that pivot rows of norm smaller that SEUIL 
C     should be postponed. 
C     Pivot rows smaller than DKEEP(1) are directly added to null space
C     and thus considered as null pivot rows. Thus we define id%DKEEP(1) 
C     relatively to SEUIL (which is based on CNTL(3))
          IF (CNTL(6).GT.0.AND.CNTL(6).LT.1) THEN
C           we want DKEEP(1) < SEUIL
            id%DKEEP(1) = SEUIL*CNTL(6)
          ELSE 
            id%DKEEP(1) = SEUIL* 0.01E0 
          ENDIF
         ELSE
C         We keep strategy currently used in MUMPS 4.10.0
          IF (CNTL3 .LT. ZERO) THEN
           id%DKEEP(1)  = abs(CNTL(3))
          ELSE IF  (CNTL3 .GT. ZERO) THEN
           id%DKEEP(1)  = CNTL3*ANORMINF
          ELSE !  (CNTL(3) .EQ. ZERO) THEN
           id%DKEEP(1)  = 1.0E-5*EPS*ANORMINF
          ENDIF
         ENDIF
         IF (PROKG) WRITE(MPG,*)
     &    ' ZERO PIVOT DETECTION ON, THRESHOLD          =',id%DKEEP(1)
         IF (CNTL5.GT.ZERO) THEN
            id%DKEEP(2) = CNTL5 * ANORMINF
            IF (PROKG) WRITE(MPG,*) 
     &    ' FIXATION FOR NULL PIVOTS                    =',id%DKEEP(2)
         ELSE
            IF (PROKG) WRITE(MPG,*) 'INFINITE FIXATION '
            IF (id%KEEP(50).EQ.0) THEN
C             Unsym
            ! the user let us choose a fixation. set in NEGATIVE
            ! to detect during facto when to set row to zero !
             id%DKEEP(2) = -max(1.0E10*ANORMINF, 
     &                sqrt(huge(ANORMINF))/1.0E8)
            ELSE
C             Sym
            id%DKEEP(2) = ZERO
            ENDIF
         ENDIF
      ENDIF
C     Find id of root node if RR is on 
      IF (KEEP(53).NE.0) THEN
        ID_ROOT =MUMPS_PROCNODE(id%PROCNODE_STEPS(id%STEP(KEEP(20))),
     &                          id%NSLAVES)
        IF ( KEEP( 46 )  .NE. 1 ) THEN
          ID_ROOT = ID_ROOT + 1
        END IF
      ENDIF
C Second pass:  set parameters for null pivot detection
C Allocate PIVNUL_LIST in case of null pivot detection
C and in case of rank revealing
      LPN_LIST = 1
      IF ( associated( id%PIVNUL_LIST) ) DEALLOCATE(id%PIVNUL_LIST)
      IF(KEEP(110) .EQ. 1) THEN
         LPN_LIST = N
      ENDIF
      IF (KEEP(19).NE.0 .AND.
     &   (ID_ROOT.EQ.id%MYID .OR. id%MYID.EQ.MASTER)) THEN
         LPN_LIST = N
      ENDIF
      ALLOCATE( id%PIVNUL_LIST(LPN_LIST),stat = IERR )
      IF ( IERR .GT. 0 ) THEN
        INFO(1)=-13
        INFO(2)=LPN_LIST
      END IF
      id%PIVNUL_LIST(1:LPN_LIST) = 0
      KEEP(109) = 0
C end set parameter for null pivot detection
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF ( INFO(1).lt.0 ) GOTO 530
C   --------------------------------------------------------------
C   STATIC PIVOTING
C     -- Static pivoting only when RR and Null pivot detection OFF
C   --------------------------------------------------------------
      IF ((KEEP(19).EQ.0).AND.(KEEP(110).EQ.0)) THEN
C       -- Set KEEP(97) and compute static pivoting threshold.
        IF (id%MYID .EQ. MASTER) CNTL4 = id%CNTL(4)
        CALL MPI_BCAST( CNTL4, 1, MPI_REAL,
     &                MASTER, id%COMM, IERR )
C 
        IF ( CNTL4 .GE. ZERO ) THEN
         KEEP(97) = 1
         IF ( CNTL4 .EQ. ZERO ) THEN
C           -- set seuil to sqrt(eps)*||A||
            IF(ANORMINF .EQ. ZERO) THEN
               CALL SMUMPS_ANORMINF(  id , ANORMINF, LSCAL )
C               WRITE(*,*) id%MYID,': ANORMINF',ANORMINF
            ENDIF
            SEUIL = sqrt(EPS) * ANORMINF
         ELSE
C           WRITE(*,*) 'id%CNTL(4)',id%CNTL(4)
            SEUIL = CNTL4
         ENDIF
         SEUIL_LDLT_NIV2 = SEUIL
C
        ELSE 
         SEUIL = ZERO
        ENDIF
      ENDIF
C     set number of tiny pivots / 2x2 pivots in types 1 /
C     2x2 pivots in types 2, to zero. This is because the
C     user can call the factorization step several times.
      KEEP(98)  = 0
      KEEP(103) = 0
      KEEP(105) = 0
      MAXS      = 1_8
C
C     The memory allowed is given by ICNTL(23) in Mbytes
C     0 means that nothing is provided.
C     Save memory available, ICNTL(23) in KEEP8(4)
C
      IF ( id%MYID.EQ.MASTER ) THEN
        ITMP = ICNTL(23)
      END IF
      CALL MPI_BCAST( ITMP, 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
C    
C     Ignore ICNTL(23) when WK_USER is provided
c     by resetting ITMP to zero on each proc where WK_USER is provided
      IF (WK_USER_PROVIDED) ITMP = 0
      ITMP8 = int(ITMP, 8)
      KEEP8(4) = ITMP8 * 1000000_8   ! convert to nb of bytes
*
*     Start allocations
*     *****************
*
C
C  The slaves can now perform the factorization
C
C
C  Allocate S on all nodes
C  or point to user provided data WK_USER when LWK_USER>0
C  =======================
C
C
      PERLU = KEEP(12)
      IF (KEEP(201) .EQ. 0) THEN
C       In-core 
        MAXS_BASE8=KEEP8(12)
       ELSE
C       OOC or no factors stored
        MAXS_BASE8=KEEP8(14)
      ENDIF
      IF (WK_USER_PROVIDED) THEN
C       -- Set MAXS to size of WK_USER_
        MAXS = KEEP8(24)
      ELSE
       IF ( MAXS_BASE8 .GT. 0_8 ) THEN
          MAXS_BASE_RELAXED8 =
     &         MAXS_BASE8 + int(PERLU,8) * ( MAXS_BASE8 / 100_8 + 1_8)
C         If PERLU < 0, we may obtain a
C         null or negative value of MAXS.
          IF (MAXS_BASE_RELAXED8 > huge(MAXS)) THEN
C           INFO(1)=-37
C           INFO(2)=int(MAXS_BASE_RELAXED8/1000000_8)
            WRITE(*,*) "Internal error: I8 overflow"
            CALL MUMPS_ABORT()
          ENDIF
          MAXS_BASE_RELAXED8 = max(MAXS_BASE_RELAXED8, 1_8)
          MAXS = MAXS_BASE_RELAXED8
C         Note that in OOC this value of MAXS will be
C         overwritten if KEEP(96) .NE. 0 or if
C         ICNTL(23) (that is, KEEP8(4)) is provided.
       ELSE
        MAXS = 1_8
        MAXS_BASE_RELAXED8 = 1_8
       END IF
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF (INFO(1) .LT. 0) THEN
        GOTO 530
      ENDIF
C
C     If KEEP(96) is provided,
C     use it without asking questions
C
      IF ((.NOT.WK_USER_PROVIDED).AND.(I_AM_SLAVE)) THEN
C
C
          IF (KEEP(96).GT.0) THEN
C           -- useful mostly for internal testing:
C           -- we can force in this way a given value
C           -- of MAXS and forget about other input values
C           -- such as ICNTL(23) (KEEP8(4)/1E6) 
C           -- that could change MAXS value.
            MAXS=int(KEEP(96),8)
          ELSE
            IF (KEEP8(4) .NE. 0_8) THEN
C             -------------------------
C             WE TRY TO USE MEM_ALLOWED (KEEP8(4)/1E6)
C             -------------------------
C             First compute what we have: TOTAL_MBYTES(PERLU)
C              and TOTAL_BYTES(PERLU)
C
              PERLU_ON = .TRUE.
              CALL SMUMPS_MAX_MEM( id%KEEP(1), id%KEEP8(1),
     &        id%MYID, id%N, id%NELT, id%NA(1), id%LNA,
     &        id%NZ, id%NA_ELT,
     &        id%NSLAVES, TOTAL_MBYTES, .FALSE., KEEP(201), 
     &        PERLU_ON, TOTAL_BYTES)
C
C             Assuming that TOTAL_BYTES is due to MAXS rather than
C             to the temporary buffers used for the distribution of
C             the matrix on the slaves (arrowheads or element distrib),
C             then we have:
C
C             KEEP8(4)-TOTAL_BYTES is the extra free space 
C
C             A simple algorithm to redistribute the extra space:
C             All extra freedom (it could be negative !) is added to MAXS:
              MAXS_BASE_RELAXED8=MAXS_BASE_RELAXED8 +
     &        (KEEP8(4)-TOTAL_BYTES)/int(KEEP(35),8) 
              IF (MAXS_BASE_RELAXED8 > int(huge(MAXS),8)) THEN
                WRITE(*,*) "Internal error: I8 overflow"
                CALL MUMPS_ABORT()
              ELSE IF (MAXS_BASE_RELAXED8 .LE. 0_8) THEN
C               We need more space in order to at least enough
                id%INFO(1)=-9
                IF ( -MAXS_BASE_RELAXED8 .GT.
     &               int(huge(id%INFO(1)),8) ) THEN
                  WRITE(*,*) "I8: OVERFLOW"
                  CALL MUMPS_ABORT()
                ENDIF
                id%INFO(2)=-int(MAXS_BASE_RELAXED8,4)
              ELSE
                MAXS=MAXS_BASE_RELAXED8
              ENDIF
            ENDIF
          ENDIF
      ENDIF ! I_AM_SLAVE
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF (INFO(1) .LT. 0) THEN
        GOTO 530
      ENDIF
      CALL SMUMPS_AVGMAX_STAT8(PROKG, MPG, MAXS, id%NSLAVES,
     & id%COMM, "effective relaxed size of S              =")
C     Next PROPINFO is there for possible negative
C     values of MAXS resulting from small MEM_ALLOWED
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF (id%INFO(1) .LT. 0) THEN
C     We jump after the call to LOAD_END and OOC_END since we didn't
C     called yet OOC_INIT and LOAD_INIT
        GOTO 530
      ENDIF
      IF ( I_AM_SLAVE ) THEN
C       ------------------
C       Dynamic scheduling
C       ------------------
        CALL SMUMPS_LOAD_SET_INICOST( dble(id%COST_SUBTREES),
     &        KEEP(64), KEEP(66),MAXS )
        K28=KEEP(28)
        MEMORY_MD_ARG = min(int(PERLU,8) * ( MAXS_BASE8 / 100_8 + 1_8 ),
C       Restrict freedom from dynamic scheduler when
C       MEM_ALLOWED=ICNTL(23) is small (case where KEEP8(4)-TOTAL_BYTES
C       is negative after call to SMUMPS_MAX_MEM)
     &                      max(0_8, MAXS-MAXS_BASE8))
        CALL SMUMPS_LOAD_INIT( id, MEMORY_MD_ARG, MAXS )
C
C       Out-Of-Core (OOC) issues. Case where we ran one factorization OOC
C       and the second one is in-core: we try to free OOC
C       related data from previous factorization.
C
        CALL SMUMPS_CLEAN_OOC_DATA(id, IERR)
        IF (IERR < 0) THEN
          INFO(1) = -90
          INFO(2) = 0
          GOTO 112
        ENDIF
        IF (KEEP(201) .GT. 0) THEN
C          -------------------
C          OOC initializations
C          -------------------
           IF (KEEP(201).EQ.1 !PANEL Version
     &         .AND.KEEP(50).EQ.0 ! Unsymmetric
     &         .AND.KEEP(251).NE.2 ! Store L to disk
     &         ) THEN
              id%OOC_NB_FILE_TYPE=2 ! declared in MUMPS_OOC_COMMON
           ELSE
              id%OOC_NB_FILE_TYPE=1 ! declared in MUMPS_OOC_COMMON
           ENDIF
C          ------------------------------
C          Dimension IO buffer, KEEP(100)
C          ------------------------------
           IF (KEEP(205) .GT. 0) THEN
             KEEP(100) = KEEP(205)
           ELSE
             IF (KEEP(201).EQ.1) THEN ! PANEL version
               I8TMP = int(id%OOC_NB_FILE_TYPE,8) *
     &               2_8 * int(KEEP(226),8)
             ELSE
               I8TMP = 2_8 * KEEP8(119)
             ENDIF
             I8TMP = I8TMP +  int(max(KEEP(12),0),8) *
     &               (I8TMP/100_8+1_8)
C            we want to avoid too large IO buffers.
C            12M corresponds to 100Mbytes given to buffers.
             I8TMP = min(I8TMP, 12000000_8)
             KEEP(100)=int(I8TMP)
           ENDIF
           IF (KEEP(201).EQ.1) THEN
C            Panel version. Force the use of a buffer. 
             IF ( KEEP(99) < 3 ) THEN
               KEEP(99) = KEEP(99) + 3
             ENDIF
           ENDIF
C          --------------------------
C          Reset KEEP(100) to 0 if no
C          buffer is used for OOC.
C          --------------------------
           IF (KEEP(99) .LT.3) KEEP(100)=0
           IF((dble(KEEP(100))*dble(KEEP(35))/dble(2)).GT.
     &       (dble(1999999999)))THEN
             IF (PROKG) THEN
               WRITE(MPG,*)id%MYID,': Warning: DIM_BUF_IO might be
     &  too big for Filesystem'
             ENDIF
           ENDIF
           ALLOCATE (id%OOC_INODE_SEQUENCE(KEEP(28),
     &          id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_INODE_SEQUENCE)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_TOTAL_NB_NODES(id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = id%OOC_NB_FILE_TYPE
              NULLIFY(id%OOC_TOTAL_NB_NODES)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_SIZE_OF_BLOCK(KEEP(28),
     &          id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_SIZE_OF_BLOCK)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_VADDR(KEEP(28),id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_VADDR)
              GOTO 112
           ENDIF
        ENDIF
      ENDIF
 112  CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF (INFO(1) < 0) THEN
C       LOAD_END must be done but not OOC_END_FACTO
        GOTO 513
      ENDIF
      IF (I_AM_SLAVE) THEN
        IF (KEEP(201) .GT. 0) THEN
           IF ((KEEP(201).EQ.1).OR.(KEEP(201).EQ.2)) THEN
             CALL SMUMPS_OOC_INIT_FACTO(id,MAXS)
           ELSE
             WRITE(*,*) "Internal error in SMUMPS_FAC_DRIVER"
             CALL MUMPS_ABORT()
           ENDIF
           IF(INFO(1).LT.0)THEN
              GOTO 111
           ENDIF
        ENDIF
#if ! defined(OLD_LOAD_MECHANISM)
C       First increment corresponds to the number of
C       floating-point operations for subtrees allocated
C       to the local processor.
        CALL SMUMPS_LOAD_UPDATE(0,.FALSE.,dble(id%COST_SUBTREES),
     &          id%KEEP(1),id%KEEP8(1))
#endif
        IF (INFO(1).LT.0) GOTO 111
      END IF
C     -----------------------
C     Manage main workarray S
C     -----------------------
      IF ( associated (id%S) ) THEN
        DEALLOCATE(id%S)
        NULLIFY(id%S)
        KEEP8(23)=0_8  ! reset space allocated to zero 
      ENDIF
#if defined (LARGEMATRICES)
      IF ( id%MYID .ne. MASTER ) THEN
#endif
      IF (.NOT.WK_USER_PROVIDED) THEN
        ALLOCATE (id%S(MAXS),stat=IERR)
        KEEP8(23) = MAXS
        IF ( IERR .GT. 0 ) THEN
          INFO(1) = -13
          CALL MUMPS_SETI8TOI4(MAXS, INFO(2))
C         On some platforms (IBM for example), an
C         allocation failure returns a non-null pointer.
C         Therefore we nullify S
          NULLIFY(id%S)
          KEEP8(23)=0_8
        ENDIF
      ELSE
       id%S => id%WK_USER(1:KEEP8(24))
      ENDIF
#if defined (LARGEMATRICES)
      END IF
#endif
C
 111  CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 514
C     --------------------------
C     Initialization of modules
C     related to data management
C     --------------------------
      NB_ACTIVE_FRONTS_ESTIM = 3
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FDM_INIT('A',NB_ACTIVE_FRONTS_ESTIM, INFO)
        IF (INFO(1) .LT. 0 ) GOTO 114
#if ! defined(NO_FDM_DESCBAND)
C         Storage of DESCBAND information
          CALL MUMPS_FDBD_INIT( NB_ACTIVE_FRONTS_ESTIM, INFO )
#endif
#if ! defined(NO_FDM_MAPROW)
C         Storage of MAPROW and ROOT2SON information
          CALL MUMPS_FMRD_INIT( NB_ACTIVE_FRONTS_ESTIM, INFO )
#endif
 114    CONTINUE
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                       id%COMM, id%MYID )
C     GOTO 500: one of the above module initializations failed
      IF ( INFO(1).LT.0 ) GOTO 500
C
C
C  Allocate space for matrix in arrowhead 
C  ======================================
C
C  CASE 1 : Matrix is assembled
C  CASE 2 : Matrix is elemental
C
      IF ( KEEP(55) .eq. 0 ) THEN
C       ------------------------------------
C       Space has been allocated already for
C       the integer part during analysis
C       Only slaves need the arrowheads.
C       ------------------------------------
        IF (associated( id%DBLARR)) THEN
          DEALLOCATE(id%DBLARR)
          NULLIFY(id%DBLARR)
        ENDIF
        IF ( I_AM_SLAVE .and. KEEP(13) .ne. 0 ) THEN
          ALLOCATE( id%DBLARR( KEEP(13) ), stat = IERR )
        ELSE
          ALLOCATE( id%DBLARR( 1 ), stat =IERR )
        END IF
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,
     &       ':Error allocating DBLARR : IERR = ', IERR
          INFO(1)=-13
          INFO(2)=KEEP(13)
          NULLIFY(id%DBLARR)
          GOTO 100
        END IF
      ELSE
C        ----------------------------------------
C        Allocate variable lists. Systematically.
C        ----------------------------------------
         IF ( associated( id%INTARR ) ) THEN
           DEALLOCATE( id%INTARR )
           NULLIFY( id%INTARR )
         END IF
         IF ( I_AM_SLAVE .and. KEEP(14) .ne. 0 ) THEN
           ALLOCATE( id%INTARR( KEEP(14) ), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = KEEP(14)
             NULLIFY(id%INTARR)
             GOTO 100
           END IF
         ELSE
           ALLOCATE( id%INTARR(1),stat=allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             NULLIFY(id%INTARR)
             GOTO 100
           END IF
         END IF
C        -----------------------------
C        Allocate real values.
C        On master, if hybrid host and
C        no scaling, avoid the copy.
C        -----------------------------
         IF (associated( id%DBLARR)) THEN
           DEALLOCATE(id%DBLARR)
           NULLIFY(id%DBLARR)
         ENDIF
         IF ( I_AM_SLAVE ) THEN
           IF (      id%MYID_NODES .eq. MASTER
     &       .AND.   KEEP(46)   .eq. 1
     &       .AND.   KEEP(52)   .eq. 0 ) THEN
C            --------------------------
C            Simple pointer association
C            --------------------------
             id%DBLARR => id%A_ELT
           ELSE
C            ----------
C            Allocation
C            ----------
             IF ( KEEP(13) .ne. 0 ) THEN
               ALLOCATE( id%DBLARR( KEEP(13) ), stat = allocok )
               IF ( allocok .GT. 0 ) THEN
                 id%INFO(1) = -13
                 id%INFO(2) = KEEP(13)
                 NULLIFY(id%DBLARR)
                 GOTO 100
               END IF
             ELSE
               ALLOCATE( id%DBLARR(1), stat = allocok )
               IF ( allocok .GT. 0 ) THEN
                 id%INFO(1) = -13
                 id%INFO(2) = 1
                 NULLIFY(id%DBLARR)
                 GOTO 100
               END IF 
             END IF
           END IF
         ELSE
           ALLOCATE( id%DBLARR(1), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             NULLIFY(id%DBLARR)
             GOTO 100
           END IF
         END IF
      END IF
C     -----------------
C     Also prepare some
C     data for the root
C     -----------------
      IF ( KEEP(38).NE.0 .AND.  I_AM_SLAVE ) THEN
         CALL SMUMPS_INIT_ROOT_FAC( id%N,
     &   id%root, id%FILS(1), KEEP(38), id%KEEP(1), id%INFO(1) )
      END IF
C
C
 100  CONTINUE
C     ----------------
C     Check for errors
C     ----------------
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
C
C       -----------------------------------
C
C       DISTRIBUTION OF THE ORIGINAL MATRIX
C
C       -----------------------------------
C
C     TIMINGS: computed (and printed) on the host
C     Next line: global time for distrib(arrowheads,elts)
C     on the host. Synchronization has been performed.
      IF (id%MYID.EQ.MASTER) CALL MUMPS_SECDEB(TIME)
C
      IF ( KEEP( 55 ) .eq. 0 ) THEN
C       ----------------------------
C       Original matrix is assembled
C       Arrowhead format to be used.
C       ----------------------------
C       KEEP(13) and KEEP(14) hold the number of entries for real/integer
C       for the matrix in arrowhead format. They have been set by the
C       analysis phase (SMUMPS_ANA_F and SMUMPS_ANA_G)
C
C       ------------------------------------------------------------------
C       Blocking is used for sending arrowhead records (I,J,VAL)
C              buffer(1) is used to store number of bytes already packed
C              buffer(2) number of records already packed
C       KEEP(39) : Number of records (blocking factor)
C       ------------------------------------------------------------------
C
C     ----------------------------------------
C     In case of parallel root compute minimum
C     size of workspace to receive arrowheads
C     of root node. Will be used to check that
C     MAXS is large enough 
C     ----------------------------------------
      IF (KEEP(38).NE.0 .AND. I_AM_SLAVE) THEN
        LWK = numroc( id%root%ROOT_SIZE, id%root%MBLOCK,
     &             id%root%MYROW, 0, id%root%NPROW ) 
        LWK = max( 1, LWK )
        LWK = LWK*
     &        numroc( id%root%ROOT_SIZE, id%root%NBLOCK,
     &        id%root%MYCOL, 0, id%root%NPCOL )
        LWK = max( 1, LWK )
      ELSE
        LWK = 1
      ENDIF
C     MAXS must be at least 1, and in case of
C     parallel root, large enough to receive
C     arrowheads of root.
      IF (MAXS .LT. int(LWK,8)) THEN
           INFO(1) = -9
           INFO(2) = LWK
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
      IF ( INFO(1).LT.0 ) GOTO 500
      IF ( KEEP(54) .eq. 0 ) THEN
C       ================================================
C       FIRST CASE : MATRIX IS NOT INITIALLY DISTRIBUTED
C       ================================================
C       A small integer workspace is needed to
C       send the arrowheads.
        IF ( id%MYID .eq. MASTER ) THEN
          ALLOCATE(IWK(id%N), stat=allocok)
          IF ( allocok .NE. 0 ) THEN
            INFO(1)=-13
            INFO(2)=id%N
          END IF
#if defined(LARGEMATRICES)
          IF ( associated (id%S) ) THEN
            DEALLOCATE(id%S)
            NULLIFY(id%S)
            KEEP8(23)=0_8
          ENDIF
          ALLOCATE (WK(LWK),stat=IERR)
          IF ( IERR .GT. 0 ) THEN
            INFO(1) = -13
            INFO(2) = LWK
            write(6,*) ' PB1 ALLOC LARGEMAT'
          ENDIF
#endif
        ENDIF
        CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
        IF ( INFO(1).LT.0 ) GOTO 500
        IF ( id%MYID .eq. MASTER ) THEN
C
C         --------------------------------
C         MASTER sends arowheads using the
C         global communicator with ranks
C         also in global communicator
C         IWK is used as temporary
C         workspace of size N.
C         --------------------------------
          IF ( .not. associated( id%INTARR ) ) THEN
            ALLOCATE( id%INTARR( 1 ) )
          ENDIF
#if defined(LARGEMATRICES)
          CALL SMUMPS_FACTO_SEND_ARROWHEADS(id%N, NZ, id%A(1),
     &      id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     &      LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     &      id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     &      min(KEEP(39),id%NZ),
     &      LP, id%COMM, id%root, KEEP,KEEP8,
     &      id%FILS(1), IWK(1), ! workspace of size N
     &
     &      id%INTARR(1), id%DBLARR(1),
     &      id%PTRAR(1), id%PTRAR(id%N+1),
     &      id%FRERE_STEPS(1), id%STEP(1), WK(1), int(LWK,8),
     &      id%ISTEP_TO_INIV2, id%I_AM_CAND,
     &      id%CANDIDATES) 
C
C         write(6,*) '!!! A,IRN,JCN are freed during factorization '
          DEALLOCATE (id%A)
          NULLIFY(id%A)
          DEALLOCATE (id%IRN)
          NULLIFY (id%IRN)
          DEALLOCATE (id%JCN)
          NULLIFY (id%JCN)
          IF (.NOT.WK_USER_PROVIDED) THEN
            ALLOCATE (id%S(MAXS),stat=IERR)
            KEEP8(23) = MAXS
            IF ( IERR .GT. 0 ) THEN
              INFO(1) = -13
              INFO(2) = MAXS
              NULLIFY(id%S)
              KEEP8(23)=0_8
              write(6,*) ' PB2 ALLOC LARGEMAT',MAXS
              CALL MUMPS_ABORT()
            ENDIF
          ELSE
            id%S => id%WK_USER(1:KEEP8(24))
          ENDIF
          id%S(MAXS-LWK+1:MAXS) = WK(1:LWK)
          DEALLOCATE (WK)
#else
          CALL SMUMPS_FACTO_SEND_ARROWHEADS(id%N, NZ, id%A(1),
     &    id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     &    LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     &    id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     &    min(KEEP(39),id%NZ),
     &    LP, id%COMM, id%root, KEEP(1),KEEP8(1),
     &    id%FILS(1), IWK(1),
     &
     &    id%INTARR(1), id%DBLARR(1),
     &    id%PTRAR(1), id%PTRAR(id%N+1),
     &    id%FRERE_STEPS(1), id%STEP(1), id%S(1), MAXS,
     &    id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     &    id%CANDIDATES(1,1) ) 
#endif
          DEALLOCATE(IWK)
        ELSE
          CALL SMUMPS_FACTO_RECV_ARROWHD2( id%N,
     &       id%DBLARR( 1 ), max(1,KEEP( 13 )),
     &       id%INTARR( 1 ), max(1,KEEP( 14 )),
     &       id%PTRAR( 1 ),
     &       id%PTRAR(id%N+1),
     &       KEEP( 1 ), KEEP8(1), id%MYID, id%COMM,
     &       min(id%KEEP(39),id%NZ),
     &
     &       id%S(1), MAXS,
     &       id%root,
     &       id%PROCNODE_STEPS(1), id%NSLAVES,
     &       id%SYM_PERM(1), id%FRERE_STEPS(1), id%STEP(1),
     &       id%INFO(1), id%INFO(2) )
        ENDIF
      ELSE
C
C     =============================================
C     SECOND CASE : MATRIX IS INITIALLY DISTRIBUTED
C     =============================================
C     Timing on master.
      IF (id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECDEB(TIME)
      END IF
      IF ( I_AM_SLAVE ) THEN
C       ---------------------------------------------------
C       In order to have possibly IRN_loc/JCN_loc/A_loc
C       of size 0, avoid to pass them inside REDISTRIBUTION
C       and pass id instead
C       NZ_locMAX gives as a maximum buffer size (send/recv) used
C        an upper bound to limit buffers on small matrices
C       ---------------------------------------------------
       NZ_locMAX = 0
       CALL MPI_ALLREDUCE(id%NZ_loc, NZ_locMAX, 1, MPI_INTEGER, 
     &                   MPI_MAX, id%COMM_NODES, IERR)
        CALL SMUMPS_REDISTRIBUTION( id%N,
     &  id%NZ_loc,
     &  id,
     &  id%DBLARR(1), KEEP(13), id%INTARR(1),
     &  KEEP(14), id%PTRAR(1), id%PTRAR(id%N+1),
     &  KEEP(1), KEEP8(1), id%MYID_NODES,
     &  id%COMM_NODES, min(id%KEEP(39),NZ_locMAX),
     &  id%S(1), MAXS, id%root, id%PROCNODE_STEPS(1),
     &  id%NSLAVES, id%SYM_PERM(1), id%STEP(1),
     &  id%ICNTL(1), id%INFO(1), NSEND, NLOCAL,
     &  id%ISTEP_TO_INIV2(1),
     &  id%CANDIDATES(1,1) )
        IF ( ( KEEP(52).EQ.7 ).OR. (KEEP(52).EQ.8) ) THEN
C         -------------------------------------------------
C         In that case, scaling arrays have been allocated
C         on all processors. They were useful for matrix
C         distribution. But we now really only need them
C         on the host. In case of distributed solution, we
C         will have to broadcast either ROWSCA or COLSCA
C         (depending on MTYPE) but this is done later.
C
C         In other words, on exit from the factorization,
C         we want to have scaling arrays available only
C         on the host.
C         -------------------------------------------------
          IF ( id%MYID > 0 ) THEN
            IF (associated(id%ROWSCA)) THEN
              DEALLOCATE(id%ROWSCA)
              NULLIFY(id%ROWSCA)
            ENDIF
            IF (associated(id%COLSCA)) THEN
              DEALLOCATE(id%COLSCA)
              NULLIFY(id%COLSCA)
            ENDIF
          ENDIF
        ENDIF
#if defined(LARGEMATRICES)
C      deallocate id%IRN_loc, id%JCN(loc) to free extra space
C      Note that in this case IRN_loc cannot be used
C      anymore during the solve phase for IR and Error analysis.
         IF (associated(id%IRN_loc)) THEN
            DEALLOCATE(id%IRN_loc)
            NULLIFY(id%IRN_loc)
         ENDIF
         IF (associated(id%JCN_loc)) THEN
            DEALLOCATE(id%JCN_loc)
            NULLIFY(id%JCN_loc)
         ENDIF
         IF (associated(id%A_loc)) THEN
            DEALLOCATE(id%A_loc)
            NULLIFY(id%A_loc)
         ENDIF
       write(6,*) ' Warning :', 
     &        ' id%A_loc, IRN_loc, JCN_loc deallocated !!! '
#endif
      IF (PROK) THEN
        WRITE(MP,120) NLOCAL, NSEND
      END IF
      END IF
      IF ( KEEP(46) .eq. 0 .AND. id%MYID.eq.MASTER ) THEN
C       ------------------------------
C       The host is not working -> had
C       no data from initial matrix
C       ------------------------------
        NSEND  = 0
        NLOCAL = 0
      END IF
C     --------------------------
C     Put into some info/infog ?
C     --------------------------
      CALL MPI_REDUCE( NSEND, NSEND_TOT, 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( NLOCAL, NLOCAL_TOT, 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      IF ( PROKG ) THEN
        WRITE(MPG,125) NLOCAL_TOT, NSEND_TOT
      END IF
C
C     -------------------------
C     Check for possible errors
C     -------------------------
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
C
      ENDIF
      ELSE
C       -------------------
C       Matrix is elemental,
C       provided on the
C       master only
C       -------------------
        IF ( id%MYID.eq.MASTER)
     &   CALL SMUMPS_MAXELT_SIZE( id%ELTPTR(1),
     &                        id%NELT,
     &                        MAXELT_SIZE )
C
C         Perform the distribution of the elements.
C         A this point,
C           PTRAIW/PTRARW have been computed.
C           INTARR/DBLARR have been allocated
C           ELTPROC gives the mapping of elements
C
        CALL SMUMPS_ELT_DISTRIB( id%N, id%NELT, id%NA_ELT,
     &     id%COMM, id%MYID,
     &     id%NSLAVES, id%PTRAR(1),
     &     id%PTRAR(id%NELT+2),
     &     id%INTARR(1), id%DBLARR(1),
     &     id%KEEP(1), id%KEEP8(1), MAXELT_SIZE,
     &     id%FRTPTR(1), id%FRTELT(1),
     &     id%S(1), MAXS, id%FILS(1),
     &     id, id%root )
C       ----------------
C       Broadcast errors
C       ----------------
        CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
        IF ( INFO( 1 ) .LT. 0 ) GOTO 500
      END IF ! Element entry
C     ------------------------
C     Time the redistribution:
C     ------------------------
      IF ( id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECFIN(TIME)
        id%DKEEP(93) = real(TIME)
        IF (PROKG) WRITE(MPG,160) TIME
      END IF
C
C     TIMINGS:
C     Next line: elapsed time for factorization
      IF (id%MYID.EQ.MASTER)  CALL MUMPS_SECDEB(TIME)
C
C  Allocate buffers on the slaves
C  ==============================
C
      IF ( I_AM_SLAVE )  THEN
        CALL SMUMPS_BUF_INI_MYID(id%MYID_NODES)
C
C  Some buffers are required to pack/unpack data and for
C  receiving MPI messages.
C  For packing/unpacking : the buffer must be large
C  enough to send several messages while receives might not
C  be posted yet.
C  It is assumed that the size of an integer is held in KEEP(34)
C  while the size of a complex is held in KEEP(35).
C  BUFR and LBUFR are declared integers, since byte is not
C  a standard datatype.
C  We now use KEEP(43) and KEEP(44) as estimated at analysis
C  to allocate appropriate buffer sizes.
C
C  Reception buffer
C  ----------------
        SMUMPS_LBUFR_BYTES8 = int(KEEP( 44 ),8) * int(KEEP( 35 ), 8)
C       -------------------
C       Ensure a reasonable
C       buffer size
C       -------------------
        SMUMPS_LBUFR_BYTES8 = max( SMUMPS_LBUFR_BYTES8,
     &                      100000_8 )
C
C  If there is pivoting, size of the message might still increase.
C  We use a relaxation (so called PERLU) to increase the estimate.
C
C  Note: PERLU is a global estimate for pivoting. 
C  It may happen that one large contribution block size is increased by more than that.
C  This is why we use an extra factor 2 relaxation coefficient for the relaxation of
C  the reception buffer in the case where pivoting is allowed.
C  A more dynamic strategy could be applied: if message to
C  be received is larger than expected, reallocate a larger
C  buffer. (But this won't work with IRECV.)
C  Finally, one may want (as we are currently doing it for moste messages)
C  to cut large messages into a series of smaller ones.
C
        PERLU = KEEP( 12 )
C       For hybrid scheduling (strategy 5), Abdou
C       wants a minimal amount of freedom even for
C       small/negative PERLU values.
C
        IF (KEEP(48).EQ.5) THEN
          MIN_PERLU=2
        ELSE
          MIN_PERLU=0
        ENDIF
C
        SMUMPS_LBUFR_BYTES8 = SMUMPS_LBUFR_BYTES8
     &        + int( 2.0E0 * real(max(PERLU,MIN_PERLU))*
     &        real(SMUMPS_LBUFR_BYTES8)/100E0, 8)
        SMUMPS_LBUFR_BYTES8 = min(SMUMPS_LBUFR_BYTES8,
     &                            int(huge (KEEP(43))-100,8))
        SMUMPS_LBUFR_BYTES  = int( SMUMPS_LBUFR_BYTES8 )
        IF (KEEP(48)==5) THEN
C          Since the buffer is allocated, use
C          it as the constraint for memory/granularity
C          in hybrid scheduler
C
           KEEP8(21) = KEEP8(22) + int( real(max(PERLU,MIN_PERLU))*
     &        real(KEEP8(22))/100E0,8)
        ENDIF
C
C  Now estimate the size for the buffer for asynchronous
C  sends of contribution blocks (so called CB). We want to be able to send at
C  least KEEP(213)/100 (two in general) messages at the
C  same time.
C
C   Send buffer
C   -----------
        SMUMPS_LBUF8 = int( real(KEEP(213)) / 100.0E0 *
     &                      real(KEEP(43)) * real(KEEP(35)), 8  )
        SMUMPS_LBUF8 = max( SMUMPS_LBUF8, 100000_8 )
        SMUMPS_LBUF8 = SMUMPS_LBUF8
     &                 + int( 2.0E0 * real(max(PERLU,MIN_PERLU))*
     &                   real(SMUMPS_LBUF8)/100E0, 8)
C       Make SMUMPS_LBUF8 small enough to be stored in a standard integer
        SMUMPS_LBUF8 = min(SMUMPS_LBUF8, int(huge (KEEP(43))-100,8))
C
C       No reason to have send buffer smaller than receive buffer.
C       This should never occur with the formulas above but just
C       in case:
        SMUMPS_LBUF8 = max(SMUMPS_LBUF8, SMUMPS_LBUFR_BYTES8+3*KEEP(34))
        SMUMPS_LBUF  = int(SMUMPS_LBUF8)
        IF(id%KEEP(48).EQ.4)THEN
           SMUMPS_LBUFR_BYTES=SMUMPS_LBUFR_BYTES*5
           SMUMPS_LBUF=SMUMPS_LBUF*5
        ENDIF
C
C  Estimate size of buffer for small messages 
C  Each node can send ( NSLAVES - 1 ) messages to (NSLAVES-1) nodes
C
C  KEEP(56) is the number of nodes of level II.
C  Messages will be sent for the symmetric case
C  for synchronisation issues.
C
C  We take an upperbound
C
        SMUMPS_LBUF_INT = ( KEEP(56) + id%NSLAVES * id%NSLAVES ) * 5
     &               * KEEP(34)
        IF ( KEEP( 38 ) .NE. 0 ) THEN
C
C
          KKKK = MUMPS_PROCNODE( id%PROCNODE_STEPS(id%STEP(KEEP(38))),
     &                           id%NSLAVES )
          IF ( KKKK .EQ. id%MYID_NODES ) THEN
             SMUMPS_LBUF_INT = SMUMPS_LBUF_INT + 
     &     10 *  ! security only
     &      2 * ( id%NE_STEPS(id%STEP(KEEP(38))) + 1 ) * id%NSLAVES
     &                      * KEEP(34)
          END IF
        END IF
        IF ( PROK ) THEN
          WRITE( MP, 9999 ) SMUMPS_LBUFR_BYTES,
     &                      SMUMPS_LBUF, SMUMPS_LBUF_INT
        END IF
 9999   FORMAT( /,' Allocated buffers',/,' ------------------',/,
     &  ' Size of reception buffer in bytes ...... = ', I10,
     &  /,
     &  ' Size of async. emission buffer (bytes).. = ', I10,/,
     &  ' Small emission buffer (bytes) .......... = ', I10)
C       ---------------------------
C       Allocate the 2 send buffers
C       ---------------------------
        CALL SMUMPS_BUF_ALLOC_SMALL_BUF( SMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,
     &   ':Error allocating small Send buffer:IERR='
     &   ,IERR
          INFO(1)= -13
C         convert to size in integer  INFO(2)= SMUMPS_LBUF_INT
          INFO(2)= (SMUMPS_LBUF_INT+KEEP(34)-1)/KEEP(34)
          GO TO 110
        END IF
        CALL SMUMPS_BUF_ALLOC_CB( SMUMPS_LBUF, IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating Send buffer:IERR='
     &   ,IERR
          INFO(1)= -13
C         convert to size in integer          INFO(2)= SMUMPS_LBUF
          INFO(2)= (SMUMPS_LBUF+KEEP(34)-1)/KEEP(34)
          GO TO 110
        END IF
C       -----------------------------
C       Allocate reception buffer and
C       keep it in the structure
C       -----------------------------
        id%LBUFR_BYTES = SMUMPS_LBUFR_BYTES
        id%LBUFR = (SMUMPS_LBUFR_BYTES+KEEP(34)-1)/KEEP(34)
        IF (associated(id%BUFR)) DEALLOCATE(id%BUFR)
        ALLOCATE( id%BUFR( id%LBUFR ),stat=IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating BUFR:IERR='
     &   ,IERR
          INFO(1)=-13
          INFO(2)=id%LBUFR
          NULLIFY(id%BUFR)
          GO TO 110
        END IF
C
C  The buffers are declared INTEGER, because BYTE is not a
C  standard data type. The sizes are in bytes, so we allocate
C  a number of INTEGERs. The allocated size in integer is the
C  size in bytes divided by KEEP(34)
C       -------------------------------
C       Allocate IS. IS will contain
C       factors and contribution blocks
C       -------------------------------
C       Relax workspace at facto now 
C       PERLU might have been modified reload initial value
        PERLU          = KEEP( 12 )
        IF (KEEP(201).GT.0) THEN
C         OOC panel or non panel (note that
C         KEEP(15)=KEEP(225) if non panel)
          MAXIS_ESTIM   = KEEP(225)
        ELSE
C         In-core or reals for factors not stored
          MAXIS_ESTIM   = KEEP(15)
        ENDIF
        MAXIS = max( 1,
     &       MAXIS_ESTIM + 2 * max(PERLU,10) * 
     &          ( MAXIS_ESTIM / 100 + 1 )
     &  )
        IF (associated(id%IS)) DEALLOCATE( id%IS )
        ALLOCATE( id%IS( MAXIS  ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
         WRITE(*,*) id%MYID,':Error allocating IS:IERR=',IERR
         INFO(1)=-13
         INFO(2)=MAXIS
         NULLIFY(id%IS)
         GO TO 110
        END IF
C       ----------------------
C       Set up subdivision of array IS. 
C       This is used by the slaves only.
C       ----------------------
        LIW = MAXIS
C       -----------------------
C       Allocate PTLUST_S. PTLUST_S
C       is used by solve later
C       -----------------------
        IF (associated( id%PTLUST_S )) DEALLOCATE(id%PTLUST_S)
        ALLOCATE( id%PTLUST_S( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocatingPTLUST:IERR = ',
     &    IERR
          INFO(1)=-13
          INFO(2)=id%KEEP(28)
          NULLIFY(id%PTLUST_S)
          GOTO 100
        END IF
        IF (associated( id%PTRFAC )) DEALLOCATE(id%PTRFAC)
        ALLOCATE( id%PTRFAC( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocatingPTRFAC:IERR = ',
     &    IERR
          INFO(1)=-13
          INFO(2)=id%KEEP(28)
          NULLIFY(id%PTRFAC)
          GOTO 100
        END IF
C       -----------------------------
C       Allocate temporary workspace :
C       IPOOL, PTRWB, ITLOC, PTRIST
C       PTRWB will be subdivided again
C       in routine SMUMPS_FAC_B
C       -----------------------------
        PTRIST = 1
        PTRWB  = PTRIST + id%KEEP(28)
        ITLOC  = PTRWB  + 3 * id%KEEP(28)
C Fwd in facto: ITLOC of size id%N + id%KEEP(253)
        IPOOL  = ITLOC  + id%N + id%KEEP(253)
C
C       --------------------------------
C       NA(1) is an upperbound for LPOOL
C       --------------------------------
C       Structure of the pool:
C     ____________________________________________________
C    | Subtrees   |         | Top nodes           | 1 2 3 |
C     ----------------------------------------------------
        LPOOL = MUMPS_GET_POOL_LENGTH(id%NA(1), id%KEEP(1),id%KEEP8(1))
        ALLOCATE( IWK(  IPOOL + LPOOL - 1 ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating IWK : IERR = ',
     &    IERR
          INFO(1)=-13
          INFO(2)=IPOOL + LPOOL - 1
          GOTO 110
        END IF
        ALLOCATE(IWK8( 2 * id%KEEP(28)), stat = IERR)
        IF ( IERR .NE. 0 ) THEN
          WRITE(*,*) id%MYID,':Error allocating IWK : IERR = ',
     &    IERR
          INFO(1)=-13
          INFO(2)=2 * id%KEEP(28)
          GOTO 110
        END IF
C
C  Return to SPMD
C
      ENDIF
C
 110  CONTINUE
C     ----------------
C     Broadcast errors
C     ----------------
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
      IF ( INFO( 1 ) .LT. 0 ) GOTO 500
C
      IF ( I_AM_SLAVE )  THEN
C
C       Store size of receive buffers in module
        CALL SMUMPS_BUF_DIST_IRECV_SIZE( id%LBUFR_BYTES )
        IF (PROK) THEN
          WRITE( MP, 170 ) MAXS, MAXIS, KEEP8(12), KEEP(15), KEEP(13),
     &    KEEP(14), KEEP8(11), KEEP(26), KEEP(27)
        ENDIF
      END IF
C
C  SPMD
C
      PERLU_ON = .TRUE.
      CALL SMUMPS_MAX_MEM( id%KEEP(1), id%KEEP8(1),
     &     id%MYID, id%N, id%NELT, id%NA(1), id%LNA, id%NZ,
     &     id%NA_ELT,
     &     id%NSLAVES, TOTAL_MBYTES, .FALSE., id%KEEP(201),
     &     PERLU_ON, TOTAL_BYTES)
      id%INFO(16) = TOTAL_MBYTES
      IF ( PROK ) THEN
          WRITE(MP,'(A,I10) ')
     &    ' ** Space in MBYTES used during factorization  :',
     &                id%INFO(16)
      END IF
C
C     ----------------------------------------------------
C     Centralize memory statistics on the host
C       INFOG(18) = size of mem in bytes for facto,
C                   for the processor using largest memory
C       INFOG(19) = size of mem in bytes for facto,
C                   sum over all processors
C     ----------------------------------------------------
C
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &                           id%INFO(16), id%INFOG(18), IRANK )
      IF ( PROKG ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** Memory relaxation parameter ( ICNTL(14)  )            :',
     &  KEEP(12)
        WRITE( MPG,'(A,I10) ')
     &  ' ** Rank of processor needing largest memory in facto     :',
     &  IRANK
        WRITE( MPG,'(A,I10) ')
     &  ' ** Space in MBYTES used by this processor for facto      :',
     &  id%INFOG(18)
        IF ( KEEP(46) .eq. 0 ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** Avg. Space in MBYTES per working proc during facto    :',
     &  ( id%INFOG(19)-id%INFO(16) ) / id%NSLAVES
        ELSE
        WRITE( MPG,'(A,I10) ')
     &  ' ** Avg. Space in MBYTES per working proc during facto    :',
     &  id%INFOG(19) / id%NSLAVES
        END IF
      END IF
C     --------------------------------------------
C     Before calling the main driver, SMUMPS_FAC_B,
C     some statistics should be initialized to 0,
C     even on the host node because they will be
C     used in REDUCE operations afterwards.
C     --------------------------------------------
C     Size of factors written. It will be set to POSFAC in
C     IC, otherwise we accumulate written factors in it.
      KEEP8(31)= 0_8
C     Number of entries in factors
      KEEP8(10) = 0_8
C     KEEP8(8) will hold the volume of extra copies due to
C              in-place stacking in fac_mem_stack.F
      KEEP8(8)=0_8
      INFO(9:14)=0
      RINFO(2:3)=ZERO
      IF ( I_AM_SLAVE ) THEN
C       ------------------------------------
C       Call effective factorization routine
C       ------------------------------------
        IF ( KEEP(55) .eq. 0 ) THEN
          LDPTRAR = id%N
        ELSE
          LDPTRAR = id%NELT + 1
        END IF
        IF ( id%KEEP(55) .NE. 0 ) THEN
          NELT = id%NELT
        ELSE
C         ------------------------------
C         Use size 1 to avoid complaints
C         when using check bound options
C         ------------------------------
          NELT = 1
        END IF
        CALL SMUMPS_FAC_B( id%N, NSTEPS,id%S(1),MAXS,id%IS(1),LIW,
     &      id%SYM_PERM(1),id%NA(1),id%LNA,id%NE_STEPS(1),
     &      id%ND_STEPS(1),id%FILS(1),id%STEP(1),id%FRERE_STEPS(1),
     &      id%DAD_STEPS(1),id%CANDIDATES(1,1),id%ISTEP_TO_INIV2(1),
     &      id%TAB_POS_IN_PERE(1,1),id%PTRAR(1),LDPTRAR,IWK(PTRIST),
     &      id%PTLUST_S(1), id%PTRFAC(1), IWK(PTRWB), IWK8, IWK(ITLOC),
     &      RHS_MUMPS(1), IWK(IPOOL), LPOOL, CNTL1, ICNTL(1), INFO(1),
     &      RINFO(1),KEEP(1),KEEP8(1), id%PROCNODE_STEPS(1), id%NSLAVES,
     &      id%COMM_NODES, id%MYID, id%MYID_NODES, id%BUFR(1),id%LBUFR,
     &      id%LBUFR_BYTES, id%INTARR(1),id%DBLARR(1), id%root, NELT,
     &      id%FRTPTR(1), id%FRTELT(1),id%COMM_LOAD, id%ASS_IRECV,SEUIL,
     &      SEUIL_LDLT_NIV2, id%MEM_DIST(0), id%DKEEP(1),
     &      id%PIVNUL_LIST(1),LPN_LIST
     &       )
        IF ( PROK .and. KEEP(38) .ne. 0 ) THEN
          WRITE( MP, 175 ) KEEP(49)
        END IF
C
C       ------------------------------
C       Deallocate temporary workspace
C       ------------------------------
        DEALLOCATE( IWK  )
        DEALLOCATE( IWK8 )
      ENDIF
C     ---------------------------------
C     Free some workspace corresponding
C     to the original matrix in
C     arrowhead or elemental format.
C                  -----
C     Note : INTARR was not allocated
C     during factorization in the case
C     of an assembled matrix.
C     ---------------------------------
        IF ( KEEP(55) .eq. 0 ) THEN
C
C         ----------------
C         Assembled matrix
C         ----------------
          IF (associated( id%DBLARR)) THEN
            DEALLOCATE(id%DBLARR)
            NULLIFY(id%DBLARR)
          ENDIF
C
        ELSE
C
C         ----------------
C         Elemental matrix
C         ----------------
          DEALLOCATE( id%INTARR)
          NULLIFY( id%INTARR )
C         ------------------------------------
C         For the master from an hybrid host
C         execution without scaling, then real
C         values have not been copied !
C         -------------------------------------
          IF (      id%MYID_NODES .eq. MASTER
     &      .AND.   KEEP(46)   .eq. 1
     &      .AND.   KEEP(52)   .eq. 0 ) THEN
            NULLIFY( id%DBLARR )
          ELSE
C           next line should be enough but ... 
C           DEALLOCATE( id%DBLARR ) 
            IF (associated( id%DBLARR)) THEN
              DEALLOCATE(id%DBLARR)
              NULLIFY(id%DBLARR)
            ENDIF
          END IF
        END IF
C     Memroy statistics
C     -----------------------------------
C     If QR (Keep(19)) is not zero, and if
C     the host does not have the information
C     (ie is not slave), send information
C     computed on the slaves during facto
C     to the host.
C     -----------------------------------
      IF ( KEEP(19) .NE. 0 ) THEN
        IF ( KEEP(46) .NE. 1 ) THEN
C         Host was not working during facto_root
C         Send him the information
          IF ( id%MYID .eq. MASTER ) THEN
            CALL MPI_RECV( KEEP(17), 1, MPI_INTEGER, 1, DEFIC_TAG,
     &                   id%COMM, STATUS, IERR )
          ELSE IF ( id%MYID .EQ. 1 ) THEN
            CALL MPI_SEND( KEEP(17), 1, MPI_INTEGER, 0, DEFIC_TAG,
     &                   id%COMM, IERR )
          END IF
        END IF
      END IF
C     ---------------------------
C     Deallocate send buffers
C     They will be reallocated
C     in the solve.
C     ------------------------
      IF (associated(id%BUFR)) THEN
        DEALLOCATE(id%BUFR)
        NULLIFY(id%BUFR)
      END IF
      CALL SMUMPS_BUF_DEALL_CB( IERR )
      CALL SMUMPS_BUF_DEALL_SMALL_BUF( IERR )
C//PIV
      IF (KEEP(219).NE.0) THEN
      CALL SMUMPS_BUF_DEALL_MAX_ARRAY()
      ENDIF
C
C     Check for errors, 
C     every slave is aware of an error. 
C     If master is included in computations, the call below should
C     not be necessary.
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &                    id%COMM, id%MYID )
C
      CALL SMUMPS_EXTRACT_SCHUR_REDRHS(id)
      IF (KEEP(201) .GT. 0) THEN
         IF ((KEEP(201).EQ.1) .OR. (KEEP(201).EQ.2)) THEN
            IF ( I_AM_SLAVE ) THEN
               CALL SMUMPS_OOC_CLEAN_PENDING(IERR)
               IF(IERR.LT.0)THEN
                  INFO(1)=IERR
                  INFO(2)=0
               ENDIF
            ENDIF
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
C             We want to collect statistics even in case of 
C             error to understand if it is due to numerical
C             issues
CC            IF ( id%INFO(1) < 0 ) GOTO 500
         END IF
      END IF
      IF (id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECFIN(TIME)
        id%DKEEP(94)=real(TIME)
        IF ( PROKG ) THEN
          IF (id%INFOG(1) .GE.0) THEN
            WRITE(MPG,180) TIME
          ELSE
            WRITE(MPG,185) TIME
          ENDIF
        ENDIF
      ENDIF
CC Made available to users on release 4.4 (April 2005)
      PERLU_ON = .TRUE.
      CALL SMUMPS_MAX_MEM( id%KEEP(1),id%KEEP8(1),
     &     id%MYID, N, id%NELT, id%NA(1), id%LNA, id%NZ,
     &     id%NA_ELT,
     &     id%NSLAVES, TOTAL_MBYTES, .TRUE., id%KEEP(201),
     &     PERLU_ON, TOTAL_BYTES)
      KEEP8(7) = TOTAL_BYTES
C     -- INFO(22) holds the effective space (in Mbytes) used by MUMPS
C     -- (it includes part of WK_USER used if provided by user)
      id%INFO(22) = TOTAL_MBYTES
      IF (PROK ) THEN
          WRITE(MP,'(A,I10) ')
     &    ' ** Effective minimum Space in MBYTES for facto  :',
     &                TOTAL_MBYTES
      ENDIF
C
      IF (I_AM_SLAVE) THEN
       K67 = KEEP8(67)
      ELSE
       K67 = 0_8
      ENDIF
C     -- Save the number of entries effectively used
C        in main working array S
      CALL MUMPS_SETI8TOI4(K67,id%INFO(21))
C
      CALL SMUMPS_AVGMAX_STAT8(PROKG, MPG, K67, id%NSLAVES,
     & id%COMM, "effective space used in S   (KEEP8(67))  =")
C
C     ----------------------------------------------------
C     Centralize memory statistics on the host
C
C       INFOG(21) = size of mem (Mbytes) for facto,
C                   for the processor using largest memory
C       INFOG(22) = size of mem (Mbytes) for facto,
C                   sum over all processors
C     ----------------------------------------------------
C
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &                    TOTAL_MBYTES, id%INFOG(21), IRANK )
      IF ( PROKG ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** EFF Min: Rank of processor needing largest memory :',
     &  IRANK
        WRITE( MPG,'(A,I10) ')
     &  ' ** EFF Min: Space in MBYTES used by this processor   :',
     &  id%INFOG(21)
        IF ( KEEP(46) .eq. 0 ) THEN
        WRITE( MPG,'(A,I10) ')
     &  ' ** EFF Min: Avg. Space in MBYTES per working proc    :',
     &  ( id%INFOG(22)-TOTAL_MBYTES ) / id%NSLAVES
        ELSE
        WRITE( MPG,'(A,I10) ')
     &  ' ** EFF Min: Avg. Space in MBYTES per working proc    :',
     &  id%INFOG(22) / id%NSLAVES
        END IF
      END IF
*     save statistics in KEEP array.
      KEEP(33) = INFO(11) ! this should be the other way round
C
C  Sum RINFO(2) : total number of flops for assemblies
C  Sum RINFO(3) : total number of flops for eliminations
C 
C  Should work even if the master does some work
C
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     &                 MPI_REAL,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
C     Reduce needed to dimension small working array
C     on all procs during SMUMPS_GATHER_SOLUTION
      KEEP(247) = 0
      CALL MPI_REDUCE( KEEP(246), KEEP(247), 1, MPI_INTEGER, 
     &                 MPI_MAX, MASTER, id%COMM, IERR)
C
C     Reduce compression times: get max compression times
      CALL MPI_REDUCE( id%DKEEP(97), id%DKEEP(98), 1,
     &     MPI_REAL,
     &     MPI_MAX, MASTER, id%COMM, IERR)
C
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     &                 MPI_REAL,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MUMPS_REDUCEI8( KEEP8(31),KEEP8(6), MPI_SUM,
     &                     MASTER, id%COMM )
      CALL MUMPS_SETI8TOI4(KEEP8(6), INFOG(9))
      CALL MPI_REDUCE( INFO(10), INFOG(10), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
C     Use MPI_MAX for this one to get largest front size
      CALL MPI_ALLREDUCE( INFO(11), INFOG(11), 1, MPI_INTEGER,
     &                 MPI_MAX, id%COMM, IERR)
C     make maximum effective frontal size available on all procs
C     for solve phase
C     (Note that INFO(11) includes root size on root master)
      KEEP(133) = INFOG(11)
      CALL MPI_REDUCE( INFO(12), INFOG(12), 3, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( KEEP(103), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(229) = INFOG(25)
      CALL MPI_REDUCE( KEEP(105), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(230) = INFOG(25)
C
      INFO(25) = KEEP(98)
      CALL MPI_ALLREDUCE( INFO(25), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, id%COMM, IERR)
C     Extra copies due to in-place stacking
      CALL MUMPS_REDUCEI8( KEEP8(8), KEEP8(108), MPI_SUM,
     &                     MASTER, id%COMM )
C     Entries in factors
      CALL MUMPS_SETI8TOI4(KEEP8(10), INFO(27))
      CALL MUMPS_REDUCEI8( KEEP8(10),KEEP8(110), MPI_SUM,
     &                     MASTER, id%COMM )
      CALL MUMPS_SETI8TOI4(KEEP8(110), INFOG(29))
C     ==============================
C     NULL PIVOTS AND RANK-REVEALING
C     ==============================
      IF(KEEP(110) .EQ. 1) THEN
C        -- make available to users the local number of null pivots detected 
C        -- with ICNTL(24) = 1.
         INFO(18) = KEEP(109)
         CALL MPI_ALLREDUCE( KEEP(109), KEEP(112), 1, MPI_INTEGER,
     &        MPI_SUM, id%COMM, IERR)
      ELSE
         INFO(18)  = 0
         KEEP(109) = 0
         KEEP(112) = 0
      ENDIF
C     INFOG(28) deficiency resulting from ICNTL(24) and ICNTL(16).
C     Note that KEEP(17) already has the same value on all procs
      INFOG(28)=KEEP(112)+KEEP(17)
C     ========================================
C     We now provide to the host the part of
C     PIVNUL_LIST resulting from the processing
C     of the root node and we update INFO(18)
C     on the processor holding the root to
C     include null pivots relative to the root
C     ========================================
      IF (KEEP(17) .NE. 0) THEN
        IF (id%MYID .EQ. ID_ROOT) THEN
C         Include in INFO(18) null pivots resulting
C         from deficiency on the root. In this way,
C         the sum of all INFO(18) is equal to INFOG(28).
          INFO(18)=INFO(18)+KEEP(17)
        ENDIF
        IF (ID_ROOT .EQ. MASTER) THEN
          IF (id%MYID.EQ.MASTER) THEN
C           --------------------------------------------------
C           Null pivots of root have been stored in
C           PIVNUL_LIST(KEEP(109)+1:KEEP(109)+KEEP(17).
C           Shift them at the end of the list because:
C           * this is what we need to build the null space
C           * we would otherwise overwrite them on the host
C             when gathering null pivots from other processors
C           --------------------------------------------------
            DO I=1, KEEP(17)
              id%PIVNUL_LIST(KEEP(112)+I)=id%PIVNUL_LIST(KEEP(109)+I)
            ENDDO
          ENDIF
        ELSE
C         ---------------------------------
C         Null pivots of root must be sent
C         from the processor responsible of
C         the root to the host (or MASTER).
C         ---------------------------------
          IF (id%MYID .EQ. ID_ROOT) THEN
            CALL MPI_SEND(id%PIVNUL_LIST(KEEP(109)+1), KEEP(17),
     &                    MPI_INTEGER, MASTER, ZERO_PIV,
     &                    id%COMM, IERR)
          ELSE IF (id%MYID .EQ. MASTER) THEN
            CALL MPI_RECV(id%PIVNUL_LIST(KEEP(112)+1), KEEP(17),
     &                    MPI_INTEGER, ID_ROOT, ZERO_PIV,
     &                    id%COMM, STATUS, IERR )
          ENDIF
        ENDIF
      ENDIF
C     ===========================
C     gather zero pivots indices
C     on the host node
C     ===========================
C     In case of non working host, the following code also
C     works considering that KEEP(109) is equal to 0 on
C     the non-working host
      IF(KEEP(110) .EQ. 1) THEN
         ALLOCATE(ITMP2(id%NPROCS),stat = IERR )  ! deallocated in 490
         IF ( IERR .GT. 0 ) THEN
            INFO(1)=-13
            INFO(2)=id%NPROCS
         END IF
         CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
         IF (INFO(1).LT.0) GOTO 490
         CALL MPI_GATHER ( KEEP(109),1, MPI_INTEGER, 
     &        ITMP2(1), 1, MPI_INTEGER, 
     &        MASTER, id%COMM, IERR)
         IF(id%MYID .EQ. MASTER) THEN
            POSBUF = ITMP2(1)+1
C           First null pivot of master is in
C           position 1 of global list
            KEEP(220)=1
            DO I = 1,id%NPROCS-1
               CALL MPI_RECV(id%PIVNUL_LIST(POSBUF), ITMP2(I+1), 
     &              MPI_INTEGER,I, 
     &              ZERO_PIV, id%COMM, STATUS, IERR)
C              Send position POSBUF of first null pivot of proc I
C              in global list. Will allow to quickly identify during
C              the solve step if one is concerned by a global position
C              K, 0 <= K <= INFOG(28).
               CALL MPI_SEND(POSBUF, 1, MPI_INTEGER, I, ZERO_PIV,
     &              id%COMM, IERR)
               POSBUF = POSBUF + ITMP2(I+1)
            ENDDO
         ELSE
            CALL MPI_SEND( id%PIVNUL_LIST(1), KEEP(109), MPI_INTEGER,
     &           MASTER,ZERO_PIV, id%COMM, IERR)
            CALL MPI_RECV( KEEP(220), 1, MPI_INTEGER, MASTER, ZERO_PIV,
     &           id%COMM, STATUS, IERR )
         ENDIF
      ENDIF
C     =====================================
C     Statistics concerning the determinant
C     =====================================
C
C     1/ on the host better take into account null pivots if scaling:
C
C     Since null pivots are excluded from the computation
C     of the determinant, we also exclude the corresponding
C     scaling entries. Since those entries have already been
C     taken into account before the factorization, we multiply
C     the determinant on the host by the scaling values corresponding
C     to pivots in PIVNUL_LIST.
      IF (id%MYID.EQ.MASTER .AND. LSCAL. AND. KEEP(258).NE.0) THEN
        DO I = 1, id%INFOG(28)
          CALL SMUMPS_UPDATEDETER(id%ROWSCA(id%PIVNUL_LIST(I)),
     &                            id%DKEEP(6), KEEP(259))
          CALL SMUMPS_UPDATEDETER(id%COLSCA(id%PIVNUL_LIST(I)),
     &                            id%DKEEP(6), KEEP(259))
        ENDDO
      ENDIF
C
C     2/ Swap signs depending on pivoting on each proc
C
      IF (KEEP(258).NE.0) THEN
C       Return the determinant in INFOG(34) and RINFOG(12/13)
C       In case of real arithmetic, initialize
C       RINFOG(13) to 0 (no imaginary part and
C       not touched by SMUMPS_DETER_REDUCTION)
        RINFOG(13)=0.0E0
        IF (KEEP(260).EQ.-1) THEN ! Local to each processor
          id%DKEEP(6)=-id%DKEEP(6)
        ENDIF
C
C       3/ Perform a reduction
C
        CALL SMUMPS_DETER_REDUCTION(
     &           id%COMM, id%DKEEP(6), KEEP(259),
     &           RINFOG(12), INFOG(34), id%NPROCS)
C
C       4/ Swap sign if needed
C
        IF (id%KEEP(50).EQ.0 .AND. id%MYID.EQ. MASTER) THEN
C         Modify sign of determinant according
C         to unsymmetric permutation (max-trans
C         of max-weighted matching)
          IF (id%KEEP(23).NE.0) THEN
            CALL SMUMPS_DETER_SIGN_PERM(
     &           RINFOG(12), id%N,
C           id%STEP: used as workspace of size N still
C                    allocated on master; restored on exit
     &           id%STEP(1),
     &           id%UNS_PERM(1) )
C           Remark that RINFOG(12/13) are modified only
C           on the host but will be broadcast on exit
C           from MUMPS (see SMUMPS_DRIVER)
          ENDIF
        ENDIF
      ENDIF
 490  IF (allocated(ITMP2)) DEALLOCATE(ITMP2)
      IF ( PROKG ) THEN
C     ---------------------------------------
C     PRINT STATISTICS  (on master)
C     -----------------------------
          WRITE(MPG,99984) RINFOG(2),RINFOG(3),KEEP8(6),INFOG(10),
     &                    INFOG(11), KEEP8(110)
          IF (id%KEEP(50) == 1 .OR. id%KEEP(50) == 2) THEN
            ! negative pivots
            WRITE(MPG, 99987) INFOG(12)
          END IF
          IF (id%KEEP(50) == 0) THEN
            ! off diag pivots
            WRITE(MPG, 99985) INFOG(12)
          END IF
          IF (id%KEEP(50) .NE. 1) THEN
            ! delayed pivots
            WRITE(MPG, 99982) INFOG(13)
          END IF
          IF (KEEP(97) .NE. 0) THEN
            ! tiny pivots
            WRITE(MPG, 99986) INFOG(25)
          ENDIF
          IF (id%KEEP(50) == 2) THEN
            !number of 2x2 pivots in type 1 nodes
             WRITE(MPG, 99988) KEEP(229)
            !number of 2x2 pivots in type 2 nodes
             WRITE(MPG, 99989) KEEP(230)
          ENDIF
          !number of zero pivots
          IF (KEEP(110) .NE.0) THEN
              WRITE(MPG, 99991) KEEP(112)
          ENDIF
          !Deficiency on root
          IF ( KEEP(17) .ne. 0 )
     &    WRITE(MPG, 99983) KEEP(17)
          !Total deficiency
          IF (KEEP(110).NE.0.OR.KEEP(17).NE.0)
     &    WRITE(MPG, 99992) KEEP(17)+KEEP(112)
          ! Memory compress
          WRITE(MPG, 99981) INFOG(14)
          ! Extra copies due to ip stack in unsym case
          ! in core case (or OLD_OOC_PANEL)
          IF (KEEP8(108) .GT. 0_8) THEN
            WRITE(MPG, 99980) KEEP8(108)
          ENDIF
          IF  ((KEEP(60).NE.0) .AND. INFOG(25).GT.0) THEN
          !  Schur on and tiny pivots set in last level 
          ! before the SChur
           WRITE(MPG, '(A)') 
     & " ** Warning Static pivoting was necessary"
           WRITE(MPG, '(A)') 
     & " ** to factor interior variables with Schur ON"
          ENDIF
          IF (KEEP(258).NE.0) THEN
            WRITE(MPG,99978) RINFOG(12)
            WRITE(MPG,99977) INFOG(34)
          ENDIF
      END IF
* ==========================================
*
*  End of Factorization Phase
*
* ==========================================
C
C  Goto 500 is done when
C  LOAD_INIT
C  OOC_INIT_FACTO
C  MUMPS_FDM_INIT
#if ! defined(NO_FDM_DESCBAND)
C  MUMPS_FDBD_INIT
#endif
#if ! defined(NO_FDM_MAPROW)
C  MUMPS_FMRD_INIT
#endif
C  are all called.
C
 500  CONTINUE
#if ! defined(NO_FDM_DESCBAND)
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FDBD_END(INFO(1))  ! INFO(1): input only
      ENDIF
#endif
#if ! defined(NO_FDM_MAPROW)
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FMRD_END(INFO(1))  ! INFO(1): input only
      ENDIF
#endif
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FDM_END('A')
      ENDIF
C
C  Goto 514 is done when an
C  error occurred in MUMPS_FDM_INIT
C  or (after FDM_INIT but before
C  OOC_INIT)
C
 514  CONTINUE
      IF ( I_AM_SLAVE ) THEN
         IF ((KEEP(201).EQ.1).OR.(KEEP(201).EQ.2)) THEN
            CALL SMUMPS_OOC_END_FACTO(id,IERR)
            IF (IERR.LT.0 .AND. INFO(1) .GE. 0) INFO(1) = IERR
         ENDIF
         IF (WK_USER_PROVIDED) THEN
C     at the end of a phase S is always freed when WK_USER provided
            NULLIFY(id%S)
         ELSE IF (KEEP(201).NE.0) THEN
C           ----------------------------------------
C           In OOC or if KEEP(201).EQ.-1 we always
C           free S at end of factorization. As id%S
C           may be unassociated in case of error
C           during or before the allocation of id%S,
C           we only free S when it was associated.
C           ----------------------------------------
            IF (associated(id%S))  DEALLOCATE(id%S)
            NULLIFY(id%S)   ! in all cases
            KEEP8(23)=0_8
         ENDIF
      ELSE  ! host not working
         IF (WK_USER_PROVIDED) THEN
C     at the end of a phase S is always freed when WK_USER provided
            NULLIFY(id%S)
         ELSE
            IF (associated(id%S))  DEALLOCATE(id%S)
            NULLIFY(id%S)   ! in all cases
            KEEP8(23)=0_8
         END IF
      END IF
C
C     Goto 513 is done in case of error where LOAD_INIT was
C     called but not OOC_INIT_FACTO.
 513  CONTINUE
      IF ( I_AM_SLAVE ) THEN
         CALL SMUMPS_LOAD_END( INFO(1), IERR )
         IF (IERR.LT.0 .AND. INFO(1) .GE. 0) INFO(1) = IERR
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), INFO(1),
     &     id%COMM, id%MYID )
C
C     Goto 530 is done when an error occurs before
C     the calls to LOAD_INIT and OOC_INIT_FACTO
 530  CONTINUE
C  Fwd in facto: free RHS_MUMPS in case
C  it was allocated.
      IF (RHS_MUMPS_ALLOCATED) DEALLOCATE(RHS_MUMPS)
      NULLIFY(RHS_MUMPS)
C
      id%KEEP(13) = KEEP13_SAVE
      RETURN
 120  FORMAT(/' LOCAL REDISTRIB: DATA LOCAL/SENT         =',I14,I14)
 125  FORMAT(/' REDISTRIB: TOTAL DATA LOCAL/SENT         =',I14,I14)
 130  FORMAT(/' ****** FACTORIZATION STEP ********'/)
 160  FORMAT(' GLOBAL TIME FOR MATRIX DISTRIBUTION       =',F12.4)
 166  FORMAT(' Convergence error after scaling for ONE-NORM',
     &       ' (option 7/8)   =',D9.2)
 170  FORMAT(/' STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     &        ' Size of internal working array S         =',I14/
     &        ' Size of internal working array IS        =',I14/
     &        ' MINIMUM (ICNTL(14)=0) size of S          =',I14/
     &        ' MINIMUM (ICNTL(14)=0) size of IS         =',I14/
     &        ' REAL SPACE FOR ORIGINAL MATRIX           =',I14/
     &        ' INTEGER SPACE FOR ORIGINAL MATRIX        =',I14/
     &        ' REAL SPACE FOR FACTORS                   =',I14/
     &        ' INTEGER SPACE FOR FACTORS                =',I14/
     &        ' MAXIMUM FRONTAL SIZE (ESTIMATED)         =',I14)
 172  FORMAT(/' GLOBAL STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     &        ' NUMBER OF WORKING PROCESSES              =',I14/
     &        ' OUT-OF-CORE OPTION (ICNTL(22))           =',I14/
     &        ' REAL SPACE FOR FACTORS                   =',I14/
     &        ' INTEGER SPACE FOR FACTORS                =',I14/
     &        ' MAXIMUM FRONTAL SIZE (ESTIMATED)         =',I14/
     &        ' NUMBER OF NODES IN THE TREE              =',I14/
     &        ' MEMORY ALLOWED (MB -- 0: N/A )           =',I14)
 173  FORMAT( ' PERFORM FORWARD DURING FACTO, NRHS       =',I14)
 175  FORMAT(/' NUMBER OF ENTRIES FOR // ROOT            =',I14)
 180  FORMAT(/' ELAPSED TIME FOR FACTORIZATION           =',F12.4)
 185  FORMAT(/' ELAPSED TIME FOR (FAILED) FACTORIZATION  =',F12.4)
99977 FORMAT( ' INFOG(34)  DETERMINANT (base 2 exponent) =',I14)
99978 FORMAT( ' RINFOG(12) DETERMINANT (real part)       =',F12.4)
99980 FORMAT( ' Extra copies due to In-Place stacking    =',I14)
99981 FORMAT( ' INFOG(14)  NUMBER OF MEMORY COMPRESS     =',I14)
99982 FORMAT( ' INFOG(13)  NUMBER OF DELAYED PIVOTS      =',I14)
99983 FORMAT( ' NB OF NULL PIVOTS DETECTED BY ICNTL(16)  =',I14)
99991 FORMAT( ' NB OF NULL PIVOTS DETECTED BY ICNTL(24)  =',I14)
99992 FORMAT( ' INFOG(28)  ESTIMATED DEFICIENCY          =',I14)
99984 FORMAT(/' GLOBAL STATISTICS '/
     &        ' RINFOG(2)  OPERATIONS IN NODE ASSEMBLY   =',1PD10.3/
     &        ' ------(3)  OPERATIONS IN NODE ELIMINATION=',1PD10.3/
     &        ' INFOG (9)  REAL SPACE FOR FACTORS        =',I14/
     &        ' INFOG(10)  INTEGER SPACE FOR FACTORS     =',I14/
     &        ' INFOG(11)  MAXIMUM FRONT SIZE            =',I14/
     &        ' INFOG(29)  NUMBER OF ENTRIES IN FACTORS  =',I14)
99985 FORMAT( ' INFOG(12)  NUMBER OF OFF DIAGONAL PIVOTS =',I14)
99986 FORMAT( ' INFOG(25)  NUMBER OF TINY PIVOTS(STATIC) =',I14)
99987 FORMAT( ' INFOG(12)  NUMBER OF NEGATIVE PIVOTS     =',I14)
99988 FORMAT( ' NUMBER OF 2x2 PIVOTS in type 1 nodes     =',I14)
99989 FORMAT( ' NUMBER OF 2x2 PIVOTS in type 2 nodes     =',I14)
      END SUBROUTINE SMUMPS_FAC_DRIVER
      SUBROUTINE SMUMPS_AVGMAX_STAT8(PROKG, MPG, VAL, NSLAVES,
     &     COMM, MSG)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      LOGICAL PROKG
      INTEGER MPG
      INTEGER(8) VAL
      INTEGER NSLAVES
      INTEGER COMM
      CHARACTER*42 MSG 
C  Local
      INTEGER(8) MAX_VAL
      INTEGER IERR, MASTER
      REAL LOC_VAL, AVG_VAL
      PARAMETER(MASTER=0)
C
      CALL MUMPS_REDUCEI8( VAL, MAX_VAL, MPI_MAX, MASTER, COMM)
      LOC_VAL = real(VAL)/real(NSLAVES)
      CALL MPI_REDUCE( LOC_VAL, AVG_VAL, 1, MPI_REAL,
     &                 MPI_SUM, MASTER, COMM, IERR )
      IF (PROKG) THEN
        WRITE(MPG,100) " Maximum ", MSG, MAX_VAL
        WRITE(MPG,100) " Average ", MSG, int(AVG_VAL,8)
      ENDIF
      RETURN
 100  FORMAT(A9,A42,I14)
      END SUBROUTINE SMUMPS_AVGMAX_STAT8
C
      SUBROUTINE SMUMPS_EXTRACT_SCHUR_REDRHS(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C     Extract the Schur and possibly also the reduced right-hand side
C     (if Fwd in facto) from the processor working on Schur and copy
C     it into the user datastructures id%SCHUR and id%REDRHS on the host.
C     This routine assumes that the integer list of the Schur has not
C     been permuted and still corresponds to LISTVAR_SCHUR.
C
C     If the Schur is centralized, the master of the Schur holds the
C     Schur and possibly also the reduced right-hand side.
C     If the Schur is distribued (already built in user's datastructure),
C     then the master of the Schur may hold the reduced right-hand side,
C     in which case it is available in root%RHS_CNTR_MASTER_ROOT.
C     
      TYPE(SMUMPS_STRUC) :: id
C
C  Local variables
C  ===============
C
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mumps_headers.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      INTEGER, PARAMETER :: MASTER = 0
      INTEGER :: ID_SCHUR, SIZE_SCHUR, LD_SCHUR, IB, BL4
      INTEGER :: ROW_LENGTH, I
      INTEGER(8) :: SURFSCHUR8, BL8, SHIFT8
      INTEGER(8) :: ISCHUR_SRC, ISCHUR_DEST, ISCHUR_SYM, ISCHUR_UNS
C
C  External functions
C  ==================
C
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
C     Quick return in case factorization did not terminate correctly
      IF (id%INFO(1) .LT. 0) RETURN
C     Quick return if Schur option off
      IF (id%KEEP(60) .EQ. 0) RETURN
C     Get Schur id
      ID_SCHUR =MUMPS_PROCNODE(
     &    id%PROCNODE_STEPS(id%STEP(max(id%KEEP(20),id%KEEP(38)))),
     &    id%NSLAVES)
      IF ( id%KEEP( 46 )  .NE. 1 ) THEN
        ID_SCHUR = ID_SCHUR + 1
      END IF
C     Get size of Schur
      IF (id%MYID.EQ.ID_SCHUR) THEN
        IF (id%KEEP(60).EQ.1) THEN
C         Sequential Schur
          LD_SCHUR =
     &    id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))+2+id%KEEP(IXSZ))
          SIZE_SCHUR = LD_SCHUR - id%KEEP(253)
        ELSE
C         Parallel Schur
          LD_SCHUR   = -999999 ! not used
          SIZE_SCHUR = id%root%TOT_ROOT_SIZE
        ENDIF
      ELSE IF (id%MYID .EQ. MASTER) THEN
        SIZE_SCHUR = id%KEEP(116)
        LD_SCHUR = -44444 ! Not used
      ELSE
C       Proc is not concerned with Schur, return
        RETURN
      ENDIF
      SURFSCHUR8 = int(SIZE_SCHUR,8)*int(SIZE_SCHUR,8)
C     =================================
C     Case of parallel Schur: if REDRHS
C     was requested, obtain it directly
C     from id%root%RHS_CNTR_MASTER_ROOT
C     =================================
      IF (id%KEEP(60) .GT. 1) THEN
        IF (id%KEEP(221).EQ.1) THEN !Implies Fwd in facto
          DO I = 1, id%KEEP(253)
            IF (ID_SCHUR.EQ.MASTER) THEN ! Necessarily = id%MYID
              CALL scopy(SIZE_SCHUR,
     &             id%root%RHS_CNTR_MASTER_ROOT((I-1)*SIZE_SCHUR+1), 1,
     &             id%REDRHS((I-1)*id%LREDRHS+1), 1)
            ELSE
              IF (id%MYID.EQ.ID_SCHUR) THEN
C               Send
                CALL MPI_SEND(
     &             id%root%RHS_CNTR_MASTER_ROOT((I-1)*SIZE_SCHUR+1),
     &             SIZE_SCHUR,
     &             MPI_REAL,
     &             MASTER, TAG_SCHUR,
     &             id%COMM, IERR )
              ELSE ! MYID.EQ.MASTER
C               Receive
                CALL MPI_RECV( id%REDRHS((I-1)*id%LREDRHS+1),
     &             SIZE_SCHUR,
     &             MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &             id%COMM, STATUS, IERR )
              ENDIF
            ENDIF
          ENDDO
C         ------------------------------
C         In case of parallel Schur, we
C         free root%RHS_CNTR_MASTER_ROOT
C         ------------------------------
          IF (id%MYID.EQ.ID_SCHUR) THEN
            DEALLOCATE(id%root%RHS_CNTR_MASTER_ROOT)
            NULLIFY   (id%root%RHS_CNTR_MASTER_ROOT)
          ENDIF
        ENDIF
C       return because this is all we need to do
C       in case of parallel Schur complement
        RETURN
      ENDIF
C     ============================
C     Centralized Schur complement
C     ============================
C     PTRAST has been freed at the moment of calling this
C     routine. Schur is available through
C     PTRFAC(IW( PTLUST_S( STEP(KEEP(20)) ) + 4 +KEEP(IXSZ) ))
      IF (id%KEEP(252).EQ.0) THEN
C       CASE 1 (ORIGINAL CODE):
C       Schur is contiguous on ID_SCHUR
        IF ( ID_SCHUR .EQ. MASTER ) THEN ! Necessarily equals id%MYID
C         ---------------------
C         Copy Schur complement
C         ---------------------
          CALL SMUMPS_COPYI8SIZE( SURFSCHUR8,
     &      id%S(id%PTRFAC(id%STEP(id%KEEP(20)))),
     &      id%SCHUR(1) )
        ELSE
C         -----------------------------------------
C         The processor responsible of the Schur
C         complement sends it to the host processor
C         -----------------------------------------
          BL8=int(huge(BL4)/id%KEEP(35)/10,8)
          DO IB=1, int((SURFSCHUR8+BL8-1_8) / BL8)
            SHIFT8 = int(IB-1,8) * BL8                ! Where to send
            BL4    = int(min(BL8,SURFSCHUR8-SHIFT8)) ! Size of block
            IF ( id%MYID .eq. ID_SCHUR ) THEN
C             Send Schur complement
              CALL MPI_SEND( id%S( SHIFT8 +
     &          id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ)))),
     &          BL4,
     &          MPI_REAL,
     &          MASTER, TAG_SCHUR,
     &          id%COMM, IERR )
            ELSE IF ( id%MYID .eq. MASTER ) THEN
C             Receive Schur complement
              CALL MPI_RECV( id%SCHUR(1_8 + SHIFT8),
     &                     BL4,
     &                     MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &                     id%COMM, STATUS, IERR )
            END IF
          ENDDO
        END IF
      ELSE
C       CASE 2 (Fwd in facto): Schur is not contiguous on ID_SCHUR,
C       process it row by row.
C
C       2.1: We first centralize Schur complement into id%SCHUR
        ISCHUR_SRC = id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &               +4+id%KEEP(IXSZ)))
        ISCHUR_DEST= 1_8
        DO I=1, SIZE_SCHUR
          ROW_LENGTH = SIZE_SCHUR
          IF (ID_SCHUR.EQ.MASTER) THEN ! Necessarily = id%MYID
            CALL scopy(ROW_LENGTH, id%S(ISCHUR_SRC), 1,
     &                 id%SCHUR(ISCHUR_DEST),1)
          ELSE
            IF (id%MYID.EQ.ID_SCHUR) THEN
C             Send
              CALL MPI_SEND( id%S(ISCHUR_SRC), ROW_LENGTH,
     &        MPI_REAL,
     &        MASTER, TAG_SCHUR,
     &        id%COMM, IERR )
            ELSE
C             Recv
              CALL MPI_RECV( id%SCHUR(ISCHUR_DEST),
     &                   ROW_LENGTH,
     &                   MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &                   id%COMM, STATUS, IERR )
            ENDIF
          ENDIF
          ISCHUR_SRC = ISCHUR_SRC+int(LD_SCHUR,8)
          ISCHUR_DEST= ISCHUR_DEST+int(SIZE_SCHUR,8)
        ENDDO
C       2.2: Get REDRHS on host
C       2.2.1: Symmetric => REDRHS is available in last KEEP(253)
C              rows of Schur structure on ID_SCHUR
C       2.2.2: Unsymmetric => REDRHS corresponds to last KEEP(253)
C              columns. However it must be transposed.
        IF (id%KEEP(221).EQ.1) THEN ! Implies Fwd in facto
          ISCHUR_SYM = id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ))) + int(SIZE_SCHUR,8) *
     &                    int(LD_SCHUR,8)
          ISCHUR_UNS =
     &                 id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ))) + int(SIZE_SCHUR,8)
          ISCHUR_DEST = 1_8
          DO I = 1, id%KEEP(253)
            IF (ID_SCHUR .EQ. MASTER) THEN ! necessarily = id%MYID
              IF (id%KEEP(50) .EQ. 0) THEN
                CALL scopy(SIZE_SCHUR, id%S(ISCHUR_UNS), LD_SCHUR,
     &                     id%REDRHS(ISCHUR_DEST), 1)
              ELSE
                CALL scopy(SIZE_SCHUR, id%S(ISCHUR_SYM), 1,
     &                     id%REDRHS(ISCHUR_DEST), 1)
              ENDIF
            ELSE
              IF (id%MYID .NE. MASTER) THEN
                IF (id%KEEP(50) .EQ. 0) THEN
C                 Use id%S(ISCHUR_SYM) as temporary contig. workspace
C                 of size SIZE_SCHUR. 
                  CALL scopy(SIZE_SCHUR, id%S(ISCHUR_UNS), LD_SCHUR,
     &            id%S(ISCHUR_SYM), 1)
                ENDIF
                CALL MPI_SEND(id%S(ISCHUR_SYM), SIZE_SCHUR,
     &          MPI_REAL, MASTER, TAG_SCHUR,
     &          id%COMM, IERR )
              ELSE
                CALL MPI_RECV(id%REDRHS(ISCHUR_DEST),
     &          SIZE_SCHUR, MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &          id%COMM, STATUS, IERR )
              ENDIF
            ENDIF
            IF (id%KEEP(50).EQ.0) THEN
              ISCHUR_UNS = ISCHUR_UNS + int(LD_SCHUR,8)
            ELSE
              ISCHUR_SYM = ISCHUR_SYM + int(LD_SCHUR,8)
            ENDIF
            ISCHUR_DEST = ISCHUR_DEST + int(id%LREDRHS,8)
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_EXTRACT_SCHUR_REDRHS
