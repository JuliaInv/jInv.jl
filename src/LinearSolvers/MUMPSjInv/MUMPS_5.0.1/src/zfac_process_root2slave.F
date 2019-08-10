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
      SUBROUTINE ZMUMPS_PROCESS_ROOT2SLAVE( TOT_ROOT_SIZE,
     &    TOT_CONT_TO_RECV, root,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM, COMM_LOAD,
     &    NBPROCFILS,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID, SLAVEF,
     &
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, PTRARW, PTRAIW,
     &    INTARR, DBLARR, ICNTL, KEEP, KEEP8, DKEEP, ND)
      USE ZMUMPS_LOAD
      USE ZMUMPS_OOC        
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'zmumps_root.h'
      TYPE (ZMUMPS_ROOT_STRUC) :: root
      INTEGER KEEP(500), ICNTL(40)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION DKEEP(130)
      INTEGER TOT_ROOT_SIZE, TOT_CONT_TO_RECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER(8) :: IPTRLU, LRLU, LRLUS, LA, POSFAC
      INTEGER(8) :: PTRFAC(KEEP(28)), PTRAST(KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      COMPLEX(kind=8) A( LA )
      INTEGER PTRIST(KEEP(28)), PTLUST(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP
      INTEGER NSTK_S( KEEP(28) ), PROCNODE_STEPS( KEEP(28) )
      INTEGER NBPROCFILS( KEEP(28) ), ND( KEEP(28) )
      INTEGER IFLAG, IERROR, COMM, COMM_LOAD
      INTEGER LPOOL, LEAF
      INTEGER IPOOL( LPOOL )
      INTEGER MYID, SLAVEF, NBFIN
      DOUBLE PRECISION OPASSW, OPELIW
      INTEGER ITLOC(N+KEEP(253)), FILS(N), PTRARW(N), PTRAIW(N)
      COMPLEX(kind=8) :: RHS_MUMPS(KEEP(255))
      INTEGER INTARR(max(1,KEEP(14)))
      COMPLEX(kind=8) DBLARR(max(1,KEEP(13)))
      INTEGER ::  allocok
      COMPLEX(kind=8), DIMENSION(:,:), POINTER :: TMP
      INTEGER NEW_LOCAL_M, NEW_LOCAL_N
      INTEGER OLD_LOCAL_M, OLD_LOCAL_N
      INTEGER I, J
      INTEGER LREQI, IROOT
      INTEGER(8) :: LREQA
      INTEGER POSHEAD, IPOS_SON,IERR
      LOGICAL MASTER_OF_ROOT
      COMPLEX(kind=8) ZERO
      PARAMETER( ZERO = (0.0D0,0.0D0) )
      INCLUDE 'mumps_headers.h'
      INTEGER numroc, MUMPS_PROCNODE
      EXTERNAL numroc, MUMPS_PROCNODE
      IROOT = KEEP( 38 )
      root%TOT_ROOT_SIZE = TOT_ROOT_SIZE
      MASTER_OF_ROOT = ( MYID .EQ. 
     &                   MUMPS_PROCNODE( PROCNODE_STEPS(STEP(IROOT)),
     &                   SLAVEF ) )
      NEW_LOCAL_M  = numroc( TOT_ROOT_SIZE, root%MBLOCK,
     &               root%MYROW, 0, root%NPROW )
      NEW_LOCAL_M  = max( 1, NEW_LOCAL_M )
      NEW_LOCAL_N  = numroc( TOT_ROOT_SIZE, root%NBLOCK,
     &               root%MYCOL, 0, root%NPCOL )
      IF ( PTRIST(STEP( IROOT )).GT.0) THEN
        OLD_LOCAL_N = -IW( PTRIST(STEP( IROOT )) + KEEP(IXSZ) )
        OLD_LOCAL_M =  IW( PTRIST(STEP( IROOT )) + 1  + KEEP(IXSZ))
      ELSE
        OLD_LOCAL_N = 0
        OLD_LOCAL_M = NEW_LOCAL_M
      ENDIF
      IF (KEEP(60) .NE. 0) THEN
        IF (root%yes) THEN
        IF ( NEW_LOCAL_M .NE. root%SCHUR_MLOC .OR.
     &       NEW_LOCAL_N .NE. root%SCHUR_NLOC ) THEN
          WRITE(*,*) "Internal error 1 in ZMUMPS_PROCESS_ROOT2SLAVE"
          CALL MUMPS_ABORT()
        ENDIF
        ENDIF
        PTLUST(STEP(IROOT)) = -4444
        PTRFAC(STEP(IROOT)) = -4445_8
        PTRIST(STEP(IROOT)) = 0
        IF ( MASTER_OF_ROOT ) THEN
          LREQI=6+2*TOT_ROOT_SIZE+KEEP(IXSZ)
          LREQA=0_8
          IF ( IWPOS + LREQI - 1. GT. IWPOSCB ) THEN
           CALL ZMUMPS_COMPRE_NEW( N, KEEP(28), IW, LIW, A, LA,
     &           LRLU, IPTRLU,
     &           IWPOS, IWPOSCB, PTRIST, PTRAST,
     &           STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &           KEEP(IXSZ),COMP,DKEEP(97),MYID )
           IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) 'PB1 compress root2slave:LRLU,LRLUS=',
     &            LRLU, LRLUS
                  IFLAG = -9
                  CALL MUMPS_SET_IERROR(LREQA-LRLUS, IERROR)
                  GOTO 700
           END IF
          ENDIF
          IF ( IWPOS + LREQI - 1. GT. IWPOSCB ) THEN
            IFLAG = -8
            IERROR = IWPOS + LREQI - 1 - IWPOSCB
            GOTO 700
          ENDIF
          PTLUST(STEP(IROOT))= IWPOS
          IWPOS = IWPOS + LREQI
          POSHEAD = PTLUST( STEP(IROOT))
          IW( POSHEAD + XXI )=LREQI
          CALL MUMPS_STOREI8( LREQA, IW(POSHEAD + XXR))
          IW( POSHEAD + XXS )=-9999
          IW(POSHEAD+XXS+1:POSHEAD+KEEP(IXSZ)-1)=-99999
          IW( POSHEAD +KEEP(IXSZ)) = 0
          IW( POSHEAD + 1 +KEEP(IXSZ)) = -1
          IW( POSHEAD + 2 +KEEP(IXSZ)) = -1
          IW( POSHEAD + 4 +KEEP(IXSZ)) = STEP(IROOT)
          IW( POSHEAD + 5 +KEEP(IXSZ)) = 0
          IW( POSHEAD + 3 +KEEP(IXSZ)) = TOT_ROOT_SIZE
        ENDIF
        GOTO 100
      ENDIF
      IF ( MASTER_OF_ROOT ) THEN
        LREQI = 6 + 2 * TOT_ROOT_SIZE+KEEP(IXSZ)
      ELSE
        LREQI = 6+KEEP(IXSZ)
      END IF
      LREQA = int(NEW_LOCAL_M, 8) * int(NEW_LOCAL_N, 8)
      IF ( LRLU . LT. LREQA .OR.
     &     IWPOS + LREQI - 1. GT. IWPOSCB )THEN
           IF ( LRLUS .LT. LREQA ) THEN
             IFLAG  = -9
             CALL MUMPS_SET_IERROR(LREQA - LRLUS, IERROR)
             GOTO 700
           END IF
           CALL ZMUMPS_COMPRE_NEW( N, KEEP(28), IW, LIW, A, LA,
     &           LRLU, IPTRLU,
     &           IWPOS, IWPOSCB, PTRIST, PTRAST,
     &           STEP, PIMASTER, PAMASTER, KEEP(216), LRLUS,
     &           KEEP(IXSZ), COMP, DKEEP(97), MYID )
           IF ( LRLU .NE. LRLUS ) THEN
                  WRITE(*,*) 'PB2 compress root2slave:LRLU,LRLUS=',
     &            LRLU, LRLUS
                  IFLAG = -9
                  CALL MUMPS_SET_IERROR(LREQA - LRLUS, IERROR)
                  GOTO 700
           END IF
           IF ( IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
              IFLAG  = -8
              IERROR = IWPOS + LREQI - 1 - IWPOSCB
              GOTO 700
           END IF
      END IF
      PTLUST(STEP( IROOT )) = IWPOS
      IWPOS           = IWPOS + LREQI
      IF (LREQA.EQ.0_8) THEN
        PTRAST (STEP(IROOT)) = max(POSFAC-1_8,1_8)
        PTRFAC (STEP(IROOT)) = max(POSFAC-1_8,1_8)
      ELSE
        PTRAST (STEP(IROOT)) = POSFAC
        PTRFAC (STEP(IROOT)) = POSFAC
      ENDIF
      POSFAC           = POSFAC + LREQA
      LRLU   = LRLU  - LREQA
      LRLUS  = LRLUS - LREQA
      KEEP8(67) = min(KEEP8(67), LRLUS)
      CALL ZMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &          LA-LRLUS,0_8,LREQA,KEEP,KEEP8,LRLUS)
      POSHEAD = PTLUST( STEP(IROOT))
      IW( POSHEAD + XXI )     = LREQI
      CALL MUMPS_STOREI8( LREQA, IW(POSHEAD + XXR))
      IW( POSHEAD + XXS ) = S_NOTFREE
      IW(POSHEAD+XXS+1:POSHEAD+KEEP(IXSZ)-1)=-99999
      IW( POSHEAD + KEEP(IXSZ) ) = 0
      IW( POSHEAD + 1 + KEEP(IXSZ) ) = NEW_LOCAL_N
      IW( POSHEAD + 2 + KEEP(IXSZ) ) = NEW_LOCAL_M
      IW( POSHEAD + 4 + KEEP(IXSZ) ) = STEP(IROOT)
      IW( POSHEAD + 5 + KEEP(IXSZ) ) = 0
      IF ( MASTER_OF_ROOT ) THEN
        IW( POSHEAD + 3 + KEEP(IXSZ) ) = TOT_ROOT_SIZE
      ELSE
        IW( POSHEAD + 3 + KEEP(IXSZ) ) = 0
      ENDIF
      IF ( PTRIST(STEP( IROOT )) .LE. 0 ) THEN
        PTRIST(STEP( IROOT ))            = 0
        PAMASTER(STEP( IROOT ))          = 0_8
        IF (LREQA.GT.0_8) A(PTRAST(STEP(IROOT)):
     &    PTRAST(STEP(IROOT))+LREQA-1_8) = ZERO
      ELSE
        OLD_LOCAL_N = -IW( PTRIST(STEP( IROOT )) + KEEP(IXSZ) )
        OLD_LOCAL_M =  IW( PTRIST(STEP( IROOT )) + 1  + KEEP(IXSZ))
        IF ( TOT_ROOT_SIZE .eq. root%ROOT_SIZE ) THEN
          IF ( LREQA .NE. int(OLD_LOCAL_M,8) * int(OLD_LOCAL_N,8) )
     &    THEN
             write(*,*) 'error 1 in PROCESS_ROOT2SLAVE',
     &       OLD_LOCAL_M, OLD_LOCAL_N
             CALL MUMPS_ABORT()
          END IF
          CALL ZMUMPS_COPYI8SIZE(LREQA,
     &                          A( PAMASTER(STEP(IROOT)) ),
     &                          A( PTRAST  (STEP(IROOT)) ) )
        ELSE
          CALL ZMUMPS_COPY_ROOT( A( PTRAST(STEP(IROOT))), 
     &        NEW_LOCAL_M,
     &        NEW_LOCAL_N, A( PAMASTER( STEP(IROOT)) ), OLD_LOCAL_M,
     &        OLD_LOCAL_N )
        END IF
        IF ( PTRIST( STEP( IROOT ) ) .GT. 0 ) THEN
           IPOS_SON= PTRIST( STEP(IROOT))
           CALL ZMUMPS_FREE_BLOCK_CB(.FALSE., MYID, N, IPOS_SON,
     &          PAMASTER(STEP(IROOT)),
     &          IW, LIW, LRLU, LRLUS, IPTRLU,
     &          IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &         )
           PTRIST(STEP( IROOT ))   = 0
           PAMASTER(STEP( IROOT )) = 0_8
        END IF
      END IF
       IF (NEW_LOCAL_M.GT.OLD_LOCAL_M) THEN
          TMP => root%RHS_ROOT
          NULLIFY(root%RHS_ROOT)
          ALLOCATE (root%RHS_ROOT(NEW_LOCAL_M, root%RHS_NLOC), 
     &                stat=allocok)
          IF ( allocok.GT.0) THEN
              IFLAG=-13
              IERROR = NEW_LOCAL_M*root%RHS_NLOC
              GOTO 700
          ENDIF
          DO J = 1, root%RHS_NLOC
            DO I = 1, OLD_LOCAL_M
              root%RHS_ROOT(I,J)=TMP(I,J)
            ENDDO
            DO I = OLD_LOCAL_M+1, NEW_LOCAL_M
              root%RHS_ROOT(I,J) = ZERO
            ENDDO
          ENDDO
          DEALLOCATE(TMP)
          NULLIFY(TMP) 
       ENDIF
 100  CONTINUE
      NBPROCFILS(STEP(IROOT))=NBPROCFILS(STEP(IROOT)) + TOT_CONT_TO_RECV
#if ! defined(NO_XXNBPR)
      KEEP(121) = KEEP(121) + TOT_CONT_TO_RECV
#endif
#if ! defined(NO_XXNBPR)
      CALL CHECK_EQUAL(NBPROCFILS(STEP(IROOT)), KEEP(121))
      IF ( KEEP(121) .eq. 0 ) THEN
#else
      IF ( NBPROCFILS(STEP(IROOT)) .eq. 0 ) THEN
#endif
         IF (KEEP(201).EQ.1) THEN 
            CALL ZMUMPS_OOC_FORCE_WRT_BUF_PANEL(IERR)
         ELSE IF (KEEP(201).EQ.2) THEN 
            CALL ZMUMPS_FORCE_WRITE_BUF(IERR)              
         ENDIF
        CALL ZMUMPS_INSERT_POOL_N( N, IPOOL, LPOOL, PROCNODE_STEPS,
     &       SLAVEF, KEEP(28), KEEP(76), KEEP(80), KEEP(47),
     &       STEP, IROOT + N )
        IF (KEEP(47) .GE. 3) THEN
           CALL ZMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &          IPOOL, LPOOL, 
     &          PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &          MYID, STEP, N, ND, FILS )
        ENDIF
      END IF
      RETURN
 700  CONTINUE
      CALL ZMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE ZMUMPS_PROCESS_ROOT2SLAVE
      SUBROUTINE ZMUMPS_COPY_ROOT
     &( NEW, M_NEW, N_NEW,OLD, M_OLD, N_OLD )
      INTEGER M_NEW, N_NEW, M_OLD, N_OLD
      COMPLEX(kind=8) NEW( M_NEW, N_NEW ), OLD( M_OLD, N_OLD )
      INTEGER J
      COMPLEX(kind=8) ZERO
      PARAMETER( ZERO = (0.0D0,0.0D0) )
      DO J = 1, N_OLD
        NEW( 1: M_OLD, J ) = OLD( 1: M_OLD, J )
        NEW( M_OLD + 1: M_NEW, J ) = ZERO
      END DO
      NEW( 1: M_NEW,N_OLD + 1: N_NEW ) = ZERO
      RETURN
      END SUBROUTINE ZMUMPS_COPY_ROOT
