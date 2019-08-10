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
      RECURSIVE SUBROUTINE ZMUMPS_TRAITER_MESSAGE_SOLVE
     &     ( BUFR, LBUFR, LBUFR_BYTES,
     &     MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     &     N, NRHS, IPOOL, LPOOL, III, LEAF,
     &     NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,
     &     PTRFAC, IWCB, LIWCB,
     &     WCB, LWCB, POSWCB,
     &     PLEFTWCB, POSIWCB,
     &     PTRICB,
     &     INFO, KEEP,KEEP8, STEP, PROCNODE_STEPS, 
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &     )
      USE ZMUMPS_OOC 
      USE ZMUMPS_COMM_BUFFER 
      IMPLICIT NONE
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIW
      INTEGER(8) :: LA
      INTEGER N, NRHS, LPOOL, III, LEAF, NBFIN, LRHSCOMP
      INTEGER LIWCB, LWCB, POSWCB, PLEFTWCB, POSIWCB
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ),  NSTK_S( N )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      INTEGER PTRICB(KEEP(28)),PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      COMPLEX(kind=8) WCB( LWCB ), A( LA )
#if defined(RHSCOMP_BYROWS)
      COMPLEX(kind=8) RHSCOMP( NRHS, LRHSCOMP )
#else
      COMPLEX(kind=8) RHSCOMP( LRHSCOMP, NRHS )
#endif
      INTEGER, intent(in) :: POSINRHSCOMP_FWD(N)
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR, K, JJ, JBDEB, JBFIN, NRHS_B
      INTEGER FINODE, FPERE, LONG, NCB, POSITION, NCV, NPIV
      INTEGER PTRX, PTRY, PDEST, I, IPOSINRHSCOMP
      INTEGER(8) :: APOS
      LOGICAL DUMMY
      LOGICAL FLAG
      EXTERNAL MUMPS_PROCNODE
      INTEGER  MUMPS_PROCNODE
      COMPLEX(kind=8) ALPHA, ONE
      PARAMETER (ONE=(1.0D0,0.0D0), ALPHA=(-1.0D0,0.0D0))
      INCLUDE 'mumps_headers.h'
      IF ( MSGTAG .EQ. RACINE_SOLVE ) THEN
         NBFIN = NBFIN - 1
         IF ( NBFIN .eq. 0 ) GOTO 270
      ELSE  IF (MSGTAG .EQ. ContVec ) THEN
         POSITION = 0
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        FINODE, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        FPERE, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        NCB, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        JBDEB, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        JBFIN, 1, MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        LONG, 1, MPI_INTEGER, COMM, IERR )
         NRHS_B = JBFIN-JBDEB+1
          IF ( NCB .eq. 0 ) THEN
             PTRICB(STEP(FINODE)) = -1
             NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
                   IPOOL( LEAF ) = FPERE
                LEAF = LEAF + 1
                IF ( LEAF > LPOOL ) THEN
                   WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
                   CALL MUMPS_ABORT()
                END IF
             END IF
          ELSE
             IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
                PTRICB(STEP(FINODE)) = NCB + 1
             END IF
             IF ( ( POSIWCB - LONG ) .LT. 0 ) THEN
                INFO( 1 ) = -14
                INFO( 2 ) = LONG
                GOTO 260
             END IF
             IF ( POSWCB - PLEFTWCB + 1 .LT. LONG * NRHS_B) THEN
                INFO( 1 ) = -11
                INFO( 2 ) = PLEFTWCB - POSWCB - 1 + LONG * NRHS_B
                GOTO 260
             END IF
             IF (LONG .GT. 0) THEN
                CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               IWCB( 1 ),
     &               LONG, MPI_INTEGER, COMM, IERR )
                DO K = 1, NRHS_B
                   CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                  WCB( PLEFTWCB ),
     &                  LONG, MPI_DOUBLE_COMPLEX, COMM, IERR )
                   DO I = 1, LONG
                   IPOSINRHSCOMP= abs(POSINRHSCOMP_FWD(IWCB(I)))
#if defined(RHSCOMP_BYROWS)
                   RHSCOMP(JBDEB+K-1,IPOSINRHSCOMP) = 
     &             RHSCOMP(JBDEB+K-1,IPOSINRHSCOMP) + 
     &                 WCB(PLEFTWCB+I-1)
#else
                   RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1) = 
     &             RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1) + 
     &                 WCB(PLEFTWCB+I-1)
#endif
                   ENDDO
                END DO
                PTRICB(STEP(FINODE)) = PTRICB(STEP(FINODE)) - LONG
             ENDIF
             IF ( PTRICB(STEP(FINODE)) == 1 ) THEN
                NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             END IF
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
                   IPOOL( LEAF ) = FPERE
                LEAF = LEAF + 1
                IF ( LEAF > LPOOL ) THEN
                   WRITE(*,*) 'Internal error 41r2 : Pool is too small.'
                   CALL MUMPS_ABORT()
                END IF
             ENDIF
          END IF
       ELSEIF ( MSGTAG .EQ. Master2Slave ) THEN
          POSITION = 0
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         FINODE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         FPERE, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         NCV, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         NPIV, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         JBDEB, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &         JBFIN, 1, MPI_INTEGER, COMM, IERR )
          NRHS_B = JBFIN-JBDEB+1
          PTRY = PLEFTWCB
          PTRX = PLEFTWCB + NCV * NRHS_B
          PLEFTWCB = PLEFTWCB + (NPIV + NCV) * NRHS_B
          IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
             INFO(1) = -11
             INFO(2) = -POSWCB + PLEFTWCB -1
             GO TO 260
          END IF
          DO K=1, NRHS_B
             CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &            WCB( PTRY + (K-1) * NCV ), NCV,
     &            MPI_DOUBLE_COMPLEX, COMM, IERR )
          ENDDO
          IF ( NPIV .GT. 0 ) THEN
             DO K=1, NRHS_B
                CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &               WCB( PTRX + (K-1)*NPIV ), NPIV,
     &               MPI_DOUBLE_COMPLEX, COMM, IERR )
             END DO
          END IF
          IF (KEEP(201).GT.0) THEN
             CALL ZMUMPS_SOLVE_GET_OOC_NODE(
     &            FINODE,PTRFAC,KEEP,A,LA,STEP,
     &            KEEP8,N,DUMMY,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
          ENDIF
          APOS = PTRFAC(STEP(FINODE))
          IF (KEEP(201).EQ.1) THEN
#if defined(MUMPS_USE_BLAS2)
             IF ( NRHS_B == 1 ) THEN
                CALL zgemv( 'N', NCV, NPIV, ALPHA, A(APOS), NCV,
     &               WCB( PTRX ), 1, ONE,
     &               WCB( PTRY ), 1 )
             ELSE
#endif
                CALL zgemm( 'N', 'N', NCV, NRHS_B, NPIV, ALPHA,
     &               A(APOS), NCV,
     &               WCB( PTRX), NPIV, ONE,
     &               WCB( PTRY), NCV )
#if defined(MUMPS_USE_BLAS2)
             ENDIF
#endif
          ELSE                  
#if defined(MUMPS_USE_BLAS2)
             IF ( NRHS_B == 1 ) THEN
                CALL zgemv( 'T', NPIV, NCV, ALPHA, A(APOS), NPIV,
     &               WCB( PTRX ), 1, ONE,
     &               WCB( PTRY ), 1 )
             ELSE
#endif
                CALL zgemm( 'T', 'N', NCV, NRHS_B, NPIV, ALPHA,
     &               A(APOS), NPIV,
     &               WCB( PTRX), NPIV, ONE,
     &               WCB( PTRY), NCV )
#if defined(MUMPS_USE_BLAS2)
             ENDIF
#endif
          ENDIF
          IF (KEEP(201).GT.0) THEN
             CALL ZMUMPS_FREE_FACTORS_FOR_SOLVE(FINODE,PTRFAC,
     &            KEEP(28),A,LA,.TRUE.,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
          ENDIF
          PLEFTWCB = PLEFTWCB - NPIV * NRHS_B
          PDEST = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(FPERE)),
     &                            SLAVEF )
          IF ( PDEST .EQ. MYID ) THEN
             IF ( PTRICB(STEP(FINODE)) .EQ. 0 ) THEN
                NCB = IW( PTRIST(STEP(FINODE)) + 2 + KEEP(IXSZ) )
                PTRICB(STEP(FINODE)) = NCB + 1
             END IF
             DO I = 1, NCV
                JJ=IW(PTRIST(STEP(FINODE))+3+I+ KEEP(IXSZ) )
                IPOSINRHSCOMP= abs(POSINRHSCOMP_FWD(JJ))
                DO K=1, NRHS_B
#if defined(RHSCOMP_BYROWS)
                   RHSCOMP(JBDEB+K-1,IPOSINRHSCOMP)= 
     &              RHSCOMP(JBDEB+K-1,IPOSINRHSCOMP) +
     &              WCB(PTRY+I-1+(K-1)*NCV)
#else
                   RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1)= 
     &              RHSCOMP(IPOSINRHSCOMP,JBDEB+K-1) +
     &              WCB(PTRY+I-1+(K-1)*NCV)
#endif
                ENDDO
             END DO
             PTRICB(STEP(FINODE)) =
     &            PTRICB(STEP(FINODE)) - NCV
             IF ( PTRICB( STEP( FINODE ) ) == 1 ) THEN
                NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
             END IF
             IF ( NSTK_S(STEP(FPERE)) .EQ. 0 ) THEN
                   IPOOL( LEAF ) = FPERE
                LEAF = LEAF + 1
                IF ( LEAF > LPOOL ) THEN
                   WRITE(*,*) 'INTERNAL Error 41r: Pool is too small.'
                   CALL MUMPS_ABORT()
                END IF
             ENDIF
          ELSE
 210         CONTINUE
             CALL ZMUMPS_BUF_SEND_VCB( NRHS_B, FINODE, FPERE,
     &            IW(PTRIST(STEP( FINODE )) + 2 + KEEP(IXSZ) ), NCV,NCV,
     &            IW(PTRIST(STEP(FINODE))+4+ KEEP(IXSZ) ),
     &            WCB( PTRY ), JBDEB, JBFIN,
     &            PDEST, ContVec, COMM, IERR )
             IF ( IERR .EQ. -1 ) THEN
                CALL ZMUMPS_SOLVE_RECV_AND_TREAT( .FALSE., FLAG,
     &               BUFR, LBUFR, LBUFR_BYTES,
     &               MYID, SLAVEF, COMM,
     &               N, NRHS, IPOOL, LPOOL, III, LEAF,
     &               NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &               IWCB, LIWCB,
     &               WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &               PTRICB, INFO, KEEP,KEEP8, STEP,
     &               PROCNODE_STEPS, 
     &               RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &               )
                IF ( INFO( 1 )  .LT. 0 )  GOTO 270
                GOTO 210
             ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     &               NCV * KEEP( 35 )
                GOTO 260
             ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = ( NCV + 4 ) * KEEP( 34 ) +
     &               NCV * KEEP( 35 )
             END IF
          END IF
          PLEFTWCB = PLEFTWCB - NCV * NRHS_B
       ELSEIF ( MSGTAG .EQ. TERREUR ) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GOTO 270
       ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     &         (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
       ELSE
          INFO(1)=-100
          INFO(2)=MSGTAG
          GO TO 260
       ENDIF
       GO TO 270
 260   CONTINUE
       CALL ZMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
 270   CONTINUE
       RETURN
       END SUBROUTINE ZMUMPS_TRAITER_MESSAGE_SOLVE
      SUBROUTINE ZMUMPS_SOLVE_NODE( INODE,
     &     BUFR, LBUFR, LBUFR_BYTES,
     &     MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     &     N, IPOOL, LPOOL, III, LEAF,
     &     NBFIN, NSTK_S,
     &     IWCB, LIWCB,
     &     WCB, LWCB, A, LA, IW, LIW,
     &     NRHS, POSWCB,
     &     PLEFTWCB, POSIWCB,
     &     PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     &     FILS, STEP, FRERE, DAD,
     &     MYROOT,
     &     INFO, KEEP,KEEP8, RHS_ROOT, MTYPE,
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD,
     &     
     &     ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &    , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE 
     &            )
      USE ZMUMPS_OOC
      USE ZMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER INODE, LBUFR, LBUFR_BYTES
      INTEGER MSGTAG, MSGSOU, MYID, SLAVEF, COMM
      INTEGER LIWCB, LWCB, LIW, POSWCB, PLEFTWCB, POSIWCB
      INTEGER(8) :: LA
      INTEGER N, LPOOL, III, LEAF, NBFIN
      INTEGER MYROOT
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      INTEGER BUFR( LBUFR )
      INTEGER IPOOL( LPOOL ), NSTK_S(KEEP(28))
      INTEGER IWCB( LIWCB ), IW( LIW )
      INTEGER NRHS
      COMPLEX(kind=8) WCB( LWCB ), A( LA )
      COMPLEX(kind=8) RHS_ROOT( * )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER FILS( N ), STEP( N ), FRERE(KEEP(28)), DAD(KEEP(28))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &     TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER POSINRHSCOMP_FWD(N), LRHSCOMP
#if defined(RHSCOMP_BYROWS)
      COMPLEX(kind=8) RHSCOMP(NRHS, LRHSCOMP)
#else
      COMPLEX(kind=8) RHSCOMP(LRHSCOMP, NRHS)
#endif
      COMPLEX(kind=8) VALPIV, A11, A22, A12, DETPIV
      LOGICAL, intent(in) :: DO_NBSPARSE
      INTEGER, intent(in) :: LRHS_BOUNDS
      INTEGER, intent(in) :: RHS_BOUNDS(LRHS_BOUNDS)
      INTEGER JBDEB, JBFIN, NRHS_B 
      EXTERNAL zgemv, ztrsv, zgemm, ztrsm, MUMPS_PROCNODE
      INTEGER MUMPS_PROCNODE
      COMPLEX(kind=8) ALPHA,ONE,ZERO
      PARAMETER (ZERO=(0.0D0,0.0D0),
     &           ONE=(1.0D0,0.0D0),
     &           ALPHA=(-1.0D0,0.0D0))
      INTEGER(8) :: APOS, APOS1, APOS2, APOSOFF
      INTEGER I, J, K, IPOS, NSLAVES, J1, J2, J3, FPERE, NPIV, NCB,
     &     IERR, IFR_ini,
     &     IFR, LIELL, JJ,
     &     NELIM, PLEFT, PCB_COURANT, PPIV_COURANT
      INTEGER IPOSINRHSCOMP
      INTEGER Effective_CB_Size, NUPDATE, ISLAVE, PDEST, FirstIndex
      LOGICAL FLAG, OMP_FLAG
      INCLUDE 'mumps_headers.h'
      INTEGER POSWCB1,POSWCB2
      INTEGER(8) :: APOSDEB
      INTEGER TempNROW, TempNCOL, PANEL_SIZE, LIWFAC, 
     &     JFIN, NBJ, NUPDATE_PANEL,
     &     PPIV_PANEL, PCB_PANEL, NBK, TYPEF
      INTEGER LD_WCBPIV         
      INTEGER LD_WCBCB          
      INTEGER LDAJ, LDAJ_FIRST_PANEL
      INTEGER TMP_NBPANELS,
     &     I_PIVRPTR, I_PIVR, IPANEL
      LOGICAL MUST_BE_PERMUTED
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY( 1 )
      DUMMY(1)=1
      IF (DO_NBSPARSE) THEN
       JBDEB= RHS_BOUNDS(2*STEP(INODE)-1)
       JBFIN= RHS_BOUNDS(2*STEP(INODE))
       NRHS_B = JBFIN-JBDEB+1
      ELSE
       JBDEB = 1
       JBFIN = NRHS
       NRHS_B = NRHS
      ENDIF
      IF (DO_NBSPARSE) THEN
       if (JBDEB.GT.JBFIN) then
         write(6,*) " Internal error 1 in nbsparse :", 
     &    JBDEB, JBFIN
         CALL MUMPS_ABORT()
       endif
       IF (JBDEB.LT.1 .OR. JBDEB.GT.NRHS .or. 
     &    JBFIN.LT.1 .OR. JBFIN.GT.NRHS ) THEN
         write(6,*) " Internal error 2 in nbsparse :", 
     &    JBDEB, JBFIN
         CALL MUMPS_ABORT()
       endif
      ENDIF
      IF ( INODE .eq. KEEP( 38 ) .OR. INODE .eq.KEEP( 20 ) ) THEN
         LIELL = IW( PTRIST( STEP(INODE)) + 3 + KEEP(IXSZ))
         NPIV  = LIELL
         NELIM = 0
         NSLAVES = 0
         IPOS = PTRIST( STEP(INODE)) + 5 + KEEP(IXSZ)
      ELSE
        IPOS = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
        LIELL = IW(IPOS-2)+IW(IPOS+1)
        NELIM = IW(IPOS-1)
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ) )
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
        IPOS = IPOS + 1
        IF (KEEP(201).GT.0) THEN
           CALL ZMUMPS_SOLVE_GET_OOC_NODE(
     &          INODE,PTRFAC,KEEP,A,LA,STEP,
     &          KEEP8,N,MUST_BE_PERMUTED,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
           IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1) THEN
           CALL ZMUMPS_OOC_PP_CHECK_PERM_FREED(
     &                 IW(IPOS+1+2*LIELL+1+NSLAVES),
     &                 MUST_BE_PERMUTED )
           ENDIF
        ENDIF                     
        NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ))
        IPOS = IPOS + 1 + NSLAVES
      END IF
      IF ( MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0 ) THEN
         J1 = IPOS + 1
         J2 = IPOS + LIELL
         J3 = IPOS + NPIV
      ELSE
         J1 = IPOS + LIELL + 1
         J2 = IPOS + 2 * LIELL
         J3 = IPOS + LIELL + NPIV
      END IF
      NCB = LIELL-NPIV
      IF ( INODE .eq. KEEP( 38 ) .OR. INODE .eq. KEEP(20) ) THEN
         IFR = 0
         IPOSINRHSCOMP = POSINRHSCOMP_FWD(IW(J1)) 
         DO JJ = J1, J3
            IFR = IFR + 1
            DO K=JBDEB,JBFIN
#if defined(RHSCOMP_BYROWS)
               RHS_ROOT(IFR+NPIV*(K-1)) = RHSCOMP
     &                                    (K,IPOSINRHSCOMP)
#else
               RHS_ROOT(IFR+NPIV*(K-1)) = RHSCOMP
     &                                    (IPOSINRHSCOMP,K) 
#endif
            END DO
            IPOSINRHSCOMP = IPOSINRHSCOMP + 1
         END DO
         IF ( NPIV .LT. LIELL ) THEN
            WRITE(*,*) ' Internal error in SOLVE_NODE for Root node'
            CALL MUMPS_ABORT()
         END IF
         MYROOT = MYROOT - 1
         IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
               CALL ZMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID,
     &              COMM, RACINE_SOLVE, SLAVEF)
            ENDIF
         END IF
         GO TO 270
      END IF
      APOS = PTRFAC(STEP(INODE))
      IF (KEEP(201).EQ.1) THEN  
        IF (MTYPE.EQ.1) THEN
            IF ((MTYPE.EQ.1).AND.NSLAVES.NE.0) THEN
              TempNROW= NPIV+NELIM
              TempNCOL= NPIV
              LDAJ_FIRST_PANEL=TempNROW
            ELSE
              TempNROW= LIELL
              TempNCOL= NPIV
              LDAJ_FIRST_PANEL=TempNROW
            ENDIF
            TYPEF=TYPEF_L
        ELSE 
            TempNCOL= LIELL
            TempNROW= NPIV
            LDAJ_FIRST_PANEL=TempNCOL
            TYPEF= TYPEF_U
        ENDIF
        LIWFAC =  IW(PTRIST(STEP(INODE))+XXI)
        PANEL_SIZE = ZMUMPS_OOC_PANEL_SIZE( LDAJ_FIRST_PANEL )
      ENDIF                     
      PLEFT    = PLEFTWCB
      PPIV_COURANT = PLEFTWCB
      PLEFTWCB = PLEFTWCB + LIELL * NRHS_B
      IF ( POSWCB - PLEFTWCB + 1 .LT. 0 ) THEN
         INFO(1) = -11
         INFO(2) = PLEFTWCB - POSWCB - 1
         GO TO 260
      END IF
      IF (KEEP(201).EQ.1) THEN  
         LD_WCBPIV = LIELL
         LD_WCBCB  = LIELL
         PCB_COURANT = PPIV_COURANT + NPIV
         DO K=JBDEB, JBFIN
            IFR = PPIV_COURANT+(K-JBDEB)*LD_WCBPIV- 1
            IPOSINRHSCOMP = POSINRHSCOMP_FWD(IW(J1)) 
            DO JJ = J1, J3
               IFR = IFR + 1
#if defined(RHSCOMP_BYROWS)
               WCB(IFR) = RHSCOMP(K,IPOSINRHSCOMP) 
#else
               WCB(IFR) = RHSCOMP(IPOSINRHSCOMP,K) 
#endif
               IPOSINRHSCOMP = IPOSINRHSCOMP + 1
            ENDDO
            IF (NCB.GT.0) THEN
               DO JJ = J3+1, J2
                  J = IW(JJ)
                  IFR = IFR + 1
                  IPOSINRHSCOMP = abs(POSINRHSCOMP_FWD(J))
#if defined(RHSCOMP_BYROWS)
                  WCB(IFR) = RHSCOMP(K,IPOSINRHSCOMP)
                  RHSCOMP (K,IPOSINRHSCOMP) = ZERO
#else
                  WCB(IFR) = RHSCOMP(IPOSINRHSCOMP,K) 
                  RHSCOMP (IPOSINRHSCOMP,K) = ZERO
#endif
               ENDDO
            ENDIF
         ENDDO
      ELSE                      
         LD_WCBPIV = NPIV
         LD_WCBCB  = NCB
         PCB_COURANT = PPIV_COURANT + NPIV*NRHS_B
         IFR = PPIV_COURANT - 1
         OMP_FLAG = NRHS_B.GT.4 .AND. .FALSE.
         IFR_ini = IFR
         IPOSINRHSCOMP = POSINRHSCOMP_FWD(IW(J1)) 
         DO 130 JJ = J1, J3
            J = IW(JJ)
            IFR = IFR_ini + (JJ-J1) + 1
            DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
               WCB(IFR+(K-JBDEB)*NPIV) = 
     &                RHSCOMP(K,IPOSINRHSCOMP+JJ-J1)
#else
               WCB(IFR+(K-JBDEB)*NPIV) = 
     &                RHSCOMP(IPOSINRHSCOMP+JJ-J1,K)
#endif
            END DO
 130     CONTINUE
         IFR = PCB_COURANT - 1
         IF (NPIV .LT. LIELL) THEN
            IFR_ini = IFR
            DO 140 JJ = J3 + 1, J2
               J = IW(JJ)
               IFR = IFR_ini + (JJ-J3)
               IPOSINRHSCOMP = abs(POSINRHSCOMP_FWD(J))
               DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
                  WCB(IFR+(K-JBDEB)*NCB) = RHSCOMP(K,IPOSINRHSCOMP)
#else
                  WCB(IFR+(K-JBDEB)*NCB) = RHSCOMP(IPOSINRHSCOMP,K)
#endif
#if defined(RHSCOMP_BYROWS)
                  RHSCOMP(K,IPOSINRHSCOMP)=ZERO
#else
                  RHSCOMP(IPOSINRHSCOMP,K)=ZERO
#endif
               ENDDO
 140        CONTINUE
         ENDIF
      ENDIF                     
      IF ( NPIV .NE. 0 ) THEN
         IF (KEEP(201).EQ.1) THEN 
        APOSDEB = APOS
        J = 1
        IPANEL = 0
  10    CONTINUE
          IPANEL = IPANEL + 1
          JFIN    = min(J+PANEL_SIZE-1, NPIV)
          IF (IW(IPOS+ LIELL + JFIN) < 0) THEN
            JFIN=JFIN+1
          ENDIF
          NBJ     = JFIN-J+1
          LDAJ    = LDAJ_FIRST_PANEL-J+1 
          IF ( (KEEP(50).NE.1).AND. MUST_BE_PERMUTED ) THEN
           CALL ZMUMPS_GET_OOC_PERM_PTR(TYPEF, TMP_NBPANELS,
     &            I_PIVRPTR, I_PIVR, IPOS+1+2*LIELL, IW, LIW)
               IF (NPIV.EQ.(IW(I_PIVRPTR+IPANEL-1)-1)) THEN
                  MUST_BE_PERMUTED=.FALSE. 
               ELSE
                  CALL ZMUMPS_PERMUTE_PANEL(
     &                 IW( I_PIVR+ IW(I_PIVRPTR+IPANEL-1)-
     &                 IW(I_PIVRPTR)),
     &                 NPIV-IW(I_PIVRPTR+IPANEL-1)+1, 
     &                 IW(I_PIVRPTR+IPANEL-1)-1,
     &                 A(APOSDEB),
     &                 LDAJ, NBJ, J-1 ) 
               ENDIF
            ENDIF 
            NUPDATE_PANEL = LDAJ - NBJ
            PPIV_PANEL = PPIV_COURANT+J-1
            PCB_PANEL  = PPIV_PANEL+NBJ
            APOS1 = APOSDEB+int(NBJ,8)
            IF  (MTYPE.EQ.1) THEN
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL ztrsv( 'L', 'N', 'U', NBJ, A(APOSDEB), LDAJ, 
     &                 WCB(PPIV_PANEL), 1 )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL zgemv('N', NUPDATE_PANEL,NBJ,ALPHA, A(APOS1),
     &                    LDAJ,  WCB(PPIV_PANEL), 1, ONE,
     &                    WCB(PCB_PANEL), 1)
                  ENDIF
               ELSE
#endif
                  CALL ztrsm( 'L','L','N','U', NBJ, NRHS_B, ONE,
     &                 A(APOSDEB), LDAJ, WCB(PPIV_PANEL),
     &                 LIELL )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL zgemm('N', 'N', NUPDATE_PANEL, NRHS_B, NBJ, 
     &                    ALPHA,
     &                    A(APOS1), LDAJ, WCB(PPIV_PANEL), LIELL, ONE,
     &                    WCB(PCB_PANEL), LIELL)
                  ENDIF
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
            ELSE
#if defined(MUMPS_USE_BLAS2)
               IF (NRHS_B == 1) THEN
                  CALL ztrsv( 'L', 'N', 'N', NBJ, A(APOSDEB), LDAJ,
     &                 WCB(PPIV_PANEL), 1 )
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL zgemv('N',NUPDATE_PANEL, NBJ, ALPHA, A(APOS1),
     &                    LDAJ, WCB(PPIV_PANEL), 1,
     &                    ONE, WCB(PCB_PANEL), 1 )
                  ENDIF
               ELSE
#endif
                  CALL ztrsm('L','L','N','N',NBJ, NRHS_B, ONE,
     &                 A(APOSDEB), LDAJ, WCB(PPIV_PANEL),
     &                 LIELL)
                  IF (NUPDATE_PANEL.GT.0) THEN
                     CALL zgemm('N', 'N', NUPDATE_PANEL, NRHS_B, NBJ, 
     &                    ALPHA,
     &                    A(APOS1), LDAJ, WCB(PPIV_PANEL), LIELL, ONE,
     &             WCB(PCB_PANEL), LIELL)
                  ENDIF
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
            ENDIF
            APOSDEB = APOSDEB+int(LDAJ,8)*int(NBJ,8)
            J=JFIN+1
            IF ( J .LE. NPIV ) GOTO 10
         ELSE                   
            IF (KEEP(50).NE.0) THEN
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL ztrsv( 'U', 'T', 'U', NPIV, A(APOS), NPIV,
     &                   WCB(PPIV_COURANT), 1 )
               ELSE
#endif
                  CALL ztrsm( 'L','U','T','U', NPIV, NRHS_B, ONE,
     &                   A(APOS), NPIV, WCB(PPIV_COURANT),
     &                   NPIV )
#if defined(MUMPS_USE_BLAS2)
               ENDIF
#endif
            ELSE
               IF ( MTYPE .eq. 1 ) THEN
#if defined(MUMPS_USE_BLAS2)
                  IF ( NRHS_B == 1)  THEN
                     CALL ztrsv( 'U', 'T', 'U', NPIV, A(APOS), LIELL, 
     &                    WCB(PPIV_COURANT), 1 )
                  ELSE
#endif
                     CALL ztrsm( 'L','U','T','U', NPIV, NRHS_B, ONE,
     &                    A(APOS), LIELL, WCB(PPIV_COURANT),
     &                    NPIV )
#if defined(MUMPS_USE_BLAS2)
                  ENDIF
#endif
               ELSE
#if defined(MUMPS_USE_BLAS2)
                  IF (NRHS_B == 1) THEN
                     CALL ztrsv( 'L', 'N', 'N', NPIV, A(APOS), LIELL,
     &                    WCB(PPIV_COURANT), 1 )
                  ELSE
#endif
                     CALL ztrsm('L','L','N','N',NPIV, NRHS_B, ONE,
     &                    A(APOS), LIELL, WCB(PPIV_COURANT),
     &                    NPIV)
#if defined(MUMPS_USE_BLAS2)
                  ENDIF
#endif
               END IF
            END IF              
         END IF                 
      END IF                    
      NCB   = LIELL - NPIV
      IF ( MTYPE .EQ. 1 ) THEN
         IF ( KEEP(50) .eq. 0 ) THEN
            APOS1 = APOS  + int(NPIV,8) * int(LIELL,8)
         ELSE
            APOS1 = APOS + int(NPIV,8) * int(NPIV,8)
         END IF
         IF ( NSLAVES .EQ. 0 .OR. NPIV .eq. 0 ) THEN
            NUPDATE = NCB
         ELSE
            NUPDATE = NELIM
         END IF
      ELSE
         APOS1 = APOS + int(NPIV,8)
         NUPDATE = NCB
      END IF
      IF (KEEP(201).NE.1) THEN  
         IF ( NPIV .NE. 0 .AND. NUPDATE.NE.0 ) THEN
            IF ( MTYPE .eq. 1 ) THEN
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL zgemv('T', NPIV, NUPDATE, ALPHA, A(APOS1),
     &            NPIV,  WCB(PPIV_COURANT), 1, ONE,
     &            WCB(PCB_COURANT), 1)
               ELSE
#endif
                  CALL zgemm('T', 'N', NUPDATE, NRHS_B, NPIV, ALPHA,
     &            A(APOS1), NPIV, WCB(PPIV_COURANT), NPIV, ONE,
     &            WCB(PCB_COURANT), NCB)
#if defined(MUMPS_USE_BLAS2)
               END IF
#endif
            ELSE                
#if defined(MUMPS_USE_BLAS2)
               IF ( NRHS_B == 1 ) THEN
                  CALL zgemv('N',NUPDATE, NPIV, ALPHA, A(APOS1),
     &                 LIELL, WCB(PPIV_COURANT), 1,
     &                 ONE, WCB(PCB_COURANT), 1 )
               ELSE
#endif
                  CALL zgemm('N', 'N', NUPDATE, NRHS_B, NPIV, ALPHA,
     &                 A(APOS1), LIELL, WCB(PPIV_COURANT), NPIV, ONE,
     &                 WCB(PCB_COURANT), NCB)
#if defined(MUMPS_USE_BLAS2)
               END IF
#endif
            END IF
         END IF
      END IF                    
      IPOSINRHSCOMP =  POSINRHSCOMP_FWD(IW(J1)) 
      IF ( KEEP(50) .eq. 0 ) THEN
#if defined(RHSCOMP_BYROWS)
         DO K=JBDEB,JBFIN
            IFR =  PPIV_COURANT + (K-JBDEB)*LD_WCBPIV
            RHSCOMP(K,IPOSINRHSCOMP:IPOSINRHSCOMP+NPIV-1) =
     &           WCB(IFR:IFR+NPIV-1)
         ENDDO
#else
         OMP_FLAG = .FALSE.
         DO K=JBDEB,JBFIN
            IFR =  PPIV_COURANT + (K-JBDEB)*LD_WCBPIV
            RHSCOMP(IPOSINRHSCOMP:IPOSINRHSCOMP+NPIV-1, K) =
     &           WCB(IFR:IFR+NPIV-1)
         ENDDO
#endif
      ELSE
         IFR = PPIV_COURANT - 1
         IF (KEEP(201).EQ.1) THEN 
            LDAJ = TempNROW  
         ELSE                
            LDAJ = NPIV 
         ENDIF
         APOS1 = APOS
         JJ    = J1
         IF (KEEP(201).EQ.1) THEN
            NBK   = 0           
         ENDIF
         DO 
            IF(JJ .GT. J3) EXIT
            IFR = IFR + 1
            IF(IW(JJ+LIELL) .GT. 0) THEN
               VALPIV  = ONE/A( APOS1 )
               DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
                  RHSCOMP(K, IPOSINRHSCOMP+JJ-J1) = 
     &                 WCB( IFR+(K-JBDEB)*LD_WCBPIV ) * VALPIV
#else
                  RHSCOMP(IPOSINRHSCOMP+JJ-J1 , K ) = 
     &                 WCB( IFR+(K-JBDEB)*LD_WCBPIV ) * VALPIV
#endif
               END DO
            IF (KEEP(201).EQ.1) THEN
              NBK = NBK+1
              IF (NBK.EQ.PANEL_SIZE) THEN
                NBK = 0
                LDAJ = LDAJ - PANEL_SIZE
              ENDIF
            ENDIF
            APOS1 = APOS1 + int(LDAJ + 1,8)
            JJ = JJ+1
         ELSE
            IF (KEEP(201).EQ.1) THEN
              NBK = NBK+1
            ENDIF
            APOS2 = APOS1+int(LDAJ+1,8)
            IF (KEEP(201).EQ.1) THEN
              APOSOFF = APOS1+int(LDAJ,8)
            ELSE
              APOSOFF=APOS1+1_8
            ENDIF
               A11 = A(APOS1)
               A22 = A(APOS2)
               A12 = A(APOSOFF)
               DETPIV = A11*A22 - A12**2
               A22 = A11/DETPIV
               A11 = A(APOS2)/DETPIV
               A12 = -A12/DETPIV
               DO K=JBDEB, JBFIN
                  POSWCB1 = IFR+(K-JBDEB)*LD_WCBPIV
                  POSWCB2 = POSWCB1+1
#if defined(RHSCOMP_BYROWS)
                  RHSCOMP(K,IPOSINRHSCOMP+JJ-J1) =
     &               WCB(POSWCB1)*A11
     &               + WCB(POSWCB2)*A12
                  RHSCOMP(K,IPOSINRHSCOMP+JJ-J1+1) = 
     &                 WCB(POSWCB1)*A12
     &                 + WCB(POSWCB2)*A22
#else
                  RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) =
     &               WCB(POSWCB1)*A11
     &               + WCB(POSWCB2)*A12
                  RHSCOMP(IPOSINRHSCOMP+JJ-J1+1,K) = 
     &                 WCB(POSWCB1)*A12
     &                 + WCB(POSWCB2)*A22
#endif
               END DO
               IF (KEEP(201).EQ.1) THEN
                  NBK = NBK+1
                  IF (NBK.GE.PANEL_SIZE) THEN
                     LDAJ = LDAJ - NBK
                     NBK = 0
                  ENDIF
               ENDIF
               APOS1 = APOS2 + int(LDAJ + 1,8)
               JJ = JJ+2
               IFR = IFR+1
            ENDIF
         ENDDO
      END IF
      IF (KEEP(201).GT.0) THEN
         CALL ZMUMPS_FREE_FACTORS_FOR_SOLVE(INODE,PTRFAC,KEEP(28),
     &        A,LA,.TRUE.,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 260
         ENDIF
      END IF
      FPERE = DAD(STEP(INODE))
      IF ( FPERE .EQ. 0 ) THEN
         MYROOT = MYROOT - 1
         PLEFTWCB = PLEFTWCB - LIELL *NRHS_B
         IF ( MYROOT .EQ. 0 ) THEN
            NBFIN = NBFIN - 1
            IF (SLAVEF .GT. 1) THEN
               CALL ZMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID,
     &             COMM, RACINE_SOLVE, SLAVEF)
            ENDIF
         END IF
         GO TO 270
      ENDIF
      IF ( NUPDATE .NE. 0 .OR. NCB.eq.0 ) THEN
         IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(FPERE)),
     &        SLAVEF) .EQ. MYID) THEN
            IF ( NCB .ne. 0 ) THEN
               PTRICB(STEP(INODE)) = NCB + 1
        OMP_FLAG = .FALSE. 
               DO 190 I = 1, NUPDATE
                  IPOSINRHSCOMP = abs(POSINRHSCOMP_FWD(IW(J3 + I)))
                  DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
                     RHSCOMP( K, IPOSINRHSCOMP ) = 
     &                    RHSCOMP( K, IPOSINRHSCOMP )
#else
                     RHSCOMP( IPOSINRHSCOMP, K ) = 
     &                    RHSCOMP( IPOSINRHSCOMP, K )
#endif
     &             + WCB(PCB_COURANT + I-1 +(K-JBDEB)*LD_WCBCB)
                  ENDDO
 190           CONTINUE
               PTRICB(STEP( INODE )) = PTRICB(STEP( INODE )) - NUPDATE
               IF ( PTRICB(STEP(INODE)) == 1 ) THEN
                  NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
                  IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN
                        IPOOL( LEAF ) = FPERE
                     LEAF = LEAF + 1
                  ENDIF
               END IF
            ELSE
               PTRICB(STEP( INODE )) = -1
               NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
               IF (NSTK_S(STEP(FPERE)) .EQ. 0) THEN 
                     IPOOL( LEAF ) = FPERE 
                  LEAF = LEAF + 1
               ENDIF            
            ENDIF
         ELSE
 210        CONTINUE
            CALL ZMUMPS_BUF_SEND_VCB( NRHS_B, INODE, FPERE, 
     &           NCB, LD_WCBCB,
     &           NUPDATE,
     &           IW( J3 + 1 ), WCB( PCB_COURANT ), JBDEB, JBFIN,
     &           MUMPS_PROCNODE(PROCNODE_STEPS(STEP(FPERE)), SLAVEF),
     &           ContVec,
     &           COMM, IERR )
            IF ( IERR .EQ. -1 ) THEN
               CALL ZMUMPS_SOLVE_RECV_AND_TREAT( .FALSE., FLAG,
     &              BUFR, LBUFR, LBUFR_BYTES,
     &              MYID, SLAVEF, COMM,
     &              N, NRHS, IPOOL, LPOOL, III, LEAF,
     &              NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &              IWCB, LIWCB,
     &              WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &              PTRICB, INFO, KEEP,KEEP8, STEP,
     &              PROCNODE_STEPS,
     &              RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &              )
               IF ( INFO( 1 )  .LT. 0 )  GOTO 270
               GOTO 210
            ELSE IF ( IERR .EQ. -2 ) THEN
               INFO( 1 ) = -1-17
               INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     &              ( NUPDATE + 3 ) * KEEP( 34 )
               GOTO 260
            ELSE IF ( IERR .EQ. -3 ) THEN
               INFO( 1 ) = -20
               INFO( 2 ) = NUPDATE * KEEP( 35 ) +
     &              ( NUPDATE + 3 ) * KEEP( 34 )
               GOTO 260
            END IF
         ENDIF
      END IF
      IF ( NSLAVES .NE. 0 .AND. MTYPE .eq. 1
     &     .and. NPIV .NE. 0 ) THEN
         DO ISLAVE = 1, NSLAVES
            PDEST = IW( PTRIST(STEP(INODE)) + 5 + ISLAVE +KEEP(IXSZ))
            CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &           KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &           ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &           ISLAVE, NCB - NELIM, 
     &           NSLAVES, 
     &           Effective_CB_Size, FirstIndex )
 222        CALL ZMUMPS_BUF_SEND_MASTER2SLAVE( NRHS_B,
     &           INODE, FPERE,
     &           Effective_CB_Size, LD_WCBCB, LD_WCBPIV, NPIV,
     &           JBDEB, JBFIN, 
     &           WCB( PCB_COURANT + NELIM + FirstIndex - 1 ),
     &           WCB( PPIV_COURANT ),
     &           PDEST, COMM, IERR )
            IF ( IERR .EQ. -1 ) THEN
               CALL ZMUMPS_SOLVE_RECV_AND_TREAT( .FALSE., FLAG,
     &              BUFR, LBUFR, LBUFR_BYTES,
     &              MYID, SLAVEF, COMM,
     &              N, NRHS, IPOOL, LPOOL, III, LEAF,
     &              NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     &              IWCB, LIWCB,
     &              WCB, LWCB, POSWCB, PLEFTWCB, POSIWCB,
     &              PTRICB, INFO, KEEP,KEEP8, STEP,
     &              PROCNODE_STEPS,
     &              RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &              )
               IF ( INFO( 1 )  .LT. 0 )  GOTO 270
               GOTO 222
            ELSE IF ( IERR .EQ. -2 ) THEN
               INFO( 1 ) = -17
               INFO( 2 ) = (NPIV+Effective_CB_Size)*NRHS_B*KEEP(35) +
     &               6 * KEEP( 34 )
               GOTO 260
            ELSE IF ( IERR .EQ. -3 ) THEN
               INFO( 1 ) = -20
               INFO( 2 ) = (NPIV+Effective_CB_Size)*NRHS_B*KEEP(35) +
     &              6 * KEEP( 34 )
               GOTO 260
            END IF
         END DO
      END IF
      PLEFTWCB = PLEFTWCB - LIELL*NRHS_B
 270  CONTINUE
      RETURN
 260  CONTINUE
      CALL ZMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE ZMUMPS_SOLVE_NODE
      RECURSIVE SUBROUTINE ZMUMPS_SOLVE_RECV_AND_TREAT( BLOQ, FLAG,
     &           BUFR, LBUFR, LBUFR_BYTES,
     &           MYID, SLAVEF, COMM,
     &           N, NRHS, IPOOL, LPOOL, III, LEAF,
     &           NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST,PTRFAC,
     &           IWCB, LIWCB,
     &           WCB, LWCB, POSWCB,
     &           PLEFTWCB, POSIWCB,
     &           PTRICB, INFO, KEEP,KEEP8, STEP, PROCNODE_STEPS,
     &           RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &            )
      IMPLICIT NONE
      LOGICAL BLOQ
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, NRHS, LPOOL, III, LEAF, NBFIN
      INTEGER LIWCB, LWCB, POSWCB, PLEFTWCB, POSIWCB
      INTEGER LIW
      INTEGER(8) :: LA
      INTEGER INFO( 40 ), KEEP( 500)
      INTEGER(8) KEEP8(150)
      INTEGER BUFR( LBUFR ), IPOOL(LPOOL)
      INTEGER NSTK_S( KEEP(28) )
      INTEGER IWCB( LIWCB )
      INTEGER IW( LIW )
      COMPLEX(kind=8) WCB( LWCB ), A( LA )
      INTEGER PTRICB(KEEP(28)), PTRIST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28))
      LOGICAL FLAG
      INTEGER LRHSCOMP, POSINRHSCOMP_FWD(N)
#if defined(RHSCOMP_BYROWS)
      COMPLEX(kind=8) RHSCOMP(NRHS,LRHSCOMP)
#else
      COMPLEX(kind=8) RHSCOMP(LRHSCOMP,NRHS)
#endif
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER MSGSOU, MSGTAG, MSGLEN
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
        CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     &                   COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     &                   FLAG, STATUS, IERR )
      END IF
      IF ( FLAG ) THEN
         MSGSOU = STATUS( MPI_SOURCE )
         MSGTAG = STATUS( MPI_TAG )
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL ZMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
         ELSE
           CALL MPI_RECV( BUFR, LBUFR_BYTES, MPI_PACKED,
     &                  MSGSOU, MSGTAG, COMM, STATUS, IERR )
           CALL ZMUMPS_TRAITER_MESSAGE_SOLVE( BUFR, LBUFR, LBUFR_BYTES,
     &          MSGTAG, MSGSOU, MYID, SLAVEF, COMM,
     &          N, NRHS, IPOOL, LPOOL, III, LEAF,
     &          NBFIN, NSTK_S, IW, LIW, A, LA, PTRIST, PTRFAC,
     &          IWCB, LIWCB,
     &          WCB, LWCB, POSWCB,
     &          PLEFTWCB, POSIWCB,
     &          PTRICB, INFO, KEEP,KEEP8, STEP,
     &          PROCNODE_STEPS, 
     &          RHSCOMP, LRHSCOMP, POSINRHSCOMP_FWD
     &          )
         END IF
      END IF
      RETURN
      END SUBROUTINE ZMUMPS_SOLVE_RECV_AND_TREAT
