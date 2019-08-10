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
      SUBROUTINE ZMUMPS_SOL_S(N, A, LA, IW, LIW, W, LWC,
     &    NRHS, 
     &    RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD,
     &    PTRICB, PTRACB, IWCB, LIWW, W2, 
     &    NE_STEPS, NA, LNA, STEP,
     &    FRERE, DAD, FILS, IPOOL, LPOOL, PTRIST, PTRFAC, 
     &    MYLEAF, ICNTL, INFO, 
     &    PROCNODE_STEPS,
     &    SLAVEF, COMM,MYID, BUFR, LBUFR, LBUFR_BYTES,
     &    KEEP,KEEP8, RHS_ROOT, LRHS_ROOT, MTYPE, 
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, PANEL_POS, LPANEL_POS
     &    , TO_PROCESS, SIZE_TO_PROCESS
     &    , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE 
     &    )
      USE ZMUMPS_OOC
      USE ZMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER(8) :: LA
      INTEGER N,LIW,LIWW,LWC,LPOOL,LNA
      INTEGER SLAVEF,MYLEAF,COMM,MYID
      INTEGER LPANEL_POS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER NA(LNA),NE_STEPS(KEEP(28))
      INTEGER IPOOL(LPOOL)
      INTEGER PANEL_POS(LPANEL_POS)
      INTEGER ICNTL(40), INFO(40)
      INTEGER PTRIST(KEEP(28)),
     &        PTRICB(KEEP(28)),PTRACB(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER NRHS
      COMPLEX(kind=8) A(LA), W(LWC)
      COMPLEX(kind=8) W2(KEEP(133))
      INTEGER IW(LIW),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)),DAD(KEEP(28)),FILS(N)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER LRHSCOMP, POSINRHSCOMP_BWD(N)
#if defined(RHSCOMP_BYROWS)
      COMPLEX(kind=8) RHSCOMP(NRHS,LRHSCOMP)
#else
      COMPLEX(kind=8) RHSCOMP(LRHSCOMP,NRHS)
#endif
      INTEGER LRHS_ROOT
      COMPLEX(kind=8) RHS_ROOT( LRHS_ROOT )
      INTEGER, intent(in)           :: SIZE_TO_PROCESS
      LOGICAL, intent(in)           :: TO_PROCESS(SIZE_TO_PROCESS)
      LOGICAL, intent(in) :: DO_NBSPARSE
      INTEGER, intent(in) :: LRHS_BOUNDS
      INTEGER, intent(in) :: RHS_BOUNDS(LRHS_BOUNDS)
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      LOGICAL FLAG
      INTEGER POSIWCB,POSWCB,K
      INTEGER(8) :: APOS, IST
      INTEGER NPIV
      INTEGER IPOS,LIELL,NELIM,IFR,JJ,I
      INTEGER J1,J2,J,NCB,NBFINF
      INTEGER NBLEAF,INODE,NBROOT,NROOT,NBFILS
      INTEGER IN,IF,LONG,POOL_FIRST_POS,TMP
      INTEGER III,IIPOOL,MYLEAFE
      INTEGER NSLAVES
      INTEGER JBDEB, JBFIN, NRHS_B
      COMPLEX(kind=8) ALPHA,ONE,ZERO
      PARAMETER (ZERO=(0.0D0,0.0D0),
     &           ONE=(1.0D0,0.0D0),
     &           ALPHA=(-1.0D0,0.0D0))
      LOGICAL BLOQ,DEBUT
      INTEGER PROCDEST, DEST
      INTEGER POSINDICES, IPOSINRHSCOMP
      INTEGER DUMMY(1)
      INTEGER PLEFTW, PTWCB
      INTEGER Offset, EffectiveSize, ISLAVE, FirstIndex
      LOGICAL LTLEVEL2, IN_SUBTREE
      INTEGER TYPENODE
      INCLUDE 'mumps_headers.h'
      LOGICAL BLOCK_SEQUENCE
      INTEGER TMP_NBPANELS, I_PIVRPTR, I_PIVR
      LOGICAL MUST_BE_PERMUTED
      LOGICAL NO_CHILDREN
      LOGICAL Exploit_Sparsity, AM1 
      LOGICAL, DIMENSION(:), ALLOCATABLE :: DEJA_SEND
      INTEGER :: allocok
      INTEGER(8) :: APOSDEB, NBENTRIES_ALLPANELS
      INTEGER LDAJ, NBJ, LIWFAC,
     &        NBJLAST, NPIV_LAST, PANEL_SIZE,
     &        PTWCB_PANEL, NCB_PANEL, TYPEF
      INTEGER BEG_PANEL
      LOGICAL TWOBYTWO
      INTEGER NPANELS, IPANEL
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR
      INTEGER MUMPS_TYPENODE
      EXTERNAL zgemv, ztrsv, ztrsm, zgemm,
     &         MUMPS_TYPENODE, 
     &         MUMPS_IN_OR_ROOT_SSARBR
      DUMMY(1)=0
      ALLOCATE(DEJA_SEND( 0:SLAVEF-1 ), stat=allocok)
      if(allocok.ne.0) then
         WRITE(6,*) ' Allocation error of DEJA_SEND in '
     &        //'routine ZMUMPS_SOL_S '
         INFO(1)=-13
         INFO(2)=SLAVEF
      endif
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM, MYID )
      IF ( INFO(1) .LT.0 ) GOTO 340
      PLEFTW = 1
      POSIWCB = LIWW
      POSWCB = LWC
      NROOT = 0
      NBLEAF = NA(1)
      NBROOT = NA(2)
      DO I = NBROOT, 1, -1
        INODE = NA(NBLEAF+I+2)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &      SLAVEF) .EQ. MYID) THEN
          NROOT = NROOT + 1
          IPOOL(NROOT) = INODE
        ENDIF
      END DO
      III = 1
      IIPOOL = NROOT + 1
      BLOCK_SEQUENCE = .FALSE.
      Exploit_Sparsity = .FALSE.
      AM1 = .FALSE.
      IF (KEEP(235).NE.0) Exploit_Sparsity = .TRUE.
      IF (KEEP(237).NE.0) AM1 = .TRUE.
      NO_CHILDREN = .FALSE.
      IF (Exploit_Sparsity .OR. AM1) MYLEAF = -1
      IF (MYLEAF .EQ. -1) THEN
        MYLEAF = 0
        DO I=1, NBLEAF
          INODE=NA(I+2)
          IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &         SLAVEF) .EQ. MYID) THEN
            MYLEAF = MYLEAF + 1
          ENDIF
        ENDDO
      ENDIF
      MYLEAFE=MYLEAF
      NBFINF = SLAVEF
      IF (MYLEAFE .EQ. 0) THEN
        CALL ZMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID, COMM, FEUILLE,
     &                  SLAVEF)
        NBFINF = NBFINF - 1
        IF (NBFINF .EQ. 0) THEN
          GOTO 340
        ENDIF
      ENDIF
 50   CONTINUE
      BLOQ = ( (  III .EQ. IIPOOL  )
     &     )
      CALL ZMUMPS_BACKSLV_RECV_AND_TREAT( BLOQ, FLAG, BUFR, LBUFR,
     &     LBUFR_BYTES, MYID, SLAVEF, COMM,
     &     N, IWCB, LIWW, POSIWCB,
     &     W, LWC, POSWCB,
     &     IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &     IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &     STEP,  FRERE, FILS, PROCNODE_STEPS,
     &     PLEFTW, KEEP,KEEP8,
     &     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &     NRHS, MTYPE, 
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &     , TO_PROCESS, SIZE_TO_PROCESS
     &     )
      IF ( INFO(1) .LT. 0 ) GOTO 340
      IF ( .NOT. FLAG ) THEN
        IF (III .NE. IIPOOL) THEN
          INODE = IPOOL(IIPOOL-1)
          IIPOOL = IIPOOL - 1
          GO TO 60
        ENDIF
      END IF                    
      IF ( NBFINF .eq. 0 ) GOTO 340
      GOTO 50
      IF (MYID.EQ.0) write(6,*) "BWD: process INODE=", INODE
   60 CONTINUE
      IF (DO_NBSPARSE) THEN
       JBDEB= RHS_BOUNDS(2*STEP(INODE)-1)
       JBFIN= RHS_BOUNDS(2*STEP(INODE))
       NRHS_B = JBFIN-JBDEB+1
      ELSE
       JBDEB = 1
       JBFIN = NRHS
       NRHS_B = NRHS
      ENDIF
      IF ( INODE .EQ. KEEP( 38 ) .OR. INODE .EQ. KEEP( 20 ) ) THEN
         IPOS = PTRIST(STEP(INODE))+KEEP(IXSZ)
          NPIV  = IW(IPOS+3)
          LIELL = IW(IPOS) + NPIV  
         IPOS =  PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ)
         IF ( MTYPE .EQ. 1 .AND. KEEP(50) .EQ. 0) THEN
            J1   = IPOS + LIELL + 1
            J2   = IPOS + LIELL + NPIV
         ELSE
            J1   = IPOS + 1
            J2   = IPOS + NPIV
         END IF
         IFR  = 0
         IPOSINRHSCOMP =  POSINRHSCOMP_BWD(IW(J1))  
         DO JJ = J1, J2
            IFR = IFR + 1
            DO K=JBDEB,JBFIN
#if defined(RHSCOMP_BYROWS)
               RHSCOMP(K,IPOSINRHSCOMP) = RHS_ROOT(IFR+NPIV*(K-1))
#else
               RHSCOMP(IPOSINRHSCOMP,K) = RHS_ROOT(IFR+NPIV*(K-1))
#endif
            END DO
            IPOSINRHSCOMP =  IPOSINRHSCOMP + 1  
         END DO 
         IN = INODE
 270     IN = FILS(IN)
         IF (IN .GT. 0) GOTO 270
         IF (IN .EQ. 0) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
               CALL ZMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &              FEUILLE, SLAVEF )
               NBFINF = NBFINF - 1
               IF (NBFINF .EQ. 0) GOTO 340
            ENDIF
            GOTO 50
         ENDIF
         IF   = -IN
         LONG = NPIV
         NBFILS = NE_STEPS(STEP(INODE))
         IF ( AM1 ) THEN
            I = NBFILS
            NBFILS = 0
            DO WHILE (I.GT.0)
               IF ( TO_PROCESS(STEP(IF)) ) NBFILS = NBFILS+1
               IF = FRERE(STEP(IF))
               I = I -1
            ENDDO
            IF (NBFILS.EQ.0) THEN
               NO_CHILDREN = .TRUE.
            ELSE
               NO_CHILDREN = .FALSE.
            ENDIF
            IF = -IN
         ENDIF
         DEBUT = .TRUE.
         DO I = 0, SLAVEF - 1
            DEJA_SEND( I ) = .FALSE.
         END DO
         POOL_FIRST_POS=IIPOOL
         DO I = 1, NBFILS
            IF ( AM1 ) THEN
 1030          IF ( .NOT.TO_PROCESS(STEP(IF)) ) THEN
                  IF = FRERE(STEP(IF))
                  GOTO 1030
               ENDIF
               NO_CHILDREN = .FALSE.
            ENDIF
            IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IF)),SLAVEF)
     &           .EQ. MYID) THEN
                  IPOOL(IIPOOL) = IF
                  IIPOOL = IIPOOL + 1
            ELSE
               PROCDEST = MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IF)),
     &              SLAVEF)
               IF (.NOT. DEJA_SEND( PROCDEST ))  THEN
 600              CALL ZMUMPS_BUF_SEND_VCB( NRHS_B, IF, 0, 0,
     &                 LONG, LONG, IW( J1 ),
     &                 RHS_ROOT( 1+NPIV*(JBDEB-1) ), 
     &                 JBDEB, JBFIN, PROCDEST,
     &                 NOEUD, COMM, IERR )
                  IF ( IERR .EQ. -1 ) THEN
                     CALL ZMUMPS_BACKSLV_RECV_AND_TREAT(
     &                    .FALSE., FLAG,
     &                    BUFR, LBUFR, LBUFR_BYTES,
     &                    MYID, SLAVEF, COMM,
     &                    N, IWCB, LIWW, POSIWCB,
     &                    W, LWC, POSWCB,
     &                    IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &                    IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &                    STEP, FRERE, FILS, PROCNODE_STEPS,
     &                    PLEFTW, KEEP,KEEP8,
     &                    PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &                    NRHS, MTYPE,
     &                    RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &                    , TO_PROCESS, SIZE_TO_PROCESS
     &                    )
                     IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                     GOTO 600
                  ELSE IF ( IERR .EQ. -2 ) THEN
                     INFO( 1 ) = -17
                     INFO( 2 ) = NRHS_B * LONG * KEEP(35) +
     &                    ( LONG + 4 ) * KEEP(34)
                     GOTO 330
                  ELSE IF ( IERR .EQ. -3 ) THEN
                     INFO( 1 ) = -20
                     INFO( 2 ) = NRHS_B * LONG * KEEP(35) +
     &                    ( LONG + 4 ) * KEEP(34)
                     GOTO 330
                  END IF
                  DEJA_SEND( PROCDEST ) = .TRUE.
               END IF
               IF ( IERR .NE. 0 ) CALL MUMPS_ABORT()
            ENDIF
            IF = FRERE(STEP(IF))
         ENDDO
         IF (AM1 .AND.NO_CHILDREN) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
               CALL ZMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &              FEUILLE, SLAVEF )
               NBFINF = NBFINF - 1
               IF (NBFINF .EQ. 0) GOTO 340
               GOTO 50
            ENDIF
         ENDIF
            IF (IIPOOL.NE.POOL_FIRST_POS) THEN
               DO I=1,(IIPOOL-POOL_FIRST_POS)/2
                  TMP=IPOOL(POOL_FIRST_POS+I-1)
                  IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
                  IPOOL(IIPOOL-I)=TMP
               ENDDO
            ENDIF
         GOTO 50
      END IF
      IN_SUBTREE = MUMPS_IN_OR_ROOT_SSARBR( 
     &               PROCNODE_STEPS(STEP(INODE)), SLAVEF ) 
      TYPENODE = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),
     &         SLAVEF)
      LTLEVEL2= ( 
     &   (TYPENODE .eq.2 ) .AND.
     &   (MTYPE.NE.1)   )
      NPIV = IW(PTRIST(STEP(INODE))+2+KEEP(IXSZ)+1)
      IF ((NPIV.NE.0).AND.(LTLEVEL2)) THEN
            IPOS  = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
            LIELL = IW(IPOS-2)+IW(IPOS+1)
            NELIM = IW(IPOS-1)
            IPOS  = IPOS + 1
            NPIV  = IW(IPOS)
            NCB   = LIELL - NPIV - NELIM
            IPOS  = IPOS + 2
            NSLAVES = IW( IPOS )
            Offset = 0  
            IPOS = IPOS + NSLAVES   
            IW(PTRIST(STEP(INODE))+XXS)= C_FINI+NSLAVES
           IF ( POSIWCB - 2 .LT. 0 .or.
     &          POSWCB - NCB*NRHS_B .LT. PLEFTW - 1 ) THEN
             CALL ZMUMPS_COMPSO( N, KEEP(28), IWCB, LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
             IF ( POSWCB - NCB*NRHS_B .LT. PLEFTW - 1 ) THEN
               INFO( 1 ) = -11
               INFO( 2 ) = NCB * NRHS_B - POSWCB - PLEFTW + 1
               GOTO 330
             END IF
             IF ( POSIWCB - 2 .LT. 0 ) THEN
               INFO( 1 ) = -14
               INFO( 2 ) = 2 - POSIWCB
               GO TO 330
             END IF
           END IF
           POSIWCB = POSIWCB - 2
           POSWCB  = POSWCB - NCB*NRHS_B
           PTRICB(STEP( INODE )) = POSIWCB + 1
           PTRACB(STEP( INODE )) = POSWCB  + 1
           IWCB( PTRICB(STEP( INODE ))     ) = NCB*NRHS_B
           IWCB( PTRICB(STEP( INODE )) + 1 ) = 1  
           IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
              POSINDICES = IPOS + LIELL + 1
           ELSE
              POSINDICES = IPOS + 1
           END IF
           IF ( NCB.EQ.0 ) THEN
             write(6,*) ' Internal Error type 2 node with no CB '
             CALL MUMPS_ABORT()
           ENDIF
           IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
               J1 = IPOS + LIELL + NPIV + NELIM +1
               J2 = IPOS + 2 * LIELL
           ELSE
               J1 = IPOS + NPIV + NELIM +1
               J2 = IPOS + LIELL
           END IF
           IFR = PTRACB(STEP( INODE )) - 1
           DO JJ = J1, J2 - KEEP(253)
               J = IW(JJ)
               IFR = IFR + 1
               IPOSINRHSCOMP =  abs(POSINRHSCOMP_BWD(J))
               DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
                 W(IFR+(K-JBDEB)*NCB) = RHSCOMP(K,IPOSINRHSCOMP)
#else
                 W(IFR+(K-JBDEB)*NCB) = RHSCOMP(IPOSINRHSCOMP,K)
#endif
               ENDDO
           ENDDO
           IF (KEEP(252).NE.0) THEN
             DO JJ = J2-KEEP(253)+1, J2
              IFR = IFR + 1
              DO K=JBDEB, JBFIN
               IF (K.EQ.JJ-J2+KEEP(253)) THEN
                 W(IFR+(K-JBDEB)*NCB) = ALPHA   
               ELSE
                 W(IFR+(K-JBDEB)*NCB) = ZERO
               ENDIF
              ENDDO
             ENDDO
           ENDIF
           DO ISLAVE = 1, NSLAVES
              CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB, 
     &                NSLAVES, 
     &                EffectiveSize,
     &                FirstIndex )
 500         DEST = IW( PTRIST(STEP(INODE))+5+ISLAVE+KEEP(IXSZ))
             CALL ZMUMPS_BUF_SEND_BACKVEC(NRHS_B, INODE,
     &             W(Offset+PTRACB(STEP(INODE))), 
     &             EffectiveSize, 
     &             NCB, DEST,
     &             BACKSLV_MASTER2SLAVE, JBDEB, JBFIN,
     &             COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                 CALL ZMUMPS_BACKSLV_RECV_AND_TREAT(
     &                .FALSE., FLAG,
     &                BUFR, LBUFR, LBUFR_BYTES,
     &                MYID, SLAVEF, COMM,
     &                N, IWCB, LIWW, POSIWCB,
     &                W, LWC, POSWCB,
     &                IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &                IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &                STEP, FRERE, FILS,
     &                PROCNODE_STEPS, PLEFTW, KEEP,KEEP8,
     &                PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &                NRHS, MTYPE,
     &                RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &                , TO_PROCESS, SIZE_TO_PROCESS
     &                )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 500
              ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = NRHS_B * EffectiveSize * KEEP(35) +
     &                            2 * KEEP(34)
                GOTO 330
              ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = NRHS_B * EffectiveSize * KEEP(35) +
     &                            2 * KEEP(34)
                GOTO 330
              END IF
              Offset = Offset + EffectiveSize
           END DO
           IWCB( PTRICB(STEP( INODE )) + 1 ) = 0
           CALL ZMUMPS_FREETOPSO(N, KEEP(28), IWCB, LIWW, W, LWC,
     &             POSWCB,POSIWCB,PTRICB,PTRACB)
           GOTO 50
      ENDIF   
      IPOS = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
      LIELL = IW(IPOS-2)+IW(IPOS+1)
      NELIM = IW(IPOS-1)
      IPOS = IPOS + 1
      NPIV = IW(IPOS)
      IPOS = IPOS + 1
      IF (KEEP(201).GT.0) THEN
         CALL ZMUMPS_SOLVE_GET_OOC_NODE(
     &        INODE,PTRFAC,KEEP,A,LA,STEP,
     &        KEEP8,N,MUST_BE_PERMUTED,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
      ENDIF                     
      APOS = PTRFAC(IW(IPOS))
      NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ) )
      IPOS = IPOS + 1 + NSLAVES
      IF (KEEP(201).EQ.1) THEN 
           LIWFAC =  IW(PTRIST(STEP(INODE))+XXI)
           IF (MTYPE.NE.1) THEN
            TYPEF = TYPEF_L
           ELSE
            TYPEF = TYPEF_U
           ENDIF
           PANEL_SIZE =  ZMUMPS_OOC_PANEL_SIZE( LIELL )
           IF (KEEP(50).NE.1) THEN
             CALL ZMUMPS_OOC_PP_CHECK_PERM_FREED(
     &                   IW(IPOS+1+2*LIELL),
     &                   MUST_BE_PERMUTED )
           ENDIF
      ENDIF  
      LONG = 0
      IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
        J1 = IPOS + LIELL + 1
        J2 = IPOS + NPIV + LIELL
      ELSE
        J1 = IPOS + 1
        J2 = IPOS + NPIV
      ENDIF
      IF (IN_SUBTREE) THEN
        PTWCB = PLEFTW
        IF ( POSWCB .LT. LIELL*NRHS_B ) THEN
          CALL ZMUMPS_COMPSO( N, KEEP(28), IWCB, LIWW, W, LWC,
     &                 POSWCB, POSIWCB, PTRICB, PTRACB)
          IF ( POSWCB .LT. LIELL*NRHS_B ) THEN
            INFO(1) = -11
            INFO(2) = LIELL*NRHS_B - POSWCB
            GOTO 330
          END IF
        END IF
      ELSE
        IF ( POSIWCB - 2 .LT. 0 .or.
     &     POSWCB - LIELL*NRHS_B .LT. PLEFTW - 1 ) THEN
          CALL ZMUMPS_COMPSO( N, KEEP(28), IWCB, LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
          IF ( POSWCB - LIELL*NRHS_B .LT. PLEFTW - 1 ) THEN
            INFO( 1 ) = -11
            INFO( 2 ) = LIELL * NRHS_B - POSWCB - PLEFTW + 1
            GOTO 330
          END IF
          IF ( POSIWCB - 2 .LT. 0 ) THEN
            INFO( 1 ) = -14
            INFO( 2 ) = 2 - POSIWCB
            GO TO 330
          END IF
        END IF
        POSIWCB = POSIWCB - 2
        POSWCB  = POSWCB - LIELL*NRHS_B
        PTRICB(STEP( INODE )) = POSIWCB + 1
        PTRACB(STEP( INODE )) = POSWCB  + 1
        IWCB( PTRICB(STEP( INODE ))     ) = LIELL*NRHS_B
        IWCB( PTRICB(STEP( INODE )) + 1 ) = 1  
        IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
           POSINDICES = IPOS + LIELL + 1
        ELSE
           POSINDICES = IPOS + 1
        END IF
        PTWCB = PTRACB(STEP( INODE )) 
      ENDIF
      IPOSINRHSCOMP = POSINRHSCOMP_BWD(IW(J1)) 
      IF (J2.GE.J1) THEN
       DO K=JBDEB, JBFIN
        IF (KEEP(252).NE.0) THEN
         DO JJ = J1, J2
          W(PTWCB+JJ-J1+(K-JBDEB)*LIELL) = ZERO
         ENDDO
       ELSE
         DO JJ = J1, J2
#if defined(RHSCOMP_BYROWS)
          W(PTWCB+JJ-J1+(K-JBDEB)*LIELL) = 
     &             RHSCOMP(K,IPOSINRHSCOMP+JJ-J1)
#else
          W(PTWCB+JJ-J1+(K-JBDEB)*LIELL) = 
     &             RHSCOMP(IPOSINRHSCOMP+JJ-J1,K)
#endif
         ENDDO
        ENDIF
       END DO
      ENDIF
      IFR   = PTWCB + NPIV - 1
      IF ( LIELL .GT. NPIV ) THEN
        IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
          J1 = IPOS + LIELL + NPIV + 1
          J2 = IPOS + 2 * LIELL
        ELSE
          J1 = IPOS + NPIV + 1
          J2 = IPOS + LIELL
        END IF
        DO JJ = J1, J2-KEEP(253)
          J = IW(JJ)
          IFR = IFR + 1
          IPOSINRHSCOMP = abs(POSINRHSCOMP_BWD(J))
          DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
            W(IFR+(K-JBDEB)*LIELL) = RHSCOMP(K,IPOSINRHSCOMP)
#else
            W(IFR+(K-JBDEB)*LIELL) = RHSCOMP(IPOSINRHSCOMP,K)
#endif
          ENDDO
        ENDDO
        IF (KEEP(252).NE.0) THEN
          DO JJ = J2-KEEP(253)+1, J2
           IFR = IFR + 1
           DO K=JBDEB, JBFIN
            IF (K.EQ.JJ-J2+KEEP(253)) THEN
              W(IFR+(K-JBDEB)*LIELL) = ALPHA   
            ELSE
              W(IFR+(K-JBDEB)*LIELL) = ZERO
            ENDIF
           ENDDO
          ENDDO
        ENDIF
        NCB = LIELL - NPIV
        IF (NPIV .EQ. 0) GOTO 160
      ENDIF
      IF (KEEP(201).EQ.1) THEN 
       J = NPIV / PANEL_SIZE 
       TWOBYTWO = KEEP(50).EQ.2 .AND.
     & ((TYPENODE.EQ.1.AND.KEEP(103).GT.0) .OR.
     &  (TYPENODE.EQ.2.AND.KEEP(105).GT.0))
       IF (TWOBYTWO) THEN 
         CALL ZMUMPS_BUILD_PANEL_POS(PANEL_SIZE, PANEL_POS, LPANEL_POS,
     &        IW(IPOS+1+LIELL), NPIV, NPANELS, LIELL,
     &        NBENTRIES_ALLPANELS)
       ELSE
         IF (NPIV.EQ.J*PANEL_SIZE) THEN
           NPIV_LAST = NPIV
           NBJLAST   = PANEL_SIZE
           NPANELS   = J
         ELSE
           NPIV_LAST = (J+1)* PANEL_SIZE
           NBJLAST   = NPIV-J*PANEL_SIZE
           NPANELS   = J+1
         ENDIF
            NBENTRIES_ALLPANELS =
     &  int(LIELL,8) * int(NPIV,8) 
     &  - int( ( J * ( J - 1 ) ) /2,8 ) 
     &    * int(PANEL_SIZE,8) * int(PANEL_SIZE,8) 
     &  - int(J,8)                       
     &    * int(mod(NPIV, PANEL_SIZE),8) 
     &    * int(PANEL_SIZE,8)    
         JJ=NPIV_LAST
       ENDIF
       APOSDEB = APOS + NBENTRIES_ALLPANELS 
       DO IPANEL = NPANELS, 1, -1
            IF (TWOBYTWO) THEN
              NBJ = PANEL_POS(IPANEL+1)-PANEL_POS(IPANEL)
              BEG_PANEL = PANEL_POS(IPANEL)
            ELSE
              IF (JJ.EQ.NPIV_LAST) THEN
                NBJ = NBJLAST
              ELSE
                NBJ = PANEL_SIZE
              ENDIF
              BEG_PANEL = JJ- PANEL_SIZE+1
            ENDIF
            LDAJ    = LIELL-BEG_PANEL+1 
            APOSDEB = APOSDEB - int(NBJ,8)*int(LDAJ,8)
            PTWCB_PANEL = PTWCB + BEG_PANEL - 1
            NCB_PANEL   = LDAJ - NBJ
            IF (KEEP(50).NE.1 .AND. MUST_BE_PERMUTED) THEN
              CALL ZMUMPS_GET_OOC_PERM_PTR(TYPEF, TMP_NBPANELS,
     &        I_PIVRPTR, I_PIVR, IPOS + 1 + 2 * LIELL, IW, LIW)
              IF (NPIV.EQ.(IW(I_PIVRPTR)-1)) THEN
                MUST_BE_PERMUTED=.FALSE. 
              ELSE
               CALL ZMUMPS_PERMUTE_PANEL(
     &         IW(I_PIVR + IW(I_PIVRPTR+IPANEL-1)-IW(I_PIVRPTR)),
     &         NPIV-IW(I_PIVRPTR+IPANEL-1)+1,
     &         IW(I_PIVRPTR+IPANEL-1)-1,
     &         A(APOSDEB),
     &         LDAJ, NBJ, BEG_PANEL-1)
              ENDIF
            ENDIF
#if defined(MUMPS_USE_BLAS2)
            IF ( NRHS_B == 1 ) THEN
              IF (NCB_PANEL.NE.0) THEN
                CALL zgemv( 'T', NCB_PANEL, NBJ, ALPHA, 
     &                A( APOSDEB + int(NBJ,8) ), LDAJ,
     &                W( NBJ + PTWCB_PANEL ),
     &                1, ONE,
     &                W(PTWCB_PANEL), 1 )
              ENDIF
              IF (MTYPE.NE.1) THEN
               CALL ztrsv('L','T','U', NBJ, A(APOSDEB), LDAJ,
     &              W(PTWCB_PANEL), 1)
              ELSE
               CALL ztrsv('L','T','N', NBJ, A(APOSDEB), LDAJ,
     &              W(PTWCB_PANEL), 1)
              ENDIF
            ELSE
#endif
              IF (NCB_PANEL.NE.0) THEN
                 CALL zgemm( 'T', 'N', NBJ, NRHS_B, NCB_PANEL, ALPHA,
     &              A(APOSDEB +int(NBJ,8)), LDAJ,
     &              W(NBJ+PTWCB_PANEL),LIELL,
     &              ONE, W(PTWCB_PANEL),LIELL)
              ENDIF
              IF (MTYPE.NE.1) THEN
               CALL ztrsm('L','L','T','U',NBJ, NRHS_B, ONE, 
     &           A(APOSDEB), 
     &           LDAJ, W(PTWCB_PANEL), LIELL)
              ELSE
               CALL ztrsm('L','L','T','N',NBJ, NRHS_B, ONE, 
     &           A(APOSDEB), 
     &           LDAJ, W(PTWCB_PANEL), LIELL)
              ENDIF
#if defined(MUMPS_USE_BLAS2)
            ENDIF
#endif
            IF (.NOT. TWOBYTWO) JJ=BEG_PANEL-1 
       ENDDO 
      ENDIF 
      IF (KEEP(201).EQ.0.OR.KEEP(201).EQ.2)THEN 
       IF ( LIELL .GT. NPIV ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          IST = APOS + int(NPIV,8)
#if defined(MUMPS_USE_BLAS2)
          IF (NRHS_B == 1) THEN
            CALL zgemv( 'T', NCB, NPIV, ALPHA, A(IST), LIELL,
     &              W(NPIV + PTWCB), 1,
     &              ONE,
     &              W(PTWCB), 1 )
          ELSE
#endif
            CALL zgemm('T','N', NPIV, NRHS_B, NCB, ALPHA, A(IST), LIELL,
     &              W(NPIV+PTWCB), LIELL, ONE,
     &              W(PTWCB), LIELL)
#if defined(MUMPS_USE_BLAS2)
          ENDIF
#endif
        ELSE
          IF ( KEEP(50) .eq. 0 ) THEN
            IST = APOS + int(NPIV,8) * int(LIELL,8)
          ELSE
            IST = APOS + int(NPIV,8) * int(NPIV,8)
          END IF
#if defined(MUMPS_USE_BLAS2)
          IF ( NRHS_B == 1 ) THEN
              CALL zgemv( 'N', NPIV, NCB, ALPHA, A( IST ), NPIV,
     &                W( NPIV + PTWCB ),
     &                1, ONE,
     &                W(PTWCB), 1 )
          ELSE
#endif
                CALL zgemm( 'N', 'N', NPIV, NRHS_B, NCB, ALPHA,
     &                A(IST), NPIV, W(NPIV+PTWCB),LIELL,
     &                ONE, W(PTWCB),LIELL)
#if defined(MUMPS_USE_BLAS2)
          END IF
#endif
        END IF 
       ENDIF  
       IF ( MTYPE .eq. 1 ) THEN
#if defined(MUMPS_USE_BLAS2)
         IF ( NRHS_B == 1 ) THEN
           CALL ztrsv('L', 'T', 'N', NPIV, A(APOS), LIELL,
     &              W(PTWCB), 1)
         ELSE
#endif
           CALL ztrsm('L','L','T','N', NPIV, NRHS_B, ONE, A(APOS),
     &              LIELL, W(PTWCB), LIELL)
#if defined(MUMPS_USE_BLAS2)
         ENDIF
#endif
       ELSE
         IF ( KEEP(50) .EQ. 0 ) THEN
           LDAJ=LIELL
         ELSE
           LDAJ=NPIV
         ENDIF
#if defined(MUMPS_USE_BLAS2)
         IF ( NRHS_B == 1 ) THEN
            CALL ztrsv('U','N','U', NPIV, A(APOS), LDAJ,
     &              W(PTWCB), 1)
         ELSE
#endif
            CALL ztrsm('L','U','N','U', NPIV, NRHS_B, ONE, A(APOS),
     &                 LDAJ,W(PTWCB),LIELL)
#if defined(MUMPS_USE_BLAS2)
         END IF
#endif
       END IF 
      ENDIF 
      IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0) THEN
        J1 = IPOS + LIELL + 1
      ELSE
        J1 = IPOS + 1
      END IF
      IPOSINRHSCOMP =  POSINRHSCOMP_BWD(IW(J1)) 
      DO 150 I = 1, NPIV
        DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
          RHSCOMP(K,IPOSINRHSCOMP) = W(PTWCB+I-1+(K-JBDEB)*LIELL)
#else
          RHSCOMP(IPOSINRHSCOMP, K) = W(PTWCB+I-1+(K-JBDEB)*LIELL)
#endif
        ENDDO
        IPOSINRHSCOMP =  IPOSINRHSCOMP + 1 
  150 CONTINUE
  160 CONTINUE
      IF (KEEP(201).GT.0) THEN
         CALL ZMUMPS_FREE_FACTORS_FOR_SOLVE(INODE,PTRFAC,KEEP(28),
     &        A,LA,.TRUE.,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
      ENDIF
      IN = INODE
  170 IN = FILS(IN)
      IF (IN .GT. 0) GOTO 170
      IF (IN .EQ. 0) THEN
        MYLEAFE = MYLEAFE - 1
        IF (MYLEAFE .EQ. 0) THEN
          CALL ZMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &                     FEUILLE, SLAVEF )
          NBFINF = NBFINF - 1
          IF (NBFINF .EQ. 0) GOTO 340
        ENDIF
        GOTO 50
      ENDIF
      IF = -IN
      NBFILS = NE_STEPS(STEP(INODE))
      IF (AM1) THEN
         I = NBFILS
         NBFILS = 0
         DO WHILE (I.GT.0)
            IF ( TO_PROCESS(STEP(IF)) ) NBFILS = NBFILS+1
            IF = FRERE(STEP(IF))
            I = I -1
         ENDDO
         IF (NBFILS.EQ.0) THEN
            NO_CHILDREN = .TRUE.
         ELSE
            NO_CHILDREN = .FALSE.
         ENDIF
         IF = -IN
      ENDIF
      IF (IN_SUBTREE) THEN
         DO I = 1, NBFILS
            IF ( AM1 ) THEN
 1010          IF ( .NOT.TO_PROCESS(STEP(IF)) )  THEN
                  IF = FRERE(STEP(IF))
                  GOTO 1010
               ENDIF
               NO_CHILDREN = .FALSE.
            ENDIF
               IPOOL((IIPOOL-I+1)+NBFILS-I) = IF
               IIPOOL = IIPOOL + 1
            IF = FRERE(STEP(IF))
         ENDDO
         IF (AM1 .AND. NO_CHILDREN) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
               CALL ZMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &              FEUILLE, SLAVEF )
               NBFINF = NBFINF - 1
               IF (NBFINF .EQ. 0) GOTO 340
               GOTO 50
            ENDIF
         ENDIF
      ELSE
        DEBUT = .TRUE.
        DO I = 0, SLAVEF - 1
          DEJA_SEND( I ) = .FALSE.
        END DO
        POOL_FIRST_POS=IIPOOL
        DO 190 I = 1, NBFILS
           IF ( AM1 ) THEN
1020      IF ( .NOT.TO_PROCESS(STEP(IF)) ) THEN
                 IF = FRERE(STEP(IF))
                 GOTO 1020
              ENDIF
              NO_CHILDREN = .FALSE.
           ENDIF
          IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IF)),
     &      SLAVEF) .EQ. MYID) THEN
                IPOOL(IIPOOL) = IF
                IIPOOL = IIPOOL + 1
            IF = FRERE(STEP(IF))
          ELSE
            PROCDEST = MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IF)),SLAVEF)
            IF (.not. DEJA_SEND( PROCDEST ))  THEN
 400          CONTINUE
              CALL ZMUMPS_BUF_SEND_VCB( NRHS_B, IF, 0, 0, LIELL,
     &         LIELL - KEEP(253),
     &         IW( POSINDICES ), 
     &         W   ( PTRACB(STEP( INODE )) ), 
     &         JBDEB, JBFIN, PROCDEST,
     &         NOEUD, COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                CALL ZMUMPS_BACKSLV_RECV_AND_TREAT(
     &          .FALSE., FLAG,
     &          BUFR, LBUFR, LBUFR_BYTES,
     &          MYID, SLAVEF, COMM,
     &          N, IWCB, LIWW, POSIWCB,
     &          W, LWC, POSWCB,
     &          IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &          IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &          STEP, FRERE, FILS, PROCNODE_STEPS,
     &          PLEFTW, KEEP,KEEP8,
     &          PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &          NRHS, MTYPE, 
     &          RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &                , TO_PROCESS, SIZE_TO_PROCESS
     &                )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 400
              ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = NRHS_B * LIELL * KEEP(35) + 4 * KEEP(34)
                GOTO 330
              ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = NRHS_B * LIELL * KEEP(35) + 4 * KEEP(34)
                GOTO 330
              END IF
              DEJA_SEND( PROCDEST ) = .TRUE.
            END IF
            IF = FRERE(STEP(IF))
          ENDIF
  190   CONTINUE
        IF (AM1 .AND. NO_CHILDREN) THEN
           MYLEAFE = MYLEAFE - 1
           IF (MYLEAFE .EQ. 0) THEN
              CALL ZMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &             FEUILLE, SLAVEF )
              NBFINF = NBFINF - 1
              IF (NBFINF .EQ. 0) GOTO 340
              GOTO 50
           ENDIF
        ENDIF
           DO I=1,(IIPOOL-POOL_FIRST_POS)/2
              TMP=IPOOL(POOL_FIRST_POS+I-1)
              IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
              IPOOL(IIPOOL-I)=TMP
           ENDDO 
        IWCB(PTRICB(STEP(INODE))+1) = IWCB(PTRICB(STEP(INODE))+1)-1
        CALL ZMUMPS_FREETOPSO(N, KEEP(28), IWCB, LIWW, 
     &     W, LWC,
     &     POSWCB,POSIWCB,PTRICB,PTRACB)
      ENDIF
      GOTO 50
  330 CONTINUE
      CALL ZMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID, COMM, TERREUR,
     & SLAVEF)
  340 CONTINUE
      CALL ZMUMPS_FINISH_RECV( MYID,COMM,BUFR,
     &                            LBUFR,LBUFR_BYTES )
      IF (ALLOCATED(DEJA_SEND)) DEALLOCATE(DEJA_SEND)
      RETURN
      END SUBROUTINE ZMUMPS_SOL_S
