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
       MODULE MUMPS_BUILD_SORT_INDEX_M
       CONTAINS
      SUBROUTINE MUMPS_BUILD_SORT_INDEX(MYID, INODE, N, IOLDPS,
     &           HF, NFRONT, NFRONT_EFF, PERM, DAD,
     &           NASS1, NASS, NUMSTK, NUMORG, IWPOSCB, 
     &           IFSON, STEP, PIMASTER, PTRIST, PTRAIW, IW, LIW, 
     &           INTARR, ITLOC, FILS, FRERE_STEPS, 
     &           SON_LEVEL2, NIV1, NBPROCFILS, KEEP,KEEP8, IFLAG,
     &           ISON_IN_PLACE, PROCNODE_STEPS, SLAVEF,
     &           SONROWS_PER_ROW, LSONROWS_PER_ROW
     & )
      IMPLICIT NONE
      INTEGER INODE, N, IOLDPS, HF, NFRONT, NASS1, LIW, NASS,
     &        NUMSTK, NUMORG, IFSON, MYID
      INTEGER, intent(in) ::  ISON_IN_PLACE
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER STEP(N), PIMASTER(KEEP(28)), PTRIST(KEEP(28)),
     &        PTRAIW(N),
     &        ITLOC(N+KEEP(253)), FILS(N), FRERE_STEPS(KEEP(28)),
     &        NBPROCFILS(KEEP(28)), PERM(N)
      INTEGER, TARGET :: IW(LIW)
      INTEGER, INTENT(IN), TARGET :: IWPOSCB
      INTEGER INTARR(max(1,KEEP(14)))
      LOGICAL, intent(in)    ::  NIV1
      INTEGER, intent(inout) :: IFLAG
      LOGICAL, intent(out)   :: SON_LEVEL2
      INTEGER, intent(out)   :: NFRONT_EFF
      INTEGER, intent(in)    :: DAD (KEEP(28))
      INTEGER, intent(in) :: PROCNODE_STEPS(KEEP(28)), SLAVEF
      INTEGER, intent(in)    :: LSONROWS_PER_ROW
      INTEGER, intent(out)   :: SONROWS_PER_ROW(LSONROWS_PER_ROW)
      INTEGER NELIM_SON_IN_PLACE 
      INTEGER NEWEL, IOLDP2, INEW, INEW1,
     &        IN, NTOTFS, ICT11, NELIM, NPIVS, NSLSON, NCOLS,
     &        ITRANS, J, JJ, J1, J2, J3, JT1, ISON, IELL, LSTK, 
     &        NROWS, HS, IP1, IP2, K1, K2, IBROT, IORG, 
     &        I, K, JDEBROW, ILOC, NEWEL_SAVE, NEWEL1_SAVE,
     &        LAST_J_ASS, JMIN, MIN_PERM
      LOGICAL LEVEL1_SON
#if ! defined(NO_XXNBPR)
      INTEGER INBPROCFILS_SON
#endif
      INTEGER TYPESPLIT
      INCLUDE 'mumps_headers.h'
      INTEGER, POINTER :: SON_IWPOSCB
      INTEGER, POINTER, DIMENSION(:) :: SON_IW
      INTEGER allocok
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PTTRI, PTLAST
      INTEGER  MUMPS_TYPESPLIT, MUMPS_TYPENODE
      EXTERNAL MUMPS_TYPESPLIT, MUMPS_TYPENODE 
#if ! defined(NO_XXNBPR)
      IW(IOLDPS+XXNBPR) = 0
#endif
      TYPESPLIT  = MUMPS_TYPESPLIT (PROCNODE_STEPS(STEP(INODE)), 
     &              SLAVEF)
      SON_LEVEL2 = .FALSE.
      IOLDP2     = IOLDPS + HF - 1
      ICT11      = IOLDP2 + NFRONT
      NTOTFS = 0
      NELIM_SON_IN_PLACE = 0
      IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6) ) THEN
        J2    = PIMASTER(STEP(IFSON))
        LSTK  = IW(J2    +KEEP(IXSZ))
        NELIM = IW(J2 + 1+KEEP(IXSZ))
        IF ( ISON_IN_PLACE > 0 ) THEN
          IF (ISON_IN_PLACE.NE.IFSON) THEN
         write(6,*) MYID, ':',
     &   ' Internal error 1 in MUMPS_BUILD_SORT_INDEX ',
     &   ' in place node is not the first son a interior split node '
         CALL MUMPS_ABORT()
          ENDIF
          NELIM_SON_IN_PLACE = NELIM
        ENDIF
        NPIVS  = IW(J2 + 3+KEEP(IXSZ))
        IF (NPIVS.LT.0) NPIVS = 0
        NSLSON = IW(J2 + 5+KEEP(IXSZ))
        IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
        LEVEL1_SON    = NSLSON.EQ.0
        NCOLS  = NPIVS + LSTK
        NROWS  = NCOLS
        ITRANS = NROWS
        IF (NIV1) THEN
          write(6,*) MYID, ':',
     &    ' Internal error 2 in MUMPS_BUILD_SORT_INDEX ',
     &    ' interior split node of type 1 '
          CALL MUMPS_ABORT()
        ENDIF
        I= MUMPS_TYPENODE(PROCNODE_STEPS(STEP(IFSON)),SLAVEF)
        J= MUMPS_TYPESPLIT(PROCNODE_STEPS(STEP(IFSON)), 
     &              SLAVEF)
        IF (LEVEL1_SON.or.J.LT.4) THEN
           write(6,*) MYID, ':',
     &     ' Internal error 3 in MUMPS_BUILD_SORT_INDEX ',
     &     ' son', IFSON, 
     &     ' of interior split node', INODE, ' of type 1 ', 
     &     ' NSLSON =', NSLSON, ' TYPE_SON=', I, 'TYPESPLIT_SON=', J
           CALL MUMPS_ABORT()
        ENDIF
#if ! defined(NO_XXNBPR)
        IF (PIMASTER(STEP(IFSON)) .GT. IWPOSCB) THEN
          INBPROCFILS_SON = PIMASTER(STEP(IFSON))+XXNBPR
        ELSE
          INBPROCFILS_SON = PTRIST(STEP(IFSON))+XXNBPR
        ENDIF
#endif
        NBPROCFILS(STEP(IFSON)) = NSLSON
        NBPROCFILS(STEP(INODE)) = NSLSON
#if ! defined(NO_XXNBPR)
        IW(IOLDPS+XXNBPR)=NSLSON
        IW(INBPROCFILS_SON) = NSLSON
        CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE)),IW(IOLDPS+XXNBPR))
#endif
        SONROWS_PER_ROW(1:NFRONT-NASS1) = 1
        IF ( J2.GT. IWPOSCB ) THEN
          NROWS = IW(J2 + 2+KEEP(IXSZ))
          ITRANS = NPIVS + NROWS
        ENDIF
        HS = NSLSON + 6 + KEEP(IXSZ)
        J1 = J2 + HS + NROWS + NPIVS
        J2 = J1 + LSTK - 1
        J3 = J1 + NELIM - 1
        IF (NELIM.GT.0) THEN
         DO JJ=J1,J3
          NTOTFS = NTOTFS + 1
          JT1 = IW(JJ)
          IW(ICT11 + NTOTFS) = JT1
          IW(JJ) = NTOTFS
          IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
         ENDDO
        ENDIF
        DO JJ =J3+1, J3+NUMORG  
         NTOTFS = NTOTFS + 1
         JT1 = IW(JJ)
         ITLOC(JT1) = NTOTFS   
         IW(JJ) = NTOTFS
         IW(ICT11 + NTOTFS) = JT1
         IW(IOLDP2 + NTOTFS) = JT1
        ENDDO
        DO JJ =J3+NUMORG+1, J2
         NTOTFS = NTOTFS + 1
         JT1 = IW(JJ)
         ITLOC(JT1) = NTOTFS 
         IW(JJ) = NTOTFS
         IW(ICT11 + NTOTFS) = JT1
         IW(IOLDP2 + NTOTFS) = JT1
        ENDDO
        NFRONT_EFF = NTOTFS
        IBROT = INODE
        DO IORG = 1, NUMORG
          K1 = PTRAIW(IBROT) + 2
          JT1 = INTARR(K1)
          INTARR(K1) = ITLOC(JT1)
          IBROT = FILS(IBROT)
         K2 = K1 + INTARR(K1 - 2) - INTARR(K1 - 1)
         K1 = K1 + 1
         IF (K1 .LE. K2) THEN
          DO JJ = K1, K2
            J = INTARR(JJ)
            INTARR(JJ) = ITLOC(J)
          ENDDO
         ENDIF
        ENDDO
        K1 = IOLDPS+HF
        DO JJ=K1+NELIM,K1+NFRONT_EFF-1
          ITLOC(IW(JJ)) = 0
        ENDDO
        RETURN   
      ENDIF
       ALLOCATE(PTTRI(NUMSTK+1), stat=allocok)
       IF (allocok .GT. 0) THEN
        IFLAG = -13
        GOTO 800
       ENDIF
       ALLOCATE(PTLAST(NUMSTK+1), stat=allocok)
       IF (allocok .GT. 0) THEN
        IFLAG = -13
        GOTO 800
       ENDIF
      NFRONT_EFF = NASS1
      IF ( ISON_IN_PLACE > 0 ) THEN
        ISON  = ISON_IN_PLACE
        J2    = PIMASTER(STEP(ISON))
        LSTK   = IW(J2    +KEEP(IXSZ))
        NELIM = IW(J2 + 1+KEEP(IXSZ))
        NPIVS  = IW(J2 + 3+KEEP(IXSZ))
        IF (NPIVS.LT.0) NPIVS = 0
        NSLSON = IW(J2 + 5+KEEP(IXSZ))
        NCOLS  = NPIVS + LSTK
        NROWS  = NCOLS
        ITRANS = NROWS
        IF ( J2.GT. IWPOSCB ) THEN
          NROWS = IW(J2 + 2+KEEP(IXSZ))
          ITRANS = NPIVS + NROWS
        ENDIF
        HS = NSLSON + 6 + KEEP(IXSZ)
        J1 = J2 + HS + NROWS + NPIVS
        J2 = J1 + LSTK - 1
        J3 = J1 + NELIM - 1
        DO JJ = J1, J3
          NTOTFS = NTOTFS + 1
          JT1 = IW(JJ)
          IW(ICT11 + NTOTFS) = JT1
          ITLOC(JT1) = NTOTFS
          IW(JJ) = NTOTFS
          IW(IOLDP2 + NTOTFS) = IW(JJ - ITRANS)
        ENDDO
        NELIM_SON_IN_PLACE = NTOTFS
      ENDIF
      IF (.NOT. NIV1) SONROWS_PER_ROW(1:NFRONT-NASS1) = 0
      IN = INODE
      INEW = IOLDPS + HF +  NTOTFS
      INEW1 = NTOTFS + 1
      JDEBROW = PTRAIW(INODE)+3
      PTTRI(NUMSTK+1)  = JDEBROW
      PTLAST(NUMSTK+1) = JDEBROW + INTARR(JDEBROW-3) - 1
   50 J1 = PTRAIW(IN) + 2
      JT1 = INTARR(J1)
      INTARR(J1) = INEW1
      ITLOC(JT1) = INEW1
      IW(INEW)         = JT1
      IW(INEW+NFRONT)  = JT1
      INEW = INEW + 1
      INEW1 = INEW1 + 1
      IN = FILS(IN)
      IF (IN .GT. 0) GOTO 50
      NTOTFS = NTOTFS + NUMORG
      IF (NUMSTK .NE. 0) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          SON_IW => IW
          SON_IWPOSCB => IWPOSCB
          LSTK   = SON_IW(J2    +KEEP(IXSZ))
          NELIM  = SON_IW(J2 + 1+KEEP(IXSZ))
          NPIVS  = SON_IW(J2 + 3+KEEP(IXSZ))
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = SON_IW(J2 + 5+KEEP(IXSZ))
          IF( NSLSON.GT.0) SON_LEVEL2 = .TRUE.
          LEVEL1_SON    = NSLSON.EQ.0
          NCOLS  = NPIVS + LSTK
          NROWS  = NCOLS
          ITRANS = NROWS
#if ! defined(NO_XXNBPR)
          IF (PIMASTER(STEP(ISON)).GT.IWPOSCB) THEN
            INBPROCFILS_SON = PIMASTER(STEP(ISON))+XXNBPR
          ELSE
            INBPROCFILS_SON = PTRIST(STEP(ISON))+XXNBPR
          ENDIF
#endif
          IF (NIV1) THEN
           NBPROCFILS(STEP(ISON)) = NSLSON
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE)) + NSLSON
#if ! defined(NO_XXNBPR)
           IW(INBPROCFILS_SON) = NSLSON
           IW(IOLDPS+XXNBPR) = IW(IOLDPS+XXNBPR) + NSLSON
           CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE)),IW(IOLDPS+XXNBPR))
           CALL CHECK_EQUAL(NBPROCFILS(STEP(ISON)),IW(INBPROCFILS_SON))
#endif
          ELSE
           IF (LEVEL1_SON) THEN
            NBPROCFILS(STEP(ISON)) = 1
#if ! defined(NO_XXNBPR)
            IW(INBPROCFILS_SON) = 1
#endif
           ELSE
            NBPROCFILS(STEP(ISON)) = NSLSON
#if ! defined(NO_XXNBPR)
            IW(INBPROCFILS_SON) = NSLSON
#endif
           ENDIF
           NBPROCFILS(STEP(INODE)) = NBPROCFILS(STEP(INODE))+
     &                               NBPROCFILS(STEP(ISON))
#if ! defined(NO_XXNBPR)
           IW(IOLDPS+XXNBPR) = IW(IOLDPS+XXNBPR) + IW(INBPROCFILS_SON)
           CALL CHECK_EQUAL(NBPROCFILS(STEP(INODE)),IW(IOLDPS+XXNBPR))
#endif
          ENDIF
          IF (J2.GT.SON_IWPOSCB) THEN
           NROWS = SON_IW(J2 + 2+KEEP(IXSZ))
           ITRANS = NPIVS + NROWS
          ENDIF
          HS = NSLSON + 6 + KEEP(IXSZ)
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1 - KEEP(253)
          J3 = J1 + NELIM - 1
          IF (NELIM .NE. 0 .AND. ISON.NE.ISON_IN_PLACE) THEN
            DO JJ = J1, J3
              NTOTFS = NTOTFS + 1
              JT1 = SON_IW(JJ)
              IW(ICT11 + NTOTFS) = JT1
              ITLOC(JT1) = NTOTFS
              SON_IW(JJ) = NTOTFS
              IW(IOLDP2 + NTOTFS) = SON_IW(JJ - ITRANS)
            ENDDO
          ENDIF
          PTTRI(IELL)  = J2+1
          PTLAST(IELL) = J2
          J1 = J3 + 1
          IF (NASS1 .NE. NFRONT - KEEP(253)) THEN
            DO JJ = J1, J2
              J = SON_IW(JJ)
              IF (ITLOC(J) .EQ. 0) THEN 
                PTTRI(IELL) = JJ
                EXIT
              ENDIF
            ENDDO
          ELSE
            DO JJ = J1, J2
              SON_IW(JJ) = ITLOC(SON_IW(JJ))
            ENDDO
            DO JJ=J2+1, J2+KEEP(253)
              SON_IW(JJ)=NFRONT-KEEP(253)+JJ-J2
           ENDDO
          ENDIF
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IF (NFRONT-KEEP(253).EQ.NASS1) GOTO 500
 199  CONTINUE
      IF ( PTTRI( NUMSTK + 1 ) .LE. PTLAST( NUMSTK + 1 ) ) THEN
      IF ( ITLOC( INTARR( PTTRI( NUMSTK + 1 ) ) ) .NE. 0 ) THEN
       PTTRI( NUMSTK + 1 ) = PTTRI( NUMSTK + 1 ) + 1
       GOTO 199
      END IF
      END IF
      MIN_PERM = N + 1
      DO IELL = 1, NUMSTK 
        SON_IW => IW
        ILOC = PTTRI( IELL )
        IF ( ILOC .LE. PTLAST( IELL ) ) THEN 
         IF ( PERM( SON_IW( ILOC ) ) .LT. MIN_PERM ) THEN
           JMIN     = SON_IW( ILOC )
           MIN_PERM = PERM( JMIN )
         END IF
        END IF
      END DO
      IELL = NUMSTK + 1
      ILOC =  PTTRI( IELL )
      IF ( ILOC .LE. PTLAST( IELL ) ) THEN
        IF ( PERM( INTARR( ILOC ) ) .LT. MIN_PERM ) THEN
         JMIN        = INTARR( ILOC )
         MIN_PERM = PERM( JMIN )
        END IF
      END IF
      NEWEL = IOLDP2 + NASS1 + NFRONT
      DO WHILE ( MIN_PERM .NE. N + 1 )
          NEWEL  = NEWEL + 1
          NFRONT_EFF = NFRONT_EFF + 1
          IW( NEWEL ) = JMIN
          ITLOC( JMIN ) = NFRONT_EFF
          LAST_J_ASS = JMIN
          MIN_PERM = N + 1
          DO IELL = 1,  NUMSTK
            SON_IW => IW
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
              IF ( SON_IW( PTTRI( IELL ) ) .eq. LAST_J_ASS )
     &        PTTRI( IELL ) = PTTRI( IELL ) + 1
            ENDIF
            IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
             IF ( PERM(SON_IW( PTTRI( IELL )) ) .LT. MIN_PERM ) THEN
                JMIN        = SON_IW( PTTRI( IELL ) )
                MIN_PERM = PERM( JMIN )
             END IF
            END IF
          END DO
          IELL = NUMSTK + 1
 145      CONTINUE
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN
            IF ( INTARR( PTTRI( IELL ) ) .eq. LAST_J_ASS ) THEN
              PTTRI( IELL ) = PTTRI( IELL ) + 1 
              GOTO 145
            END IF
          END IF
          IF ( PTTRI( IELL ) .LE. PTLAST( IELL ) ) THEN 
            IF (PERM(INTARR( PTTRI(IELL) )) .LT. MIN_PERM) THEN
              JMIN        = INTARR( PTTRI(IELL) )
              MIN_PERM = PERM( JMIN )
            END IF
          END IF
      END DO
      NEWEL_SAVE  = NEWEL
      NEWEL1_SAVE = NFRONT_EFF
      IF (NEWEL1_SAVE.LT.NFRONT - KEEP(253)) THEN 
      IBROT = INODE
      DO IORG = 1, NUMORG
         J1    = PTRAIW(IBROT) + 2
         J2    = J1 + INTARR(J1 - 2) - INTARR(J1-1)
         IBROT = FILS( IBROT )
         IF ( IORG.EQ. 1) THEN
           IF ( KEEP(50).NE.0 ) CYCLE
           J1 = J1 + 1 + INTARR(J1-2)
         ELSE
           J1 = J1 + 1
         ENDIF
         DO JJ = J1, J2
           J     = INTARR( JJ )
           IF ( ITLOC( J ) .eq. 0 ) THEN
            NEWEL  = NEWEL + 1
            NFRONT_EFF = NFRONT_EFF + 1
            IW( NEWEL ) = J
            ITLOC( J ) = NFRONT_EFF
           END IF
         ENDDO
      ENDDO
       IF ( (TYPESPLIT.EQ.4).AND.
     &      (NFRONT_EFF.LT.NFRONT-KEEP(253)) ) THEN
         IBROT = INODE
         DO WHILE
     &      (
     &        ( MUMPS_TYPESPLIT 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.5 
     &        )
     &        .OR.
     &        ( MUMPS_TYPESPLIT 
     &           (PROCNODE_STEPS(STEP(DAD(STEP(IBROT)))),SLAVEF)
     &           .EQ.6  
     &        )
     &      )
          IBROT = DAD(STEP(IBROT))
          IN = IBROT
          DO WHILE (IN.GT.0.AND.NFRONT_EFF.LT.NFRONT-KEEP(253))
            J1    = PTRAIW(IN) + 2
            J2    = J1 + INTARR(J1 - 2) - INTARR(J1-1)
            IN = FILS( IN )
            DO JJ = J1, J2
              J     = INTARR( JJ )
              IF ( ITLOC( J ) .eq. 0 ) THEN
                NEWEL  = NEWEL + 1
                NFRONT_EFF = NFRONT_EFF + 1
                IW( NEWEL ) = J
                ITLOC( J ) = NFRONT_EFF
              END IF
            ENDDO
          ENDDO
          IF (NFRONT_EFF.EQ.NFRONT-KEEP(253)) EXIT
        ENDDO
       ENDIF
      ENDIF
      IF ( NEWEL1_SAVE .eq. NFRONT_EFF ) THEN
         DO JJ=NASS1+1, NFRONT_EFF
           IW( IOLDP2+JJ ) = IW( ICT11+JJ )
         ENDDO
      ELSE
        CALL MUMPS_SORT( N, PERM, 
     &           IW( NEWEL_SAVE + 1 ), NFRONT_EFF - NEWEL1_SAVE )
        CALL MUMPS_SORTED_MERGE( N, NASS1, PERM, ITLOC,
     &    IW( NEWEL_SAVE + 1), NFRONT_EFF - NEWEL1_SAVE,
     &    IW( ICT11  + NASS1 + 1 ), NEWEL1_SAVE - NASS1,
     &    IW( IOLDP2 + NASS1 + 1 ), NFRONT_EFF - NASS1 )
        DO JJ = NASS1+1, NFRONT_EFF
          IW(ICT11 + JJ) = IW(IOLDP2+JJ)
        ENDDO
      END IF
  500 CONTINUE
      IF ( KEEP(253).GT.0) THEN
        IP1 = IOLDPS +  HF + NFRONT_EFF  
        IP2 = IOLDPS + HF + NFRONT + NFRONT_EFF 
        DO I= 1, KEEP(253)
          IW(IP1+I-1) = N+I
          IW(IP2+I-1) = N+I
          ITLOC(N+I)  = NFRONT_EFF + I
        ENDDO
        NFRONT_EFF = NFRONT_EFF + KEEP(253)
      ENDIF
      IF (NFRONT.GT.NFRONT_EFF) THEN
        IP1 = IOLDPS + NFRONT + HF
        IP2 = IOLDPS + NFRONT_EFF + HF
        DO I=1, NFRONT_EFF
          IW(IP2+I-1)=IW(IP1+I-1)
        ENDDO
      ELSE IF (NFRONT .LT. NFRONT_EFF) THEN
        WRITE(*,*) "Internal error in MUMPS_ELT_BUILD_SORT",
     &             NFRONT, NFRONT_EFF
        CALL MUMPS_ABORT()
      ENDIF
      IF ( NUMSTK .NE. 0  
     &    .AND. (NFRONT-KEEP(253).GT.NASS1)
     &  ) THEN
        ISON = IFSON
        DO IELL = 1, NUMSTK
          J2 = PIMASTER(STEP(ISON))
          SON_IW => IW
          SON_IWPOSCB => IWPOSCB
          LSTK = SON_IW(J2+KEEP(IXSZ))
          NELIM = SON_IW(J2 + 1 +KEEP(IXSZ))
          NPIVS = SON_IW(J2 + 3 +KEEP(IXSZ))
          IF (NPIVS.LT.0) NPIVS = 0
          NSLSON = SON_IW(J2 + 5 +KEEP(IXSZ))
          LEVEL1_SON = (NSLSON .EQ. 0)
          NCOLS = NPIVS + LSTK
          NROWS = NCOLS
          IF (J2.GT.SON_IWPOSCB) THEN
           NROWS = SON_IW(J2 + 2+KEEP(IXSZ))
          ENDIF
          HS = NSLSON + 6 +KEEP(IXSZ)
          J1 = J2 + HS + NROWS + NPIVS
          J2 = J1 + LSTK - 1
          J3 = J1 + NELIM - 1
          J1 = J3 + 1
          IF (NFRONT-KEEP(253).GT.NASS1) THEN
            DO JJ = J1, J2
              J = SON_IW(JJ)
              SON_IW(JJ) = ITLOC(J)
              IF (NIV1 .AND. NSLSON.EQ.0) THEN
              ELSE
                IF (SON_IW(JJ) .LE. NASS1 .OR. NIV1) THEN
                ELSE
                  SONROWS_PER_ROW(SON_IW(JJ)-NASS1) =
     &                         SONROWS_PER_ROW(SON_IW(JJ)-NASS1) + 1
                ENDIF
              ENDIF
            ENDDO
          ELSE
              IF (.not. NIV1) THEN
                WRITE(*,*) "Internal error 1 in CMUMPS_BUILD_SORT_INDEX"
                CALL MUMPS_ABORT() 
              ENDIF
              IF (.not.LEVEL1_SON) THEN
              ENDIF
          ENDIF
          ISON = FRERE_STEPS(STEP(ISON))
        ENDDO
      ENDIF
      IBROT = INODE
      DO IORG = 1, NUMORG
        J1 = PTRAIW(IBROT) + 2
        IBROT = FILS(IBROT)
        J2 = J1 + INTARR(J1 - 2) - INTARR(J1 - 1)
        J1 = J1 + 1
        DO JJ = J1, J2
          J = INTARR(JJ)
            INTARR(JJ) = ITLOC(J)
        ENDDO
      ENDDO
        K1 = IOLDPS + HF
        K2 = K1 + NFRONT_EFF -1
        IF (KEEP(50).EQ.0) K2 = K2 + NELIM_SON_IN_PLACE
        DO K = K1, K2
          I = IW(K)
          ITLOC(I) = 0
        ENDDO
        IF (KEEP(50).EQ.0) THEN
          K1 = IOLDPS+HF+NFRONT_EFF+NELIM_SON_IN_PLACE+NUMORG
          K2 = K1 + NASS -NELIM_SON_IN_PLACE - 1
          DO K = K1, K2
            I = IW(K)
            ITLOC(I) = 0
          ENDDO
        ENDIF
  800 CONTINUE
      IF (allocated(PTTRI)) DEALLOCATE(PTTRI)
      IF (allocated(PTLAST)) DEALLOCATE(PTLAST)
      RETURN
      END SUBROUTINE MUMPS_BUILD_SORT_INDEX
      END MODULE MUMPS_BUILD_SORT_INDEX_M
      SUBROUTINE MUMPS_SORT( N, PERM, IW, LIW )
      IMPLICIT NONE
      INTEGER N, LIW
      INTEGER PERM( N ), IW( LIW )
      INTEGER I, SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, LIW - 1
          IF ( PERM( IW( I ) ) .GT. PERM( IW( I + 1 ) ) ) THEN
            DONE = .FALSE.
            SWAP  = IW( I + 1 )
            IW( I + 1 ) = IW( I )
            IW( I ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_SORT
      SUBROUTINE MUMPS_SORTED_MERGE( N, NASS1, PERM, ITLOC,
     &                             SMALL, LSMALL,
     &                             LARGE, LLARGE,
     &                             MERGE, LMERGE )
      IMPLICIT NONE
      INTEGER N, NASS1, LSMALL, LLARGE, LMERGE
      INTEGER PERM( N ), ITLOC( N ) 
      INTEGER SMALL(LSMALL), LARGE(LLARGE), MERGE(LMERGE)
      INTEGER PSMALL, PLARGE, PMERGE, VSMALL, VLARGE, VMERGE
      PSMALL = 1
      PLARGE = 1
      PMERGE = 1
      DO WHILE ( PSMALL .LE. LSMALL .or. PLARGE.LE. LLARGE )
        IF ( PSMALL .GT. LSMALL ) THEN
          VMERGE = LARGE( PLARGE )
          PLARGE = PLARGE + 1
        ELSE IF ( PLARGE .GT. LLARGE ) THEN
          VMERGE = SMALL( PSMALL )
          PSMALL = PSMALL + 1
        ELSE
          VSMALL = SMALL( PSMALL )
          VLARGE = LARGE( PLARGE )
          IF ( PERM( VSMALL ) .LT. PERM( VLARGE ) ) THEN
            VMERGE = VSMALL
            PSMALL   = PSMALL + 1
          ELSE
            VMERGE = VLARGE
            PLARGE   = PLARGE + 1
          END IF
        END IF
        MERGE( PMERGE ) = VMERGE
        ITLOC( VMERGE ) = PMERGE + NASS1
        PMERGE = PMERGE + 1
      END DO
      PMERGE = PMERGE - 1
      RETURN
      END SUBROUTINE MUMPS_SORTED_MERGE
