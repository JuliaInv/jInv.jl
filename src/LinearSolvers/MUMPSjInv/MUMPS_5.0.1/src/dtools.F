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
      SUBROUTINE DMUMPS_COMPRESS_LU(SIZE_INPLACE,
     &MYID,N,IOLDPS,TYPE,IW, LIW, A, LA,
     &POSFAC, LRLU, LRLUS, IWPOS, PTRAST, PTRFAC, STEP, KEEP,KEEP8,
     &SSARBR,INODE,IERR)
      USE DMUMPS_LOAD
      USE DMUMPS_OOC
      IMPLICIT NONE
      INTEGER MYID
      INTEGER IOLDPS, TYPE, LIW, N, KEEP(500)
      INTEGER(8) :: SIZE_INPLACE, LA, POSFAC, LRLU, LRLUS
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER(8) KEEP8(150)
      INTEGER IW( LIW )
      DOUBLE PRECISION A( LA )
      INTEGER IWPOS, LDLT
      INTEGER STEP( N )
      INTEGER (8) :: PTRFAC(KEEP(28))
      LOGICAL SSARBR
      INTEGER IOLDSHIFT, IPSSHIFT
      INCLUDE 'mumps_headers.h'
      INTEGER LCONT, NELIM, NROW, NPIV, INTSIZ
      INTEGER NFRONT, NSLAVES
      INTEGER IPS, IPSIZE
      INTEGER(8) :: SIZELU, SIZECB, IAPOS, I
      LOGICAL MOVEPTRAST
      INTEGER INODE
      INTEGER IERR
      IERR=0
      LDLT = KEEP(50)
      IOLDSHIFT = IOLDPS + KEEP(IXSZ)
      IF ( IW( IOLDSHIFT ) < 0 ) THEN
        write(*,*) ' ERROR 1 compressLU:Should not point to a band.'
        CALL MUMPS_ABORT()
      ELSE IF ( IW( IOLDSHIFT + 2 ) < 0 ) THEN
        write(*,*) ' ERROR 2 compressLU:Stack not performed yet',
     &  IW(IOLDSHIFT + 2)
        CALL MUMPS_ABORT()
      ENDIF
      LCONT  = IW( IOLDSHIFT )
      NELIM  = IW( IOLDSHIFT + 1 )
      NROW   = IW( IOLDSHIFT + 2 )
      NPIV   = IW( IOLDSHIFT + 3 )
      IAPOS  = PTRFAC(IW( IOLDSHIFT + 4 ))
      NSLAVES= IW( IOLDSHIFT + 5 )
      NFRONT = LCONT + NPIV
      INTSIZ = IW(IOLDPS+XXI)
      IF ( (NSLAVES > 0  .AND. TYPE .NE. 2) .OR. 
     &   (NSLAVES .eq. 0 .AND. TYPE .EQ. 2 ) ) THEN
          WRITE(*,*) ' ERROR 3 compressLU: problem with level of inode'
          CALL MUMPS_ABORT()
      END IF
      IF (LDLT.EQ.0) THEN
        SIZELU = int(LCONT + NROW, 8) * int(NPIV,8)
      ELSE
        SIZELU =   int(NROW,8) * int(NPIV,8)
      ENDIF
      IF ( TYPE .EQ. 2 ) THEN
        IF (LDLT.EQ.0) THEN
          SIZECB = int(NELIM,8) * int(LCONT,8)
        ELSE
          IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
            SIZECB = int(NELIM+1,8) * int(NELIM + NPIV,8)
          ELSE
            SIZECB = int(NELIM,8) * int(NELIM + NPIV,8)
          ENDIF
        ENDIF
      ELSE
        IF (LDLT.EQ.0) THEN
         SIZECB = int(LCONT,8) * int(LCONT,8)
        ELSE
         SIZECB = int(NROW,8) * int(LCONT,8)
        ENDIF
      END IF
      CALL MUMPS_SUBTRI8TOARRAY( IW(IOLDPS+XXR), SIZECB )
      IF ((SIZECB.EQ.0_8).AND.(KEEP(201).EQ.0)) THEN
         GOTO 500
      ENDIF
      IF (KEEP(201).EQ.2) THEN
         KEEP8(31)=KEEP8(31)+SIZELU
         CALL DMUMPS_NEW_FACTOR(INODE,PTRFAC,KEEP,KEEP8,
     &        A,LA,SIZELU, IERR)
         IF(IERR.LT.0)THEN
            WRITE(*,*)MYID,': Internal error in DMUMPS_NEW_FACTOR'
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
      IF ( IOLDPS + INTSIZ .NE. IWPOS ) THEN
         IPS = IOLDPS + INTSIZ
         MOVEPTRAST = .FALSE.
         DO WHILE ( IPS .NE. IWPOS )
           IPSIZE = IW(IPS+XXI)
           IPSSHIFT = IPS + KEEP(IXSZ)
           IF ( IW( IPSSHIFT + 2 ) < 0 ) THEN
             NFRONT = IW( IPSSHIFT )
             IF(KEEP(201).EQ.0)THEN
               PTRFAC(IW( IPSSHIFT + 4 )) = 
     &                      PTRFAC(IW( IPSSHIFT + 4 )) - SIZECB
             ELSE
               PTRFAC(IW(IPSSHIFT+4))=PTRFAC(IW(IPSSHIFT+4)) -
     &               SIZECB - SIZELU
             ENDIF
             MOVEPTRAST = .TRUE.
             IF(KEEP(201).EQ.0)THEN
               PTRAST(IW(IPSSHIFT+4))=PTRAST(IW(IPSSHIFT+4))-SIZECB
             ELSE
               PTRAST(IW(IPSSHIFT+4))=PTRAST(IW(IPSSHIFT+4))-SIZECB
     &               - SIZELU
             ENDIF
           ELSE IF ( IW( IPSSHIFT ) < 0 ) THEN
             IF(KEEP(201).EQ.0)THEN
               PTRFAC(IW(IPSSHIFT+3)) = PTRFAC(IW(IPSSHIFT+3))-SIZECB
             ELSE
               PTRFAC(IW(IPSSHIFT+3)) = PTRFAC(IW(IPSSHIFT+3))
     &                                  -SIZECB-SIZELU
             ENDIF
           ELSE
             NFRONT = IW( IPSSHIFT ) + IW( IPSSHIFT + 3 )
             IF(KEEP(201).EQ.0)THEN
                PTRFAC(IW( IPSSHIFT + 4 )) = 
     &                    PTRFAC(IW( IPSSHIFT + 4 )) - SIZECB
             ELSE
                PTRFAC(IW( IPSSHIFT + 4 )) = 
     &               PTRFAC(IW( IPSSHIFT + 4 )) - SIZECB
     &               - SIZELU
             ENDIF
           END IF
           IPS = IPS + IPSIZE
         END DO
         IF ((SIZECB .NE. 0_8).OR.(KEEP(201).NE.0)) THEN
            IF (KEEP(201).NE.0) THEN
               DO I=IAPOS, POSFAC - SIZECB - SIZELU - 1_8
                  A( I ) = A( I + SIZECB + SIZELU)
               END DO
            ELSE
               DO I=IAPOS + SIZELU, POSFAC - SIZECB - 1_8
                  A( I ) = A( I + SIZECB )
               END DO
            ENDIF
         END IF
      ENDIF
      IF (KEEP(201).NE.0) THEN
        POSFAC = POSFAC  - (SIZECB+SIZELU)
        LRLU   = LRLU    + (SIZECB+SIZELU)
        LRLUS  = LRLUS   + (SIZECB+SIZELU) - SIZE_INPLACE
      ELSE
        POSFAC = POSFAC - SIZECB
        LRLU   = LRLU   + SIZECB
        LRLUS  = LRLUS  + SIZECB - SIZE_INPLACE
      ENDIF
 500  CONTINUE
      CALL DMUMPS_LOAD_MEM_UPDATE(SSARBR,.FALSE.,
     &     LA-LRLUS,SIZELU,-SIZECB+SIZE_INPLACE,KEEP,KEEP8,LRLUS)
      RETURN
      END SUBROUTINE DMUMPS_COMPRESS_LU
      SUBROUTINE DMUMPS_STACK_BAND( N, ISON, 
     &    PTRIST, PTRAST, PTLUST_S, PTRFAC, IW, LIW, A, LA, 
     &    LRLU, LRLUS, IWPOS, IWPOSCB, POSFAC, COMP, 
     &    IPTRLU, OPELIW, STEP, PIMASTER, PAMASTER,
     &    IFLAG, IERROR, SLAVEF, MYID, COMM,
     &    KEEP, KEEP8, DKEEP, TYPE_SON
     &     )
      USE DMUMPS_OOC
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER(8) :: LA, LRLU, LRLUS, POSFAC, IPTRLU
      INTEGER N, ISON, LIW, IWPOS, IWPOSCB,
     &        COMP, IFLAG, IERROR, SLAVEF, MYID, COMM,
     &        TYPE_SON
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION DKEEP(130)
      INTEGER(8) :: PTRAST(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), STEP(N), 
     & PIMASTER(KEEP(28)), IW(LIW)
      INTEGER PTLUST_S(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      DOUBLE PRECISION OPELIW
      DOUBLE PRECISION FLOP1, FLOP1_EFFECTIVE
      DOUBLE PRECISION A( LA )
      INTEGER(8) :: LREQA, POSA, POSALOC, OLDPOS, JJ
      INTEGER  NFRONT, NCOL_L, NROW_L, LREQI, NSLAVES_L,
     &         POSI, I, IROW_L, ICOL_L, LDA_BAND, NASS
      LOGICAL NONEED_TO_COPY_FACTORS
      INTEGER(8) :: LAFAC, LREQA_HEADER
      INTEGER LIWFAC, STRAT, TYPEFile, NextPivDummy,
     &        IOLDPS_CB
      LOGICAL LAST_CALL
      TYPE(IO_BLOCK) :: MonBloc 
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0d0)
      FLOP1 = ZERO
      NCOL_L = IW( PTRIST(STEP( ISON )) + 3 + KEEP(IXSZ) )
      NROW_L = IW( PTRIST(STEP( ISON )) + 2 + KEEP(IXSZ) )
      NSLAVES_L = IW( PTRIST(STEP( ISON )) + 5 + KEEP(IXSZ) )
      LDA_BAND = NCOL_L + IW( PTRIST(STEP( ISON )) + KEEP(IXSZ) )
      IF  ( KEEP(50) .eq. 0 ) THEN
        NFRONT = LDA_BAND
      ELSE
        NFRONT = IW( PTRIST(STEP( ISON )) + 7 + KEEP(IXSZ) )
      END IF
      IF (KEEP(201).EQ.1) THEN 
          IOLDPS_CB = PTRIST(STEP( ISON ))
          CALL MUMPS_GETI8(LAFAC, IW(IOLDPS_CB+XXR))
          LIWFAC    = IW(IOLDPS_CB+XXI)
          TYPEFile  = TYPEF_L
          NextPivDummy      = -8888 
          MonBloc%INODE    = ISON
          MonBloc%MASTER   = .FALSE.   
          MonBloc%Typenode =  2        
          MonBloc%NROW     = NROW_L
          MonBloc%NCOL     = LDA_BAND
          MonBloc%NFS      = IW(IOLDPS_CB+1+KEEP(IXSZ))
          MonBloc%LastPiv  = NCOL_L    
          MonBloc%LastPanelWritten_L=-9999 
          MonBloc%LastPanelWritten_U=-9999 
          NULLIFY(MonBloc%INDICES)
          STRAT        = STRAT_WRITE_MAX
          LAST_CALL    = .TRUE.
          MonBloc%Last = .TRUE.
          CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEFile, 
     &           A(PTRAST(STEP(ISON))), LAFAC, MonBloc,
     &           NextPivDummy, NextPivDummy,
     &           IW(IOLDPS_CB), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG,LAST_CALL )
          IF ((NCOL_L.EQ.0).OR.(NROW_L.EQ.0)) THEN 
          ENDIF
      ENDIF  
      NONEED_TO_COPY_FACTORS = (KEEP(201).EQ.1) .OR. (KEEP(201).EQ.-1)
      IF ((NCOL_L.EQ.0).OR.(NROW_L.EQ.0)) THEN 
        GOTO 80
      ENDIF
      LREQI   = 4 + NCOL_L + NROW_L + KEEP(IXSZ)
      LREQA_HEADER =  int(NCOL_L,8) * int(NROW_L,8)
      IF (NONEED_TO_COPY_FACTORS) THEN 
        LREQA = 0_8
      ELSE
        LREQA   = LREQA_HEADER
      ENDIF
      IF ( LRLU .LT. LREQA .OR.
     &  IWPOS + LREQI - 1 .GT. IWPOSCB ) THEN
        IF ( LRLUS .LT. LREQA ) THEN
          IFLAG  = -9
          CALL MUMPS_SET_IERROR(LREQA - LRLUS, IERROR)
          GO TO 700
        END IF
        CALL DMUMPS_COMPRE_NEW( N,KEEP(28), IW, LIW, A, LA,
     &        LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, KEEP(216),LRLUS,
     &        KEEP(IXSZ), COMP, DKEEP(97), MYID )
        IF ( LRLU .NE. LRLUS ) THEN
               WRITE(*,*) 'PB compress DMUMPS_STACK_BAND:LRLU,LRLUS=',
     &         LRLU, LRLUS
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
      IF (.NOT. NONEED_TO_COPY_FACTORS) THEN
        POSA = POSFAC
        POSFAC = POSFAC + LREQA
        LRLU = LRLU - LREQA
        LRLUS = LRLUS - LREQA
        KEEP8(67) = min(LRLUS, KEEP8(67))
        IF(KEEP(201).NE.2)THEN
           CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &          LA-LRLUS,LREQA,LREQA,KEEP,KEEP8,LRLUS)
        ELSE
           CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &          LA-LRLUS,0_8,LREQA,KEEP,KEEP8,LRLUS)
        ENDIF
      ENDIF
      POSI = IWPOS
      IWPOS = IWPOS + LREQI
      PTLUST_S(STEP( ISON )) = POSI
      IW(POSI+XXI)=LREQI
      CALL MUMPS_STOREI8(LREQA, IW(POSI+XXR))
      CALL MUMPS_STOREI8(LREQA_HEADER, IW(POSI+XXR))
      IW(POSI+XXS)=-9999
      IW(POSI+XXS+1:POSI+KEEP(IXSZ)-1)=-99999
      POSI=POSI+KEEP(IXSZ)
      IW( POSI     ) = - NCOL_L
      IW( POSI + 1 ) =   NROW_L
      IW( POSI + 2 ) =   NFRONT - NCOL_L
      IW( POSI + 3 ) =   STEP(ISON)
      IF (.NOT. NONEED_TO_COPY_FACTORS) THEN
        PTRFAC(STEP(ISON)) = POSA
      ELSE
        PTRFAC(STEP(ISON)) = -77777_8
      ENDIF
      IROW_L = PTRIST(STEP(ISON)) + 6 + NSLAVES_L + KEEP(IXSZ)
      ICOL_L = PTRIST(STEP(ISON)) + 6 + NROW_L + NSLAVES_L + KEEP(IXSZ)
      DO I = 1, NROW_L
        IW( POSI+3+I ) = IW( IROW_L+I-1 )
      ENDDO
      DO I = 1, NCOL_L
        IW( POSI+NROW_L+3+I) = IW( ICOL_L+I-1 )
      ENDDO
      IF (.NOT.NONEED_TO_COPY_FACTORS) THEN
        POSALOC = POSA
        DO I = 1, NROW_L
          OLDPOS =  PTRAST( STEP(ISON)) + int(I-1,8)*int(LDA_BAND,8)
          DO JJ = 0_8, int(NCOL_L-1,8)
            A( POSALOC+JJ ) = A( OLDPOS+JJ )
          ENDDO
          POSALOC = POSALOC + int(NCOL_L,8)
        END DO
      ENDIF
      IF (KEEP(201).EQ.2) THEN
       KEEP8(31)=KEEP8(31)+LREQA
      ENDIF
      KEEP8(10) = KEEP8(10) + int(NCOL_L,8) * int(NROW_L,8)
      IF (KEEP(201).EQ.2) THEN 
        CALL DMUMPS_NEW_FACTOR(ISON,PTRFAC,KEEP,KEEP8,A,LA,LREQA,IFLAG)
        IF(IFLAG.LT.0)THEN
          WRITE(*,*)MYID,': Internal error in DMUMPS_NEW_FACTOR'
          IERROR=0
          GOTO 700
        ENDIF
        POSFAC = POSFAC - LREQA
        LRLU = LRLU + LREQA
        LRLUS = LRLUS + LREQA
        CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &            LA-LRLUS,LREQA,0_8,KEEP,KEEP8,LRLUS)
      ENDIF
  80  CONTINUE
      IF (TYPE_SON == 1) THEN
         GOTO 90
      ENDIF
      IF ( KEEP(50) .eq. 0 ) THEN
         FLOP1 = dble( NCOL_L * NROW_L) +
     &     dble(NROW_L*NCOL_L)*dble(2*NFRONT-NCOL_L-1)
      ELSE
         FLOP1 = dble( NCOL_L ) * dble( NROW_L )
     &         * dble( 2 * LDA_BAND - NROW_L - NCOL_L + 1)
      END IF
      OPELIW = OPELIW + FLOP1
      FLOP1_EFFECTIVE = FLOP1
      NASS = IW( PTRIST(STEP( ISON )) + 4 + KEEP(IXSZ) )
      IF ( NCOL_L .NE. NASS ) THEN
        IF ( KEEP(50).eq.0 ) THEN
           FLOP1 = dble( NASS * NROW_L) +
     &     dble(NROW_L*NASS)*dble(2*NFRONT-NASS-1)
        ELSE
           FLOP1 = dble( NASS ) * dble( NROW_L ) *
     &     dble( 2 * LDA_BAND - NROW_L - NASS + 1)
        END IF
      END IF
      CALL DMUMPS_LOAD_UPDATE(1,.FALSE.,FLOP1_EFFECTIVE-FLOP1,
     &                        KEEP,KEEP8)
      CALL DMUMPS_LOAD_UPDATE(2,.FALSE.,-FLOP1,KEEP,KEEP8)
 90   CONTINUE
      RETURN
 700  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
      RETURN
      END SUBROUTINE DMUMPS_STACK_BAND
      SUBROUTINE DMUMPS_FREE_BAND( N, ISON, 
     &    PTRIST, PTRAST, IW, LIW, A, LA, 
     &    LRLU, LRLUS, IWPOSCB,
     &    IPTRLU, STEP, MYID, KEEP, TYPE_SON
     &     )
      IMPLICIT NONE
      include 'mumps_headers.h'
      INTEGER(8) :: LRLU, LRLUS, IPTRLU, LA
      INTEGER ISON, MYID, N, IWPOSCB, TYPE_SON
      INTEGER KEEP(500), STEP(N)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: PTRAST(KEEP(28))
      INTEGER PTRIST(KEEP(28))
      INTEGER LIW
      INTEGER IW(LIW)
      DOUBLE PRECISION A(LA)
      INTEGER ISTCHK
      ISTCHK = PTRIST(STEP(ISON))
      CALL DMUMPS_FREE_BLOCK_CB(.FALSE.,MYID, N, ISTCHK,
     &     PTRAST(STEP(ISON)),
     &     IW, LIW, LRLU, LRLUS, IPTRLU,
     &     IWPOSCB, LA, KEEP,KEEP8, .FALSE.
     &     )
      PTRIST(STEP( ISON )) = -9999888
      PTRAST(STEP( ISON )) = -9999888_8
      RETURN
      END SUBROUTINE DMUMPS_FREE_BAND
      SUBROUTINE DMUMPS_MAX_MEM( KEEP,KEEP8,
     &           MYID, N, NELT, NA, LNA, NZ, NA_ELT, NSLAVES,
     &           MEMORY_MBYTES, EFF, OOC_STRAT, PERLU_ON,
     &           MEMORY_BYTES )
      IMPLICIT NONE
      LOGICAL,   INTENT(IN)  :: EFF, PERLU_ON
      INTEGER,   INTENT(IN)  :: OOC_STRAT
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER MYID, N, NELT, NSLAVES, LNA, NZ, NA_ELT
      INTEGER(8), INTENT(IN) :: NA(LNA)
      INTEGER(8), INTENT(OUT) :: MEMORY_BYTES
      INTEGER,   INTENT(OUT) :: MEMORY_MBYTES
      INTEGER  :: MUMPS_GET_POOL_LENGTH
      EXTERNAL :: MUMPS_GET_POOL_LENGTH
      LOGICAL    :: I_AM_SLAVE, I_AM_MASTER
      INTEGER    :: PERLU, NBRECORDS
      INTEGER(8) :: NB_REAL, MAXS_MIN
      INTEGER(8) :: TEMP, NB_BYTES, NB_INT
      INTEGER    :: DMUMPS_LBUF_INT
      INTEGER(8) :: DMUMPS_LBUFR_BYTES8, DMUMPS_LBUF8
      INTEGER    :: NBUFS
      INTEGER(8) :: TEMPI
      INTEGER(8) :: TEMPR
      INTEGER    :: MIN_PERLU
      INTEGER(8) :: BUF_OOC, BUF_OOC_PANEL, BUF_OOC_NOPANEL
      INTEGER(8) :: OOC_NB_FILE_TYPE
      INTEGER(8) :: NSTEPS8, N8, NELT8
      INTEGER(8) :: I8OVERI
      I8OVERI   = int(KEEP(10),8)
      PERLU     = KEEP(12)
      NSTEPS8   = int(KEEP(28),8)
      N8        = int(N,8)
      NELT8     = int(NELT,8)
      IF (.NOT.PERLU_ON) PERLU = 0
      I_AM_MASTER = ( MYID .eq. 0 )
      I_AM_SLAVE  = ( KEEP(46).eq. 1 .or. MYID .ne. 0 )
      TEMP    = 0_8
      NB_REAL = 0_8
      NB_BYTES = 0_8
      NB_INT  = 0_8
      IF (KEEP(235) .NE. 0 .OR. KEEP(237) .NE. 0) THEN
         NB_INT  = NB_INT + NSTEPS8
      ENDIF
      NB_INT = NB_INT + 5_8 * NSTEPS8
      NB_INT = NB_INT + NSTEPS8 + int(KEEP(56),8)*int(NSLAVES+2,8)
      NB_INT = NB_INT + 3_8 * N8
      IF (KEEP(23).ne.0 .and. I_AM_MASTER) NB_INT=NB_INT + N8
      IF (KEEP(55).eq.0) THEN
        NB_INT = NB_INT + 2_8 * N8
      ELSE
        NB_INT = NB_INT + 2_8 * ( NELT8 + 1_8 )
      ENDIF
      IF (KEEP(55) .ne. 0 ) THEN
        NB_INT = NB_INT + N8 + 1_8 + NELT8
      END IF
      NB_INT = NB_INT + int(LNA,8)
      IF ( OOC_STRAT .GT. 0 .OR. OOC_STRAT .EQ. -1 ) THEN
        MAXS_MIN = KEEP8(14)
      ELSE
        MAXS_MIN = KEEP8(12)
      ENDIF
      IF ( .NOT. EFF ) THEN
        IF ( KEEP8(24).EQ.0_8 ) THEN
         NB_REAL = NB_REAL + MAXS_MIN +
     &             int(PERLU,8)*(MAXS_MIN / 100_8 + 1_8 )
        ENDIF
      ELSE
        NB_REAL = NB_REAL + KEEP8(67)
      ENDIF
      IF ( OOC_STRAT .GT. 0 .AND. I_AM_SLAVE ) THEN
        BUF_OOC_NOPANEL = 2_8 * KEEP8(119)
        IF (KEEP(50).EQ.0)THEN
          BUF_OOC_PANEL = 8_8 * int(KEEP(226),8)
        ELSE
          BUF_OOC_PANEL = 4_8 * int(KEEP(226),8)
        ENDIF
        IF (OOC_STRAT .EQ. 2) THEN
          BUF_OOC = BUF_OOC_NOPANEL
        ELSE
          BUF_OOC = BUF_OOC_PANEL
        ENDIF
        NB_REAL = NB_REAL + min(BUF_OOC + int(max(PERLU,0),8) *
     &          (BUF_OOC/100_8+1_8),12000000_8)
        IF (OOC_STRAT .EQ. 2) THEN
          OOC_NB_FILE_TYPE = 1_8
        ELSE
          IF (KEEP(50).EQ.0) THEN
            OOC_NB_FILE_TYPE = 2_8
          ELSE
            OOC_NB_FILE_TYPE = 1_8
          ENDIF
        ENDIF
        NB_INT = NB_INT + OOC_NB_FILE_TYPE * NSTEPS8 * I8OVERI
        NB_INT = NB_INT + OOC_NB_FILE_TYPE * NSTEPS8 * I8OVERI
        NB_INT = NB_INT + OOC_NB_FILE_TYPE * NSTEPS8
      ENDIF
      NB_REAL = NB_REAL + int(KEEP(13),8)
      IF (KEEP(252).EQ.1 .AND. .NOT. I_AM_MASTER) THEN
        NB_REAL = NB_REAL + N8
      ENDIF
      IF ( .not. ( I_AM_SLAVE .and. I_AM_MASTER .and. KEEP(52) .eq. 0
     &         .and. KEEP(55) .ne. 0 ) ) THEN
        NB_INT  = NB_INT  + int(KEEP(14),8)
      END IF
      IF ( I_AM_SLAVE .and. KEEP(38) .ne. 0 ) THEN
        NB_INT = NB_INT + 2_8 * N8
      END IF
      TEMPI= 0_8
      TEMPR = 0_8
      NBRECORDS = KEEP(39)
      IF (KEEP(55).eq.0) THEN
        NBRECORDS = min(KEEP(39), NZ)
      ELSE
        NBRECORDS = min(KEEP(39), NA_ELT)
      ENDIF
      IF ( KEEP(54) .eq. 0 ) THEN
        IF ( I_AM_MASTER ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            NBUFS = NSLAVES 
          ELSE
            NBUFS = NSLAVES - 1
            IF (KEEP(55) .eq. 0 )
     &      TEMPI = TEMPI + 2_8 * N8
          END IF
          TEMPI = TEMPI + 2_8 * int(NBRECORDS,8) * int(NBUFS,8)
          TEMPR = TEMPR + int(NBRECORDS,8) * int(NBUFS,8)
        ELSE
          IF ( KEEP(55) .eq. 0 )THEN
            TEMPI = TEMPI + 2_8 * int(NBRECORDS,8)
            TEMPR = TEMPR + int(NBRECORDS,8)
          END IF
        END IF
      ELSE
        IF ( I_AM_SLAVE ) THEN
          TEMPI = TEMPI + int(1+4*NSLAVES,8) * int(NBRECORDS,8)
          TEMPR = TEMPR + int(1+2*NSLAVES,8) * int(NBRECORDS,8)
        END IF
      END IF
      TEMP = max( NB_BYTES + (NB_INT + TEMPI) * int(KEEP(34),8)
     &           + (NB_REAL+TEMPR) * int(KEEP(35),8)
     &            , TEMP )
      IF ( I_AM_SLAVE ) THEN
        DMUMPS_LBUFR_BYTES8 = int(KEEP(44),8) * int(KEEP(35),8)
        DMUMPS_LBUFR_BYTES8 = max( DMUMPS_LBUFR_BYTES8,
     &                      100000_8 )
        IF (KEEP(48).EQ.5) THEN
          MIN_PERLU=2
        ELSE
          MIN_PERLU=0
        ENDIF
        DMUMPS_LBUFR_BYTES8 = DMUMPS_LBUFR_BYTES8
     &        + int( 2.0D0 * dble(max(PERLU,MIN_PERLU))*
     &        dble(DMUMPS_LBUFR_BYTES8)/100D0)
        DMUMPS_LBUFR_BYTES8 = min(DMUMPS_LBUFR_BYTES8,
     &                            int(huge (KEEP(43))-100,8))
        NB_BYTES = NB_BYTES + DMUMPS_LBUFR_BYTES8
        DMUMPS_LBUF8 = int( dble(KEEP(213)) / 100.0D0
     &                     * dble(KEEP( 43 ) * KEEP( 35 )), 8 )
        DMUMPS_LBUF8 = max( DMUMPS_LBUF8, 100000_8 )
        DMUMPS_LBUF8 = DMUMPS_LBUF8
     &                 + int( 2.0D0 * dble(max(PERLU,MIN_PERLU))*
     &                   dble(DMUMPS_LBUF8)/100D0, 8)
        DMUMPS_LBUF8 = min(DMUMPS_LBUF8, int(huge (KEEP(43)-100),8))
        DMUMPS_LBUF8 = max(DMUMPS_LBUF8, DMUMPS_LBUFR_BYTES8+
     &                 3_8*int(KEEP(34),8))
        NB_BYTES = NB_BYTES + DMUMPS_LBUF8
        DMUMPS_LBUF_INT = ( KEEP(56) + 
     &         NSLAVES * NSLAVES ) * 5
     &               * KEEP(34)
        NB_BYTES = NB_BYTES + int(DMUMPS_LBUF_INT,8)
        IF ( EFF ) THEN
          IF (OOC_STRAT .GT. 0) THEN
            NB_INT = NB_INT + int(KEEP(225),8)
          ELSE
            NB_INT = NB_INT + int(KEEP(15),8)
          ENDIF
        ELSE
          IF (OOC_STRAT .GT. 0) THEN
            NB_INT = NB_INT +  int(
     &           KEEP(225) + 2 * max(PERLU,10) *
     &           ( KEEP(225) / 100 + 1 )
     &                              ,8)
          ELSE
            NB_INT = NB_INT +  int(
     &           KEEP(15) + 2 * max(PERLU,10) *
     &           ( KEEP(15) / 100 + 1 )
     &                              ,8)
          ENDIF
        ENDIF
        NB_INT = NB_INT + NSTEPS8
        NB_INT = NB_INT + NSTEPS8 * I8OVERI
        NB_INT = NB_INT + N8 + 4_8 * NSTEPS8 +
     &           MUMPS_GET_POOL_LENGTH(NA(1), KEEP, KEEP8)
        NB_INT = NB_INT + 2_8 * NSTEPS8 * I8OVERI
      END IF
      MEMORY_BYTES = NB_BYTES + NB_INT * int(KEEP(34),8) +
     &               NB_REAL * int(KEEP(35),8)
      MEMORY_BYTES = max( MEMORY_BYTES, TEMP )
      MEMORY_MBYTES = int( MEMORY_BYTES / 1000000_8 )  + 1
      RETURN
      END SUBROUTINE DMUMPS_MAX_MEM
      SUBROUTINE DMUMPS_SETMAXTOZERO(M_ARRAY, M_SIZE)
      IMPLICIT NONE
      INTEGER M_SIZE
      DOUBLE PRECISION M_ARRAY(M_SIZE)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      M_ARRAY=ZERO
      RETURN
      END SUBROUTINE DMUMPS_SETMAXTOZERO
      SUBROUTINE DMUMPS_COMPUTE_MAXPERCOL(
     &     A,ASIZE,NCOL,NROW,
     &     M_ARRAY,NMAX,COMPRESSCB,LROW1)
      IMPLICIT NONE
      INTEGER(8) :: ASIZE
      INTEGER NROW,NCOL,NMAX,LROW1
      LOGICAL COMPRESSCB
      DOUBLE PRECISION A(ASIZE)
      DOUBLE PRECISION M_ARRAY(NMAX)
      INTEGER I
      INTEGER(8):: APOS, J, LROW
      DOUBLE PRECISION ZERO,TMP
      PARAMETER (ZERO=0.0D0)
      M_ARRAY(1:NMAX) = ZERO
      APOS = 0_8
      IF (COMPRESSCB) THEN
        LROW=int(LROW1,8)
      ELSE
        LROW=int(NCOL,8)
      ENDIF
      DO I=1,NROW
         DO J=1_8,int(NMAX,8)
            TMP = abs(A(APOS+J))
            IF(TMP.GT.M_ARRAY(J)) M_ARRAY(J) = TMP
         ENDDO
         APOS = APOS + LROW
         IF (COMPRESSCB) LROW=LROW+1_8
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_COMPUTE_MAXPERCOL
      SUBROUTINE DMUMPS_SIZE_IN_STRUCT (id, NB_INT,NB_CMPLX )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE(DMUMPS_STRUC) :: id
      INTEGER(8) NB_INT, NB_CMPLX 
      INTEGER(8) NB_REAL
      NB_INT        = 0_8
      NB_CMPLX      = 0_8
      NB_REAL       = 0_8
      IF (associated(id%IS))          NB_INT=NB_INT+size(id%IS)
      IF (associated(id%IS1))         NB_INT=NB_INT+size(id%IS1)
      NB_INT=NB_INT+size(id%KEEP)
      NB_INT=NB_INT+size(id%ICNTL)
      NB_INT=NB_INT+size(id%INFO)
      NB_INT=NB_INT+size(id%INFOG)
      IF (associated(id%MAPPING))     NB_INT=NB_INT+size(id%MAPPING)
      IF (associated(id%POIDS))       NB_INT=NB_INT+size(id%POIDS)
      IF (associated(id%BUFR))        NB_INT=NB_INT+size(id%BUFR)
      IF (associated(id%STEP))        NB_INT=NB_INT+size(id%STEP)
      IF (associated(id%NE_STEPS  ))  NB_INT=NB_INT+size(id%NE_STEPS  )
      IF (associated(id%ND_STEPS))    NB_INT=NB_INT+size(id%ND_STEPS)
      IF (associated(id%Step2node))   NB_INT=NB_INT+size(id%Step2node)
      IF (associated(id%FRERE_STEPS)) NB_INT=NB_INT+size(id%FRERE_STEPS)
      IF (associated(id%DAD_STEPS))   NB_INT=NB_INT+size(id%DAD_STEPS)
      IF (associated(id%FILS))        NB_INT=NB_INT+size(id%FILS)
      IF (associated(id%PTRAR))       NB_INT=NB_INT+size(id%PTRAR)
      IF (associated(id%FRTPTR))      NB_INT=NB_INT+size(id%FRTPTR)
      NB_INT=NB_INT+size(id%KEEP8) * id%KEEP(10)
      IF (associated(id%PTRFAC)) NB_INT=NB_INT+size(id%PTRFAC) *
     &                                         id%KEEP(10)
      IF (associated(id%FRTELT))      NB_INT=NB_INT+size(id%FRTELT)
      IF (associated(id%NA))          NB_INT=NB_INT+size(id%NA)
      IF       (associated(id%PROCNODE_STEPS))
     &  NB_INT=NB_INT+size(id%PROCNODE_STEPS)
      IF (associated(id%PTLUST_S)) NB_INT=NB_INT+size(id%PTLUST_S)
      IF (associated(id%PROCNODE)) NB_INT=NB_INT+size(id%PROCNODE)
      IF (associated(id%INTARR)) NB_INT=NB_INT+size(id%INTARR)
      IF (associated(id%ELTPROC))  NB_INT=NB_INT+size(id%ELTPROC)
      IF (associated(id%CANDIDATES))
     &  NB_INT=NB_INT+size(id%CANDIDATES)
      IF       (associated(id%ISTEP_TO_INIV2))
     &  NB_INT=NB_INT+size(id%ISTEP_TO_INIV2)
      IF       (associated(id%FUTURE_NIV2))
     &  NB_INT=NB_INT+size(id%FUTURE_NIV2)
      IF (associated(id%TAB_POS_IN_PERE))
     &  NB_INT=NB_INT+size(id%TAB_POS_IN_PERE)
      IF (associated(id%I_AM_CAND))
     &  NB_INT=NB_INT+size(id%I_AM_CAND)
      IF (associated(id%MEM_DIST)) 
     &  NB_INT=NB_INT+size(id%MEM_DIST)
      IF (associated(id%POSINRHSCOMP_ROW))
     &  NB_INT=NB_INT+size(id%POSINRHSCOMP_ROW)
      IF       (associated(id%MEM_SUBTREE))
     &  NB_INT=NB_INT+size(id%MEM_SUBTREE)
      IF       (associated(id%MY_ROOT_SBTR))
     &  NB_INT=NB_INT+size(id%MY_ROOT_SBTR)
      IF       (associated(id%MY_FIRST_LEAF))
     &  NB_INT=NB_INT+size(id%MY_FIRST_LEAF)
      IF (associated(id%MY_NB_LEAF)) NB_INT=NB_INT+size(id%MY_NB_LEAF)
      IF (associated(id%DEPTH_FIRST)) NB_INT=NB_INT+size(id%DEPTH_FIRST)
      IF (associated(id%COST_TRAV)) NB_INT=NB_INT+size(id%COST_TRAV)
      IF (associated(id%CB_SON_SIZE)) NB_INT=NB_INT+size(id%CB_SON_SIZE)
      IF       (associated(id%OOC_INODE_SEQUENCE))
     &  NB_INT=NB_INT+size(id%OOC_INODE_SEQUENCE)
      IF       (associated(id%OOC_SIZE_OF_BLOCK))
     &  NB_INT=NB_INT+size(id%OOC_SIZE_OF_BLOCK)
      IF       (associated(id%OOC_VADDR)) 
     &  NB_INT=NB_INT+size(id%OOC_VADDR)
      IF       (associated(id%OOC_TOTAL_NB_NODES))
     &  NB_INT=NB_INT+size(id%OOC_TOTAL_NB_NODES)
      IF       (associated(id%OOC_NB_FILES))
     &  NB_INT=NB_INT+size(id%OOC_NB_FILES)
      IF       (associated(id%OOC_FILE_NAME_LENGTH))
     &  NB_INT=NB_INT+size(id%OOC_FILE_NAME_LENGTH)
      IF (associated(id%PIVNUL_LIST)) NB_INT=NB_INT+size(id%PIVNUL_LIST)
      IF (associated(id%SUP_PROC))    NB_INT=NB_INT+size(id%SUP_PROC)
      IF (associated(id%DBLARR))  NB_CMPLX=NB_CMPLX+size(id%DBLARR)
      IF (associated(id%RHSCOMP)) NB_CMPLX=NB_CMPLX+size(id%RHSCOMP)
      IF (associated(id%S))       NB_CMPLX=NB_CMPLX+id%KEEP8(23)
      IF (associated(id%COLSCA))  NB_REAL=NB_REAL+size(id%COLSCA)
      IF (associated(id%ROWSCA))  NB_REAL=NB_REAL+size(id%ROWSCA)
      NB_REAL=NB_REAL+size(id%CNTL)
      NB_REAL=NB_REAL+size(id%RINFO)
      NB_REAL=NB_REAL+size(id%RINFOG)
      NB_REAL=NB_REAL+size(id%DKEEP)
      NB_CMPLX = NB_CMPLX + NB_REAL
      RETURN
      END SUBROUTINE DMUMPS_SIZE_IN_STRUCT 
      SUBROUTINE DMUMPS_COPYI8SIZE(N8,SRC,DEST)
      IMPLICIT NONE
      INTEGER(8) :: N8
      DOUBLE PRECISION, intent(in)  :: SRC(N8)
      DOUBLE PRECISION, intent(out) :: DEST(N8)
      INTEGER(8) :: SHIFT8, HUG8
      INTEGER    :: I, I4SIZE
      IF(int(huge(I4SIZE),8) .EQ. int(huge(HUG8),8)) THEN
         CALL dcopy(N8, SRC(1), 1, DEST(1), 1)
      ELSE
         HUG8=int(huge(I4SIZE),8)
         DO I = 1, int(( N8 + HUG8 - 1_8 ) / HUG8)
            SHIFT8 = 1_8 + int(I-1,8) * HUG8
            I4SIZE = int(min(HUG8, N8-SHIFT8+1_8))
            CALL dcopy(I4SIZE, SRC(SHIFT8), 1, DEST(SHIFT8), 1)
         ENDDO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_COPYI8SIZE
      SUBROUTINE DMUMPS_SET_TMP_PTR( THE_ADDRESS, THE_SIZE )
      USE DMUMPS_STATIC_PTR_M
      INTEGER, INTENT(IN) :: THE_SIZE
      DOUBLE PRECISION, INTENT(IN) :: THE_ADDRESS(THE_SIZE)
      CALL DMUMPS_SET_STATIC_PTR(THE_ADDRESS(1:THE_SIZE)) 
      RETURN
      END SUBROUTINE DMUMPS_SET_TMP_PTR
