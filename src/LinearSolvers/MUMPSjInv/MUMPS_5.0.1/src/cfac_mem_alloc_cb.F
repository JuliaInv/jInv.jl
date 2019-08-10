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
      SUBROUTINE CMUMPS_ALLOC_CB( INPLACE, MIN_SPACE_IN_PLACE,
     &   SSARBR, PROCESS_BANDE,
     &   MYID,N, KEEP,KEEP8,DKEEP,
     &   IW, LIW, A, LA,
     &   LRLU, IPTRLU,IWPOS,IWPOSCB,
     &   PTRIST,PTRAST,STEP,PIMASTER,PAMASTER,
     &   LREQ, LREQCB, NODE_ARG, STATE_ARG, SET_HEADER,
     &   COMP, LRLUS, IFLAG, IERROR )
      USE CMUMPS_LOAD
      IMPLICIT NONE
      INTEGER N,LIW, KEEP(500)
      INTEGER(8) LA, LRLU, IPTRLU, LRLUS, LREQCB
      INTEGER(8) PAMASTER(KEEP(28)), PTRAST(KEEP(28))
      INTEGER IWPOS,IWPOSCB
      INTEGER(8) :: MIN_SPACE_IN_PLACE
      INTEGER NODE_ARG, STATE_ARG
      INTEGER(8) KEEP8(150)
      REAL DKEEP(130)
      INTEGER IW(LIW),PTRIST(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER MYID, IXXP
      COMPLEX A(LA)
      LOGICAL INPLACE, PROCESS_BANDE, SSARBR, SET_HEADER
      INTEGER COMP, LREQ, IFLAG, IERROR
      INCLUDE 'mumps_headers.h'
      INTEGER INODE_LOC,NPIV,NASS,NROW,NCB
      INTEGER ISIZEHOLE
      INTEGER(8) :: MEM_GAIN, RSIZEHOLE, LREQCB_EFF, LREQCB_WISHED
      LOGICAL DONE
      IF ( INPLACE ) THEN
        LREQCB_EFF = MIN_SPACE_IN_PLACE
        IF ( MIN_SPACE_IN_PLACE > 0_8 ) THEN
          LREQCB_WISHED = LREQCB
        ELSE
          LREQCB_WISHED = 0_8
        ENDIF
      ELSE
        LREQCB_EFF = LREQCB
        LREQCB_WISHED = LREQCB
      ENDIF
      IF (IWPOSCB.EQ.LIW) THEN
        IF (LREQ.NE.KEEP(IXSZ).OR.LREQCB.NE.0_8
     &      .OR. .NOT. SET_HEADER) THEN
          WRITE(*,*) "Internal error in CMUMPS_ALLOC_CB",
     &      SET_HEADER, LREQ, LREQCB
          CALL MUMPS_ABORT()
        ENDIF
        IF (IWPOSCB-IWPOS+1 .LT. KEEP(IXSZ)) THEN
          WRITE(*,*) "Problem with integer stack size",IWPOSCB,
     &               IWPOS, KEEP(IXSZ)
          IFLAG  = -8
          IERROR = LREQ
          RETURN
        ENDIF
        IWPOSCB=IWPOSCB-KEEP(IXSZ)
        IW(IWPOSCB+1+XXI)=KEEP(IXSZ)
        CALL MUMPS_STOREI8(0_8,IW(IWPOSCB+1+XXR))
        IW(IWPOSCB+1+XXN)=-919191
        IW(IWPOSCB+1+XXS)=S_NOTFREE
        IW(IWPOSCB+1+XXP)=TOP_OF_STACK
        RETURN
      ENDIF
      IF (KEEP(214).EQ.1.AND.
     &    KEEP(216).EQ.1.AND.
     &    IWPOSCB.NE.LIW) THEN
       IF (IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG.OR.
     &     IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG38) THEN
        NCB  = IW( IWPOSCB+1 + KEEP(IXSZ) )
        NROW = IW( IWPOSCB+1 + KEEP(IXSZ) + 2)
        NPIV = IW( IWPOSCB+1 + KEEP(IXSZ) + 3)
        INODE_LOC= IW( IWPOSCB+1 + XXN)
        CALL CMUMPS_GET_SIZEHOLE(IWPOSCB+1,IW,LIW,
     &                          ISIZEHOLE,RSIZEHOLE)
        IF (IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG) THEN
          CALL CMUMPS_MAKECBCONTIG(A,LA,IPTRLU+1_8,
     &                           NROW,NCB,NPIV+NCB,0,
     &                           IW(IWPOSCB+1 + XXS),RSIZEHOLE)
          IW(IWPOSCB+1 + XXS) =S_NOLCLEANED
          MEM_GAIN            = int(NROW,8)*int(NPIV,8)
        ENDIF
        IF (IW(IWPOSCB+1 + XXS).EQ.S_NOLCBNOCONTIG38) THEN
          NASS = IW( IWPOSCB+1 + KEEP(IXSZ) + 4)
          CALL CMUMPS_MAKECBCONTIG(A,LA,IPTRLU+1_8,
     &                           NROW,NCB,NPIV+NCB,NASS-NPIV,
     &                           IW(IWPOSCB+1 + XXS),RSIZEHOLE)
          IW(IWPOSCB+1 + XXS) =S_NOLCLEANED38
          MEM_GAIN = int(NROW,8)*int(NPIV+NCB-(NASS-NPIV),8)
        ENDIF
        IF (ISIZEHOLE.NE.0) THEN
          CALL CMUMPS_ISHIFT( IW,LIW,IWPOSCB+1,
     &                       IWPOSCB+IW(IWPOSCB+1+XXI),
     &                       ISIZEHOLE )
          IWPOSCB=IWPOSCB+ISIZEHOLE
          IW(IWPOSCB+1+XXP+IW(IWPOSCB+1+XXI))=IWPOSCB+1
          PTRIST(STEP(INODE_LOC))=PTRIST(STEP(INODE_LOC))+
     &    ISIZEHOLE
        ENDIF
        CALL MUMPS_SUBTRI8TOARRAY(IW(IWPOSCB+1+XXR), MEM_GAIN)
        IPTRLU              = IPTRLU+MEM_GAIN+RSIZEHOLE
        LRLU                = LRLU+MEM_GAIN+RSIZEHOLE
        PTRAST(STEP(INODE_LOC))=
     &  PTRAST(STEP(INODE_LOC))+MEM_GAIN+RSIZEHOLE
       ENDIF
      ENDIF
      DONE =.FALSE.
      IF ((IPTRLU.LT.LREQCB_WISHED).OR.(LRLU.LT.LREQCB_WISHED)) THEN
        IF (LRLUS.LT.LREQCB_EFF) THEN
          GOTO 620
        ELSE
          CALL CMUMPS_COMPRE_NEW(N,KEEP(28),IW,LIW,A,LA,
     &                    LRLU,IPTRLU,IWPOS,IWPOSCB,
     &                    PTRIST,PTRAST,
     &                    STEP, PIMASTER,PAMASTER,KEEP(216),LRLUS,
     &                    KEEP(IXSZ), COMP, DKEEP(97), MYID)
          IF ( LRLU .NE. LRLUS ) THEN
            WRITE(*,*) 'PB compress... CMUMPS_ALLOC_CB',
     &      'LRLU,LRLUS=',LRLU,LRLUS
            GOTO 620
          END IF
          DONE = .TRUE.
        ENDIF
      ENDIF
      IF (IWPOSCB-IWPOS+1 .LT. LREQ) THEN
       IF (DONE) GOTO 600
                 CALL CMUMPS_COMPRE_NEW(N,KEEP(28),IW,LIW,A,LA,
     &                    LRLU,IPTRLU,IWPOS,IWPOSCB,
     &                    PTRIST,PTRAST,
     &                    STEP, PIMASTER,PAMASTER,KEEP(216),LRLUS,
     &                    KEEP(IXSZ), COMP, DKEEP(97), MYID)
                 IF ( LRLU .NE. LRLUS ) THEN
                   WRITE(*,*) 'PB compress... CMUMPS_ALLOC_CB',
     &             'LRLU,LRLUS=',LRLU,LRLUS
                   GOTO 620
                 END IF
          IF (IWPOSCB-IWPOS+1 .LT. LREQ) GOTO 600
      ENDIF
      IXXP=IWPOSCB+XXP+1
      IF (IXXP.GT.LIW) THEN
        WRITE(*,*) "Internal error 3 in CMUMPS_ALLOC_CB",IXXP
      ENDIF
      IF (IW(IXXP).GT.0) THEN
        WRITE(*,*) "Internal error 2 in CMUMPS_ALLOC_CB",IW(IXXP),IXXP
      ENDIF
      IWPOSCB = IWPOSCB - LREQ
      IF (SET_HEADER) THEN
        IW(IXXP)= IWPOSCB + 1
        IW(IWPOSCB+1+XXI)=LREQ
        CALL MUMPS_STOREI8(LREQCB, IW(IWPOSCB+1+XXR))
        IW(IWPOSCB+1+XXS)=STATE_ARG
        IW(IWPOSCB+1+XXN)=NODE_ARG
        IW(IWPOSCB+1+XXP)=TOP_OF_STACK
        IW(IWPOSCB+1+XXP+1:IWPOSCB+1+KEEP(IXSZ))=-99999 
#if ! defined(NO_XXNBPR)
        IW(IWPOSCB+1+XXNBPR)=0 
#endif
      ENDIF
      IPTRLU = IPTRLU - LREQCB
      LRLU   = LRLU - LREQCB
      LRLUS  = LRLUS - LREQCB_EFF
      KEEP8(67) = min(LRLUS, KEEP8(67))
#if ! defined(OLD_LOAD_MECHANISM)
      CALL CMUMPS_LOAD_MEM_UPDATE(SSARBR,PROCESS_BANDE,
     &              LA-LRLUS,0_8,LREQCB_EFF,KEEP,KEEP8,LRLUS)
#else
#if defined (CHECK_COHERENCE)
      CALL CMUMPS_LOAD_MEM_UPDATE(SSARBR,PROCESS_BANDE,
     &              LA-LRLUS,0_8,LREQCB_EFF,KEEP,KEEP8,LRLUS)
#else
      CALL CMUMPS_LOAD_MEM_UPDATE(SSARBR,.FALSE.,
     &              LA-LRLUS,0_8,LREQCB_EFF,KEEP,KEEP8,LRLUS)
#endif
#endif
      RETURN
 600  IFLAG  = -8
      IERROR = LREQ
      RETURN
 620  IFLAG  = -9
      CALL MUMPS_SET_IERROR(LREQCB_EFF - LRLUS, IERROR)
      RETURN
      END SUBROUTINE CMUMPS_ALLOC_CB
