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
      SUBROUTINE SMUMPS_PROCESS_NODE( MYID,KEEP,KEEP8,DKEEP,
     &           BUFR, LBUFR, LBUFR_BYTES,
     &           IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS,
     &           N, IW, LIW, A, LA,
     &           PTRIST, PTRAST, STEP, PIMASTER, PAMASTER,
     &           NSTK_S, COMP,
     &           FPERE, FLAG, IFLAG, IERROR, COMM,
     &           ITLOC, RHS_MUMPS )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      INTEGER MYID
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER KEEP(500), BUFR( LBUFR )
      INTEGER(8) KEEP8(150)
      REAL DKEEP(130)
      INTEGER(8) :: LA, IPTRLU, LRLU, LRLUS
      INTEGER IWPOS, IWPOSCB
      INTEGER N, LIW
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER(8) :: PTRAST  (KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST( KEEP(28) )
      INTEGER STEP(N), PIMASTER(KEEP(28))
      INTEGER COMP, FPERE
      LOGICAL FLAG
      INTEGER NSTK_S( KEEP(28) ), ITLOC( N + KEEP(253) )
      REAL :: RHS_MUMPS(KEEP(255))
      INTEGER IFLAG, IERROR, COMM
      INTEGER POSITION, FINODE, FLCONT, LREQ
      INTEGER(8) :: LREQCB
      INTEGER(8) :: IPOS_NODE, ISHIFT_PACKET
      INTEGER SIZE_PACKET
      INTEGER NBROWS_ALREADY_SENT, NBROWS_PACKET
      INCLUDE 'mumps_headers.h'
      LOGICAL COMPRESSCB
      FLAG = .FALSE.
      POSITION = 0
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &                FINODE, 1, MPI_INTEGER,
     &                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &                FPERE, 1, MPI_INTEGER, 
     &                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &                FLCONT, 1, MPI_INTEGER,
     &                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &                NBROWS_ALREADY_SENT, 1, MPI_INTEGER,
     &                COMM, IERR)
      CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &                NBROWS_PACKET, 1, MPI_INTEGER,
     &                COMM, IERR)
      COMPRESSCB = (FLCONT.LT.0) 
      IF (COMPRESSCB) THEN
        FLCONT   = -FLCONT
        LREQCB  = (int(FLCONT,8) * int(FLCONT+1,8)) / 2_8
      ELSE
        LREQCB  = int(FLCONT,8) * int(FLCONT,8)
      ENDIF
      IF (NBROWS_ALREADY_SENT == 0) THEN
        LREQ    = 2 * FLCONT + 6 + KEEP(IXSZ)
        IF (IPTRLU.LT.0_8) WRITE(*,*) 'before alloc_cb:IPTRLU = ',IPTRLU
        CALL SMUMPS_ALLOC_CB( .FALSE., 0_8, .FALSE.,.FALSE.,
     &  MYID,N, KEEP,KEEP8, DKEEP, IW, LIW, A, LA,
     &  LRLU, IPTRLU,IWPOS,IWPOSCB,
     &  PTRIST,PTRAST,STEP, PIMASTER, PAMASTER,
     &  LREQ, LREQCB, FINODE, S_NOTFREE, .TRUE.,
     &  COMP, LRLUS, IFLAG, IERROR
     &     )
        IF (IPTRLU.LT.0_8) WRITE(*,*) 'after alloc_cb:IPTRLU = ',IPTRLU
        IF ( IFLAG .LT. 0 ) RETURN
        PIMASTER(STEP( FINODE )) = IWPOSCB + 1
        PAMASTER(STEP( FINODE )) = IPTRLU  + 1_8
        IF (COMPRESSCB)  IW(IWPOSCB + 1 + XXS ) = S_CB1COMP
        CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &        IW(IWPOSCB + 1+KEEP(IXSZ)), LREQ-KEEP(IXSZ),
     &        MPI_INTEGER, COMM, IERR)
      ENDIF
      IF (COMPRESSCB) THEN
        ISHIFT_PACKET = int(NBROWS_ALREADY_SENT,8) *
     &                  int(NBROWS_ALREADY_SENT+1,8) / 2_8
        SIZE_PACKET = (NBROWS_PACKET * (NBROWS_PACKET+1))/2 +
     &                 NBROWS_ALREADY_SENT * NBROWS_PACKET
      ELSE
        ISHIFT_PACKET = int(NBROWS_ALREADY_SENT,8) * int(FLCONT,8)
        SIZE_PACKET = NBROWS_PACKET * FLCONT
      ENDIF
      IF (NBROWS_PACKET.NE.0) THEN
        IF ( LREQCB .ne. 0_8 ) THEN
        IPOS_NODE = PAMASTER(STEP(FINODE))-1_8
        CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &        A(IPOS_NODE + 1_8 + ISHIFT_PACKET),
     &        SIZE_PACKET, MPI_REAL, COMM, IERR)
        END IF
      ENDIF
      IF (NBROWS_ALREADY_SENT+NBROWS_PACKET == FLCONT) THEN
        NSTK_S(STEP(FPERE)) = NSTK_S(STEP(FPERE)) - 1
        IF ( NSTK_S(STEP(FPERE)).EQ.0 ) THEN
          FLAG = . TRUE.
        END IF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_PROCESS_NODE
