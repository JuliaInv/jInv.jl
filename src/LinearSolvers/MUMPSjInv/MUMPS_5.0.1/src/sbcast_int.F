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
      SUBROUTINE SMUMPS_MCAST2(DATA, LDATA, MPITYPE, ROOT, COMMW, TAG,
     &SLAVEF)
      USE SMUMPS_COMM_BUFFER
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR
      INTEGER LDATA, ROOT, COMMW, TAG, MPITYPE, SLAVEF
      INTEGER DEST
      INTEGER DATA(LDATA)
      DO 10 DEST = 0, SLAVEF - 1
        IF (DEST .NE. ROOT) THEN
          IF ( LDATA .EQ. 1 .and. MPITYPE .EQ. MPI_INTEGER ) THEN
            CALL SMUMPS_BUF_SEND_1INT( DATA(1), DEST, TAG, 
     &                                COMMW, IERR )
          ELSE
            WRITE(*,*) 'Error : bad argument to SMUMPS_MCAST2'
            CALL MUMPS_ABORT()
          END IF
        ENDIF
   10 CONTINUE
      RETURN
      END SUBROUTINE SMUMPS_MCAST2
      SUBROUTINE SMUMPS_BDC_ERROR( MYID, SLAVEF, COMM )
      INTEGER MYID, SLAVEF, COMM
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER DUMMY (1)
      CALL SMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID,
     &                 COMM, TERREUR, SLAVEF )
      RETURN
      END SUBROUTINE SMUMPS_BDC_ERROR
