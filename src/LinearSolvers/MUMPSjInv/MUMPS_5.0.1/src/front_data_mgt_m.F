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
      MODULE MUMPS_FRONT_DATA_MGT_M
      IMPLICIT NONE
      PRIVATE
C     --------------------------------------------
C     This module contains routines to manage
C     handlers of various data associated to
C     active fronts *during the factorization*.
C
C     It should be initialized at the beginning
C     of the factorization and terminated at the
C     end of the factorization.
C
C     There are two types of data, see below.
C
C     'A' is for active type 2 fronts: list must
C         be empty at the end of the factorization
C
C     'F' will be for general fronts -- not used currently
C
C     Only handlers are managed in this module.
C     The data itself is in the module above using it.
C     For example, FAC_MAPROW_DATA_M manages MAPROW
C     messages that arrive too early. It handles an
C     array that contains all early MAPROW messages
C     and that is indexed with the handlers managed
C     by MUMPS_FRONT_DATA_MGT_M.
C
C     --------------------------------------------
C
C     ===============
C     Public routines
C     ===============
      PUBLIC :: MUMPS_FDM_INIT,
     &          MUMPS_FDM_END,
     &          MUMPS_FDM_START_IDX,
     &          MUMPS_FDM_END_IDX
C     STACK_FREE_IDX(1:NB_FREE_IDX) holds the NB_FREE_IDX indices
C                                   of free handlers
C     STACK_FREE_IDX(NB_FREE_IDX+1:size(STACK_FREE_IDX)) is trash data
      TYPE FDM_STRUC_T
        INTEGER :: NB_FREE_IDX
        INTEGER, DIMENSION(:), POINTER :: STACK_FREE_IDX => null()
        INTEGER, DIMENSION(:), POINTER :: COUNT_ACCESS   => null()
      END TYPE FDM_STRUC_T
      TYPE (FDM_STRUC_T), TARGET, SAVE :: FDM_A, FDM_F
      CONTAINS
C
      SUBROUTINE MUMPS_FDM_INIT(WHAT, INITIAL_SIZE, INFO)
C
C     Purpose:
C     =======
C
C     Initialize handler data ('A' or 'F')
C
C     Arguments:
C     =========
C
      INTEGER, INTENT(IN) :: INITIAL_SIZE
      CHARACTER, INTENT(IN) :: WHAT  ! 'A' or 'F'
      INTEGER, INTENT(INOUT) :: INFO(2)
C
C     Local variables:
C     ===============
C
      INTEGER :: IERR
      TYPE (FDM_STRUC_T), POINTER :: FDM_PTR
C
      CALL MUMPS_FDM_SET_PTR(WHAT, FDM_PTR)
      ALLOCATE( FDM_PTR%STACK_FREE_IDX(INITIAL_SIZE),
     &          FDM_PTR%COUNT_ACCESS  (INITIAL_SIZE), stat=IERR )
      IF (IERR < 0) THEN
        INFO(1) = -13
        INFO(2) = INITIAL_SIZE * 2
        RETURN
      ENDIF
      CALL MUMPS_FDM_SET_ALL_FREE(FDM_PTR)
      RETURN
      END SUBROUTINE MUMPS_FDM_INIT
C
      SUBROUTINE MUMPS_FDM_END(WHAT)
C
C     Purpose:
C     =======
C     Free module datastructures associated to "WHAT" at
C     the end of a phase (typically factorization).
C
      CHARACTER, INTENT(IN) :: WHAT
C
C     Local variables
C     ===============
C
      TYPE (FDM_STRUC_T), POINTER :: FDM_PTR
C
      CALL MUMPS_FDM_SET_PTR(WHAT, FDM_PTR)
      IF (associated(FDM_PTR%STACK_FREE_IDX)) THEN
          DEALLOCATE(FDM_PTR%STACK_FREE_IDX)
          NULLIFY(FDM_PTR%STACK_FREE_IDX)
          FDM_PTR%NB_FREE_IDX=0
      ELSE
C         Should not be called twice or when array is unassociated
          WRITE(*,*) "Internal error 1 in MUMPS_FDM_END", WHAT
          CALL MUMPS_ABORT()
      ENDIF
      IF (associated(FDM_PTR%COUNT_ACCESS)) THEN
          DEALLOCATE(FDM_PTR%COUNT_ACCESS)
          NULLIFY(FDM_PTR%COUNT_ACCESS)
      ELSE
C         Should not be called twice or when array is unassociated
          WRITE(*,*) "Internal error 1 in MUMPS_FDM_END", WHAT
          CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_FDM_END
C
      SUBROUTINE MUMPS_FDM_START_IDX(WHAT, FROM, IWHANDLER, INFO)
C
C     Purpose:
C     =======
C
C     Return a new free index/handler
C     (typically stored in IW)
C
      CHARACTER, INTENT(IN)  :: WHAT
      CHARACTER(LEN=*), INTENT(IN)  :: FROM !For debugging purposes only
      INTEGER, INTENT(INOUT) :: IWHANDLER
      INTEGER, INTENT(INOUT) :: INFO(2)
C
C     Local variables
C     ===============
C
      INTEGER :: OLD_SIZE, NEW_SIZE, IERR
      INTEGER :: I
      INTEGER, DIMENSION(:), POINTER :: TMP_COUNT_ACCESS
      TYPE(FDM_STRUC_T), POINTER :: FDM_PTR
      CALL MUMPS_FDM_SET_PTR(WHAT, FDM_PTR)
C
      IF (IWHANDLER .GT. 0) THEN
C       Already started, counter should at least be 1
        IF (FDM_PTR%COUNT_ACCESS(IWHANDLER) .LT. 1) THEN
          WRITE(*,*) "Internal error 1 in MUMPS_FDM_START_IDX",
     &    FDM_PTR%COUNT_ACCESS(IWHANDLER)
          CALL MUMPS_ABORT()
        ENDIF
        GOTO 100
        RETURN
      ENDIF
C
      IF (FDM_PTR%NB_FREE_IDX .EQ. 0) THEN
        OLD_SIZE = size(FDM_PTR%STACK_FREE_IDX)
        NEW_SIZE = (OLD_SIZE * 3) / 2 + 1 ! or something else
        FDM_PTR%NB_FREE_IDX = NEW_SIZE - OLD_SIZE
        DEALLOCATE(FDM_PTR%STACK_FREE_IDX)
        ALLOCATE(FDM_PTR%STACK_FREE_IDX(NEW_SIZE),
     &           TMP_COUNT_ACCESS(NEW_SIZE), stat=IERR)
        IF (IERR < 0) THEN
          INFO(1) = -13
          INFO(2) = NEW_SIZE
          RETURN
        ENDIF
C       All new handlers indices are created 
        DO I=1, FDM_PTR%NB_FREE_IDX
          FDM_PTR%STACK_FREE_IDX(I)=NEW_SIZE-I+1
        ENDDO
C       Count access: copy old ones
        DO I=1, OLD_SIZE
          TMP_COUNT_ACCESS(I)=FDM_PTR%COUNT_ACCESS(I)
        ENDDO
        DO I=OLD_SIZE+1, NEW_SIZE
          TMP_COUNT_ACCESS(I)=0
        ENDDO
        DEALLOCATE(FDM_PTR%COUNT_ACCESS)
        FDM_PTR%COUNT_ACCESS=>TMP_COUNT_ACCESS
      ENDIF
C
      IWHANDLER = FDM_PTR%STACK_FREE_IDX(FDM_PTR%NB_FREE_IDX)
      FDM_PTR%NB_FREE_IDX = FDM_PTR%NB_FREE_IDX - 1
 100  CONTINUE
C     Number of modules accessing this handler
      FDM_PTR%COUNT_ACCESS(IWHANDLER)=FDM_PTR%COUNT_ACCESS(IWHANDLER)+1
#if defined(DBG_FDM)
      WRITE(*,*) "DBG_FDM: IWHANDLER=",IWHANDLER, "Starting FROM=",FROM
#endif
      RETURN
      END SUBROUTINE MUMPS_FDM_START_IDX
C
      SUBROUTINE MUMPS_FDM_END_IDX(WHAT, FROM, IWHANDLER)
C
C     Purpose:
C     =======
C
C     Notify than an index/handler has been freed.
C     Mark it free for future reuse.
C
      CHARACTER, INTENT(IN) :: WHAT
      CHARACTER(LEN=*), INTENT(IN) :: FROM ! for debug purposes only
      INTEGER, INTENT(INOUT) :: IWHANDLER
      TYPE(FDM_STRUC_T), POINTER :: FDM_PTR
C
      CALL MUMPS_FDM_SET_PTR(WHAT, FDM_PTR)
      IF (IWHANDLER .LE.0) THEN
C       Already ended
        WRITE(*,*) "Internal error 1 in MUMPS_FDM_END_IDX",IWHANDLER
        CALL MUMPS_ABORT()
      ENDIF
#if defined(DBG_FDM)
      WRITE(*,*) "DBG_FDM: IWHANDLER=",IWHANDLER, "Ending FROM=",FROM
#endif
      FDM_PTR%COUNT_ACCESS(IWHANDLER)=FDM_PTR%COUNT_ACCESS(IWHANDLER)-1
      IF (FDM_PTR%COUNT_ACCESS(IWHANDLER) .LT. 0) THEN
C       Negative counter!
        WRITE(*,*) "Internal error 2 in MUMPS_FDM_END_IDX",
     &  IWHANDLER, FDM_PTR%COUNT_ACCESS(IWHANDLER)
        CALL MUMPS_ABORT()
      ENDIF
      IF (FDM_PTR%COUNT_ACCESS(IWHANDLER) .EQ.0 ) THEN
        IF (FDM_PTR%NB_FREE_IDX .GE. size(FDM_PTR%STACK_FREE_IDX)) THEN
          WRITE(*,*) "Internal error 3 in MUMPS_FDM_END_IDX"
          CALL MUMPS_ABORT()
        ENDIF
        FDM_PTR%NB_FREE_IDX = FDM_PTR%NB_FREE_IDX + 1
C       Having incremented the nb of free handlers we
C       store the index (IWHANDLER) that has been
C       effectively released for future reuse.
        FDM_PTR%STACK_FREE_IDX(FDM_PTR%NB_FREE_IDX) = IWHANDLER
        IWHANDLER = -8888 ! has been used and is now free
      ENDIF
C
      RETURN
      END SUBROUTINE MUMPS_FDM_END_IDX
C     ===================
C     Private subroutines
C     ===================
      SUBROUTINE MUMPS_FDM_SET_PTR(WHAT, FDM_PTR)
      CHARACTER, INTENT(IN) :: WHAT
#if defined(MUMPS_F2003)
      TYPE(FDM_STRUC_T), POINTER, INTENT(OUT) :: FDM_PTR
#else
      TYPE(FDM_STRUC_T), POINTER :: FDM_PTR
#endif
C
      IF ( WHAT .EQ. 'A' ) THEN
        FDM_PTR => FDM_A
      ELSE IF ( WHAT .EQ. 'F' ) THEN
        FDM_PTR => FDM_F
      ELSE
C       Should be called with either A or F
        WRITE(*,*) "Internal error 1 in MUMPS_FDM_INIT"
        WRITE(*,*) "Allowed arguments for WHAT are A or F"
        CALL MUMPS_ABORT()
      ENDIF
      END SUBROUTINE MUMPS_FDM_SET_PTR
      SUBROUTINE MUMPS_FDM_SET_ALL_FREE(FDM_PTR)
C
C     Purpose:
C     =======
C     Initialize the stack of free elements for the first time
C
      TYPE(FDM_STRUC_T), POINTER :: FDM_PTR
      INTEGER :: I
      FDM_PTR%NB_FREE_IDX = size(FDM_PTR%STACK_FREE_IDX)
      DO I = 1, FDM_PTR%NB_FREE_IDX
        FDM_PTR%STACK_FREE_IDX(I)=FDM_PTR%NB_FREE_IDX-I+1
        FDM_PTR%COUNT_ACCESS  (I)=0
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_FDM_SET_ALL_FREE
C
      END MODULE MUMPS_FRONT_DATA_MGT_M
