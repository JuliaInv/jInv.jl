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
      MODULE MUMPS_L0_OMP_M
        LOGICAL, DIMENSION(:), POINTER :: NB_CORE_PER_THREAD_CHANGED
        INTEGER, DIMENSION(:), POINTER :: NB_CORE_PER_THREAD
        INTEGER :: THREAD_ID
        LOGICAL :: IS_ROOT_OF_L0_OMP
!$OMP   THREADPRIVATE ( THREAD_ID , IS_ROOT_OF_L0_OMP )
      END MODULE MUMPS_L0_OMP_M
