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
      MODULE ZMUMPS_STATIC_PTR_M
      PUBLIC :: ZMUMPS_TMP_PTR, ZMUMPS_GET_TMP_PTR
      COMPLEX(kind=8), DIMENSION(:), POINTER, SAVE :: ZMUMPS_TMP_PTR
      CONTAINS
      SUBROUTINE ZMUMPS_SET_STATIC_PTR(ARRAY)
      COMPLEX(kind=8), DIMENSION(:), TARGET :: ARRAY
      ZMUMPS_TMP_PTR => ARRAY
      RETURN
      END SUBROUTINE ZMUMPS_SET_STATIC_PTR
      SUBROUTINE ZMUMPS_GET_TMP_PTR(PTR)
#if defined(MUMPS_F2003)
      COMPLEX(kind=8), DIMENSION(:), POINTER, INTENT(OUT) :: PTR
#else
      COMPLEX(kind=8), DIMENSION(:), POINTER :: PTR
#endif
      PTR => ZMUMPS_TMP_PTR
      RETURN
      END SUBROUTINE ZMUMPS_GET_TMP_PTR
      END MODULE ZMUMPS_STATIC_PTR_M
