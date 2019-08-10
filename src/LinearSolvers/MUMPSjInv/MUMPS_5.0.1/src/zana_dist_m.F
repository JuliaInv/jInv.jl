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
      SUBROUTINE ZMUMPS_ANA_DISTM(MYID, N, STEP, FRERE, FILS,
     &     NA, LNA, NE, DAD, ND, PROCNODE, SLAVEF,
     &     NRLADU, NIRADU, NIRNEC, NRLNEC,
     &     NRLNEC_ACTIVE, 
     &     NIRADU_OOC, NIRNEC_OOC,
     &     MAXFR, OPSA,
     &     KEEP,KEEP8, LOCAL_M, LOCAL_N, SBUF_RECOLD,
     &     SBUF_SEND, SBUF_REC, OPS_SUBTREE, NSTEPS,
     &     I_AM_CAND,NMB_PAR2, ISTEP_TO_INIV2, CANDIDATES, 
     &     IFLAG, IERROR
     &     ,MAX_FRONT_SURFACE_LOCAL
     &     ,MAX_SIZE_FACTOR, ENTRIES_IN_FACTORS_LOC
     &     ,ENTRIES_IN_FACTORS_LOC_MASTERS, ROOT_yes
     &     ,ROOT_NPROW, ROOT_NPCOL
     &     )
      IMPLICIT NONE
      LOGICAL, intent(in) :: ROOT_yes
      INTEGER, intent(in) :: ROOT_NPROW, ROOT_NPCOL
      INTEGER  MYID, N, LNA, IFLAG, IERROR
      INTEGER  NIRADU, NIRNEC
      INTEGER(8) NRLADU, NRLNEC, NRLNEC_ACTIVE
      INTEGER(8) NRLADU_CURRENT, NRLADU_ROOT_3
      INTEGER NIRADU_OOC, NIRNEC_OOC
      INTEGER MAXFR, NSTEPS
      INTEGER(8) MAX_FRONT_SURFACE_LOCAL
      INTEGER STEP(N)
      INTEGER FRERE(NSTEPS), FILS(N), NA(LNA), NE(NSTEPS),
     &        ND(NSTEPS), PROCNODE(NSTEPS), DAD(NSTEPS)
      INTEGER  SLAVEF, KEEP(500), LOCAL_M, LOCAL_N
      INTEGER(8) KEEP8(150)
      INTEGER(8) ENTRIES_IN_FACTORS_LOC,
     &          ENTRIES_IN_FACTORS_LOC_MASTERS
      INTEGER  SBUF_SEND, SBUF_REC
      INTEGER(8) SBUF_RECOLD
      INTEGER  NMB_PAR2
      INTEGER  ISTEP_TO_INIV2( KEEP(71) )
      LOGICAL  I_AM_CAND(NMB_PAR2)
      INTEGER  CANDIDATES( SLAVEF+1, NMB_PAR2 )
      DOUBLE PRECISION OPSA
      DOUBLE PRECISION OPSA_LOC 
      INTEGER(8) MAX_SIZE_FACTOR
      DOUBLE PRECISION OPS_SUBTREE
      DOUBLE PRECISION OPS_SBTR_LOC 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TNSTK, IPOOL, LSTKI 
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: LSTKR
      LOGICAL  OUTER_SENDS  
      INTEGER(8) SBUFS_CB, SBUFR_CB
      INTEGER SBUFR, SBUFS
      INTEGER BLOCKING_RHS
      INTEGER ITOP,NELIM,NFR
      INTEGER(8) ISTKR, LSTK
      INTEGER ISTKI,  STKI, ISTKI_OOC
      INTEGER K,NSTK, IFATH
      INTEGER INODE, LEAF, NBROOT, IN
      INTEGER LEVEL, MAXITEMPCB
      INTEGER(8) CURRENT_ACTIVE_MEM, MAXTEMPCB
      LOGICAL UPDATE, UPDATEF, MASTER, MASTERF, INSSARBR
      INTEGER LEVELF, NCB, SIZECBI
      INTEGER(8) NCB8
      INTEGER(8) NFR8, NELIM8
      INTEGER(8) SIZECB, SIZECBINFR, SIZECB_SLAVE
      INTEGER SIZEHEADER, SIZEHEADER_OOC, XSIZE_OOC
      INTEGER EXTRA_PERM_INFO_OOC
      INTEGER NBROWMAX, NSLAVES_LOC, NSLAVES_PASSED,
     &         NELIMF, NFRF, NCBF,
     &         NBROWMAXF, LKJIB,
     &         LKJIBT, NBR, NBCOLFAC
      INTEGER(8) LEV3MAXREC, CBMAXR, CBMAXS
      INTEGER ALLOCOK
      INTEGER PANEL_SIZE
      LOGICAL COMPRESSCB
      DOUBLE PRECISION OPS_NODE, OPS_NODE_MASTER, OPS_NODE_SLAVE
      INTEGER(8) ENTRIES_NODE_UPPER_PART, ENTRIES_NODE_LOWER_PART
      INCLUDE 'mumps_headers.h'
      INTEGER WHAT
      INTEGER(8) IDUMMY8
      INTRINSIC min, int, real
      INTEGER ZMUMPS_OOC_GET_PANEL_SIZE
      EXTERNAL ZMUMPS_OOC_GET_PANEL_SIZE
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR
      INTEGER MUMPS_BLOC2_GET_NSLAVESMAX
      EXTERNAL MUMPS_MAX_SURFCB_NBROWS, MUMPS_BLOC2_GET_NSLAVESMAX
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, 
     &         MUMPS_IN_OR_ROOT_SSARBR
      logical :: FORCE_CAND, CONCERNED, UPDATES, STACKCB, MASTERSON
      integer :: IFSON, LEVELSON
      IF (KEEP(50).eq.2) THEN
        EXTRA_PERM_INFO_OOC = 1
      ELSE IF (KEEP(50).eq.0) THEN
        EXTRA_PERM_INFO_OOC = 2
      ELSE
        EXTRA_PERM_INFO_OOC = 0
      ENDIF
      COMPRESSCB=( KEEP(215).EQ.0 .AND. KEEP(50).NE.0 )
      MAX_FRONT_SURFACE_LOCAL=0_8
      MAX_SIZE_FACTOR=0_8
      ALLOCATE( LSTKR(NSTEPS), TNSTK(NSTEPS), IPOOL(NSTEPS),
     &          LSTKI(NSTEPS) , stat=ALLOCOK)
      if (ALLOCOK .GT. 0) THEN
        IFLAG  =-7
        IERROR = 4*NSTEPS
        RETURN
      endif
      LKJIB = max(KEEP(5),KEEP(6))
      OUTER_SENDS = (KEEP(263).NE.0)
      IF (OUTER_SENDS) LKJIB = max(LKJIB,KEEP(488))
      TNSTK = NE
      LEAF = NA(1)+1
      IPOOL(1:LEAF-1) = NA(3:3+LEAF-2)
      NBROOT = NA(2)
#if defined(OLD_OOC_NOPANEL)
      XSIZE_OOC=XSIZE_OOC_NOPANEL
#else
      IF (KEEP(50).EQ.0) THEN
              XSIZE_OOC=XSIZE_OOC_UNSYM
      ELSE
              XSIZE_OOC=XSIZE_OOC_SYM
      ENDIF
#endif
      SIZEHEADER_OOC = XSIZE_OOC+6 
      SIZEHEADER = XSIZE_IC + 6  
      ISTKR      = 0_8
      ISTKI      = 0
      ISTKI_OOC  = 0
      OPSA_LOC   = 0.0D0
      ENTRIES_IN_FACTORS_LOC = 0_8
      ENTRIES_IN_FACTORS_LOC_MASTERS = 0_8
      OPS_SBTR_LOC = 0.0D0
      NRLADU     = 0_8
      NIRADU     = 0
      NIRADU_OOC = 0
      NRLADU_CURRENT = 0_8
      NRLADU_ROOT_3 = 0_8
      NRLNEC_ACTIVE = 0_8
      NRLNEC     = 0_8
      NIRNEC     = 0
      NIRNEC_OOC = 0
      MAXFR      = 0
      ITOP       = 0
      MAXTEMPCB  = 0_8
      MAXITEMPCB = 0
      SBUFS_CB   = 1_8
      SBUFS      = 1
      SBUFR_CB   = 1_8
      SBUFR      = 1
      IF (KEEP(38) .NE. 0 .AND. KEEP(60).EQ.0) THEN
        INODE  = KEEP(38)
        NRLADU_ROOT_3 = int(LOCAL_M,8)*int(LOCAL_N,8)
        NRLADU = NRLADU_ROOT_3
        NRLNEC_ACTIVE = NRLADU_CURRENT
        MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_ROOT_3)
        NRLNEC = NRLADU
        IF (MUMPS_PROCNODE(PROCNODE(STEP(INODE)),SLAVEF)
     &                                       .EQ. MYID) THEN
          NIRADU     = SIZEHEADER+2*(ND(STEP(INODE))+KEEP(253))
          NIRADU_OOC = SIZEHEADER_OOC+2*(ND(STEP(INODE))+KEEP(253))
        ELSE
          NIRADU     = SIZEHEADER
          NIRADU_OOC = SIZEHEADER_OOC
        ENDIF
        NIRNEC     = NIRADU
        NIRNEC_OOC = NIRADU_OOC
      ENDIF
      IF((KEEP(24).eq.0).OR.(KEEP(24).eq.1)) THEN
         FORCE_CAND=.FALSE.           
      ELSE
         FORCE_CAND=(mod(KEEP(24),2).eq.0)
      END IF
 90   CONTINUE
      IF (LEAF.NE.1) THEN
         LEAF = LEAF - 1
         INODE = IPOOL(LEAF)
      ELSE 
         WRITE(MYID+6,*) ' ERROR 1 in ZMUMPS_ANA_DISTM '
         CALL MUMPS_ABORT()
      ENDIF
 95   CONTINUE 
      NFR    = ND(STEP(INODE))+KEEP(253)
      NFR8   = int(NFR,8)
      NSTK   = NE(STEP(INODE))
      NELIM = 0 
        IN = INODE
 100    NELIM = NELIM + 1 
      NELIM8=int(NELIM,8)
        IN = FILS(IN)
        IF (IN .GT. 0 ) GOTO 100
      IFSON = -IN
      IFATH = DAD(STEP(INODE))
      MASTER = MUMPS_PROCNODE(PROCNODE(STEP(INODE)),SLAVEF)
     &           .EQ. MYID
      LEVEL  = MUMPS_TYPENODE(PROCNODE(STEP(INODE)),SLAVEF)
      INSSARBR = MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &        SLAVEF)
      UPDATE=.FALSE.
       if(.NOT.FORCE_CAND) then
         UPDATE = ( (MASTER.AND.(LEVEL.NE.3) ).OR. LEVEL.EQ.2 )
       else
         if(MASTER.and.(LEVEL.ne.3)) then
            UPDATE = .TRUE.
         else if(LEVEL.eq.2) then
            if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(INODE)))) THEN
              UPDATE = .TRUE.
            end if
         end if
       end if
      NCB      = NFR-NELIM
      NCB8     = int(NCB,8)
      SIZECBINFR = NCB8*NCB8
      IF (KEEP(50).EQ.0) THEN
        SIZECB = SIZECBINFR
      ELSE
        IFATH = DAD(STEP(INODE))
        IF ( IFATH.NE.KEEP(38) .AND. COMPRESSCB ) THEN
          SIZECB    = (NCB8*(NCB8+1_8))/2_8
        ELSE
          SIZECB    = SIZECBINFR
        ENDIF
      ENDIF
      SIZECBI      = 2* NCB  + SIZEHEADER
      IF (LEVEL.NE.2) THEN
        NSLAVES_LOC     = -99999999
        SIZECB_SLAVE = -99999997_8
        NBROWMAX        = NCB
      ELSE
        IF (KEEP(48) .EQ. 5) THEN
          WHAT = 5 
          IF (FORCE_CAND) THEN
            NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     &                    ISTEP_TO_INIV2(STEP(INODE)))
          ELSE
            NSLAVES_LOC=SLAVEF-1
          ENDIF
          NSLAVES_PASSED=NSLAVES_LOC
        ELSE
          WHAT = 2 
          NSLAVES_PASSED=SLAVEF
          NSLAVES_LOC   =SLAVEF-1
        ENDIF
         CALL MUMPS_MAX_SURFCB_NBROWS(WHAT, KEEP,KEEP8,
     &     NCB, NFR, NSLAVES_PASSED, NBROWMAX, SIZECB_SLAVE
     &    )
      ENDIF
      IF (KEEP(60).GT.1) THEN
         IF (MASTER .AND. INODE.EQ.KEEP(38)) THEN
          NIRADU     = NIRADU+SIZEHEADER+2*(ND(STEP(INODE))+KEEP(253))
          NIRADU_OOC = NIRADU_OOC+SIZEHEADER_OOC+
     &                 2*(ND(STEP(INODE))+KEEP(253))
         ENDIF
      ENDIF
      IF (LEVEL.EQ.3) THEN
         IF ( 
     &     KEEP(60).LE.1
     &      ) THEN
           NRLNEC = max(NRLNEC,NRLADU+ISTKR+
     &                 int(LOCAL_M,8)*int(LOCAL_N,8))
           NRLADU_CURRENT = int(LOCAL_M,8)*int(LOCAL_N,8)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_ROOT_3 + 
     &                        NRLADU_CURRENT+ISTKR)
         ENDIF
         IF (MASTER) THEN 
            IF (NFR.GT.MAXFR) MAXFR = NFR
         ENDIF
      ENDIF
      IF(KEEP(86).EQ.1)THEN
         IF(MASTER.AND.(.NOT.MUMPS_IN_OR_ROOT_SSARBR(
     &        PROCNODE(STEP(INODE)), SLAVEF))
     &     )THEN
            IF(LEVEL.EQ.1)THEN
               MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &              NFR8*NFR8)
            ELSEIF(LEVEL.EQ.2)THEN
               IF(KEEP(50).EQ.0)THEN
                 MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                 NFR8*NELIM8)
               ELSE
                 MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                 NELIM8*NELIM8)
                 IF (KEEP(219).NE.0.AND.KEEP(50).EQ.2) THEN
                  MAX_FRONT_SURFACE_LOCAL=max(MAX_FRONT_SURFACE_LOCAL,
     &                  NELIM8*(NELIM8+1_8))
                 ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF (LEVEL.EQ.2) THEN
        IF (MASTER) THEN
          IF (KEEP(50).EQ.0) THEN
             SBUFS = max(SBUFS, NFR*LKJIB+LKJIB+4)
          ELSE
             SBUFS = max(SBUFS, NELIM*LKJIB+NELIM+6)
          ENDIF
        ELSEIF (UPDATE) THEN
            if (KEEP(50).EQ.0) THEN
              SBUFR   = max(SBUFR, NFR*LKJIB+LKJIB+4)
            else
              SBUFR = max( SBUFR, NELIM*LKJIB+NELIM+6 )
              IF (KEEP(50).EQ.1) THEN
                LKJIBT  = LKJIB
              ELSE
                LKJIBT  = min( NELIM, LKJIB * 2 )
              ENDIF
              SBUFS = max(SBUFS,
     &                        LKJIBT*NBROWMAX+6)
              SBUFR = max( SBUFR, NBROWMAX*LKJIBT+6 )
            endif
        ENDIF
      ENDIF
      IF ( UPDATE ) THEN
          IF ( (MASTER) .AND. (LEVEL.EQ.1) ) THEN
            NIRADU     = NIRADU + 2*NFR + SIZEHEADER
            NIRADU_OOC = NIRADU_OOC + 2*NFR + SIZEHEADER_OOC
            PANEL_SIZE = ZMUMPS_OOC_GET_PANEL_SIZE(
     &      2_8*int(KEEP(226),8), NFR, KEEP(227), KEEP(50))
            NIRADU_OOC = NIRADU_OOC +
     &      EXTRA_PERM_INFO_OOC*(2+NELIM + NELIM/PANEL_SIZE+1)
            IF (KEEP(50).EQ.0) THEN
             NRLADU_CURRENT = int(NELIM,8)*int(2*NFR-NELIM,8)
             NRLADU = NRLADU + NRLADU_CURRENT
             MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
            ELSE
             NRLADU_CURRENT = int(NELIM,8)*int(NFR,8)
             NRLADU = NRLADU + NRLADU_CURRENT
             MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
            ENDIF
            SIZECBI        = 2* NCB  + 6 + 3
          ELSEIF (LEVEL.EQ.2) THEN
            IF (MASTER) THEN
              NIRADU     = NIRADU+SIZEHEADER +SLAVEF-1+2*NFR 
              NIRADU_OOC = NIRADU_OOC+SIZEHEADER_OOC +SLAVEF-1+2*NFR 
              IF (KEEP(50).EQ.0) THEN
                NBCOLFAC=NFR
              ELSE
                NBCOLFAC=NELIM
              ENDIF
              PANEL_SIZE = ZMUMPS_OOC_GET_PANEL_SIZE(
     &        2_8*int(KEEP(226),8), NBCOLFAC, KEEP(227), KEEP(50))
              NIRADU_OOC = NIRADU_OOC +
     &        EXTRA_PERM_INFO_OOC*(2+NELIM + NELIM/PANEL_SIZE+1)
              NRLADU_CURRENT = int(NBCOLFAC,8)*int(NELIM,8)
              NRLADU = NRLADU + NRLADU_CURRENT
              MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
               SIZECB     = 0_8
               SIZECBINFR = 0_8
               SIZECBI    = NCB + 5 +  SLAVEF - 1
            ELSE
             SIZECB=SIZECB_SLAVE
             SIZECBINFR = SIZECB
             NIRADU       = NIRADU+4+NELIM+NBROWMAX
             NIRADU_OOC   = NIRADU_OOC+4+NELIM+NBROWMAX
             IF (KEEP(50).EQ.0) THEN
               NRLADU_CURRENT   =  int(NELIM,8)*int(NBROWMAX,8)
               NRLADU   = NRLADU + NRLADU_CURRENT
             ELSE
               NRLADU_CURRENT   = int(NELIM,8)*int(NCB/NSLAVES_LOC,8)
               NRLADU   = NRLADU + NRLADU_CURRENT
             ENDIF
             MAX_SIZE_FACTOR=max(MAX_SIZE_FACTOR,NRLADU_CURRENT)
             SIZECBI       = 4 + NBROWMAX + NCB
             IF (KEEP(50).NE.0) THEN 
                     SIZECBI=SIZECBI+NSLAVES_LOC+
     &                                  XTRA_SLAVES_SYM
             ELSE
                     SIZECBI=SIZECBI+NSLAVES_LOC+
     &                                  XTRA_SLAVES_UNSYM 
             ENDIF
            ENDIF
         ENDIF
         NIRNEC = max0(NIRNEC,
     &             NIRADU+ISTKI+SIZECBI+MAXITEMPCB)
         NIRNEC_OOC = max0(NIRNEC_OOC,
     &             NIRADU_OOC+ISTKI_OOC+SIZECBI+MAXITEMPCB +
     &             (XSIZE_OOC-XSIZE_IC) ) 
         CURRENT_ACTIVE_MEM = ISTKR+SIZECBINFR
         IF (NSTK .NE. 0 .AND. INSSARBR .AND.
     &     KEEP(234).NE.0 .AND. KEEP(55).EQ.0) THEN
           CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM - LSTKR(ITOP)
         ENDIF
         IF (KEEP(50).NE.0.AND.UPDATE.AND.LEVEL.EQ.1) THEN
             CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM +
     &              int(NELIM,8)*int(NCB,8)
         ENDIF
         IF (MASTER .AND.  KEEP(219).NE.0.AND.
     &       KEEP(50).EQ.2.AND.LEVEL.EQ.2) THEN
             CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM + int(NELIM,8)
         ENDIF
         IF (SLAVEF.EQ.1) THEN
           NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
         ELSE
           NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+MAXTEMPCB)
           NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &             NRLADU_ROOT_3+CURRENT_ACTIVE_MEM+MAXTEMPCB)
         ENDIF
         IF (NFR.GT.MAXFR) MAXFR = NFR
         IF (NSTK.GT.0) THEN
            DO 70 K=1,NSTK
               LSTK = LSTKR(ITOP)
               ISTKR = ISTKR - LSTK
               IF (K==1 .AND. INSSARBR.AND.KEEP(234).NE.0
     &            .AND.KEEP(55).EQ.0) THEN
               ELSE
                 CURRENT_ACTIVE_MEM = CURRENT_ACTIVE_MEM - LSTK
               ENDIF
               STKI = LSTKI( ITOP )
               ISTKI = ISTKI - STKI
               ISTKI_OOC = ISTKI_OOC - STKI - (XSIZE_OOC-XSIZE_IC)
               ITOP = ITOP - 1
               IF (ITOP.LT.0) THEN
                  write(*,*) MYID,
     &            ': ERROR 2 in ZMUMPS_ANA_DISTM. ITOP = ',ITOP
                  CALL MUMPS_ABORT()
               ENDIF
 70         CONTINUE
         ENDIF 
      ELSE IF (LEVEL.NE.3) THEN
         DO WHILE (IFSON.GT.0) 
            UPDATES=.FALSE.
            MASTERSON = MUMPS_PROCNODE(PROCNODE(STEP(IFSON)),SLAVEF)
     &                  .EQ.MYID
            LEVELSON  = MUMPS_TYPENODE(PROCNODE(STEP(IFSON)),SLAVEF)
            if(.NOT.FORCE_CAND) then
               UPDATES =((MASTERSON.AND.(LEVELSON.NE.3)).OR. 
     &                   LEVELSON.EQ.2)
            else
               if(MASTERSON.and.(LEVELSON.ne.3)) then
                  UPDATES = .TRUE.
               else if(LEVELSON.eq.2) then
                  if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(IFSON)))) then
                    UPDATES = .TRUE.
                  end if
               end if
            end if
            IF (UPDATES) THEN
              LSTK = LSTKR(ITOP)
              ISTKR = ISTKR - LSTK
              STKI = LSTKI( ITOP )
              ISTKI = ISTKI - STKI
              ISTKI_OOC = ISTKI_OOC - STKI - (XSIZE_OOC-XSIZE_IC)
              ITOP = ITOP - 1
              IF (ITOP.LT.0) THEN
                write(*,*) MYID,
     &          ': ERROR 2 in ZMUMPS_ANA_DISTM. ITOP = ',ITOP
                CALL MUMPS_ABORT()
              ENDIF
            ENDIF
            IFSON = FRERE(STEP(IFSON)) 
         END DO
      ENDIF
      IF (
     &        ( (INODE.NE.KEEP(20)).OR.(KEEP(60).EQ.0) ) 
     &       .AND.
     &        ( (INODE.NE.KEEP(38)).OR.(KEEP(60).LE.1) ) 
     &      )
     &  THEN
            ENTRIES_NODE_LOWER_PART = int(NFR-NELIM,8) * int(NELIM,8)
            IF ( KEEP(50).EQ.0 ) THEN
              ENTRIES_NODE_UPPER_PART = int(NFR,8) * int(NELIM,8)
            ELSE
              ENTRIES_NODE_UPPER_PART =
     &        (int(NELIM,8)*int(NELIM+1,8))/2_8
            ENDIF
            IF (KEEP(50).EQ.2 .AND. LEVEL.EQ.3) THEN
              CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM, 0,
     &           1,OPS_NODE)
            ELSE
              CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           1,OPS_NODE)
            ENDIF
            IF (LEVEL.EQ.2) THEN
              CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           2,OPS_NODE_MASTER)
              OPS_NODE_SLAVE=OPS_NODE-OPS_NODE_MASTER
            ENDIF
      ELSE
           OPS_NODE = 0.0D0
           ENTRIES_NODE_UPPER_PART = 0_8
           ENTRIES_NODE_LOWER_PART = 0_8
      ENDIF
      IF ( MASTER )  THEN
        ENTRIES_IN_FACTORS_LOC_MASTERS = 
     &                     ENTRIES_IN_FACTORS_LOC_MASTERS +
     &                            ENTRIES_NODE_UPPER_PART +
     &                            ENTRIES_NODE_LOWER_PART
      ENDIF
      IF (UPDATE.OR.LEVEL.EQ.3) THEN
         IF ( LEVEL .EQ. 3 ) THEN
            IF (ROOT_yes) THEN
              CALL MUMPS_UPDATE_FLOPS_ROOT( OPSA_LOC, KEEP(50), NFR,
     &             NFR, ROOT_NPROW, ROOT_NPCOL, MYID )
              ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                            ENTRIES_NODE_UPPER_PART /
     &                            int(ROOT_NPROW*ROOT_NPCOL,8)
              IF (MASTER) THEN
                ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                                 mod(ENTRIES_NODE_UPPER_PART,
     &                                 int(SLAVEF,8))
              ENDIF
            ENDIF
         ELSE IF (MASTER .AND. LEVEL.EQ.2) THEN
            OPSA_LOC = OPSA_LOC + OPS_NODE_MASTER
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                      ENTRIES_NODE_UPPER_PART +
     &                      mod(ENTRIES_NODE_LOWER_PART,
     &                          int(NSLAVES_LOC,8))
         ELSE IF (MASTER .AND. LEVEL.EQ.1) THEN
            OPSA_LOC = OPSA_LOC + dble(OPS_NODE)
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC +
     &                               ENTRIES_NODE_UPPER_PART +
     &                               ENTRIES_NODE_LOWER_PART
         ELSE IF (UPDATE) THEN 
            OPSA_LOC = OPSA_LOC + 
     &            dble(OPS_NODE_SLAVE)/dble(NSLAVES_LOC)
            ENTRIES_IN_FACTORS_LOC = ENTRIES_IN_FACTORS_LOC
     &                 + ENTRIES_NODE_LOWER_PART / 
     &                 int(NSLAVES_LOC,8)
         ENDIF
         IF (MUMPS_IN_OR_ROOT_SSARBR(PROCNODE(STEP(INODE)),
     &        SLAVEF) .OR. NE(STEP(INODE))==0) THEN
           IF (LEVEL == 1) THEN
             OPS_SBTR_LOC = OPS_SBTR_LOC + OPS_NODE
           ELSE
             CALL MUMPS_GET_FLOPS_COST(NFR, 
     &           NELIM, NELIM,KEEP(50),
     &           1,OPS_NODE)
             OPS_SBTR_LOC = OPS_SBTR_LOC + OPS_NODE
           ENDIF
         ENDIF
        ENDIF
      IF (IFATH .EQ. 0) THEN
         NBROOT = NBROOT - 1
         IF (NBROOT.EQ.0) GOTO 115
         GOTO 90
      ELSE
         NFRF = ND(STEP(IFATH))+KEEP(253)
         IF (DAD(STEP(IFATH)).EQ.0) THEN
           NELIMF = NFRF
         ELSE
           NELIMF = 0
           IN = IFATH
           DO WHILE (IN.GT.0)
              IN = FILS(IN)
              NELIMF = NELIMF+1
           ENDDO
         ENDIF
         NCBF = NFRF - NELIMF
         LEVELF = MUMPS_TYPENODE(PROCNODE(STEP(IFATH)),SLAVEF)
         MASTERF= MUMPS_PROCNODE(PROCNODE(STEP(IFATH)),SLAVEF).EQ.MYID
         UPDATEF= .FALSE.
         if(.NOT.FORCE_CAND) then
            UPDATEF= ((MASTERF.AND.(LEVELF.NE.3)).OR.LEVELF.EQ.2)
         else
            if(MASTERF.and.(LEVELF.ne.3)) then
               UPDATEF = .TRUE.
            else if (LEVELF.eq.2) then
               if ( I_AM_CAND(ISTEP_TO_INIV2(STEP(IFATH)))) THEN
                 UPDATEF = .TRUE.
               end if
            end if
         end if
         CONCERNED  = UPDATEF .OR. UPDATE
         IF (LEVELF .NE. 2) THEN
           NBROWMAXF = -999999
         ELSE
           IF (KEEP(48) .EQ. 5) THEN
               WHAT = 4
               IF (FORCE_CAND) THEN
                 NSLAVES_LOC=CANDIDATES(SLAVEF+1,
     &               ISTEP_TO_INIV2(STEP(IFATH)))
               ELSE
                 NSLAVES_LOC=SLAVEF-1
               ENDIF
           ELSE
               WHAT = 1 
               NSLAVES_LOC=SLAVEF
           ENDIF
           CALL MUMPS_MAX_SURFCB_NBROWS( WHAT, KEEP, KEEP8,
     &     NCBF, NFRF, NSLAVES_LOC, NBROWMAXF, IDUMMY8 
     &          )
         ENDIF
         IF(LEVEL.EQ.1.AND.UPDATE.AND.
     &      (UPDATEF.OR.LEVELF.EQ.2)
     &      .AND.LEVELF.NE.3) THEN
             IF ( INSSARBR .AND. KEEP(234).NE.0) THEN
               NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &           NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
               NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM)
             ELSE
               NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,NRLADU_CURRENT+
     &           NRLADU_ROOT_3+CURRENT_ACTIVE_MEM+SIZECB)
               NRLNEC = max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+SIZECB)
             ENDIF
         ENDIF
         IF (UPDATE .AND. LEVEL.EQ.2 .AND. .NOT. MASTER) THEN
             NRLNEC =
     &         max(NRLNEC,NRLADU+CURRENT_ACTIVE_MEM+NRLADU_CURRENT)
             NRLNEC_ACTIVE = max(NRLNEC_ACTIVE,2_8*NRLADU_CURRENT+
     &         NRLADU_ROOT_3+CURRENT_ACTIVE_MEM)
         ENDIF
        IF (LEVELF.EQ.3) THEN
          IF (LEVEL.EQ.1) THEN
            LEV3MAXREC = int(min(NCB,LOCAL_M),8) *
     &                   int(min(NCB,LOCAL_N),8)
          ELSE
            LEV3MAXREC = min(SIZECB,
     &                 int(min(NBROWMAX,LOCAL_M),8)
     &                *int(min(NCB,LOCAL_N),8))
          ENDIF
          MAXTEMPCB  = max(MAXTEMPCB, LEV3MAXREC)
          MAXITEMPCB = max(MAXITEMPCB,SIZECBI+SIZEHEADER)
          SBUFR_CB   = max(SBUFR_CB, LEV3MAXREC+int(SIZECBI,8))
          NIRNEC = max(NIRNEC,NIRADU+ISTKI+
     &    min(NCB,LOCAL_M)+ min(NCB,LOCAL_N)+SIZEHEADER)
          NIRNEC_OOC = max(NIRNEC_OOC,NIRADU_OOC+ISTKI_OOC+
     &    min(NCB,LOCAL_M)+ min(NCB,LOCAL_N)+SIZEHEADER)
        ENDIF
        IF (CONCERNED) THEN
         IF (LEVELF.EQ.2) THEN
           IF (UPDATE.AND.(LEVEL.NE.2.OR..NOT.MASTER)) THEN
             IF(MASTERF)THEN
                 NBR = min(NBROWMAXF,NBROWMAX)
             ELSE
                 NBR = min(max(NELIMF,NBROWMAXF),NBROWMAX)
             ENDIF
             IF (KEEP(50).EQ.0) THEN
               CBMAXS = int(NBR,8)*int(NCB,8)
             ELSE
               CBMAXS = int(NBR,8)*int(NCB,8) -
     &                  (int(NBR,8)*int(NBR-1,8))/2_8
             ENDIF
           ELSE
              CBMAXS = 0_8
           END IF
           IF (MASTERF) THEN
             IF (LEVEL.EQ.1) THEN
                IF (.NOT.UPDATE) THEN
                  NBR = min(NELIMF, NCB)
                ELSE
                  NBR = 0
                ENDIF
             ELSE
                NBR = min(NELIMF, NBROWMAX)
             ENDIF
             IF (KEEP(50).EQ.0) THEN
                CBMAXR = int(NBR,8)*NCB8
             ELSE
                CBMAXR = int(NBR,8)*int(min(NCB,NELIMF),8)-
     &                   (int(NBR,8)*int(NBR-1,8))/2_8
                CBMAXR = min(CBMAXR, int(NELIMF,8)*int(NELIMF+1,8)/2_8)
                CBMAXR = min(CBMAXR, SIZECB)
                IF ((LEVEL.EQ.1).AND.(.NOT. COMPRESSCB)) THEN
                  CBMAXR = min(CBMAXR,(NCB8*(NCB8+1_8))/2_8)
                ENDIF
             ENDIF
           ELSE IF (UPDATEF) THEN
              NBR = min(NBROWMAXF,NBROWMAX)
              CBMAXR = int(NBR,8) * NCB8
              IF (KEEP(50).NE.0) THEN
                CBMAXR = CBMAXR - (int(NBR,8)*(int(NBR-1,8)))/2_8
              ENDIF
           ELSE
              CBMAXR = 0_8
           ENDIF
         ELSEIF (LEVELF.EQ.3) THEN
           CBMAXR = LEV3MAXREC
           IF (UPDATE.AND. .NOT. (MASTER.AND.LEVEL.EQ.2)) THEN
             CBMAXS = LEV3MAXREC
           ELSE
             CBMAXS = 0_8
           ENDIF
         ELSE
           IF (MASTERF) THEN
             CBMAXS = 0_8
             NBR = min(NFRF,NBROWMAX)
             IF ((LEVEL.EQ.1).AND.UPDATE) THEN
                NBR = 0
             ENDIF
             CBMAXR = int(NBR,8)*int(min(NFRF,NCB),8)
             IF (LEVEL.EQ.2)
     &       CBMAXR = min(CBMAXR, SIZECB_SLAVE)
             IF ( KEEP(50).NE.0 )  THEN
              CBMAXR = min(CBMAXR,(int(NFRF,8)*int(NFRF+1,8))/2_8)
             ELSE
              CBMAXR = min(CBMAXR,int(NFRF,8)*int(NFRF,8))
             ENDIF
           ELSE
             CBMAXR = 0_8
             CBMAXS = SIZECB
           ENDIF
         ENDIF
         IF (UPDATE) THEN
           CBMAXS = min(CBMAXS, SIZECB)
           IF ( .not. ( LEVELF .eq. 1 .AND. UPDATEF ) )THEN
              SBUFS_CB = max(SBUFS_CB, CBMAXS+int(SIZECBI,8))
           ENDIF
         ENDIF
         STACKCB = .FALSE.
         IF (UPDATEF) THEN
          STACKCB = .TRUE.
          SIZECBI     = 2 * NFR + SIZEHEADER
          IF (LEVEL.EQ.1) THEN
             IF (KEEP(50).NE.0.AND.LEVELF.NE.3
     &           .AND.COMPRESSCB) THEN
                 SIZECB = (NCB8*(NCB8+1_8))/2_8
             ELSE
                 SIZECB = NCB8*NCB8
             ENDIF
             IF (MASTER) THEN
               SIZECBI     = 2+ XSIZE_IC
             ELSE IF (LEVELF.EQ.1) THEN
               SIZECB  = min(CBMAXR,SIZECB)
               SIZECBI    = 2 * NCB +  9 
               SBUFR_CB   = max(SBUFR_CB, int(SIZECBI,8)+SIZECB)
               SIZECBI    =  2 * NCB + SIZEHEADER     
             ELSE 
               SIZECBI    = 2 * NCB +  9 
               SBUFR_CB   = max(SBUFR_CB, 
     &                      min(SIZECB,CBMAXR) + int(SIZECBI,8))
               MAXTEMPCB  = max(MAXTEMPCB, min(SIZECB,CBMAXR)) 
               SIZECBI    =  2 * NCB + SIZEHEADER 
               MAXITEMPCB = max(MAXITEMPCB, SIZECBI)
               SIZECBI     = 0
               SIZECB      = 0_8
             ENDIF
          ELSE 
             SIZECB = SIZECB_SLAVE
             MAXTEMPCB  = max(MAXTEMPCB, min(CBMAXR,SIZECB) )
             MAXITEMPCB = max(MAXITEMPCB,NBROWMAX+NCB+SIZEHEADER)
             IF (.NOT. 
     &        (UPDATE.AND.(.NOT.MASTER).AND.(NSLAVES_LOC.EQ.1))
     &          ) 
     &       SBUFR_CB = max(SBUFR_CB, 
     &            min(CBMAXR,SIZECB) + int(NBROWMAX + NCB + 6,8))
             IF (MASTER) THEN
              SIZECBI     =  NCB + 5 +  SLAVEF - 1 + XSIZE_IC
              SIZECB  = 0_8
             ELSE IF (UPDATE) THEN
              SIZECBI      =  NFR + 6 + SLAVEF - 1 + XSIZE_IC
              IF (KEEP(50).EQ.0) THEN
                SIZECBI = SIZECBI + NBROWMAX + NFR + 
     &                    SIZEHEADER
              ELSE
                SIZECBI = SIZECBI + NBROWMAX + NFR +
     &                    SIZEHEADER+ NSLAVES_LOC
              ENDIF
             ELSE
              SIZECB      = 0_8
              SIZECBI     = 0
             ENDIF
          ENDIF
         ELSE
           IF (LEVELF.NE.3) THEN
               STACKCB     = .TRUE.
               SIZECB      = 0_8
               SIZECBI     = 0
               IF ( (LEVEL.EQ.1) .AND. (LEVELF.NE.1) ) THEN
                  IF (COMPRESSCB) THEN 
                      SIZECB  = (NCB8*(NCB8+1_8))/2_8
                  ELSE
                      SIZECB  = NCB8*NCB8
                  ENDIF
                  SIZECBI     = 2 * NCB + SIZEHEADER
               ELSE IF (LEVEL.EQ.2) THEN
                 IF (MASTER) THEN
                   SIZECBI     =  NCB + 5 +  SLAVEF - 1 + XSIZE_IC
                 ELSE 
                   SIZECB  = SIZECB_SLAVE
                   SIZECBI = SIZECBI + NBROWMAX + NFR + SIZEHEADER
                 ENDIF 
               ENDIF
           ENDIF
         ENDIF
         IF (STACKCB) THEN
           IF (FRERE(STEP(INODE)).EQ.0) THEN
                  write(*,*) ' ERROR 3 in ZMUMPS_ANA_DISTM'
                  CALL MUMPS_ABORT()
           ENDIF
           ITOP = ITOP + 1
           IF ( ITOP .GT. NSTEPS ) THEN
             WRITE(*,*) 'ERROR 4 in ZMUMPS_ANA_DISTM '
           ENDIF
           LSTKI(ITOP) = SIZECBI
           ISTKI=ISTKI + SIZECBI
           ISTKI_OOC = ISTKI_OOC + SIZECBI + (XSIZE_OOC-XSIZE_IC)
           LSTKR(ITOP) = SIZECB
           ISTKR = ISTKR + LSTKR(ITOP)
           NRLNEC = max(NRLNEC,NRLADU+ISTKR+MAXTEMPCB)
           NIRNEC = max0(NIRNEC,NIRADU+ISTKI+MAXITEMPCB)
           NIRNEC_OOC = max0(NIRNEC_OOC,NIRADU_OOC+ISTKI_OOC+
     &          MAXITEMPCB + 
     &          (XSIZE_OOC-XSIZE_IC) ) 
         ENDIF 
        ENDIF 
         TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
         IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN
            INODE = IFATH 
            GOTO 95
         ELSE
            GOTO 90
         ENDIF
      ENDIF 
 115  CONTINUE
      BLOCKING_RHS = KEEP(84)
      IF (KEEP(84).EQ.0) BLOCKING_RHS=1
      NRLNEC = max(NRLNEC, 
     &         NRLADU+int(4*KEEP(127)*abs(BLOCKING_RHS),8))
      IF (BLOCKING_RHS .LT. 0) THEN
        BLOCKING_RHS = - 2 * BLOCKING_RHS
      ENDIF
      NRLNEC_ACTIVE = max(NRLNEC_ACTIVE, MAX_SIZE_FACTOR+
     &                    int(4*KEEP(127)*BLOCKING_RHS,8))
      SBUF_RECOLD = max(int(SBUFR,8),SBUFR_CB)
      SBUF_RECOLD = max(SBUF_RECOLD,
     &        MAXTEMPCB+int(MAXITEMPCB,8)) + 10_8
      SBUF_REC = max(SBUFR, int(min(100000_8,SBUFR_CB)))
      SBUF_REC = SBUF_REC   + 17
      SBUF_REC = SBUF_REC + 2 * KEEP(127) + SLAVEF - 1 + 7
      SBUF_SEND = max(SBUFS, int(min(100000_8,SBUFR_CB)))
      SBUF_SEND = SBUF_SEND + 17 
      IF(KEEP(219).NE.0.AND.KEEP(50) .EQ. 2) THEN
         SBUF_RECOLD = SBUF_RECOLD+int(KEEP(108)+1,8)
         SBUF_REC = SBUF_REC+KEEP(108)+1
         SBUF_SEND = SBUF_SEND+KEEP(108)+1
      ENDIF
      IF (SLAVEF.EQ.1) THEN 
         SBUF_RECOLD = 1_8
         SBUF_REC = 1
         SBUF_SEND= 1
      ENDIF
      DEALLOCATE( LSTKR, TNSTK, IPOOL,
     &          LSTKI )
      OPS_SUBTREE = dble(OPS_SBTR_LOC)
      OPSA        = dble(OPSA_LOC)
      KEEP(66)    = int(OPSA_LOC/1000000.d0)
      RETURN
      END SUBROUTINE ZMUMPS_ANA_DISTM
