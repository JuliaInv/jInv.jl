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
      SUBROUTINE SMUMPS_STOP(id,OUTFILE)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      CHARACTER(len=*), intent(in) :: OUTFILE
      INTEGER::i1,i2,NBVARIABLES,NBVARIABLES_ROOT
      INTEGER::OUT,err,OUTROOT
      CHARACTER(len=100), allocatable, dimension(:)::VARIABLES
      CHARACTER(len=100), allocatable, dimension(:)::VARIABLES_ROOT
      CHARACTER(len=100):: OUTFILEMAIN,OUTFILEROOT
      CHARACTER(len=3) :: STRING_ID
      LOGICAL :: error
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      TYPE (SMUMPS_STRUC) :: id
      NBVARIABLES=172
      allocate(VARIABLES(NBVARIABLES))
      VARIABLES(172)="NB_SINGULAR_VALUES"
      VARIABLES(171)="SINGULAR_VALUES"
      VARIABLES(170)="IF_RESTARTING"
      VARIABLES(169)="L0_OMP_MAPPING"
      VARIABLES(168)="PTR_LEAFS_L0_OMP"
      VARIABLES(167)="PERM_L0_OMP"
      VARIABLES(166)="VIRT_L0_OMP"
      VARIABLES(165)="PHYS_L0_OMP"
      VARIABLES(164)="IPOOL_AFTER_L0_OMP"
      VARIABLES(163)="IPOOL_BEFORE_L0_OMP"
      VARIABLES(162)="THREAD_LA"
      VARIABLES(161)="LL0_OMP_MAPPING"
      VARIABLES(160)="L_VIRT_L0_OMP"
      VARIABLES(159)="L_PHYS_L0_OMP"
      VARIABLES(158)="LPOOL_BEFORE_L0_OMP"
      VARIABLES(157)="LPOOL_AFTER_L0_OMP"
      VARIABLES(156)="NBGRP"
      VARIABLES(155)="LRGROUPS"
      VARIABLES(154)="root"
      VARIABLES(153)="WORKING"
      VARIABLES(152)="IPTR_WORKING"
      VARIABLES(151)="pad14"
      VARIABLES(150)="SUP_PROC"
      VARIABLES(149)="PIVNUL_LIST"
      VARIABLES(148)="OOC_FILE_NAME_LENGTH"
      VARIABLES(147)="OOC_FILE_NAMES"
      VARIABLES(146)="OOC_NB_FILE_TYPE"
      VARIABLES(145)="OOC_NB_FILES"
      VARIABLES(144)="OOC_TOTAL_NB_NODES"
      VARIABLES(143)="OOC_VADDR"
      VARIABLES(142)="OOC_SIZE_OF_BLOCK"
      VARIABLES(141)="pad13"
      VARIABLES(140)="OOC_INODE_SEQUENCE"
      VARIABLES(139)="OOC_MAX_NB_NODES_FOR_ZONE"
      VARIABLES(138)="INSTANCE_NUMBER"
      VARIABLES(137)="pad12"
      VARIABLES(136)="CB_SON_SIZE"
      VARIABLES(135)="DKEEP"
      VARIABLES(134)="LWK_USER"
      VARIABLES(133)="NBSA_LOCAL"
      VARIABLES(132)="WK_USER"
      VARIABLES(131)="CROIX_MANU"
      VARIABLES(130)="SCHED_SBTR"
      VARIABLES(129)="SCHED_GRP"
      VARIABLES(128)="SCHED_DEP"
      VARIABLES(127)="SBTR_ID"
      VARIABLES(126)="DEPTH_FIRST_SEQ"
      VARIABLES(125)="DEPTH_FIRST"
      VARIABLES(124)="MY_NB_LEAF"
      VARIABLES(123)="MY_FIRST_LEAF"
      VARIABLES(122)="MY_ROOT_SBTR"
      VARIABLES(121)="COST_TRAV"
      VARIABLES(120)="MEM_SUBTREE"
      VARIABLES(119)="RHSCOMP"
      VARIABLES(118)="pad111"
      VARIABLES(117)="POSINRHSCOMP_COL_ALLOC"
      VARIABLES(116)="POSINRHSCOMP_COL"
      VARIABLES(115)="POSINRHSCOMP_ROW"
      VARIABLES(114)="MEM_DIST"
      VARIABLES(113)="I_AM_CAND"
      VARIABLES(112)="TAB_POS_IN_PERE"
      VARIABLES(111)="FUTURE_NIV2"
      VARIABLES(110)="ISTEP_TO_INIV2"
      VARIABLES(109)="CANDIDATES"
      VARIABLES(108)="ELTPROC"
      VARIABLES(107)="pad11"
      VARIABLES(106)="NA_ELT"
      VARIABLES(105)="LELTVAR"
      VARIABLES(104)="NELT_loc"
      VARIABLES(103)="DBLARR"
      VARIABLES(102)="INTARR"
      VARIABLES(101)="PROCNODE"
      VARIABLES(100)="S"
      VARIABLES(99)="PTRFAC"
      VARIABLES(98)="PTLUST_S"
      VARIABLES(97)="PROCNODE_STEPS"
      VARIABLES(96)="NA"
      VARIABLES(95)="FRTELT"
      VARIABLES(94)="FRTPTR"
      VARIABLES(93)="PTRAR"
      VARIABLES(92)="FILS"
      VARIABLES(91)="DAD_STEPS"
      VARIABLES(90)="FRERE_STEPS"
      VARIABLES(89)="Step2node"
      VARIABLES(88)="ND_STEPS"
      VARIABLES(87)="NE_STEPS"
      VARIABLES(86)="STEP"
      VARIABLES(85)="NBSA"
      VARIABLES(84)="LNA"
      VARIABLES(83)="KEEP"
      VARIABLES(82)="Deficiency"
      VARIABLES(81)="MAXIS1"
      VARIABLES(80)="IS1"
      VARIABLES(79)="IS"
      VARIABLES(78)="BUFR"
      VARIABLES(77)="POIDS"
      VARIABLES(76)="LBUFR_BYTES"
      VARIABLES(75)="LBUFR"
      VARIABLES(74)="ASS_IRECV"
      VARIABLES(73)="NSLAVES"
      VARIABLES(72)="NPROCS"
      VARIABLES(71)="MYID"
      VARIABLES(70)="COMM_LOAD"
      VARIABLES(69)="MYID_NODES"
      VARIABLES(68)="COMM_NODES"
      VARIABLES(67)="INST_Number"
      VARIABLES(66)="MAX_SURF_MASTER"
      VARIABLES(65)="KEEP8"
      VARIABLES(64)="pad8"
      VARIABLES(63)="WRITE_PROBLEM"
      VARIABLES(62)="OOC_PREFIX"
      VARIABLES(61)="OOC_TMPDIR"
      VARIABLES(60)="VERSION_NUMBER"
      VARIABLES(59)="MAPPING"
      VARIABLES(58)="LISTVAR_SCHUR"
      VARIABLES(57)="SCHUR_CINTERFACE"
      VARIABLES(56)="SCHUR"
      VARIABLES(55)="SIZE_SCHUR"
      VARIABLES(54)="SCHUR_LLD"
      VARIABLES(53)="SCHUR_NLOC"
      VARIABLES(52)="SCHUR_MLOC"
      VARIABLES(51)="NBLOCK"
      VARIABLES(50)="MBLOCK"
      VARIABLES(49)="NPCOL"
      VARIABLES(48)="NPROW"
      VARIABLES(47)="UNS_PERM"
      VARIABLES(46)="SYM_PERM"
      VARIABLES(45)="RINFOG"
      VARIABLES(44)="RINFO"
      VARIABLES(43)="CNTL"
      VARIABLES(42)="COST_SUBTREES"
      VARIABLES(41)="INFOG"
      VARIABLES(40)="INFO"
      VARIABLES(39)="ICNTL"
      VARIABLES(38)="pad5"
      VARIABLES(37)="LREDRHS"
      VARIABLES(36)="LSOL_loc"
      VARIABLES(35)="NZ_RHS"
      VARIABLES(34)="NRHS"
      VARIABLES(33)="LRHS"
      VARIABLES(32)="ISOL_loc"
      VARIABLES(31)="IRHS_PTR"
      VARIABLES(30)="IRHS_SPARSE"
      VARIABLES(29)="SOL_loc"
      VARIABLES(28)="RHS_SPARSE"
      VARIABLES(27)="REDRHS"
      VARIABLES(26)="RHS"
      VARIABLES(25)="PERM_IN"
      VARIABLES(24)="pad4"
      VARIABLES(23)="A_ELT"
      VARIABLES(22)="ELTVAR"
      VARIABLES(21)="ELTPTR"
      VARIABLES(20)="pad3"
      VARIABLES(19)="NELT"
      VARIABLES(18)="pad2"
      VARIABLES(17)="A_loc"
      VARIABLES(16)="JCN_loc"
      VARIABLES(15)="IRN_loc"
      VARIABLES(14)="pad1"
      VARIABLES(13)="NZ_loc"
      VARIABLES(12)="pad0"
      VARIABLES(11)="ROWSCA"
      VARIABLES(10)="COLSCA"
      VARIABLES(9)="JCN"
      VARIABLES(8)="IRN"
      VARIABLES(7)="A"
      VARIABLES(6)="NZ"
      VARIABLES(5)="N"
      VARIABLES(4)="JOB"
      VARIABLES(3)="PAR"
      VARIABLES(2)="SYM"
      VARIABLES(1)="COMM"
      NBVARIABLES_ROOT=34
      allocate(VARIABLES_ROOT(NBVARIABLES_ROOT))
      VARIABLES_ROOT(34)="NB_SINGULAR_VALUES"
      VARIABLES_ROOT(33)="SINGULAR_VALUES"
      VARIABLES_ROOT(32)="SVD_VT"
      VARIABLES_ROOT(31)="SVD_U"
      VARIABLES_ROOT(30)="gridinit_done"
      VARIABLES_ROOT(29)="yes"
      VARIABLES_ROOT(28)="rootpad3"
      VARIABLES_ROOT(27)="QR_RCOND"
      VARIABLES_ROOT(26)="rootpad"
      VARIABLES_ROOT(25)="RHS_ROOT"
      VARIABLES_ROOT(24)="rootpad2"
      VARIABLES_ROOT(23)="QR_TAU"
      VARIABLES_ROOT(22)="SCHUR_POINTER"
      VARIABLES_ROOT(21)="RHS_CNTR_MASTER_ROOT"
      VARIABLES_ROOT(20)="rootpad1"
      VARIABLES_ROOT(19)="IPIV"
      VARIABLES_ROOT(18)="RG2L_COL"
      VARIABLES_ROOT(17)="RG2L_ROW"
      VARIABLES_ROOT(16)="rootpad0"
      VARIABLES_ROOT(15)="LPIV"
      VARIABLES_ROOT(14)="CNTXT_BLACS"
      VARIABLES_ROOT(13)="DESCRIPTOR"
      VARIABLES_ROOT(12)="TOT_ROOT_SIZE"
      VARIABLES_ROOT(11)="ROOT_SIZE"
      VARIABLES_ROOT(10)="RHS_NLOC"
      VARIABLES_ROOT(9)="SCHUR_LLD"
      VARIABLES_ROOT(8)="SCHUR_NLOC"
      VARIABLES_ROOT(7)="SCHUR_MLOC"
      VARIABLES_ROOT(6)="MYCOL"
      VARIABLES_ROOT(5)="MYROW"
      VARIABLES_ROOT(4)="NPCOL"
      VARIABLES_ROOT(3)="NPROW"
      VARIABLES_ROOT(2)="NBLOCK"
      VARIABLES_ROOT(1)="MBLOCK"
      if(((id%ICNTL(3).GT.0).AND.(id%MYID .EQ. MASTER))) then
         write(*,*) "DUMPING MUMPS STRUCTURE IN FILE:",OUTFILE
      endif
      error=.false.
      write (STRING_ID, '(i3)') id%MYID
      OUT=42+id%MYID
      OUTFILEMAIN=trim(adjustl(OUTFILE)) // trim(adjustl(STRING_ID))
      open(UNIT=OUT,FILE=OUTFILEMAIN,STATUS='replace',
     &     form='unformatted',iostat=err)
      if(err.ne.0) THEN
         id%INFOG(1)=-91
         write(*,*) "IN SMUMPS_STOP CANNOT OPEN FILE: "
     &        //trim(adjustl(OUTFILEMAIN))
         error=.true.
         goto 100
      endif      
      DO i1=1,NBVARIABLES 
         SELECT CASE(trim(adjustl(VARIABLES(i1))))
         CASE("COMM") 
         CASE("SYM")
         CASE("PAR")
         CASE("JOB")
            write(OUT) id%JOB
         CASE("N")
            write(OUT) id%N
         CASE("ICNTL")
            write(OUT) id%ICNTL
         CASE("INFO")
            write(OUT) id%INFO
         CASE("INFOG")
            write(OUT) id%INFOG
         CASE("COST_SUBTREES")
            write(OUT) id%COST_SUBTREES
         CASE("CNTL")
            write(OUT) id%CNTL
         CASE("RINFO")
            write(OUT) id%RINFO
         CASE("RINFOG")
            write(OUT) id%RINFOG
         CASE("KEEP8")
            write(OUT) id%KEEP8
         CASE("KEEP")
            write(OUT) id%KEEP
         CASE("DKEEP")
            write(OUT) id%DKEEP    
         CASE("NZ")
            write(OUT) id%NZ
         CASE("A")
            IF(associated(id%A)) THEN
               write(OUT) size(id%A,1)
               write(OUT) id%A
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("IRN")
            IF(associated(id%IRN)) THEN
               write(OUT) size(id%IRN,1)
               write(OUT) id%IRN
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("JCN")
            IF(associated(id%JCN)) THEN
               write(OUT) size(id%JCN,1)
               write(OUT) id%JCN
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("COLSCA")
            IF(associated(id%COLSCA)) THEN
               write(OUT) size(id%COLSCA,1)
               write(OUT) id%COLSCA
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF 
         CASE("ROWSCA")
            IF(associated(id%ROWSCA)) THEN
               write(OUT) size(id%ROWSCA,1)
               write(OUT) id%ROWSCA
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("NZ_loc")
            write(OUT) id%NZ_loc
         CASE("IRN_loc")
            IF(associated(id%IRN_loc)) THEN
               write(OUT) size(id%IRN_loc,1)
               write(OUT) id%IRN_loc
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("JCN_loc")
            IF(associated(id%JCN_loc)) THEN
               write(OUT) size(id%JCN_loc,1)
               write(OUT) id%JCN_loc
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("A_loc")
            IF(associated(id%A_loc)) THEN
               write(OUT) size(id%A_loc,1)
               write(OUT) id%A_loc
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("NELT")
            write(OUT) id%NELT
         CASE("ELTPTR")
            IF(associated(id%ELTPTR)) THEN
               write(OUT) size(id%ELTPTR,1)
               write(OUT) id%ELTPTR
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("ELTVAR")
            IF(associated(id%ELTVAR)) THEN
               write(OUT) size(id%ELTVAR,1)
               write(OUT) id%ELTVAR
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("A_ELT")
            IF(associated(id%A_ELT)) THEN
               write(OUT) size(id%A_ELT,1)
               write(OUT) id%A_ELT
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("PERM_IN")
            IF(associated(id%PERM_IN)) THEN
               write(OUT) size(id%PERM_IN,1)
               write(OUT) id%PERM_IN
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("RHS")
            IF(associated(id%RHS)) THEN
               write(OUT) size(id%RHS,1)
               write(OUT) id%RHS
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("REDRHS")
            IF(associated(id%REDRHS)) THEN
               write(OUT) size(id%REDRHS,1)
               write(OUT) id%REDRHS
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("RHS_SPARSE")
            IF(associated(id%RHS_SPARSE)) THEN
               write(OUT) size(id%RHS_SPARSE,1)
               write(OUT) id%RHS_SPARSE
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("SOL_loc")
            IF(associated(id%SOL_loc)) THEN
               write(OUT) size(id%SOL_loc,1)
               write(OUT) id%SOL_loc
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("IRHS_SPARSE")
            IF(associated(id%IRHS_SPARSE)) THEN
               write(OUT) size(id%IRHS_SPARSE,1)
               write(OUT) id%IRHS_SPARSE
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("IRHS_PTR")
            IF(associated(id%IRHS_PTR)) THEN
               write(OUT) size(id%IRHS_PTR,1)
               write(OUT) id%IRHS_PTR
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("ISOL_loc")
            IF(associated(id%ISOL_loc)) THEN
               write(OUT) size(id%ISOL_loc,1)
               write(OUT) id%ISOL_loc
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("LRHS")
            write(OUT) id%LRHS
         CASE("NRHS")
            write(OUT) id%NRHS
         CASE("NZ_RHS")
            write(OUT) id%NZ_RHS
         CASE("LSOL_loc")
            write(OUT) id%LSOL_loc
         CASE("LREDRHS")
            write(OUT) id%LREDRHS
         CASE("SYM_PERM")
            IF(associated(id%SYM_PERM)) THEN
               write(OUT) size(id%SYM_PERM,1)
               write(OUT) id%SYM_PERM
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("UNS_PERM")
            IF(associated(id%UNS_PERM)) THEN
               write(OUT) size(id%UNS_PERM,1)
               write(OUT) id%UNS_PERM
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("NPROW")
            write(OUT) id%NPROW
         CASE("NPCOL")
            write(OUT) id%NPCOL
         CASE("MBLOCK")
            write(OUT) id%MBLOCK
         CASE("NBLOCK")
            write(OUT) id%NBLOCK
         CASE("SCHUR_MLOC")
            write(OUT) id%SCHUR_MLOC
         CASE("SCHUR_NLOC")
            write(OUT) id%SCHUR_NLOC
         CASE("SCHUR_LLD")
            write(OUT) id%SCHUR_LLD
         CASE("SIZE_SCHUR")
            write(OUT) id%SIZE_SCHUR
         CASE("SCHUR")
            IF(associated(id%SCHUR)) THEN
               write(OUT) size(id%SCHUR,1)
               write(OUT) id%SCHUR
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("SCHUR_CINTERFACE")
            IF(associated(id%SCHUR_CINTERFACE)) THEN
               write(OUT) size(id%SCHUR_CINTERFACE,1)
               write(OUT) id%SCHUR_CINTERFACE
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("LISTVAR_SCHUR")
            IF(associated(id%LISTVAR_SCHUR)) THEN
               write(OUT) size(id%LISTVAR_SCHUR,1)
               write(OUT) id%LISTVAR_SCHUR
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("MAPPING")
            IF(associated(id%MAPPING)) THEN
               write(OUT) size(id%MAPPING,1)
               write(OUT) id%MAPPING
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("VERSION_NUMBER")
            write(OUT) id%VERSION_NUMBER
         CASE("OOC_TMPDIR")
            write(OUT) id%OOC_TMPDIR
         CASE("OOC_PREFIX")
            write(OUT) id%OOC_PREFIX
         CASE("WRITE_PROBLEM")
            write(OUT) id%WRITE_PROBLEM
         CASE("MAX_SURF_MASTER")
            write(OUT) id%MAX_SURF_MASTER
         CASE("INST_Number")
            write(OUT) id%INST_Number
         CASE("COMM_NODES")
            write(OUT) id%COMM_NODES
         CASE("MYID_NODES")
            write(OUT) id%MYID_NODES
         CASE("COMM_LOAD")
            write(OUT) id%COMM_LOAD
         CASE("MYID")
            write(OUT) id%MYID
         CASE("NPROCS")
            write(OUT) id%NPROCS
         CASE("NSLAVES")
            write(OUT) id%NSLAVES
         CASE("ASS_IRECV")
            write(OUT) id%ASS_IRECV
         CASE("LBUFR")
            write(OUT) id%LBUFR
         CASE("LBUFR_BYTES")
            write(OUT) id%LBUFR_BYTES
         CASE("POIDS")
            IF(associated(id%POIDS)) THEN 
               write(OUT) size(id%POIDS,1)
               write(OUT) id%POIDS
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("BUFR")
            IF(associated(id%BUFR)) THEN 
               write(OUT) size(id%BUFR,1)
               write(OUT) id%BUFR
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("IS")
            IF(associated(id%IS)) THEN 
               write(OUT) size(id%IS,1)
               write(OUT) id%IS
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("IS1")
            IF(associated(id%IS1)) THEN 
               write(OUT) size(id%IS1,1)
               write(OUT) id%IS1
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("MAXIS1")
            write(OUT) id%MAXIS1
         CASE("Deficiency")
            write(OUT) id%Deficiency
         CASE("LNA")
            write(OUT) id%LNA
         CASE("NBSA")
            write(OUT) id%NBSA
         CASE("STEP")
            IF(associated(id%STEP)) THEN
               write(OUT) size(id%STEP,1)
               write(OUT) id%STEP
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("NE_STEPS")
            IF(associated(id%NE_STEPS)) THEN
               write(OUT) size(id%NE_STEPS,1)
               write(OUT) id%NE_STEPS
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("ND_STEPS")
            IF(associated(id%ND_STEPS)) THEN
               write(OUT) size(id%ND_STEPS,1)
               write(OUT) id%ND_STEPS
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("Step2node")
            IF(associated(id%Step2node)) THEN
               write(OUT) size(id%Step2node,1)
               write(OUT) id%Step2node
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("FRERE_STEPS")
            IF(associated(id%FRERE_STEPS)) THEN
               write(OUT) size(id%FRERE_STEPS,1)
               write(OUT) id%FRERE_STEPS
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("DAD_STEPS")
            IF(associated(id%DAD_STEPS)) THEN
               write(OUT) size(id%DAD_STEPS,1)
               write(OUT) id%DAD_STEPS
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("FILS")
            IF(associated(id%FILS)) THEN
               write(OUT) size(id%FILS,1)
               write(OUT) id%FILS
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("PTRAR")
            IF(associated(id%PTRAR)) THEN
               write(OUT) size(id%PTRAR,1)
               write(OUT) id%PTRAR
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("FRTPTR")
            IF(associated(id%FRTPTR)) THEN
               write(OUT) size(id%FRTPTR,1)
               write(OUT) id%FRTPTR
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("FRTELT")
            IF(associated(id%FRTELT)) THEN
               write(OUT) size(id%FRTELT,1)
               write(OUT) id%FRTELT
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("NA")
            IF(associated(id%NA)) THEN
               write(OUT) size(id%NA,1)
               write(OUT) id%NA
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("PROCNODE_STEPS")
            IF(associated(id%PROCNODE_STEPS)) THEN
               write(OUT) size(id%PROCNODE_STEPS,1)
               write(OUT) id%PROCNODE_STEPS
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("PTLUST_S")
            IF(associated(id%PTLUST_S)) THEN
               write(OUT) size(id%PTLUST_S,1)
               write(OUT) id%PTLUST_S
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("PTRFAC")
            IF(associated(id%PTRFAC)) THEN
               write(OUT) size(id%PTRFAC,1)
               write(OUT) id%PTRFAC
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("S")
            IF(associated(id%S)) THEN
               write(OUT) id%KEEP8(23)
               write(OUT) id%S
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("PROCNODE")
            IF(associated(id%PROCNODE)) THEN
               write(OUT) size(id%PROCNODE,1)
               write(OUT) id%PROCNODE
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("INTARR")
            IF(associated(id%INTARR)) THEN
               write(OUT) size(id%INTARR,1)
               write(OUT) id%INTARR
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("DBLARR")
            IF(associated(id%DBLARR)) THEN
               write(OUT) size(id%DBLARR,1)
               write(OUT) id%DBLARR
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("NELT_loc")
            write(OUT) id%NELT_loc
         CASE("LELTVAR")
            write(OUT) id%LELTVAR
         CASE("NA_ELT")
            write(OUT) id%NA_ELT
         CASE("ELTPROC")
            IF(associated(id%ELTPROC)) THEN
               write(OUT) size(id%ELTPROC,1)
               write(OUT) id%ELTPROC
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("CANDIDATES")
            IF(associated(id%CANDIDATES)) THEN
               write(OUT) size(id%CANDIDATES,1)
     &              ,size(id%CANDIDATES,2)
               write(OUT) id%CANDIDATES
            ELSE
               write(OUT) -999,-998
               write(OUT) -999
            ENDIF
         CASE("ISTEP_TO_INIV2")
            IF(associated(id%ISTEP_TO_INIV2)) THEN
               write(OUT) size(id%ISTEP_TO_INIV2,1)
               write(OUT) id%ISTEP_TO_INIV2
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("FUTURE_NIV2")
            IF(associated(id%FUTURE_NIV2)) THEN
               write(OUT) size(id%FUTURE_NIV2,1)
               write(OUT) id%FUTURE_NIV2
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("TAB_POS_IN_PERE")
            IF(associated(id%TAB_POS_IN_PERE)) THEN
               write(OUT) size(id%TAB_POS_IN_PERE,1)
     &              ,size(id%TAB_POS_IN_PERE,2)
               write(OUT) id%TAB_POS_IN_PERE
            ELSE
               write(OUT) -999,-998
               write(OUT) -999
            ENDIF
         CASE("I_AM_CAND")
            IF(associated(id%I_AM_CAND)) THEN
               write(OUT) size(id%I_AM_CAND,1)
               write(OUT) id%I_AM_CAND
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("MEM_DIST")
            IF(associated(id%MEM_DIST)) THEN
               write(OUT) size(id%MEM_DIST,1)
               write(OUT) id%MEM_DIST
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("POSINRHSCOMP_ROW")
            IF(associated(id%POSINRHSCOMP_ROW)) THEN 
               write(OUT) size(id%POSINRHSCOMP_ROW,1)
               write(OUT) id%POSINRHSCOMP_ROW
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("POSINRHSCOMP_COL")
            IF(associated(id%POSINRHSCOMP_COL)) THEN
               write(OUT) size(id%POSINRHSCOMP_COL,1)
               write(OUT) id%POSINRHSCOMP_COL
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("POSINRHSCOMP_COL_ALLOC")
            write(OUT) id%POSINRHSCOMP_COL_ALLOC
         CASE("RHSCOMP")
            IF(associated(id%RHSCOMP)) THEN
               write(OUT) size(id%RHSCOMP,1)
               write(OUT) id%RHSCOMP
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("MEM_SUBTREE")
            IF(associated(id%MEM_SUBTREE)) THEN
               write(OUT) size(id%MEM_SUBTREE,1)
               write(OUT) id%MEM_SUBTREE
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("COST_TRAV")
            IF(associated(id%COST_TRAV)) THEN
               write(OUT) size(id%COST_TRAV,1)
               write(OUT) id%COST_TRAV
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("MY_ROOT_SBTR")
            IF(associated(id%MY_ROOT_SBTR)) THEN
               write(OUT) size(id%MY_ROOT_SBTR,1)
               write(OUT) id%MY_ROOT_SBTR
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("MY_FIRST_LEAF")
            IF(associated(id%MY_FIRST_LEAF)) THEN
               write(OUT) size(id%MY_FIRST_LEAF,1)
               write(OUT) id%MY_FIRST_LEAF
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("MY_NB_LEAF")
            IF(associated(id%MY_NB_LEAF)) THEN
               write(OUT) size(id%MY_NB_LEAF,1)
               write(OUT) id%MY_NB_LEAF
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("DEPTH_FIRST")
            IF(associated(id%DEPTH_FIRST)) THEN
               write(OUT) size(id%DEPTH_FIRST,1)
               write(OUT) id%DEPTH_FIRST
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("DEPTH_FIRST_SEQ")
            IF(associated(id%DEPTH_FIRST_SEQ)) THEN
               write(OUT) size(id%DEPTH_FIRST_SEQ,1)
               write(OUT) id%DEPTH_FIRST_SEQ
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("SBTR_ID")
            IF(associated(id%SBTR_ID)) THEN
               write(OUT) size(id%SBTR_ID,1)
               write(OUT) id%SBTR_ID
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("SCHED_DEP")
            IF(associated(id%SCHED_DEP)) THEN
               write(OUT) size(id%SCHED_DEP,1)
               write(OUT) id%SCHED_DEP
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("SCHED_GRP")
            IF(associated(id%SCHED_GRP)) THEN
               write(OUT) size(id%SCHED_GRP,1)
               write(OUT) id%SCHED_GRP
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("CROIX_MANU")
            IF(associated(id%CROIX_MANU)) THEN
               write(OUT) size(id%CROIX_MANU,1)
               write(OUT) id%CROIX_MANU
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("WK_USER")
            IF(associated(id%WK_USER)) THEN
               write(OUT) id%KEEP8(24)
               write(OUT) id%WK_USER
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("NBSA_LOCAL")
            write(OUT) id%NBSA_LOCAL
         CASE("LWK_USER")
            write(OUT) id%LWK_USER
         CASE("CB_SON_SIZE")
            IF(associated(id%CB_SON_SIZE)) THEN
               write(OUT) size(id%CB_SON_SIZE,1)
               write(OUT) id%CB_SON_SIZE
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("INSTANCE_NUMBER")
            write(OUT) id%INSTANCE_NUMBER
         CASE("OOC_MAX_NB_NODES_FOR_ZONE")
            write(OUT) id%OOC_MAX_NB_NODES_FOR_ZONE
         CASE("OOC_INODE_SEQUENCE")
            IF(associated(id%OOC_INODE_SEQUENCE)) THEN
               write(OUT) size(id%OOC_INODE_SEQUENCE,1)
     &              ,size(id%OOC_INODE_SEQUENCE,2)
               write(OUT) id%OOC_INODE_SEQUENCE
            ELSE
               write(OUT) -999,-998
               write(OUT) -999
            ENDIF
         CASE("OOC_SIZE_OF_BLOCK")
            IF(associated(id%OOC_SIZE_OF_BLOCK)) THEN
               write(OUT) size(id%OOC_SIZE_OF_BLOCK,1)
     &              ,size(id%OOC_SIZE_OF_BLOCK,2)  
               write(OUT) id%OOC_SIZE_OF_BLOCK
            ELSE
               write(OUT) -999,-998
               write(OUT) -999
            ENDIF
         CASE("OOC_VADDR")
            IF(associated(id%OOC_VADDR)) THEN
               write(OUT) size(id%OOC_VADDR,1),size(id%OOC_VADDR,2)
               write(OUT) id%OOC_VADDR
            ELSE
               write(OUT) -999,-998
               write(OUT) -999
            ENDIF
         CASE("OOC_TOTAL_NB_NODES")
            IF(associated(id%OOC_TOTAL_NB_NODES)) THEN
               write(OUT) size(id%OOC_TOTAL_NB_NODES,1)
               write(OUT) id%OOC_TOTAL_NB_NODES
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("OOC_NB_FILES")
            IF(associated(id%OOC_NB_FILES)) THEN
               write(OUT) size(id%OOC_NB_FILES,1)
               write(OUT) id%OOC_NB_FILES
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("OOC_NB_FILE_TYPE")
            write(OUT) id%OOC_NB_FILE_TYPE
         CASE("OOC_FILE_NAMES")
            IF(associated(id%OOC_FILE_NAMES)) THEN
               write(OUT) size(id%OOC_FILE_NAMES,1)
     &              ,size(id%OOC_FILE_NAMES,2)
               write(OUT) id%OOC_FILE_NAMES
            ELSE
               write(OUT) -999,-998
               write(OUT) -999
            ENDIF
         CASE("OOC_FILE_NAME_LENGTH")
            IF(associated(id%OOC_FILE_NAME_LENGTH)) THEN
               write(OUT) size(id%OOC_FILE_NAME_LENGTH,1)
               write(OUT) id%OOC_FILE_NAME_LENGTH
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("PIVNUL_LIST")
            IF(associated(id%PIVNUL_LIST)) THEN
               write(OUT) size(id%PIVNUL_LIST,1)
               write(OUT) id%PIVNUL_LIST
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("SUP_PROC")
            IF(associated(id%SUP_PROC)) THEN
               write(OUT) size(id%SUP_PROC,1),size(id%SUP_PROC,2)
               write(OUT) id%SUP_PROC
            ELSE
               write(OUT) -999,-998
               write(OUT) -999
            ENDIF
         CASE("IPTR_WORKING")
            IF(associated(id%IPTR_WORKING)) THEN
               write(OUT) size(id%IPTR_WORKING,1)
               write(OUT) id%IPTR_WORKING
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("WORKING")
            IF(associated(id%WORKING)) THEN
               write(OUT) size(id%WORKING,1)
               write(OUT) id%WORKING
            ELSE
               write(OUT) -999
               write(OUT) -999
            ENDIF
         CASE("root")
         CASE("NBGRP")
           write(OUT) id%NBGRP
        CASE("LRGROUPS")
           IF(associated(id%LRGROUPS)) THEN
              write(OUT) size(id%LRGROUPS,1)
              write(OUT) id%LRGROUPS
           ELSE
              write(OUT) -999
              write(OUT) -999
           ENDIF
        CASE("SCHED_SBTR")
           IF(associated(id%SCHED_SBTR)) THEN
              write(OUT) size(id%SCHED_SBTR,1)
              write(OUT) id%SCHED_SBTR
           ELSE
              write(OUT) -999
              write(OUT) -999
           ENDIF
        CASE("LPOOL_AFTER_L0_OMP")
           write(OUT) id%LPOOL_AFTER_L0_OMP
        CASE("LPOOL_BEFORE_L0_OMP")
           write(OUT) id%LPOOL_BEFORE_L0_OMP
        CASE("L_PHYS_L0_OMP")
           write(OUT) id%L_PHYS_L0_OMP
        CASE("L_VIRT_L0_OMP")        
           write(OUT) id%L_VIRT_L0_OMP
        CASE("LL0_OMP_MAPPING")
           write(OUT) id%LL0_OMP_MAPPING
        CASE("THREAD_LA")
           write(OUT) id%THREAD_LA
        CASE("IPOOL_AFTER_L0_OMP")
           IF(associated(id%IPOOL_AFTER_L0_OMP)) THEN
              write(OUT) size(id%IPOOL_AFTER_L0_OMP,1)
              write(OUT) id%IPOOL_AFTER_L0_OMP
           ELSE
              write(OUT) -999
              write(OUT) -999
           ENDIF
        CASE("IPOOL_BEFORE_L0_OMP")
           IF(associated(id%IPOOL_BEFORE_L0_OMP)) THEN
              write(OUT) size(id%IPOOL_BEFORE_L0_OMP,1)
              write(OUT) id%IPOOL_BEFORE_L0_OMP
           ELSE
              write(OUT) -999
              write(OUT) -999
           ENDIF
        CASE("PHYS_L0_OMP")
           IF(associated(id%PHYS_L0_OMP)) THEN
              write(OUT) size(id%PHYS_L0_OMP,1)
              write(OUT) id%PHYS_L0_OMP
           ELSE
              write(OUT) -999
              write(OUT) -999
           ENDIF
        CASE("VIRT_L0_OMP")
           IF(associated(id%VIRT_L0_OMP)) THEN
              write(OUT) size(id%VIRT_L0_OMP,1)
              write(OUT) id%VIRT_L0_OMP
           ELSE
              write(OUT) -999
              write(OUT) -999
           ENDIF
        CASE("PERM_L0_OMP")
           IF(associated(id%PERM_L0_OMP)) THEN
              write(OUT) size(id%PERM_L0_OMP,1)
              write(OUT) id%PERM_L0_OMP
           ELSE
              write(OUT) -999
              write(OUT) -999
           ENDIF
        CASE("PTR_LEAFS_L0_OMP")
           IF(associated(id%PTR_LEAFS_L0_OMP)) THEN
              write(OUT) size(id%PTR_LEAFS_L0_OMP,1)
              write(OUT) id%PTR_LEAFS_L0_OMP
           ELSE
              write(OUT) -999
              write(OUT) -999
           ENDIF
        CASE("L0_OMP_MAPPING")
          IF(associated(id%L0_OMP_MAPPING)) THEN
              write(OUT) size(id%L0_OMP_MAPPING,1)
              write(OUT) id%L0_OMP_MAPPING
           ELSE
              write(OUT) -999
              write(OUT) -999
           ENDIF
        CASE("IF_RESTARTING")
           write(OUT) id%IF_RESTARTING
      CASE("SINGULAR_VALUES")
         IF(associated(id%SINGULAR_VALUES)) THEN
            write(OUT) size(id%SINGULAR_VALUES,1)
            write(OUT) id%SINGULAR_VALUES
         ELSE
            write(OUT) -999,-998
            write(OUT) -999
         ENDIF 
      CASE("NB_SINGULAR_VALUES")
        write(OUT) id%NB_SINGULAR_VALUES
        CASE("pad0","pad1","pad2","pad3","pad4","pad5","pad8","pad11",
     &          "pad111", "pad12","pad13","pad14")
        CASE DEFAULT
           id%INFOG(1)=-92
           write(*,*) "IN SMUMPS_STOP UNKNOWN VARIABLE IN STRUCTURE: "
     &          // trim(adjustl(VARIABLES(i1)))
           error=.true.
           goto 100
        END SELECT
      ENDDO
 100  continue
      CLOSE(OUT)
      if(error) goto 101
      OUTROOT=OUT+1000
      OUTFILEROOT=trim(adjustl(OUTFILE)) // trim(adjustl("ROOT")) 
     &     // trim(adjustl(STRING_ID))
      open(UNIT=OUTROOT,FILE=OUTFILEROOT,STATUS='replace',
     &     form='unformatted',iostat=err)
      if(err.ne.0) THEN
         id%INFOG(1)=-91
         write(*,*) "IN SMUMPS_STOP CANNOT OPEN FILE: "
     &        //trim(adjustl(OUTFILEROOT))
         goto 100
      endif
      DO i2=1,NBVARIABLES_ROOT
         SELECT CASE(VARIABLES_ROOT(i2))
      CASE("MBLOCK")
         write(OUTROOT) id%root%MBLOCK
      CASE("NBLOCK")
         write(OUTROOT) id%root%NBLOCK
      CASE("NPROW")
         write(OUTROOT) id%root%NPROW
      CASE("NPCOL")
         write(OUTROOT) id%root%NPCOL
      CASE("MYROW")
         write(OUTROOT) id%root%MYROW
      CASE("MYCOL")
         write(OUTROOT) id%root%MYCOL
      CASE("SCHUR_MLOC")
         write(OUTROOT) id%root%SCHUR_MLOC
      CASE("SCHUR_NLOC")
         write(OUTROOT) id%root%SCHUR_NLOC
      CASE("SCHUR_LLD")
         write(OUTROOT) id%root%SCHUR_LLD
      CASE("RHS_NLOC")
         write(OUTROOT) id%root%RHS_NLOC
      CASE("ROOT_SIZE")
         write(OUTROOT) id%root%ROOT_SIZE
      CASE("TOT_ROOT_SIZE")
         write(OUTROOT) id%root%TOT_ROOT_SIZE
      CASE("DESCRIPTOR")
         write(OUTROOT) id%root%DESCRIPTOR
      CASE("CNTXT_BLACS")
         write(OUTROOT) id%root%CNTXT_BLACS
      CASE("LPIV")
         write(OUTROOT) id%root%LPIV
      CASE("RG2L_ROW")
         IF(associated(id%root%RG2L_ROW)) THEN
            write(OUTROOT) size(id%root%RG2L_ROW,1)
            write(OUTROOT) id%root%RG2L_ROW
         ELSE
            write(OUTROOT) -999
            write(OUTROOT) -999
         ENDIF
      CASE("RG2L_COL")
         IF(associated(id%root%RG2L_COL)) THEN
            write(OUTROOT) size(id%root%RG2L_COL,1)
            write(OUTROOT) id%root%RG2L_COL
         ELSE
            write(OUTROOT) -999
            write(OUTROOT) -999
         ENDIF
      CASE("IPIV")
         IF(associated(id%root%IPIV)) THEN
            write(OUTROOT) size(id%root%IPIV,1)
            write(OUTROOT) id%root%IPIV
         ELSE
            write(OUTROOT) -999
            write(OUTROOT) -999
         ENDIF
      CASE("RHS_CNTR_MASTER_ROOT")
         IF(associated(id%root%RHS_CNTR_MASTER_ROOT)) THEN
            write(OUTROOT) size(id%root%RHS_CNTR_MASTER_ROOT,1)
            write(OUTROOT) id%root%RHS_CNTR_MASTER_ROOT
         ELSE
            write(OUTROOT) -999
            write(OUTROOT) -999
         ENDIF
      CASE("SCHUR_POINTER")
         IF(associated(id%root%SCHUR_POINTER)) THEN
            write(OUTROOT) size(id%root%SCHUR_POINTER,1)
            write(OUTROOT) id%root%SCHUR_POINTER
         ELSE
            write(OUTROOT) -999
            write(OUTROOT) -999
         ENDIF
      CASE("QR_TAU")
         IF(associated(id%root%QR_TAU)) THEN
            write(OUTROOT) size(id%root%QR_TAU,1)
            write(OUTROOT) id%root%QR_TAU
         ELSE
            write(OUTROOT) -999
            write(OUTROOT) -999
         ENDIF
      CASE("RHS_ROOT")
         IF(associated(id%root%RHS_ROOT)) THEN
            write(OUTROOT) size(id%root%RHS_ROOT,1)
     &           ,size(id%root%RHS_ROOT,2)
            write(OUTROOT) id%root%RHS_ROOT
         ELSE
            write(OUTROOT) -999,-998
            write(OUTROOT) -999
         ENDIF
      CASE("QR_RCOND")
         write(OUTROOT) id%root%QR_RCOND
      CASE("yes")
         write(OUTROOT) id%root%yes
      CASE("gridinit_done")
         write(OUTROOT) id%root%gridinit_done
      CASE("SVD_U")
         IF(associated(id%root%SVD_U)) THEN
            write(OUTROOT) size(id%root%SVD_U,1)
     &           ,size(id%root%SVD_U,2)
            write(OUTROOT) id%root%SVD_U
         ELSE
            write(OUTROOT) -999,-998
            write(OUTROOT) -999
         ENDIF   
      CASE("SVD_VT")
         IF(associated(id%root%SVD_VT)) THEN
            write(OUTROOT) size(id%root%SVD_VT,1)
     &           ,size(id%root%SVD_VT,2)
            write(OUTROOT) id%root%SVD_VT
         ELSE
            write(OUTROOT) -999,-998
            write(OUTROOT) -999
         ENDIF  
      CASE("SINGULAR_VALUES")
         IF(associated(id%root%SINGULAR_VALUES)) THEN
            write(OUTROOT) size(id%root%SINGULAR_VALUES,1)
            write(OUTROOT) id%root%SINGULAR_VALUES
         ELSE
            write(OUTROOT) -999,-998
            write(OUTROOT) -999
         ENDIF 
      CASE("NB_SINGULAR_VALUES")
         write(OUTROOT) id%root%NB_SINGULAR_VALUES
      CASE("rootpad0","rootpad1","rootpad2","rootpad",
     &        "rootpad3")
      CASE DEFAULT
         id%INFOG(1)=-92
         write(*,*) "IN SMUMPS_STOP UNKNOWN VARIABLE"
     &        // "IN ROOT: " // trim(adjustl(VARIABLES_ROOT(i2)))
         goto 101
      END SELECT  
      ENDDO
 101  continue
      CLOSE(OUTROOT)
      id%IF_RESTARTING=.TRUE.
      END SUBROUTINE SMUMPS_STOP
      SUBROUTINE SMUMPS_RESTART(id,INFILE)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER::IN,err,NBVARIABLES,i1,NBVARIABLES_ROOT,dummy
      INTEGER:: i2,size_array1,size_array2,INROOT
      CHARACTER(len=100),allocatable,DIMENSION(:)::VARIABLES
      CHARACTER(len=100),allocatable,DIMENSION(:)::VARIABLES_ROOT
      CHARACTER(len=*), intent(in):: INFILE
      CHARACTER(len=100):: INFILEMAIN,INFILEROOT
      CHARACTER(len=3) :: STRING_ID
      LOGICAL :: error
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      TYPE (SMUMPS_STRUC) :: id
      error=.false.
      write (STRING_ID, '(i3)') id%MYID
      NBVARIABLES=172
      allocate(VARIABLES(NBVARIABLES))
      VARIABLES(172)="NB_SINGULAR_VALUES"
      VARIABLES(171)="SINGULAR_VALUES"
      VARIABLES(170)="IF_RESTARTING"
      VARIABLES(169)="L0_OMP_MAPPING"
      VARIABLES(168)="PTR_LEAFS_L0_OMP"
      VARIABLES(167)="PERM_L0_OMP"
      VARIABLES(166)="VIRT_L0_OMP"
      VARIABLES(165)="PHYS_L0_OMP"
      VARIABLES(164)="IPOOL_AFTER_L0_OMP"
      VARIABLES(163)="IPOOL_BEFORE_L0_OMP"
      VARIABLES(162)="THREAD_LA"
      VARIABLES(161)="LL0_OMP_MAPPING"
      VARIABLES(160)="L_VIRT_L0_OMP"
      VARIABLES(159)="L_PHYS_L0_OMP"
      VARIABLES(158)="LPOOL_BEFORE_L0_OMP"
      VARIABLES(157)="LPOOL_AFTER_L0_OMP"
      VARIABLES(156)="NBGRP"
      VARIABLES(155)="LRGROUPS"
      VARIABLES(154)="root"
      VARIABLES(153)="WORKING"
      VARIABLES(152)="IPTR_WORKING"
      VARIABLES(151)="pad14"
      VARIABLES(150)="SUP_PROC"
      VARIABLES(149)="PIVNUL_LIST"
      VARIABLES(148)="OOC_FILE_NAME_LENGTH"
      VARIABLES(147)="OOC_FILE_NAMES"
      VARIABLES(146)="OOC_NB_FILE_TYPE"
      VARIABLES(145)="OOC_NB_FILES"
      VARIABLES(144)="OOC_TOTAL_NB_NODES"
      VARIABLES(143)="OOC_VADDR"
      VARIABLES(142)="OOC_SIZE_OF_BLOCK"
      VARIABLES(141)="pad13"
      VARIABLES(140)="OOC_INODE_SEQUENCE"
      VARIABLES(139)="OOC_MAX_NB_NODES_FOR_ZONE"
      VARIABLES(138)="INSTANCE_NUMBER"
      VARIABLES(137)="pad12"
      VARIABLES(136)="CB_SON_SIZE"
      VARIABLES(135)="DKEEP"
      VARIABLES(134)="LWK_USER"
      VARIABLES(133)="NBSA_LOCAL"
      VARIABLES(132)="WK_USER"
      VARIABLES(131)="CROIX_MANU"
      VARIABLES(130)="SCHED_SBTR"
      VARIABLES(129)="SCHED_GRP"
      VARIABLES(128)="SCHED_DEP"
      VARIABLES(127)="SBTR_ID"
      VARIABLES(126)="DEPTH_FIRST_SEQ"
      VARIABLES(125)="DEPTH_FIRST"
      VARIABLES(124)="MY_NB_LEAF"
      VARIABLES(123)="MY_FIRST_LEAF"
      VARIABLES(122)="MY_ROOT_SBTR"
      VARIABLES(121)="COST_TRAV"
      VARIABLES(120)="MEM_SUBTREE"
      VARIABLES(119)="RHSCOMP"
      VARIABLES(118)="pad111"
      VARIABLES(117)="POSINRHSCOMP_COL_ALLOC"
      VARIABLES(116)="POSINRHSCOMP_COL"
      VARIABLES(115)="POSINRHSCOMP_ROW"
      VARIABLES(114)="MEM_DIST"
      VARIABLES(113)="I_AM_CAND"
      VARIABLES(112)="TAB_POS_IN_PERE"
      VARIABLES(111)="FUTURE_NIV2"
      VARIABLES(110)="ISTEP_TO_INIV2"
      VARIABLES(109)="CANDIDATES"
      VARIABLES(108)="ELTPROC"
      VARIABLES(107)="pad11"
      VARIABLES(106)="NA_ELT"
      VARIABLES(105)="LELTVAR"
      VARIABLES(104)="NELT_loc"
      VARIABLES(103)="DBLARR"
      VARIABLES(102)="INTARR"
      VARIABLES(101)="PROCNODE"
      VARIABLES(100)="S"
      VARIABLES(99)="PTRFAC"
      VARIABLES(98)="PTLUST_S"
      VARIABLES(97)="PROCNODE_STEPS"
      VARIABLES(96)="NA"
      VARIABLES(95)="FRTELT"
      VARIABLES(94)="FRTPTR"
      VARIABLES(93)="PTRAR"
      VARIABLES(92)="FILS"
      VARIABLES(91)="DAD_STEPS"
      VARIABLES(90)="FRERE_STEPS"
      VARIABLES(89)="Step2node"
      VARIABLES(88)="ND_STEPS"
      VARIABLES(87)="NE_STEPS"
      VARIABLES(86)="STEP"
      VARIABLES(85)="NBSA"
      VARIABLES(84)="LNA"
      VARIABLES(83)="KEEP"
      VARIABLES(82)="Deficiency"
      VARIABLES(81)="MAXIS1"
      VARIABLES(80)="IS1"
      VARIABLES(79)="IS"
      VARIABLES(78)="BUFR"
      VARIABLES(77)="POIDS"
      VARIABLES(76)="LBUFR_BYTES"
      VARIABLES(75)="LBUFR"
      VARIABLES(74)="ASS_IRECV"
      VARIABLES(73)="NSLAVES"
      VARIABLES(72)="NPROCS"
      VARIABLES(71)="MYID"
      VARIABLES(70)="COMM_LOAD"
      VARIABLES(69)="MYID_NODES"
      VARIABLES(68)="COMM_NODES"
      VARIABLES(67)="INST_Number"
      VARIABLES(66)="MAX_SURF_MASTER"
      VARIABLES(65)="KEEP8"
      VARIABLES(64)="pad8"
      VARIABLES(63)="WRITE_PROBLEM"
      VARIABLES(62)="OOC_PREFIX"
      VARIABLES(61)="OOC_TMPDIR"
      VARIABLES(60)="VERSION_NUMBER"
      VARIABLES(59)="MAPPING"
      VARIABLES(58)="LISTVAR_SCHUR"
      VARIABLES(57)="SCHUR_CINTERFACE"
      VARIABLES(56)="SCHUR"
      VARIABLES(55)="SIZE_SCHUR"
      VARIABLES(54)="SCHUR_LLD"
      VARIABLES(53)="SCHUR_NLOC"
      VARIABLES(52)="SCHUR_MLOC"
      VARIABLES(51)="NBLOCK"
      VARIABLES(50)="MBLOCK"
      VARIABLES(49)="NPCOL"
      VARIABLES(48)="NPROW"
      VARIABLES(47)="UNS_PERM"
      VARIABLES(46)="SYM_PERM"
      VARIABLES(45)="RINFOG"
      VARIABLES(44)="RINFO"
      VARIABLES(43)="CNTL"
      VARIABLES(42)="COST_SUBTREES"
      VARIABLES(41)="INFOG"
      VARIABLES(40)="INFO"
      VARIABLES(39)="ICNTL"
      VARIABLES(38)="pad5"
      VARIABLES(37)="LREDRHS"
      VARIABLES(36)="LSOL_loc"
      VARIABLES(35)="NZ_RHS"
      VARIABLES(34)="NRHS"
      VARIABLES(33)="LRHS"
      VARIABLES(32)="ISOL_loc"
      VARIABLES(31)="IRHS_PTR"
      VARIABLES(30)="IRHS_SPARSE"
      VARIABLES(29)="SOL_loc"
      VARIABLES(28)="RHS_SPARSE"
      VARIABLES(27)="REDRHS"
      VARIABLES(26)="RHS"
      VARIABLES(25)="PERM_IN"
      VARIABLES(24)="pad4"
      VARIABLES(23)="A_ELT"
      VARIABLES(22)="ELTVAR"
      VARIABLES(21)="ELTPTR"
      VARIABLES(20)="pad3"
      VARIABLES(19)="NELT"
      VARIABLES(18)="pad2"
      VARIABLES(17)="A_loc"
      VARIABLES(16)="JCN_loc"
      VARIABLES(15)="IRN_loc"
      VARIABLES(14)="pad1"
      VARIABLES(13)="NZ_loc"
      VARIABLES(12)="pad0"
      VARIABLES(11)="ROWSCA"
      VARIABLES(10)="COLSCA"
      VARIABLES(9)="JCN"
      VARIABLES(8)="IRN"
      VARIABLES(7)="A"
      VARIABLES(6)="NZ"
      VARIABLES(5)="N"
      VARIABLES(4)="JOB"
      VARIABLES(3)="PAR"
      VARIABLES(2)="SYM"
      VARIABLES(1)="COMM"
      NBVARIABLES_ROOT=34
      allocate(VARIABLES_ROOT(NBVARIABLES_ROOT))
      VARIABLES_ROOT(34)="NB_SINGULAR_VALUES"
      VARIABLES_ROOT(33)="SINGULAR_VALUES"
      VARIABLES_ROOT(32)="SVD_VT"
      VARIABLES_ROOT(31)="SVD_U"
      VARIABLES_ROOT(30)="gridinit_done"
      VARIABLES_ROOT(29)="yes"
      VARIABLES_ROOT(28)="rootpad3"
      VARIABLES_ROOT(27)="QR_RCOND"
      VARIABLES_ROOT(26)="rootpad"
      VARIABLES_ROOT(25)="RHS_ROOT"
      VARIABLES_ROOT(24)="rootpad2"
      VARIABLES_ROOT(23)="QR_TAU"
      VARIABLES_ROOT(22)="SCHUR_POINTER"
      VARIABLES_ROOT(21)="RHS_CNTR_MASTER_ROOT"
      VARIABLES_ROOT(20)="rootpad1"
      VARIABLES_ROOT(19)="IPIV"
      VARIABLES_ROOT(18)="RG2L_COL"
      VARIABLES_ROOT(17)="RG2L_ROW"
      VARIABLES_ROOT(16)="rootpad0"
      VARIABLES_ROOT(15)="LPIV"
      VARIABLES_ROOT(14)="CNTXT_BLACS"
      VARIABLES_ROOT(13)="DESCRIPTOR"
      VARIABLES_ROOT(12)="TOT_ROOT_SIZE"
      VARIABLES_ROOT(11)="ROOT_SIZE"
      VARIABLES_ROOT(10)="RHS_NLOC"
      VARIABLES_ROOT(9)="SCHUR_LLD"
      VARIABLES_ROOT(8)="SCHUR_NLOC"
      VARIABLES_ROOT(7)="SCHUR_MLOC"
      VARIABLES_ROOT(6)="MYCOL"
      VARIABLES_ROOT(5)="MYROW"
      VARIABLES_ROOT(4)="NPCOL"
      VARIABLES_ROOT(3)="NPROW"
      VARIABLES_ROOT(2)="NBLOCK"
      VARIABLES_ROOT(1)="MBLOCK"
      IN=52+id%MYID
      INFILEMAIN=trim(adjustl(INFILE)) // trim(adjustl(STRING_ID))
      open(UNIT=IN,FILE=INFILEMAIN, STATUS='old',FORM='unformatted'
     &     ,iostat=err)
      if(err.ne.0) THEN
         id%INFOG(1)=-91
         write(*,*) "IN SMUMPS_RESTART CANNOT OPEN FILE: "
     &        //trim(adjustl(INFILEMAIN))
         error=.true.
         goto 101
      endif
      DO i1=4,NBVARIABLES
         size_array1=0
         size_array2=0
         SELECT CASE(VARIABLES(i1))
         CASE("JOB")
            read(IN) id%JOB
         CASE("N")
            read(IN) id%N
         CASE("ICNTL")
            read(IN) id%ICNTL
         CASE("INFO")
            read(IN) id%INFO
         CASE("INFOG")
            read(IN) id%INFOG
         CASE("COST_SUBTREES")
            read(IN) id%COST_SUBTREES
         CASE("CNTL")
            read(IN) id%CNTL
         CASE("RINFO")
            read(IN) id%RINFO
         CASE("RINFOG")
            read(IN) id%RINFOG
         CASE("KEEP8")
            read(IN) id%KEEP8
         CASE("KEEP")
            read(IN) id%KEEP
         CASE("DKEEP")
            read(IN) id%DKEEP    
         CASE("NZ")
            read(IN) id%NZ
         CASE("A")
            nullify(id%A)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%A(size_array1))
               read(IN) id%A
            endif
         CASE("IRN")
            nullify(id%IRN)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%IRN(size_array1))
               read(IN) id%IRN
            endif
         CASE("JCN")
            nullify(id%JCN)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%JCN(size_array1))
               read(IN) id%JCN
            endif
         CASE("COLSCA")
            nullify(id%COLSCA)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%COLSCA(size_array1))
               read(IN) id%COLSCA
            endif
         CASE("ROWSCA")
            nullify(id%ROWSCA)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%ROWSCA(size_array1))
               read(IN) id%ROWSCA
            endif
         CASE("NZ_loc")
            read(IN) id%NZ_loc
         CASE("IRN_loc")
            nullify(id%IRN_loc)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%IRN_loc(size_array1))
               read(IN) id%IRN_loc
            endif
         CASE("JCN_loc")
            nullify(id%JCN_loc)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%JCN_loc(size_array1))
               read(IN) id%JCN_loc
            endif
         CASE("A_loc")
            nullify(id%A_loc)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%A_loc(size_array1))
               read(IN) id%A_loc
            endif
         CASE("NELT")
            read(IN) id%NELT
         CASE("ELTPTR")
            nullify(id%ELTPTR)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%ELTPTR(size_array1))
               read(IN) id%ELTPTR
            endif
         CASE("ELTVAR")
            nullify(id%ELTVAR)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%ELTVAR(size_array1))
               read(IN) id%ELTVAR
            endif
         CASE("A_ELT")
            nullify(id%A_ELT)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%A_ELT(size_array1))
               read(IN) id%A_ELT
            endif
         CASE("PERM_IN")
            nullify(id%PERM_IN)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%PERM_IN(size_array1))
               read(IN) id%PERM_IN
            endif
         CASE("RHS")
            nullify(id%RHS)
            read(IN) size_array1
            allocate(id%RHS(size_array1))
            read(IN) id%RHS
         CASE("REDRHS")
            nullify(id%REDRHS)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%REDRHS(size_array1))
               read(IN) id%REDRHS
            endif
         CASE("RHS_SPARSE")
            nullify(id%RHS_SPARSE)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%RHS_SPARSE(size_array1))
               read(IN) id%RHS_SPARSE
            endif
         CASE("SOL_loc")
            nullify(id%SOL_loc)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%SOL_loc(size_array1))
               read(IN) id%SOL_loc
            endif
         CASE("IRHS_SPARSE")
            nullify(id%IRHS_SPARSE)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%IRHS_SPARSE(size_array1))
               read(IN) id%IRHS_SPARSE
            endif
         CASE("IRHS_PTR")
            nullify(id%IRHS_PTR)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%IRHS_PTR(size_array1))
               read(IN) id%IRHS_PTR
            endif
         CASE("ISOL_loc")
            nullify(id%ISOL_loc)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%ISOL_loc(size_array1))
               read(IN) id%ISOL_loc
            endif
         CASE("LRHS")
            read(IN) id%LRHS
         CASE("NRHS")
            read(IN) id%NRHS
         CASE("NZ_RHS")
            read(IN) id%NZ_RHS
         CASE("LSOL_loc")
            read(IN) id%LSOL_loc
         CASE("LREDRHS")
            read(IN) id%LREDRHS
         CASE("SYM_PERM")
            nullify(id%SYM_PERM)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%SYM_PERM(size_array1))
               read(IN) id%SYM_PERM
            endif
         CASE("UNS_PERM")
            nullify(id%UNS_PERM)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%UNS_PERM(size_array1))
               read(IN) id%UNS_PERM
            endif
         CASE("NPROW")
            read(IN) id%NPROW
         CASE("NPCOL")
            read(IN) id%NPCOL
         CASE("MBLOCK")
            read(IN) id%MBLOCK
         CASE("NBLOCK")
            read(IN) id%NBLOCK
         CASE("SCHUR_MLOC")
            read(IN) id%SCHUR_MLOC
         CASE("SCHUR_NLOC")
            read(IN) id%SCHUR_NLOC
         CASE("SCHUR_LLD")
            read(IN) id%SCHUR_LLD
         CASE("SIZE_SCHUR")
            read(IN) id%SIZE_SCHUR
         CASE("SCHUR")
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%SCHUR(size_array1))
               read(IN) id%SCHUR
            endif
         CASE("SCHUR_CINTERFACE")
            nullify(id%SCHUR_CINTERFACE)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%SCHUR_CINTERFACE(size_array1))
               read(IN) id%SCHUR_CINTERFACE
            endif
         CASE("LISTVAR_SCHUR")
            nullify(id%LISTVAR_SCHUR)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%LISTVAR_SCHUR(size_array1))
               read(IN) id%LISTVAR_SCHUR
            endif
         CASE("MAPPING")
            nullify(id%MAPPING)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%MAPPING(size_array1))
               read(IN) id%MAPPING
            endif
         CASE("VERSION_NUMBER")
            read(IN) id%VERSION_NUMBER
         CASE("OOC_TMPDIR")
            read(IN) id%OOC_TMPDIR
         CASE("OOC_PREFIX")
            read(IN) id%OOC_PREFIX
         CASE("WRITE_PROBLEM")
            read(IN) id%WRITE_PROBLEM
         CASE("MAX_SURF_MASTER")
            read(IN) id%MAX_SURF_MASTER
         CASE("INST_Number")
            read(IN) id%INST_Number
         CASE("COMM_NODES")
            read(IN) id%COMM_NODES
         CASE("MYID_NODES")
            read(IN) id%MYID_NODES
         CASE("COMM_LOAD")
            read(IN) id%COMM_LOAD
         CASE("MYID")
            read(IN) id%MYID
         CASE("NPROCS")
            read(IN) id%NPROCS
         CASE("NSLAVES")
            read(IN) id%NSLAVES
         CASE("ASS_IRECV")
            read(IN) id%ASS_IRECV
         CASE("LBUFR")
            read(IN) id%LBUFR
         CASE("LBUFR_BYTES")
            read(IN) id%LBUFR_BYTES
         CASE("POIDS")
            nullify(id%POIDS)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%POIDS(size_array1))
               read(IN) id%POIDS
            endif
         CASE("BUFR")
            nullify(id%BUFR)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%BUFR(size_array1))
               read(IN) id%BUFR
            endif
         CASE("IS")
            nullify(id%IS)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%IS(size_array1))
               read(IN) id%IS
            endif
         CASE("IS1")
            nullify(id%IS1)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%IS1(size_array1))
               read(IN) id%IS1
            endif
         CASE("MAXIS1")
            read(IN) id%MAXIS1
         CASE("Deficiency")
            read(IN) id%Deficiency
         CASE("LNA")
            read(IN) id%LNA
         CASE("NBSA")
            read(IN) id%NBSA
         CASE("STEP")
            nullify(id%STEP)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%STEP(size_array1))
               read(IN) id%STEP
            endif
         CASE("NE_STEPS")
            nullify(id%NE_STEPS)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%NE_STEPS(size_array1))
               read(IN) id%NE_STEPS
            endif
         CASE("ND_STEPS")
            nullify(id%ND_STEPS)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%ND_STEPS(size_array1))
               read(IN) id%ND_STEPS
            endif
         CASE("Step2node")
            nullify(id%Step2node)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%Step2node(size_array1))
               read(IN) id%Step2node
            endif
         CASE("FRERE_STEPS")
            nullify(id%FRERE_STEPS)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%FRERE_STEPS(size_array1))
               read(IN) id%FRERE_STEPS
            endif
         CASE("DAD_STEPS")
            nullify(id%DAD_STEPS)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%DAD_STEPS(size_array1))
               read(IN) id%DAD_STEPS
            endif
         CASE("FILS")
            nullify(id%FILS)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%FILS(size_array1))
               read(IN) id%FILS
            endif
         CASE("PTRAR")
            nullify(id%PTRAR)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%PTRAR(size_array1))
               read(IN) id%PTRAR
            endif
         CASE("FRTPTR")
            nullify(id%FRTPTR)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%FRTPTR(size_array1))
               read(IN) id%FRTPTR
            endif
         CASE("FRTELT")
            nullify(id%FRTELT)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%FRTELT(size_array1))
               read(IN) id%FRTELT
            endif
         CASE("NA")
            nullify(id%NA)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%NA(size_array1))
               read(IN) id%NA
            endif
         CASE("PROCNODE_STEPS")
            nullify(id%PROCNODE_STEPS)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%PROCNODE_STEPS(size_array1))
               read(IN) id%PROCNODE_STEPS
            endif
         CASE("PTLUST_S")
            nullify(id%PTLUST_S)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%PTLUST_S(size_array1))
               read(IN) id%PTLUST_S
            endif
         CASE("PTRFAC")
            nullify(id%PTRFAC)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%PTRFAC(size_array1))
               read(IN) id%PTRFAC
            endif
         CASE("S")
            nullify(id%S)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%S(size_array1))
               read(IN) id%S
            endif
         CASE("PROCNODE")
            nullify(id%PROCNODE)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%PROCNODE(size_array1))
               read(IN) id%PROCNODE
            endif
         CASE("INTARR")
            nullify(id%INTARR)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%INTARR(size_array1))
               read(IN) id%INTARR
            endif
         CASE("DBLARR")
            nullify(id%DBLARR)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%DBLARR(size_array1))
               read(IN) id%DBLARR
            endif
         CASE("NELT_loc")
            read(IN) id%NELT_loc
         CASE("LELTVAR")
            read(IN) id%LELTVAR
         CASE("NA_ELT")
            read(IN) id%NA_ELT
         CASE("ELTPROC")
            nullify(id%ELTPROC)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%ELTPROC(size_array1))
               read(IN) id%ELTPROC
            endif
         CASE("CANDIDATES")
            nullify(id%CANDIDATES)
            read(IN) size_array1,size_array2
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%CANDIDATES(size_array1,size_array2))
               read(IN) id%CANDIDATES
            endif
         CASE("ISTEP_TO_INIV2")
            nullify(id%ISTEP_TO_INIV2)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%ISTEP_TO_INIV2(size_array1))
               read(IN) id%ISTEP_TO_INIV2
            endif
         CASE("FUTURE_NIV2")
            nullify(id%FUTURE_NIV2)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%FUTURE_NIV2(size_array1))
               read(IN) id%FUTURE_NIV2
            endif
         CASE("TAB_POS_IN_PERE")
            nullify(id%TAB_POS_IN_PERE)
            read(IN) size_array1,size_array2
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%TAB_POS_IN_PERE(size_array1,size_array2))
               read(IN) id%TAB_POS_IN_PERE
            endif
         CASE("I_AM_CAND")
            nullify(id%I_AM_CAND)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%I_AM_CAND(size_array1))
               read(IN) id%I_AM_CAND
            endif
         CASE("MEM_DIST")
            nullify(id%MEM_DIST)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%MEM_DIST(0:size_array1-1))
               read(IN) id%MEM_DIST
            endif
         CASE("POSINRHSCOMP_ROW")
            nullify(id%POSINRHSCOMP_ROW)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%POSINRHSCOMP_ROW(size_array1))
               read(IN) id%POSINRHSCOMP_ROW
            endif
         CASE("POSINRHSCOMP_COL")
            nullify(id%POSINRHSCOMP_COL)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%POSINRHSCOMP_COL(size_array1))
               read(IN) id%POSINRHSCOMP_COL
            endif
         CASE("POSINRHSCOMP_COL_ALLOC")
            read(IN) id%POSINRHSCOMP_COL_ALLOC
         CASE("RHSCOMP")
            nullify(id%RHSCOMP)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%RHSCOMP(size_array1))
               read(IN) id%RHSCOMP
            endif
         CASE("MEM_SUBTREE")
            nullify(id%MEM_SUBTREE)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%MEM_SUBTREE(size_array1))
               read(IN) id%MEM_SUBTREE
            endif
         CASE("COST_TRAV")
            nullify(id%COST_TRAV)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%COST_TRAV(size_array1))
               read(IN) id%COST_TRAV
            endif
         CASE("MY_ROOT_SBTR")
            nullify(id%MY_ROOT_SBTR)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%MY_ROOT_SBTR(size_array1))
               read(IN) id%MY_ROOT_SBTR
            endif
         CASE("MY_FIRST_LEAF")
            nullify(id%MY_FIRST_LEAF)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%MY_FIRST_LEAF(size_array1))
               read(IN) id%MY_FIRST_LEAF
            endif
         CASE("MY_NB_LEAF")
            nullify(id%MY_NB_LEAF)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%MY_NB_LEAF(size_array1))
               read(IN) id%MY_NB_LEAF
            endif
         CASE("DEPTH_FIRST")
            nullify(id%DEPTH_FIRST)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%DEPTH_FIRST(size_array1))
               read(IN) id%DEPTH_FIRST
            endif
         CASE("DEPTH_FIRST_SEQ")
            nullify(id%DEPTH_FIRST_SEQ)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%DEPTH_FIRST_SEQ(size_array1))
               read(IN) id%DEPTH_FIRST_SEQ
            endif
         CASE("SBTR_ID")
            nullify(id%SBTR_ID)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%SBTR_ID(size_array1))
               read(IN) id%SBTR_ID
            endif
         CASE("SCHED_DEP")
            nullify(id%SCHED_DEP)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%SCHED_DEP(size_array1))
               read(IN) id%SCHED_DEP
            endif
         CASE("SCHED_GRP")
            nullify(id%SCHED_GRP)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%SCHED_GRP(size_array1))
               read(IN) id%SCHED_GRP
            endif
         CASE("CROIX_MANU")
            nullify(id%CROIX_MANU)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%CROIX_MANU(size_array1))
               read(IN) id%CROIX_MANU
            endif
         CASE("WK_USER")
            nullify(id%WK_USER)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%WK_USER(size_array1))
               read(IN) id%WK_USER
            endif
         CASE("NBSA_LOCAL")
            read(IN) id%NBSA_LOCAL
         CASE("LWK_USER")
            read(IN) id%LWK_USER
         CASE("CB_SON_SIZE")
            nullify(id%CB_SON_SIZE)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%CB_SON_SIZE(size_array1))
               read(IN) id%CB_SON_SIZE
            endif
         CASE("INSTANCE_NUMBER")
            read(IN) id%INSTANCE_NUMBER
         CASE("OOC_MAX_NB_NODES_FOR_ZONE")
            read(IN) id%OOC_MAX_NB_NODES_FOR_ZONE
         CASE("OOC_INODE_SEQUENCE")
            nullify(id%OOC_INODE_SEQUENCE)
            read(IN) size_array1,size_array2
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%OOC_INODE_SEQUENCE(size_array1,size_array2))
               read(IN) id%OOC_INODE_SEQUENCE
            endif
         CASE("OOC_SIZE_OF_BLOCK")
            nullify(id%OOC_SIZE_OF_BLOCK)
            read(IN) size_array1,size_array2
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%OOC_SIZE_OF_BLOCK(size_array1,size_array2))
               read(IN) id%OOC_SIZE_OF_BLOCK
            endif
         CASE("OOC_VADDR")
            nullify(id%OOC_VADDR)
            read(IN) size_array1,size_array2
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%OOC_VADDR(size_array1,size_array2))
               read(IN) id%OOC_VADDR
            endif
         CASE("OOC_TOTAL_NB_NODES")
            nullify(id%OOC_TOTAL_NB_NODES)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%OOC_TOTAL_NB_NODES(size_array1))
               read(IN) id%OOC_TOTAL_NB_NODES
            endif
         CASE("OOC_NB_FILES")
            nullify(id%OOC_NB_FILES)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%OOC_NB_FILES(size_array1))
               read(IN) id%OOC_NB_FILES
            endif
         CASE("OOC_NB_FILE_TYPE")
            read(IN) id%OOC_NB_FILE_TYPE
         CASE("OOC_FILE_NAMES")
            nullify(id%OOC_FILE_NAMES)
            read(IN) size_array1,size_array2
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%OOC_FILE_NAMES(size_array1,size_array2))
               read(IN) id%OOC_FILE_NAMES
            endif
         CASE("OOC_FILE_NAME_LENGTH")
            nullify(id%OOC_FILE_NAME_LENGTH)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%OOC_FILE_NAME_LENGTH(size_array1))
               read(IN) id%OOC_FILE_NAME_LENGTH
            endif
         CASE("PIVNUL_LIST")
            nullify(id%PIVNUL_LIST)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%PIVNUL_LIST(size_array1))
               read(IN) id%PIVNUL_LIST
            endif
         CASE("SUP_PROC")
            nullify(id%SUP_PROC)
            read(IN) size_array1,size_array2
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%SUP_PROC(size_array1,size_array2))
               read(IN) id%SUP_PROC
            endif
         CASE("IPTR_WORKING")
            nullify(id%IPTR_WORKING)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%IPTR_WORKING(size_array1))
               read(IN) id%IPTR_WORKING
            endif
         CASE("WORKING")
            nullify(id%WORKING)
            read(IN) size_array1
            if(size_array1.EQ.-999) then
               read(IN) dummy
            else
               allocate(id%WORKING(size_array1))
               read(IN) id%WORKING
            endif
         CASE("root")
        CASE("NBGRP")
           read(IN) id%NBGRP
        CASE("LRGROUPS")
           nullify(id%LRGROUPS)
           read(IN) size_array1
           if(size_array1.EQ.-999) then
              read(IN) dummy
           else
              allocate(id%LRGROUPS(size_array1))
              read(IN) id%LRGROUPS
           endif
         CASE("SCHED_SBTR")
            nullify(id%SCHED_SBTR)
            read(IN) size_array1
            IF(size_array1.EQ.-999) THEN
               read(IN) dummy
            ELSE
               allocate(id%SCHED_SBTR(size_array1))
               read(IN) id%SCHED_SBTR
            ENDIF 
        CASE("LPOOL_AFTER_L0_OMP")
           read(IN) id%LPOOL_AFTER_L0_OMP
        CASE("LPOOL_BEFORE_L0_OMP")
           read(IN) id%LPOOL_BEFORE_L0_OMP
        CASE("L_PHYS_L0_OMP")
           read(IN) id%L_PHYS_L0_OMP
        CASE("L_VIRT_L0_OMP")  
           read(IN) id%L_VIRT_L0_OMP
        CASE("LL0_OMP_MAPPING")
           read(IN) id%LL0_OMP_MAPPING
        CASE("THREAD_LA")
           read(IN) id%THREAD_LA
        CASE("IPOOL_AFTER_L0_OMP")
           nullify(id%IPOOL_AFTER_L0_OMP)
           read(IN) size_array1
           IF(size_array1.EQ.-999) THEN
              read(IN) dummy
           ELSE
              allocate(id%IPOOL_AFTER_L0_OMP(size_array1))
              read(IN) id%IPOOL_AFTER_L0_OMP
           ENDIF 
        CASE("IPOOL_BEFORE_L0_OMP")
           nullify(id%IPOOL_BEFORE_L0_OMP)
           read(IN) size_array1
           IF(size_array1.EQ.-999) THEN
              read(IN) dummy
           ELSE
              allocate(id%IPOOL_BEFORE_L0_OMP(size_array1))
              read(IN) id%IPOOL_BEFORE_L0_OMP
           ENDIF 
        CASE("PHYS_L0_OMP")
           nullify(id%PHYS_L0_OMP)
           read(IN) size_array1
           IF(size_array1.EQ.-999) THEN
              read(IN) dummy
           ELSE
              allocate(id%PHYS_L0_OMP(size_array1))
              read(IN) id%PHYS_L0_OMP
           ENDIF 
        CASE("VIRT_L0_OMP")
           nullify(id%VIRT_L0_OMP)
           read(IN) size_array1
           IF(size_array1.EQ.-999) THEN
              read(IN) dummy
           ELSE
              allocate(id%VIRT_L0_OMP(size_array1))
              read(IN) id%VIRT_L0_OMP
           ENDIF
        CASE("PERM_L0_OMP")
           nullify(id%PERM_L0_OMP)
           read(IN) size_array1
           IF(size_array1.EQ.-999) THEN
              read(IN) dummy
           ELSE
              allocate(id%PERM_L0_OMP(size_array1))
              read(IN) id%PERM_L0_OMP
           ENDIF
        CASE("PTR_LEAFS_L0_OMP")
           nullify(id%PTR_LEAFS_L0_OMP)
           read(IN) size_array1
           IF(size_array1.EQ.-999) THEN
              read(IN) dummy
           ELSE
              allocate(id%PTR_LEAFS_L0_OMP(size_array1))
              read(IN) id%PTR_LEAFS_L0_OMP
           ENDIF
        CASE("L0_OMP_MAPPING")
           nullify(id%L0_OMP_MAPPING)
           read(IN) size_array1
           IF(size_array1.EQ.-999) THEN
              read(IN) dummy
           ELSE
              allocate(id%L0_OMP_MAPPING(size_array1))
              read(IN) id%L0_OMP_MAPPING
           ENDIF 
        CASE("IF_RESTARTING")
           read(IN) id%IF_RESTARTING 
        CASE("SINGULAR_VALUES")
         read(IN) size_array1
         if(size_array1.EQ.-999) then
            read(IN) dummy
         else
            allocate(id%SINGULAR_VALUES(size_array1))
            read(IN) id%SINGULAR_VALUES
         endif  
      CASE("NB_SINGULAR_VALUES")
        write(IN) id%NB_SINGULAR_VALUES
        CASE("pad0","pad1","pad2","pad3","pad4","pad5","pad8","pad11",
     &          "pad111", "pad12","pad13","pad14")
        CASE DEFAULT
           id%INFOG(1)=-92
           write(*,*) "IN SMUMPS_RESTART UNKNOWN VARIABLE IN "
     &       // "STRUCTURE: "//  trim(adjustl(VARIABLES(i1)))
           goto 101
        END SELECT
      ENDDO
 101  continue
      close(IN)
      if(error) goto 102
      INROOT=1000+IN
      INFILEROOT=trim(adjustl(INFILE)) // trim(adjustl("ROOT")) 
     &     // trim(adjustl(STRING_ID)) 
      open(UNIT=INROOT,FILE=INFILEROOT, STATUS='old',FORM='unformatted'
     &     ,iostat=err)
      if(err.ne.0) THEN
         id%INFOG(1)=-91
         write(*,*) "IN SMUMPS_RESTART CANNOT OPEN FILE: "
     &        //trim(adjustl(INFILEROOT))
         goto 102
      endif
      DO i2=1,NBVARIABLES_ROOT
         SELECT CASE(VARIABLES_ROOT(i2))
      CASE("MBLOCK")
         read(INROOT) id%root%MBLOCK
      CASE("NBLOCK")
         read(INROOT) id%root%NBLOCK
      CASE("NPROW")
         read(INROOT) id%root%NPROW
      CASE("NPCOL")
         read(INROOT) id%root%NPCOL
      CASE("MYROW")
         read(INROOT) id%root%MYROW
      CASE("MYCOL")
         read(INROOT) id%root%MYCOL
      CASE("SCHUR_MLOC")
         read(INROOT) id%root%SCHUR_MLOC
      CASE("SCHUR_NLOC")
         read(INROOT) id%root%SCHUR_NLOC
      CASE("SCHUR_LLD")
         read(INROOT) id%root%SCHUR_LLD
      CASE("RHS_NLOC")
         read(INROOT) id%root%RHS_NLOC
      CASE("ROOT_SIZE")
         read(INROOT) id%root%ROOT_SIZE
      CASE("TOT_ROOT_SIZE")
         read(INROOT) id%root%TOT_ROOT_SIZE
      CASE("DESCRIPTOR")
         read(INROOT) id%root%DESCRIPTOR
      CASE("CNTXT_BLACS")
         read(INROOT) id%root%CNTXT_BLACS
      CASE("LPIV")
         read(INROOT) id%root%LPIV
      CASE("RG2L_ROW")
         nullify(id%root%RG2L_ROW)
         read(INROOT) size_array1
         if(size_array1.EQ.-999) then
            read(INROOT) dummy
         else
            allocate(id%root%RG2L_ROW(size_array1))
            read(INROOT) id%root%RG2L_ROW
         endif
      CASE("RG2L_COL")
         nullify(id%root%RG2L_COL)
         read(INROOT) size_array1
         if(size_array1.EQ.-999) then
            read(INROOT) dummy
         else
            allocate(id%root%RG2L_COL(size_array1))
            read(INROOT) id%root%RG2L_COL
         endif
      CASE("IPIV")
         nullify(id%root%IPIV)
         read(INROOT) size_array1
         if(size_array1.EQ.-999) then
            read(INROOT) dummy
         else
            allocate(id%root%IPIV(size_array1))
            read(INROOT) id%root%IPIV
         endif
      CASE("RHS_CNTR_MASTER_ROOT")
         nullify(id%root%RHS_CNTR_MASTER_ROOT)
         read(INROOT) size_array1
         if(size_array1.EQ.-999) then
            read(INROOT) dummy
         else
            allocate(id%root%RHS_CNTR_MASTER_ROOT(size_array1))
            read(INROOT) id%root%RHS_CNTR_MASTER_ROOT
         endif
      CASE("SCHUR_POINTER")
         nullify(id%root%SCHUR_POINTER)
         read(INROOT) size_array1
         if(size_array1.EQ.-999) then
            read(INROOT) dummy
         else
            allocate(id%root%SCHUR_POINTER(size_array1))
            read(INROOT) id%root%SCHUR_POINTER
         endif
      CASE("QR_TAU")
         nullify(id%root%QR_TAU)
         read(INROOT) size_array1
         if(size_array1.EQ.-999) then
            read(INROOT) dummy
         else
            allocate(id%root%QR_TAU(size_array1))
            read(INROOT) id%root%QR_TAU
         endif
      CASE("RHS_ROOT")
         nullify(id%root%RHS_ROOT)
         read(INROOT) size_array1,size_array2
         if(size_array1.EQ.-999) then
            read(INROOT) dummy
         else
            allocate(id%root%RHS_ROOT(size_array1,size_array2))
            read(INROOT) id%root%RHS_ROOT
         endif
      CASE("QR_RCOND")
         read(INROOT) id%root%QR_RCOND
      CASE("yes")
         read(INROOT) id%root%yes
      CASE("gridinit_done")
         read(INROOT) id%root%gridinit_done
      CASE("SVD_U")
         read(INROOT) size_array1,size_array2
         if(size_array1.EQ.-999) then
            read(INROOT) dummy
         else
            allocate(id%root%SVD_U(size_array1,size_array2))
            read(INROOT) id%root%SVD_U
         endif
      CASE("SVD_VT")
         read(INROOT) size_array1,size_array2
         if(size_array1.EQ.-999) then
            read(INROOT) dummy
         else
            allocate(id%root%SVD_VT(size_array1,size_array2))
            read(INROOT) id%root%SVD_VT
         endif
      CASE("SINGULAR_VALUES")
         read(INROOT) size_array1
         if(size_array1.EQ.-999) then
            read(INROOT) dummy
         else
            allocate(id%root%SINGULAR_VALUES(size_array1))
            read(INROOT) id%root%SINGULAR_VALUES
         endif  
      CASE("NB_SINGULAR_VALUES")
        write(INROOT) id%root%NB_SINGULAR_VALUES
         CASE("rootpad0","rootpad1","rootpad2","rootpad",
     &        "rootpad3")
      CASE DEFAULT
         id%INFOG(1)=-92
         write(*,*) "IN SMUMPS_RESTART UNKNOWN VARIABLE "
     &        // "IN ROOT: "// trim(adjustl(VARIABLES_ROOT(i2)))
         goto 102
      END SELECT  
      ENDDO
      if(id%root%gridinit_done) then
         id%root%CNTXT_BLACS = id%COMM_NODES
         CALL blacs_gridinit( id%root%CNTXT_BLACS, 'R',
     &        id%root%NPROW, id%root%NPCOL )
         id%root%gridinit_done = .TRUE.
      endif
 102  continue
      close(INROOT)
      END SUBROUTINE SMUMPS_RESTART
