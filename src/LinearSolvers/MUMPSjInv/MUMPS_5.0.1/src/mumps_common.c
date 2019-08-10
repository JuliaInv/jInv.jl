/*
 *
 *  This file is part of MUMPS 5.0.1, released
 *  on Thu Jul 23 17:08:29 UTC 2015
 *
 *
 *  Copyright 1991-2015 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license:
 *  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
 *
 */
#include "mumps_common.h"
/* Special case of mapping and pivnul_list -- allocated from MUMPS */
static MUMPS_INT * MUMPS_MAPPING;
static MUMPS_INT * MUMPS_PIVNUL_LIST;
/* as uns_perm and sym_perm */
static MUMPS_INT * MUMPS_SYM_PERM;
static MUMPS_INT * MUMPS_UNS_PERM;
MUMPS_INT*
mumps_get_mapping()
{
    return MUMPS_MAPPING;
}
void MUMPS_CALL
MUMPS_ASSIGN_MAPPING(MUMPS_INT * f77mapping)
{
    MUMPS_MAPPING = f77mapping;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_MAPPING()
{
    MUMPS_MAPPING = 0;
}
MUMPS_INT*
mumps_get_pivnul_list()
{
    return MUMPS_PIVNUL_LIST;
}
void MUMPS_CALL
MUMPS_ASSIGN_PIVNUL_LIST(MUMPS_INT * f77pivnul_list)
{
    MUMPS_PIVNUL_LIST = f77pivnul_list;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_PIVNUL_LIST()
{
    MUMPS_PIVNUL_LIST = 0;
}
MUMPS_INT*
mumps_get_sym_perm()
{
    return MUMPS_SYM_PERM;
}
void MUMPS_CALL
MUMPS_ASSIGN_SYM_PERM(MUMPS_INT * f77sym_perm)
{
    MUMPS_SYM_PERM = f77sym_perm;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_SYM_PERM()
{
    MUMPS_SYM_PERM = 0;
}
MUMPS_INT*
mumps_get_uns_perm()
{
    return MUMPS_UNS_PERM;
}
void MUMPS_CALL
MUMPS_ASSIGN_UNS_PERM(MUMPS_INT * f77uns_perm)
{
    MUMPS_UNS_PERM = f77uns_perm;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_UNS_PERM()
{
    MUMPS_UNS_PERM = 0;
}
