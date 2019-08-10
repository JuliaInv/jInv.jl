

!-------------------------------------------------

function factor_mumps_cmplx( n, sym, ooc, A, jA, iA, ierr )  result(pm_out)
!DIR$ ATTRIBUTES DLLEXPORT :: factor_mumps_cmplx
!DIR$ ATTRIBUTES ALIAS: 'factor_mumps_cmplx_':: factor_mumps_cmplx

use mumps_cmplx_mod, only: init, convert_to_mumps_format, factor_matrix, destroy
implicit none

INCLUDE 'zmumps_struc.h'

integer:: pm_out   ! mumps pointer
integer,intent(in):: n  ! # of rows in A
integer,intent(in):: sym  ! =0 unsymmetric, =1 symm. pos def, =2 general symm.
integer,intent(in):: ooc  ! = 0 in-core factorization, = 1 out-of-core factorization
complex(kind=8),intent(in):: A(*)
integer,intent(in):: jA(*), iA(n+1)

integer,intent(out):: ierr
! integer,intent(out):: ierr  ! =0 no error, < 0 error (memory, singular, etc.)

TYPE(ZMUMPS_STRUC),pointer:: pmumps_par
TYPE(ZMUMPS_STRUC):: mumps_par

pointer ( pm, mumps_par )

allocate(pmumps_par)
pm = loc(pmumps_par)

call init( mumps_par, sym)
call convert_to_mumps_format(mumps_par, n, A,jA,iA, ierr )
if (ierr < 0) return  ! allocation error

call factor_matrix(mumps_par, ooc, ierr)

if ( ierr < 0 ) then
   ! Error in factorization.
   call destroy(mumps_par)
   pm_out = 0
else
   pm_out = pm
end if

return
end  function factor_mumps_cmplx

!-------------------------------------------------

subroutine solve_mumps_cmplx( pm_in, nrhs, rhs, x, transpose )
! Solve A*x = rhs

!DIR$ ATTRIBUTES DLLEXPORT :: solve_mumps_cmplx
!DIR$ ATTRIBUTES ALIAS: 'solve_mumps_cmplx_':: solve_mumps_cmplx
use mumps_cmplx_mod, only: solve

implicit none

INCLUDE 'zmumps_struc.h'

integer,intent(in):: pm_in  ! mumps pointer
integer,intent(in):: nrhs  ! # of right-hand-sides
complex(kind=8),intent(in):: rhs(*)   ! right-hand-side
complex(kind=8),intent(out):: x(*)    ! solution
integer,intent(in):: transpose   ! =1 for transpose

TYPE(ZMUMPS_STRUC):: mumps_par
pointer ( pm, mumps_par )
pm = pm_in

call solve(mumps_par, nrhs, rhs, x, (transpose==1) )

return
end subroutine solve_mumps_cmplx

!-------------------------------------------------

subroutine solve_mumps_cmplx_sparse_rhs( pm_in, nzrhs, nrhs, rhsval, irhsrow, irhscolptr, x, transpose )
! Solve A*x = rhs where rhs is input as a CSC sparse matrix. MUMPS automatically chooses whether or
! not to exploit rhs sparsity in the solution procedure.
! See pg. 36 of MUMPS 5.0.0 users' guide for description of sparse rhs input.

use mumps_cmplx_mod, only: solve_sparse_rhs

implicit none

INCLUDE 'zmumps_struc.h'

integer,intent(in):: pm_in  ! mumps pointer
integer,intent(in):: nzrhs,nrhs  ! total # of non-zeros across all rhs,# of right-hand-sides
integer,intent(in),target:: irhsrow(*), irhscolptr(*)  !Right hand side indexing pointers
complex(kind=8),intent(in),target:: rhsval(*)   ! right-hand-side values
complex(kind=8),intent(out):: x(*)    ! solution
integer,intent(in):: transpose   ! =1 for transpose

TYPE(ZMUMPS_STRUC):: mumps_par
pointer ( pm, mumps_par )
pm = pm_in

! Setup sparse rhs 
mumps_par%NZ_RHS      = nzrhs
mumps_par%NRHS        = nrhs
mumps_par%LRHS        = mumps_par%N  ! size of system
mumps_par%RHS_SPARSE  => rhsval(1:nzrhs)
mumps_par%IRHS_SPARSE => irhsrow(1:nzrhs)
mumps_par%IRHS_PTR    => irhscolptr(1:nrhs+1)

call solve_sparse_rhs(mumps_par, nrhs, x, (transpose==1) )

return
end subroutine solve_mumps_cmplx_sparse_rhs

!-------------------------------------------------

subroutine destroy_mumps_cmplx( pm_in )
!  Destroy the instance (deallocate internal data structures)

!DIR$ ATTRIBUTES DLLEXPORT :: destroy_mumps_cmplx
!DIR$ ATTRIBUTES ALIAS: 'destroy_mumps_cmplx_':: destroy_mumps_cmplx
use mumps_cmplx_mod, only: destroy

implicit none
INCLUDE 'zmumps_struc.h'

integer,intent(in):: pm_in  ! mumps pointer
TYPE(ZMUMPS_STRUC):: mumps_par
pointer ( pm, mumps_par )
pm = pm_in

call destroy(mumps_par)

return
end subroutine destroy_mumps_cmplx

!-------------------------------------------------

