


module mumps_mod

   save
  ! integer,parameter:: ui_out = 10  ! unit for mumps log file
   !integer,parameter:: st_unit=11  ! unit for output statistics file


   INCLUDE 'dmumps_struc.h'

   !logical,private:: mumps_initialized = .false.


contains

!---------------------------------------------------------------

subroutine init( mumps_par, sym )

implicit none

TYPE (DMUMPS_STRUC),intent(inout):: mumps_par
integer,intent(in):: sym  ! =0 unsymmetric, =1 symm. pos def, =2 general symm.x

   !mumps_par%COMM = MPI_COMM_WORLD
   mumps_par%SYM = sym  ! =0 unsymmetric, =1 symm. pos def, =2 general symm.
   mumps_par%PAR = 1

   mumps_par%JOB = -1  ! initialize
   CALL DMUMPS(mumps_par)


   !if (mumps_par%MYID == 0) then
    !  open(unit=ui_out, file='mumps.log', action='write')
   !end if

   !mumps_initialized = .true.

return
end subroutine init

!---------------------------------------------------------------

subroutine convert_to_mumps_format( mumps_par, n, A,jA,iA, ierr )
! A is actually a transpose.
implicit none

TYPE (DMUMPS_STRUC),intent(inout):: mumps_par
integer,intent(in):: n  ! # of rows in A
real(kind=8),intent(in):: A(*)
integer,intent(in):: jA(*), iA(n+1)
integer,intent(out):: ierr

integer nonz, i,j, j1,j2, ind, jcol, istat


nonz = iA(n+1) - 1
!n = size(iA) - 1

mumps_par%N = n

if (mumps_par%SYM == 0) then  ! not symmetric matrix
   mumps_par%NZ = nonz
   allocate( mumps_par%IRN(nonz), mumps_par%JCN(nonz), mumps_par%A(nonz), stat=istat)
   if ( istat /= 0 ) then
      ierr = -13  ! allocation error
      return
   end if


   do i = 1, n
      mumps_par%JCN( iA(i) : iA(i+1)-1 )  =  i
   end do  ! i

   mumps_par%IRN = jA(1:nonz)
   mumps_par%A    = A(1:nonz)

else  ! symmetric matrix
   ! Keep only the lower half of the matrix (row >= column).

   nonz = nonz/2 + n  ! should be +n/2, but I'm using n just in case.
   allocate( mumps_par%IRN(nonz), mumps_par%JCN(nonz), mumps_par%A(nonz), stat=istat)
   if ( istat /= 0 ) then
      ierr = -13  ! allocation error
      return
   end if

   ind = 0
   do i = 1, n

      j1 = iA(i)
      j2 = iA(i+1) - 1
      do j = j1, j2
         jcol = jA(j)

         if (i >= jcol) then
            ind = ind + 1
            mumps_par%A(ind) = A(j)
            mumps_par%JCN(ind) = jcol
            mumps_par%IRN(ind) = i
         end if

      end do  ! j
   end do  ! i

   mumps_par%NZ = ind
   !if (nonz < ind) call errorstop('nonz < ind')  ! debug
end if

ierr = 0
return
end subroutine convert_to_mumps_format

!---------------------------------------------------------------

subroutine factor_matrix( mumps_par, ooc, ierr )

implicit none

TYPE (DMUMPS_STRUC),intent(inout):: mumps_par
integer,intent(in):: ooc ! = 0 in-core factorization, = 1 out-of-core factorization
integer,intent(out):: ierr

mumps_par%icntl(2) = -1  ! output stream for diagnostics
mumps_par%icntl(4) = -1 ! 1  ! amount of output


mumps_par%icntl(2) = 0 ! ui_out  ! output stream for diagnostic printing
mumps_par%icntl(3) = 0 ! ui_out  ! output stream for global information

!mumps_par%icntl(14) = 40 ! % increase in working space

mumps_par%icntl(22) = ooc ! out-of-core factorization

!mumps_par%icntl(11) = 0 ! return statistics
! mumps_par%icntl(7) = 5  ! ordering type
!mumps_par%cntl(1) = 0.d0  ! rel. threshold for numerical pivoting


mumps_factorization: do

   mumps_par%JOB = 1  !  analysis
   CALL DMUMPS(mumps_par)

   mumps_par%JOB = 2  !  factorization
   CALL DMUMPS(mumps_par)

   ierr = mumps_par%INFOG(1)


   if ( mumps_par%INFOG(1) == -9 .or. mumps_par%INFOG(1) == -8 ) then
      ! Main internal real/complex workarray S too small.
      mumps_par%icntl(14) = mumps_par%icntl(14) + 10
      if ( mumps_par%MYID == 0 ) then
         write(*,30) mumps_par%icntl(14)
         30 format(/'MUMPS percentage increase in the estimated working space',  &
                   /'increased to',i4)
      end if

   else
      exit mumps_factorization ! factorization successful
   end if

end do mumps_factorization


if ( mumps_par%MYID == 0 ) then
   !flush(ui_out)
   ! Turn off mumps output.
   mumps_par%icntl(2) = 0  ! output stream for diagnostic printing
   mumps_par%icntl(3) = 0  ! output stream for global information
end if


return
end subroutine factor_matrix

!--------------------------------------------------------------

subroutine solve( mumps_par, nrhs, rhs, x, transpose )
! Solve A*x = rhs

implicit none

TYPE (DMUMPS_STRUC),intent(inout):: mumps_par
integer,intent(in):: nrhs   ! # of right-hand-sides
real(kind=8),intent(in):: rhs(nrhs * mumps_par%N)
real(kind=8),intent(out),target:: x(nrhs * mumps_par%N)  ! solution
logical,intent(in):: transpose  ! if .true. take the transpose

   x = rhs

   mumps_par%RHS => x

   ! The following is significant only on the host cpu.
   mumps_par%NRHS = nrhs  ! # of right-hand-sides
   mumps_par%LRHS = mumps_par%N  ! size of system
   mumps_par%icntl(20) = 0  !Use dense rhs

   if (transpose) then
      mumps_par%icntl(9) = 0  ! for solving A'x = b
   else
      mumps_par%icntl(9) = 1  ! for solving Ax = b
   end if


   mumps_par%JOB = 3  ! Solve the system.
   CALL DMUMPS(mumps_par)
   ! At this point mumps_par%RHS (rhs) contains the solution.


return
end subroutine solve

!--------------------------------------------------------------

subroutine solve_sparse_rhs( mumps_par, nrhs, x, transpose )
! Solve A*x = rhs

implicit none

TYPE (DMUMPS_STRUC),intent(inout):: mumps_par
integer,intent(in):: nrhs   ! # of right-hand-sides
real(kind=8),intent(out),target:: x(nrhs * mumps_par%N)  ! solution
logical,intent(in):: transpose  ! if .true. take the transpose

   mumps_par%RHS       => x
   mumps_par%icntl(20) = 1 !Use sparse rhs and let MUMPS choose whether
                           !or not to exploit it during solution

   if (transpose) then
      mumps_par%icntl(9) = 0  ! for solving A'x = b
   else
      mumps_par%icntl(9) = 1  ! for solving Ax = b
   end if


   mumps_par%JOB = 3  ! Solve the system.
   CALL DMUMPS(mumps_par)
   ! At this point mumps_par%RHS (rhs) contains the solution.


return
end subroutine solve_sparse_rhs

!--------------------------------------------------------------

subroutine destroy(mumps_par)

implicit none
TYPE (DMUMPS_STRUC),intent(inout):: mumps_par

!if (.not. mumps_initialized)  return

!  Destroy the instance (deallocate internal data structures)
mumps_par%JOB = -2
CALL DMUMPS(mumps_par)

if (associated(mumps_par%A)) then
   deallocate(mumps_par%IRN, mumps_par%JCN, mumps_par%A)
end if

!if (mumps_par%MYID == 0) then
!   close(ui_out)
!end if

return
end subroutine destroy


end module mumps_mod
