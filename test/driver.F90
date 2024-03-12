program driver

  use MPI
  use OMP_LIB

  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: error_type
  !Tests:
  use Functionality_Suite, only: get_quasienergies_of_BC2N
  use Randomized_Suite, only: randomized_input_parameters

  implicit none

  integer :: rank, ierror

  type(error_type), allocatable :: error

  call MPI_INIT(ierror)

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  if (rank == 0) write (error_unit, "(A)") "Suite: Functionality Suite:"

  if (rank == 0) write (error_unit, "(A)") &
    "Test (1/1): Testing quasienergies of BC2N. Comparing with reference."
  call get_quasienergies_of_BC2N(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (1/1): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (1/1): PASSED."

  if (rank == 0) write (error_unit, "(A)") "Suite: Randomized Suite:"

  if (rank == 0) write (error_unit, "(A)") &
    "Test (1/1): Testing for runtime errors and NaN/Infinity. Random driving parameters."
  call randomized_input_parameters(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (1/1): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (1/1): PASSED."

  call MPI_FINALIZE(ierror)

end program driver
