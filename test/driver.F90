program driver

  use MPI
  use OMP_LIB

  use, intrinsic :: iso_fortran_env, only: error_unit, input_unit
  use testdrive, only: error_type
  !Tests:
  use Q_Functionality_Suite, only: Q_fnc => get_quasienergies_of_BC2N
  use Q_Randomized_Suite, only: Q_rnd => randomized_input_parameters
  use C_Functionality_Suite, only: C_fnc => get_current_of_BC2N
  use C_Randomized_Suite, only: C_rnd => randomized_input_parameters

  implicit none

  integer :: rank, ierror

  type(error_type), allocatable :: error

  call MPI_INIT(ierror)

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  if (rank == 0) write (error_unit, "(A)") "Testing: Quasienergy module:"

  if (rank == 0) write (error_unit, "(A)") "Suite: Functionality Suite:"

  if (rank == 0) write (error_unit, "(A)") &
    "Test (1/1): Testing quasienergies of BC2N. Comparing with reference."
  call Q_fnc(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (1/1): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (1/1): PASSED."

  if (rank == 0) write (error_unit, "(A)") "Suite: Randomized Suite:"

  if (rank == 0) write (error_unit, "(A)") &
    "Test (1/1): Testing for runtime errors and NaN/Infinity. Random driving parameters."
  call Q_rnd(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (1/1): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (1/1): PASSED."
  if (rank == 0) write (error_unit, "(A)") "Testing: Currents module:"

  if (rank == 0) write (error_unit, "(A)") "Suite: Functionality Suite:"

  if (rank == 0) write (error_unit, "(A)") &
    "Test (1/1): Testing current of BC2N. Comparing with reference."
  call C_fnc(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (1/1): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (1/1): PASSED."

  if (rank == 0) write (error_unit, "(A)") "Suite: Randomized Suite:"

  if (rank == 0) write (error_unit, "(A)") &
    "Test (1/1): Testing for runtime errors and NaN/Infinity. Random driving parameters."
  call C_rnd(error)
  if (rank == 0) then
    if (allocated(error)) error stop "Test (1/1): FAILED."
  endif
  if (rank == 0) write (error_unit, "(A)") "Test (1/1): PASSED."

  call MPI_FINALIZE(ierror)

end program driver
