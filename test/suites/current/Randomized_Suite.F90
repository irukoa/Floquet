module C_Randomized_Suite
  use, intrinsic :: iso_fortran_env, only: error_unit
  use MPI
  use OMP_LIB

  use Floquet_kinds, only: dp
  use Floquet_defs, only: pi
  use Currents, only: currents_calc_tsk
  use WannInt, only: crystal

  use testdrive, only: error_type

  implicit none
  private

  real(dp) :: inf = huge(1.0_dp)

  public :: randomized_input_parameters

contains

  subroutine randomized_input_parameters(error)

    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: BC2N
    type(currents_calc_tsk) :: tsk

    real(dp) :: klist(3, 1), rNH, omgst, omgnd, t0st, t0nd, rnt, lmdst, lmdnd, rns, rsmr
    complex(dp), allocatable :: store_at(:, :, :)
    integer :: ic_way, i, j, NH, NT, NS

    real(dp), allocatable :: rndff(:, :), rrndstps(:, :)
    integer, allocatable :: rndstps(:, :)

    integer :: rank, ierror

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

    call BC2N%construct(name="BC2N", &
                        from_file="./material_data/BC2N_tb.dat", &
                        fermi_energy=7.8461_dp)

    call random_seed()

    call random_number(klist)
    klist = klist - 0.5_dp

    call random_number(rNH)
    NH = nint(1.0_dp + 0.0_dp*rNH)
    !NH = 1 !DBG.
    call MPI_BCAST(NH, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    allocate (rndff(12, NH), rrndstps(9, NH), rndstps(9, NH))
    call random_number(rndff)
    rndff = 10.0_dp**(5.0_dp + 7.0_dp*rndff)
    call MPI_BCAST(rndff, size(rndff), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call random_number(rrndstps)
    rndstps = nint(1.0_dp + 1.0_dp*rrndstps)
    !rndstps = 2 !DBG.
    call MPI_BCAST(rndstps, size(rndstps), MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call random_number(omgst)
    omgst = 0.05_dp + 1.15_dp*omgst
    call MPI_BCAST(omgst, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call random_number(omgnd)
    omgnd = 0.05_dp + 1.15_dp*omgnd
    call MPI_BCAST(omgnd, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call random_number(t0st)
    t0st = (2*pi/omgst)*t0st
    call MPI_BCAST(t0st, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call random_number(t0nd)
    t0nd = (2*pi/omgst)*t0nd
    call MPI_BCAST(t0nd, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call random_number(lmdst)
    lmdst = 0.05_dp + 3.15_dp*lmdst
    call MPI_BCAST(omgst, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call random_number(lmdnd)
    lmdnd = 0.05_dp + 3.15_dp*lmdnd
    call MPI_BCAST(lmdnd, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call random_number(rnt)
    NT = 2**(nint(1.0_dp + 9.0_dp*rnt)) + 1
    !NT = 2**(9) + 1 !DBG.
    call MPI_BCAST(NT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call random_number(rns)
    NS = nint(1.0_dp + 14.0_dp*rns)
    !NS = 15 !DBG.
    call MPI_BCAST(NS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    call random_number(rsmr)
    call MPI_BCAST(rsmr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

    if (rank == 0) write (error_unit, "(A)") "Info:"
    if (rank == 0) write (error_unit, "(A, i0, A)") "  Number of harmonics = ", NH, "."
    if (rank == 0) write (error_unit, "(A, i0, A)") "  Number of t-points = ", NT, "."
    if (rank == 0) write (error_unit, "(A, i0, A)") "  Number of omega-harmonics = ", NS, "."

    do ic_way = -2, 4

      call tsk%build_floquet_task(crys=BC2N, &
                                  Nharm=NH, &
                                  axstart=rndff(1, :), axend=rndff(2, :), axsteps=rndstps(1, :), &
                                  pxstart=rndff(3, :), pxend=rndff(4, :), pxsteps=rndstps(2, :), &
                                  aystart=rndff(5, :), ayend=rndff(6, :), aysteps=rndstps(3, :), &
                                  pystart=rndff(7, :), pyend=rndff(8, :), pysteps=rndstps(4, :), &
                                  azstart=rndff(9, :), azend=rndff(10, :), azsteps=rndstps(5, :), &
                                  pzstart=rndff(11, :), pzend=rndff(12, :), pzsteps=rndstps(6, :), &
                                  omegastart=omgst, omegaend=omgnd, omegasteps=rndstps(7, 1), &
                                  t0start=t0st, t0end=t0nd, t0steps=rndstps(8, 1), &
                                  lambdastart=lmdst, lambdaend=lmdnd, lambdasteps=rndstps(9, 1), &
                                  delta_smr=rsmr, &
                                  Nt=NT, Ns=NS, htk_calc_method=ic_way)

      if ((rank == 0) .and. (ic_way == -2)) write (error_unit, "(A, i0, A)") "  Number of cont. variables = ", tsk%cdims%rank(), "."
      if ((rank == 0) .and. (ic_way == -2)) write (error_unit, "(A)") "  Shape = "
      if ((rank == 0) .and. (ic_way == -2)) write (error_unit, *) tsk%cdims%shape()
      if ((rank == 0) .and. (ic_way == -2)) write (error_unit, "(A, i0, A)") "  Number of cont. variable permutations = ", &
        tsk%cdims%size(), "."

      call tsk%sample(BC2N, klist, store_at, parallelization="none")

      if ((rank == 0) .and. (ic_way == -2)) write (error_unit, "(A, i0, A)") "  Size of store_at = ", size(store_at), "."

      if (size(store_at) /= tsk%idims%size()*tsk%cdims%size()) then
        allocate (error)
        return
      endif

      do i = 1, tsk%idims%size()
        do j = 1, tsk%cdims%size()
          if ((abs(real(store_at(1, i, j), dp)) > inf) .or. &
              (real(store_at(1, i, j), dp) /= real(store_at(1, i, j), dp))) then
            allocate (error)
            return
          endif
        enddo
      enddo

    enddo

    deallocate (rndff, rrndstps, rndstps, store_at)

  end subroutine randomized_input_parameters

end module C_Randomized_Suite
