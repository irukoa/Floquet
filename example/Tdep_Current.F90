program Tdep_Current

  use MPI
  use OMP_LIB

  use Floquet_kinds, only: dp
  use Tdep_Currents, only: tdep_currents_calc_tsk
  use WannInt, only: crystal

  implicit none

  type(crystal) :: CRYS

  type(tdep_currents_calc_tsk)  :: tsk

  complex(dp), allocatable :: store_at(:, :)

  !Declarations for MPI parallelization.
  integer :: nProcs, rank, ierror

  !Auxiliary
  integer :: out, i

  !Simulation parameters.
  real(dp), parameter :: polarization_vec(3) = [1.0_dp, 0.0_dp, 0.0_dp]
  real(dp), parameter :: omega = 3.5_dp !Energy corresponding to the driving field frequency in eV.
  real(dp), parameter :: E = 5.2E8_dp !Electric field amplitude in V/m.
  integer, parameter :: n_periods = 2 !Number of periods to consider in time discretization.
  !Constant 2*pi*(\hbar/e)*10^{15}, to pass from frequency in eV to period in fs.
  real(dp), parameter :: factor = 4.135667694_dp

  !Time: start, end, steps
  real(dp) :: time
  !Starting and ending times in fs.
  real(dp) :: start = 0.0_dp, &
              end = factor*n_periods
  integer :: steps = 257

  !Sampling points discretization.
  integer :: kpart(3) = [101, 101, 101]

  call MPI_INIT(ierror)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  !DEFINE TIGHT-BINDING CRYSTAL
  call CRYS%construct(name="CRYS", &
                      from_file="./material_data/BC2N_tb.dat", &
                      fermi_energy=7.8461_dp)
  !END DEFINE TIGHT-BINDING CRYSTAL.

  call tsk%build_floquet_task(crys=CRYS, &
                              Nharm=1, &
                              axstart=[E*polarization_vec(1)], axend=[E*polarization_vec(1)], axsteps=[1], &
                              pxstart=[0.0_dp], pxend=[0.0_dp], pxsteps=[1], &
                              aystart=[E*polarization_vec(2)], ayend=[E*polarization_vec(2)], aysteps=[1], &
                              pystart=[0.0_dp], pyend=[0.0_dp], pysteps=[1], &
                              azstart=[E*polarization_vec(3)], azend=[E*polarization_vec(3)], azsteps=[1], &
                              pzstart=[0.0_dp], pzend=[0.0_dp], pzsteps=[1], &
                              omegastart=omega, omegaend=omega, omegasteps=1, &
                              tstart=start, tend=end, tsteps=steps, &
                              t0start=0.0_dp, t0end=0.0_dp, t0steps=1, &
                              Nt=257, Ns=24, htk_calc_method=3)

  call tsk%sample(crys=CRYS, &
                  kpart=kpart, &
                  store_at=store_at, &
                  parallelization="MPI+OMP")
  !Right now, store_at holds the current, j^l(t),
  !in units of eV*A (those of the momentum dH/dk - i[\xi, H]).
  !We have to divide by the cell volume in A^3 and
  !the number of sampling points, so the units are eV/A^2.
  store_at = store_at/(product(kpart)*CRYS%cell_volume())

  !We have to multiply by q, and pass from eV/(A^2) to J/(A^2). We also pass to
  !J/(m^2) by multiplying by 10^20.
  store_at = store_at*16.02176634_dp
  !Lastly, we multiply by q and divide by \hbar to pass to C*J/(J*s*m^2) = Amp/m^2.
  store_at = store_at*1.602176634_dp/1.05457182_dp
  store_at = store_at*(10.0_dp**(15.0_dp))

  if (rank == 0) open (newunit=out, action="write", file="CRYS_tdcurr_x.dat", status="unknown")
  do i = 1, steps
    time = start + (end - start)*real(i - 1, dp)/real(steps - 1, dp)
    if (rank == 0) write (unit=out, fmt="(F15.8, 2(2xES15.8))") time, &
      real(store_at(1, i), dp), aimag(store_at(1, i))
  enddo
  if (rank == 0) close (unit=out)

  if (rank == 0) open (newunit=out, action="write", file="CRYS_tdcurr_y.dat", status="unknown")
  do i = 1, steps
    time = start + (end - start)*real(i - 1, dp)/real(steps - 1, dp)
    if (rank == 0) write (unit=out, fmt="(F15.8, 2(2xES15.8))") time, &
      real(store_at(2, i), dp), aimag(store_at(2, i))
  enddo
  if (rank == 0) close (unit=out)

  if (rank == 0) open (newunit=out, action="write", file="CRYS_tdcurr_z.dat", status="unknown")
  do i = 1, steps
    time = start + (end - start)*real(i - 1, dp)/real(steps - 1, dp)
    if (rank == 0) write (unit=out, fmt="(F15.8, 2(2xES15.8))") time, &
      real(store_at(3, i), dp), aimag(store_at(3, i))
  enddo
  if (rank == 0) close (unit=out)

  call MPI_FINALIZE(ierror)

end program Tdep_Current
