program BC2N_Current

  use MPI
  use OMP_LIB

  use Floquet_kinds, only: dp
  use Currents, only: currents_calc_tsk
  use WannInt, only: crystal

  implicit none

  !==========================================!
  !                                          !
  ! In this example we calculate the         !
  ! Fourier transform of the current,        !
  ! j^l(\lambda) generated in BC2N if the    !
  ! system is subjected to a periodic        !
  ! driving \bm{E}(t) = E_x cos(w*t) e_x,    !
  ! where E_x = 5.2 \times 10^8 V/m and      !
  ! \hbar w = 2.5 eV.                        !
  ! We consider method #2 to calculate       !
  ! H(k, t):                                 !
  ! Velocity gauge: A*p + A*A*P2W            !
  ! interaction. The smearing in the         !
  ! resonant delta function is set at        !
  ! 0.1 eV.                                  !
  !                                          !
  !==========================================!

  type(crystal) :: BC2N

  type(currents_calc_tsk)  :: tsk

  complex(dp), allocatable :: store_at(:, :)

  !Declarations for MPI parallelization.
  integer :: nProcs, rank, ierror

  !Auxiliary
  integer :: out, i

  !Lambda: start, end, steps
  real(dp) :: lambda
  real(dp), parameter :: start = 0.0_dp, end = 10.0_dp
  integer, parameter :: steps = 257

  !Sampling points discretization.
  integer :: kpart(3) = [101, 101, 101]

  call MPI_INIT(ierror)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  !DEFINE TIGHT-BINDING CRYSTAL
  !OF BC2N.
  call BC2N%construct(name="BC2N", &
                      from_file="./material_data/BC2N_tb.dat", &
                      fermi_energy=7.8461_dp)
  !END DEFINE TIGHT-BINDING CRYSTAL.

  call tsk%build_floquet_task(crys=BC2N, &
                              Nharm=1, &
                              axstart=[5.2E8_dp], axend=[5.2E8_dp], axsteps=[1], &
                              pxstart=[0.0_dp], pxend=[0.0_dp], pxsteps=[1], &
                              aystart=[0.0_dp], ayend=[0.0_dp], aysteps=[1], &
                              pystart=[0.0_dp], pyend=[0.0_dp], pysteps=[1], &
                              azstart=[0.0_dp], azend=[0.0_dp], azsteps=[1], &
                              pzstart=[0.0_dp], pzend=[0.0_dp], pzsteps=[1], &
                              omegastart=2.5_dp, omegaend=2.5_dp, omegasteps=1, &
                              lambdastart=start, lambdaend=end, lambdasteps=steps, &
                              t0start=0.0_dp, t0end=0.0_dp, t0steps=1, &
                              delta_smr=0.01_dp, &
                              Nt=257, Ns=24, htk_calc_method=3)

  call tsk%sample(crys=BC2N, &
                  kpart=kpart, &
                  store_at=store_at, &
                  parallelization="MPI+OMP")
  !Right now, store_at holds the Fourier transform of the current, j^l(\lambda),
  !in units of A (those of the momentum dH/dk - i[\xi, H] (eV*A)) times a delta function
  !with its argument having units of energy (1/eV).
  !We have to divide by the cell volume in A^3 and
  !the number of sampling points, so the units are 1/A^2.
  store_at = store_at/(product(kpart)*BC2N%cell_volume())

  !Finally, we multiply by q and pass from C/A^2 to C/m^2 by multiplying
  !by 10^20.
  store_at = store_at*16.02176634_dp

  if (rank == 0) open (newunit=out, action="write", file="BC2N_curr_x.dat", status="unknown")
  do i = 1, steps
    lambda = start + (end - start)*real(i - 1, dp)/real(steps - 1, dp)
    if (rank == 0) write (unit=out, fmt="(F15.8, 2(2xES15.8))") lambda, &
      real(store_at(1, i), dp), aimag(store_at(1, i))
  enddo
  if (rank == 0) close (unit=out)

  if (rank == 0) open (newunit=out, action="write", file="BC2N_curr_y.dat", status="unknown")
  do i = 1, steps
    lambda = start + (end - start)*real(i - 1, dp)/real(steps - 1, dp)
    if (rank == 0) write (unit=out, fmt="(F15.8, 2(2xES15.8))") lambda, &
      real(store_at(2, i), dp), aimag(store_at(2, i))
  enddo
  if (rank == 0) close (unit=out)

  if (rank == 0) open (newunit=out, action="write", file="BC2N_curr_z.dat", status="unknown")
  do i = 1, steps
    lambda = start + (end - start)*real(i - 1, dp)/real(steps - 1, dp)
    if (rank == 0) write (unit=out, fmt="(F15.8, 2(2xES15.8))") lambda, &
      real(store_at(3, i), dp), aimag(store_at(3, i))
  enddo
  if (rank == 0) close (unit=out)

  !j^l(t) is related to j^l(\labmda) by
  !j^l(t) = \int_{-\infty}^{\infty}d\lambda j^l(\lambda) e^{-i\lambda t}.

  call MPI_FINALIZE(ierror)

end program BC2N_Current
