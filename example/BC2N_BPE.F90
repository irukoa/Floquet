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
  ! Fourier series component of the current, !
  ! j^l_{\lambda=0} (BPE) generated in BC2N  !
  ! if the system is subjected to a periodic !
  ! driving \bm{E}(t) = E_x cos(w*t) e_x,    !
  ! where E_x = 5.2 \times 10^8 V/m and      !
  ! \hbar w = 2.5 eV.                        !
  ! We consider method #2 to calculate       !
  ! H(k, t):                                 !
  ! Velocity gauge: A*p + A*A*P2W            !
  ! interaction.                             !
  !                                          !
  !==========================================!

  type(crystal) :: BC2N

  type(currents_calc_tsk)  :: tsk

  complex(dp), allocatable :: store_at(:, :)

  !Declarations for MPI parallelization.
  integer :: nProcs, rank, ierror

  !Auxiliary
  integer :: out, i

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
                              lambdastart=0.0_dp, lambdaend=0.0_dp, lambdasteps=1, &
                              t0start=0.0_dp, t0end=0.0_dp, t0steps=1, &
                              FS_component_calc=.true., FS_kpt_tolerance=0.001_dp, &
                              Nt=257, Ns=24, htk_calc_method=3)

  call tsk%sample(crys=BC2N, &
                  kpart=kpart, &
                  store_at=store_at, &
                  parallelization="MPI+OMP")
  !Right now, store_at holds the Fourier series component of the current, j^l_{0},
  !in units of eV*A (those of the momentum dH/dk - i[\xi, H]).
  !We have to divide by the cell volume in A^3 and
  !the number of sampling points, so the units are eV/A^2.
  store_at = store_at/(product(kpart)*BC2N%cell_volume())

  !We have to multiply by q, and pass from eV/(A^2) to J/(A^2). We also pass to
  !J/(m^2) by multiplying by 10^20.
  store_at = store_at*16.02176634_dp
  !Lastly, we multiply by q and divide by \hbar to pass to C*J/(J*s*m^2) = Amp/m^2.
  store_at = store_at*1.602176634_dp/1.05457182_dp
  store_at = store_at*(10.0_dp**(15.0_dp))

  if (rank == 0) open (newunit=out, action="write", file="BC2N_BPE.dat", status="unknown")
  if (rank == 0) write (unit=out, fmt="(A, 2(2xES15.8))") "x component: Re, Im: ", &
    real(store_at(1, 1), dp), aimag(store_at(1, 1))
  if (rank == 0) write (unit=out, fmt="(A, 2(2xES15.8))") "y component: Re, Im: ", &
    real(store_at(2, 1), dp), aimag(store_at(2, 1))
  if (rank == 0) write (unit=out, fmt="(A, 2(2xES15.8))") "z component: Re, Im: ", &
    real(store_at(3, 1), dp), aimag(store_at(3, 1))
  if (rank == 0) close (unit=out)

  !j^l(t) is related to j^l_{\labmda} by
  !j^l(t) = \sum_{\lambda = -\infty}^{\infty} j^l_{\lambda} e^{-i\lambda t}.

  call MPI_FINALIZE(ierror)

end program BC2N_Current
