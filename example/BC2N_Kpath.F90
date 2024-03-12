program BC2N_Kpath

  use MPI
  use OMP_LIB

  use Floquet_kinds, only: dp
  use Quasienergies, only: quasienergies_calc_tsk
  use WannInt, only: crystal
  use SsTC_driver, only: kpath

  implicit none

  !==========================================!
  !                                          !
  ! In this example we calculate the         !
  ! quasienergy spectrum of BC2N in the      !
  ! path Gamma - X - S - Y - Gamma when the  !
  ! system is subjected to a periodic        !
  ! driving \bm{E}(t) = E_x cos(w*t) e_x,    !
  ! where E_x = 5.2 \times 10^8 V/m and      !
  ! \hbar w = 0.5 eV.                        !
  ! We consider method #2 to calculate       !
  ! H(k, t):                                 !
  ! Velocity gauge: A*p + A*A*P2W            !
  ! interaction.                             !
  !                                          !
  !==========================================!

  !Declarations for MPI parallelization.
  integer :: nProcs, rank, ierror

  !Declarations to define tight-binding crystal.
  type(crystal) :: BC2N

  !Declarations to define quasienergy calculation task.
  integer, parameter :: Nt = 513 !Number of discretization points of [t_0, t_0 + T].
  real(dp) :: vecs(5, 3)
  real(dp), allocatable :: path(:, :)

  type(quasienergies_calc_tsk) :: tsk
  real(dp) :: omega, forc
  complex(dp), allocatable :: store_at(:, :, :)

  !Auxiliary.
  integer :: i, j
  integer :: out

  call MPI_INIT(ierror)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  !DEFINE TIGHT-BINDING CRYSTAL
  !OF BC2N.
  call BC2N%construct(name="BC2N", &
                      from_file="./material_data/BC2N_tb.dat", &
                      fermi_energy=7.8461_dp)
  !END DEFINE TIGHT-BINDING CRYSTAL.

  !Choose omega, kpath vectors and field amplitude.
  omega = 0.5_dp !In eV.

  vecs(1, :) = [0.0_dp, 0.0_dp, 0.0_dp] !G
  vecs(2, :) = [0.5_dp, 0.0_dp, 0.0_dp] !X
  vecs(3, :) = [0.5_dp, 0.5_dp, 0.0_dp] !S
  vecs(4, :) = [0.0_dp, 0.5_dp, 0.0_dp] !Y
  vecs(5, :) = [0.0_dp, 0.0_dp, 0.0_dp] !G

  forc = 5.2E8_dp !In V/m.

  !Define kpath.
  path = kpath(vecs=vecs, &
               nkpts=[100, 100, 100, 100])

  !GET QUASIENERGY SPECTRUM.
  call tsk%build_floquet_task(crys=BC2N, &
                              Nharm=1, &
                              axstart=[forc], axend=[forc], axsteps=[1], &
                              pxstart=[0.0_dp], pxend=[0.0_dp], pxsteps=[1], &
                              aystart=[0.0_dp], ayend=[0.0_dp], aysteps=[1], &
                              pystart=[0.0_dp], pyend=[0.0_dp], pysteps=[1], &
                              azstart=[0.0_dp], azend=[0.0_dp], azsteps=[1], &
                              pzstart=[0.0_dp], pzend=[0.0_dp], pzsteps=[1], &
                              omegastart=omega, omegaend=omega, omegasteps=1, &
                              t0start=0.0_dp, t0end=0.0_dp, t0steps=1, &
                              Nt=Nt, htk_calc_method=2)

  call tsk%sample(crys=BC2N, &
                  klist=path, &
                  store_at=store_at, &
                  parallelization="MPI+OMP")

  !Print to files.
  if (rank == 0) open (newunit=out, action="write", file="BC2N_QS-kpath.dat", status="unknown")
  do j = 1, BC2N%num_bands()
    do i = 1, size(path(1, :))
      if (rank == 0) write (unit=out, fmt="(i0, 4(2xF15.8))") i, path(1, i), path(2, i), path(3, i), &
        real(store_at(i, j, 1), dp)/omega
    enddo
    if (rank == 0) write (unit=out, fmt="(a)") ""
  enddo
  if (rank == 0) close (unit=out)

  !END GET QUASIENERGY SPECTRUM.

  call MPI_FINALIZE(ierror)

end program BC2N_Kpath
