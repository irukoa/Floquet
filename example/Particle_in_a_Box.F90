program Particle_in_a_Box

  use MPI
  use OMP_LIB

  use Floquet_kinds, only: dp
  use Floquet_defs, only: cmplx_0, cmplx_1, cmplx_i, pi
  use Quasienergies, only: quasienergies_calc_tsk
  use WannInt, only: crystal

  use MPI_Partition, only: MPI_task_partition

  implicit none

  !==========================================!
  !                                          !
  ! In this example we calculate the         !
  ! quasienergy spectrum of a particle in a  !
  ! box subjected to a periodic driving      !
  ! E(t) = E cos(w*t), where                 !
  ! \hbar w = 95(e_2 - e_1)/100 and we scan  !
  ! the range E_x\in[0, -10\hbar w/(q * a)]. !
  ! We consider different methods to         !
  ! calculate H(k, t):                       !
  ! #-1: Length gauge: E*r interaction.      !
  ! #0: Velocity gauge: exponential form.    !
  ! #1: Velocity gauge: A*p interaction.     !
  !                                          !
  !==========================================!

  !Declarations for MPI parallelization.
  integer :: nProcs, rank, ierror
  integer, allocatable :: counts(:), displs(:)
  complex(dp), allocatable :: local_results(:, :)

  !Declarations to define tight-binding crystal.
  type(crystal) :: PiaB

  real(dp), parameter :: a = 1.0_dp, & !Lattice constant in A.
                         const = 9.40075413240204_dp !\hbar^2 * pi^2 / (8 * M) In eV*A^2.

  integer, parameter :: nbnds = 20
  real(dp) :: dlb(3, 3), eig(nbnds)
  integer :: rpts(1, 3)
  complex(dp) :: hr(1, nbnds, nbnds), &
                 rr(1, nbnds, nbnds, 3)

  !Declarations to define quasienergy calculation task.
  integer, parameter :: Ns = 1000, & !Number of discretization points of driving field.
                        Nt = 513     !Number of discretization points of [t_0, t_0 + T].

  real(dp) :: klist(3, 1)

  type(quasienergies_calc_tsk) :: tsk
  real(dp) :: omega, forc, max_forc, fzp
  complex(dp), allocatable :: store_at(:, :, :)
  complex(dp) :: results(NS, nbnds)

  !Auxiliary.
  integer :: i, j, n, m
  integer :: out
  integer :: method
  character(len=400) :: num

  call MPI_INIT(ierror)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  !DEFINE TIGHT-BINDING CRYSTAL
  !OF THE PARTICLE IN A BOX.
  dlb(1, :) = [1.0_dp, 0.0_dp, 0.0_dp]*a
  dlb(2, :) = [0.0_dp, 1.0_dp, 0.0_dp]*a
  dlb(3, :) = [0.0_dp, 0.0_dp, 1.0_dp]*a

  rpts(1, :) = [0, 0, 0]

  hr = cmplx_0
  rr = cmplx_0

  do n = 1, nbnds
    eig(n) = (const/(a**2))*(real(n, dp))**2
    hr(1, n, n) = cmplx_1*eig(n)
    do m = 1, nbnds
      if (mod(n + m, 2) == 1) &
        rr(1, n, m, 1) = cmplx_1* &
                         (-16.0_dp*a/(pi)**2)*real(m, dp)*real(n, dp)/ &
                         ((real(m, dp))**2 - (real(n, dp))**2)**2
    enddo
  enddo

  call PiaB%construct(name="Particle_in_a_Box", &
                      direct_lattice_basis=dlb, &
                      num_bands=nbnds, &
                      R_points=rpts, &
                      tunnellings=hr, &
                      dipoles=rr, &
                      fermi_energy=0.0_dp)
  !END DEFINE TIGHT-BINDING CRYSTAL.

  !Choose omega and sampling range.
  omega = 0.95_dp*(eig(2) - eig(1))
  max_forc = 10.0_dp
  klist(:, 1) = [0.0_dp, 0.0_dp, 0.0_dp]

  !GET QUASIENERGY SPECTRUM
  !FOR METHODS #-1, #0, #1.
  do method = -1, 1

    call MPI_task_partition(task_size=Ns, nodes=nProcs, &
                            counts=counts, displs=displs)

    allocate (local_results(displs(rank) + 1:displs(rank) + counts(rank), &
                            nbnds))

    local_results = cmplx_0

    !Distribute the driving field amplitude sampling between MPI ranks and OpenMP threads.
!$OMP     PARALLEL PRIVATE(i, forc, tsk, store_at) SHARED(local_results, klist, displs, counts, rank)
!$OMP     DO
    do i = displs(rank) + 1, displs(rank) + counts(rank)

      !Sampling.
      forc = (max_forc*omega/(1.0E-10_dp*a))*(real(i - 1, dp)/real(Ns - 1, dp)) !In V/m.

      call tsk%build_floquet_task(crys=PiaB, &
                                  Nharm=1, &
                                  axstart=[forc], axend=[forc], axsteps=[1], &
                                  pxstart=[0.0_dp], pxend=[0.0_dp], pxsteps=[1], &
                                  aystart=[0.0_dp], ayend=[0.0_dp], aysteps=[1], &
                                  pystart=[0.0_dp], pyend=[0.0_dp], pysteps=[1], &
                                  azstart=[0.0_dp], azend=[0.0_dp], azsteps=[1], &
                                  pzstart=[0.0_dp], pzend=[0.0_dp], pzsteps=[1], &
                                  omegastart=omega, omegaend=omega, omegasteps=1, &
                                  t0start=0.0_dp, t0end=0.0_dp, t0steps=1, &
                                  Nt=Nt, htk_calc_method=method)

      call tsk%sample(PiaB, klist, store_at, parallelization="none")

      local_results(i, :) = store_at(1, :, 1) !Retrieve quasienergy eigenvalues.

    enddo
!$OMP     END DO
!$OMP     END PARALLEL

    !Gather to root node.
    do i = 1, nbnds
      call MPI_ALLGATHERV(local_results(:, i), &
                          size(local_results(:, i)), &
                          MPI_COMPLEX16, &
                          results(:, i), &
                          counts, &
                          displs, &
                          MPI_COMPLEX16, &
                          MPI_COMM_WORLD, &
                          ierror)
    enddo

    deallocate (counts, displs, local_results)

    !Print to files.
    write (num, fmt="(i0)") method
    if (rank == 0) open (newunit=out, action="write", file="PiaB_QS-method_"//trim(adjustl(num))//".dat", status="unknown")
    do j = 1, nbnds
      do i = 1, Ns
        if (rank == 0) write (unit=out, fmt="(F15.7, 2x, F15.7)") max_forc*(real(i - 1, dp)/real(Ns - 1, dp)), &
          real(results(i, j), dp)/omega
      enddo
      if (rank == 0) write (unit=out, fmt="(a)") ""
    enddo
    if (rank == 0) close (unit=out)

  enddo
  !END GET QUASIENERGY SPECTRUM.

  !GET REPRESENTATIVE TRACEBACK.
  if (rank == 0) open (newunit=out, action="write", file="PiaB_QS-traceback.dat", status="unknown")
  do i = 1, nbnds
    fzp = real(cmplx_i* &
               log(exp(-2*pi*cmplx_i*eig(i)/omega)) &
               /(2*pi), dp)
    if (rank == 0) write (unit=out, fmt="(A, F10.7, A, i0, A, i0)") &
      "E(t) = 0 quasienergy: ", fzp, &
      " Band label: ", i, &
      " Representative: ", nint(eig(i)/omega - fzp)
  enddo
  if (rank == 0) close (unit=out)
  !END GET REPRESENTATIVE TRACEBACK.

  call MPI_FINALIZE(ierror)

end program Particle_in_a_Box
