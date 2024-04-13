program BC2N_Driving_Scan

  use MPI
  use OMP_LIB

  use Floquet_kinds, only: dp
  use Floquet_defs, only: cmplx_0, cmplx_i, pi
  use Quasienergies, only: quasienergies_calc_tsk
  use WannInt, only: crystal, diagonalize

  use MPI_Partition, only: MPI_task_partition

  implicit none

  !==========================================!
  !                                          !
  ! In this example we calculate the         !
  ! quasienergy spectrum of BC2N in the      !
  ! point k = (1/2, 1/4, 0) when the system  !
  ! is subjected to a periodic driving       !
  ! \bm{E}(t) = E_x cos(w*t) e_x, where      !
  ! \hbar w = 0.5 eV and we scan the         !
  ! logarithmically spaced range             !
  ! E_x\in[10^7, 10^11].                     !
  ! We consider different methods to         !
  ! calculate H(k, t):                       !
  ! #-1: Length gauge: "no intraband" appr.  !
  ! #0: Velocity gauge: "no curvature" appr. !
  ! #1: Velocity gauge: A*p interaction.     !
  ! #2: Velocity gauge: A*p + A*A*P2W        !
  ! interaction.                             !
  ! #3: Velocity gauge: A*p + A*A*P2W +      !
  !                     A*A*A*P3W            !
  ! #4: Velocity gauge: A*p + A*A*P2W +      !
  !                     A*A*A*P3W +          !
  !                     A*A*A*A*P4W          !
  ! interaction.                             !
  !                                          !
  !==========================================!

  !Declarations for MPI parallelization.
  integer :: nProcs, rank, ierror
  integer, allocatable :: counts(:), displs(:)
  complex(dp), allocatable :: local_results(:, :)

  !Declarations to define tight-binding crystal.
  type(crystal) :: BC2N
  complex(dp), allocatable :: HW(:, :), rot(:, :)
  real(dp), allocatable :: eig(:)

  !Declarations to define quasienergy calculation task.
  integer, parameter :: Ns = 1000, & !Number of discretization points of driving field.
                        Nt = 513     !Number of discretization points of [t_0, t_0 + T].

  real(dp) :: klist(3, 1)

  type(quasienergies_calc_tsk) :: tsk
  real(dp) :: omega, forc, min, max, fzp
  complex(dp), allocatable :: store_at(:, :, :)
  complex(dp), allocatable :: results(:, :)

  !Auxiliary.
  integer :: i, j
  integer :: out
  integer :: method
  character(len=400) :: num

  call MPI_INIT(ierror)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)

  !DEFINE TIGHT-BINDING CRYSTAL
  !OF BC2N.
  call BC2N%construct(name="BC2N", &
                      from_file="./material_data/BC2N_tb.dat", &
                      fermi_energy=7.8461_dp)
  !END DEFINE TIGHT-BINDING CRYSTAL.

  !Choose omega, kpt and sampling range.
  omega = 0.5_dp !In eV.
  klist(:, 1) = [0.5_dp, 0.25_dp, 0.0_dp] !A point in the S-X line.
  min = 7.0_dp   !10E7 V/m.
  max = 11.0_dp  !10E11 V/m.

  !GET QUASIENERGY SPECTRUM
  !FOR METHODS #-1, #0, #1, #2, #3, #4.
  do method = -1, 4

    call MPI_task_partition(task_size=Ns, nodes=nProcs, &
                            counts=counts, displs=displs)

    allocate (local_results(displs(rank) + 1:displs(rank) + counts(rank), &
                            BC2N%num_bands()))
    allocate (results(Ns, BC2N%num_bands()))

    local_results = cmplx_0

    !Distribute the driving field amplitude sampling between MPI ranks and OpenMP threads.
!$OMP     PARALLEL PRIVATE(i, forc, tsk, store_at) SHARED(local_results, klist, displs, counts, rank)
!$OMP     DO
    do i = displs(rank) + 1, displs(rank) + counts(rank)

      !Sampling.
      forc = 10.0_dp**(min + (max - min)*real(i - 1, dp)/real(Ns - 1, dp)) !In V/m.

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
                                  Nt=Nt, htk_calc_method=method)

      call tsk%sample(BC2N, klist, store_at, parallelization="none")

      local_results(i, :) = store_at(1, :, 1) !Retrieve quasienergy eigenvalues.

    enddo
!$OMP     END DO
!$OMP     END PARALLEL

    !Gather to root node.
    do i = 1, BC2N%num_bands()
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
    if (rank == 0) open (newunit=out, action="write", file="BC2N_QS-method_"//trim(adjustl(num))//".dat", status="unknown")
    do j = 1, BC2N%num_bands()
      do i = 1, Ns
        forc = 10.0_dp**(min + (max - min)*real(i - 1, dp)/real(Ns - 1, dp)) !In V/m.
        if (rank == 0) write (unit=out, fmt="(E15.7, 2x, F15.7)") forc, &
          real(results(i, j), dp)/omega
      enddo
      if (rank == 0) write (unit=out, fmt="(a)") ""
    enddo
    if (rank == 0) close (unit=out)

    deallocate (results)

  enddo
  !END GET QUASIENERGY SPECTRUM.

  !GET REPRESENTATIVE TRACEBACK.
  allocate (HW(BC2N%num_bands(), BC2N%num_bands()), &
            rot(BC2N%num_bands(), BC2N%num_bands()), &
            eig(BC2N%num_bands()))
  HW = BC2N%hamiltonian(kpt=klist(:, 1))
  call diagonalize(matrix=HW, P=rot, eig=eig)
  if (rank == 0) open (newunit=out, action="write", file="BC2N_QS-traceback.dat", status="unknown")
  do i = 1, BC2N%num_bands()
    fzp = real(cmplx_i* &
               log(exp(-2*pi*cmplx_i*eig(i)/omega)) &
               /(2*pi), dp)
    if (rank == 0) write (unit=out, fmt="(A, F10.7, A, i0, A, i0)") &
      "E(t) = 0 quasienergy: ", fzp, &
      " Band label: ", i, &
      " Representative: ", nint(eig(i)/omega - fzp)
  enddo
  if (rank == 0) close (unit=out)
  deallocate (HW, rot, eig)
  !END GET REPRESENTATIVE TRACEBACK.

  call MPI_FINALIZE(ierror)

end program BC2N_Driving_Scan
