module Functionality_Suite
  use, intrinsic :: iso_fortran_env, only: error_unit
  use MPI
  use OMP_LIB

  use Floquet_kinds, only: dp
  use Quasienergies, only: quasienergies_calc_tsk
  use WannInt, only: crystal

  use testdrive, only: error_type

  implicit none
  private

  real(dp) :: tol = 1.0E3_dp

  public :: get_quasienergies_of_BC2N

contains

  subroutine get_quasienergies_of_BC2N(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: BC2N
    type(quasienergies_calc_tsk) :: tsk

    real(dp) :: test(8), reference(-1:3, 8)
    real(dp) :: klist(3, 1)
    complex(dp), allocatable :: store_at(:, :, :)
    integer :: ic_way

    reference(-1, :) = [-0.246711257250613_dp, -0.135638065337363_dp, &
                        -7.816303141478401E-002_dp, -7.561611766076404E-002_dp, &
                        -3.947301771345653E-002_dp, -3.598925512842017E-004_dp, &
                        9.195174330975003E-002_dp, 0.153647684457003_dp]

    reference(0, :) = [-0.24671277831681396_dp, -0.13563858599075562_dp, &
                       -7.8162623918914806E-002_dp, -7.5615168338419955E-002_dp, &
                       -3.9473864282245932E-002_dp, -3.5803774317198485E-004_dp, &
                       9.1952130474646732E-002_dp, 0.15364697395417892_dp]

    reference(1, :) = [-0.168218482258199_dp, -0.148901734999523_dp, &
                       -0.106991792827316_dp, -0.105154763117103_dp, &
                       -9.854146120104357E-002_dp, 2.675293389500047E-002_dp, &
                       6.021747588107929E-002_dp, 0.210475870465593_dp]

    reference(2, :) = [-0.161914419546202_dp, -0.146236499523516_dp, &
                       -9.621702730945489E-002_dp, -9.164518499000121E-002_dp, &
                       -8.895083591555945E-002_dp, 4.270507188051712E-002_dp, &
                       5.524258070670914E-002_dp, 0.233900225862371_dp]

    reference(3, :) = [-0.16138564054336371_dp, -0.14262274413577328_dp, &
                       -9.5433263875139160E-002_dp, -9.0818138744458168E-002_dp, &
                       -9.0231177095024151E-002_dp, 3.9328208981623931E-002_dp, &
                       5.7088077179957196E-002_dp, 0.23095858939704494_dp]

    klist(:, 1) = [0.5_dp, 0.25_dp, 0.0_dp]

    call BC2N%construct(name="BC2N", &
                        from_file="./material_data/BC2N_tb.dat", &
                        fermi_energy=7.8461_dp)

    do ic_way = -1, 3

      call tsk%build_floquet_task(crys=BC2N, &
                                  Nharm=1, &
                                  axstart=[5.2E8_dp], axend=[5.2E8_dp], axsteps=[1], &
                                  pxstart=[0.0_dp], pxend=[0.0_dp], pxsteps=[1], &
                                  aystart=[0.0_dp], ayend=[0.0_dp], aysteps=[1], &
                                  pystart=[0.0_dp], pyend=[0.0_dp], pysteps=[1], &
                                  azstart=[0.0_dp], azend=[0.0_dp], azsteps=[1], &
                                  pzstart=[0.0_dp], pzend=[0.0_dp], pzsteps=[1], &
                                  omegastart=0.5_dp, omegaend=0.5_dp, omegasteps=1, &
                                  t0start=0.0_dp, t0end=0.0_dp, t0steps=1, &
                                  Nt=513, htk_calc_method=ic_way)

      call tsk%sample(BC2N, klist, store_at, parallelization="none")

      test = real(store_at(1, :, 1), dp)

      if (sum(abs(test - reference(ic_way, :))) > tol*epsilon(1.0_dp)) then
        allocate (error)
        return
      endif

    enddo

  end subroutine get_quasienergies_of_BC2N

end module Functionality_Suite
