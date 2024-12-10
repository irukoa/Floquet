module C_Functionality_Suite
  use, intrinsic :: iso_fortran_env, only: error_unit
  use MPI
  use OMP_LIB

  use Floquet_kinds, only: dp
  use Currents, only: currents_calc_tsk
  use WannInt, only: crystal

  use testdrive, only: error_type

  implicit none
  private

  real(dp) :: tol = 1.0E5_dp

  public :: get_current_of_BC2N

contains

  subroutine get_current_of_BC2N(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: BC2N
    type(currents_calc_tsk) :: tsk

    complex(dp) :: test(3), reference(-2:4, 3)
    integer :: kpart(3)
    complex(dp), allocatable :: store_at(:, :)
    integer :: ic_way

    reference(-2, :) = [cmplx(-0.68103541231949760_dp, 46.611180003265154_dp, dp), &
                        cmplx(0.10437311866695376_dp, -6.1668012171582856_dp, dp), &
                        cmplx(0.23722429791929187_dp, 5.7509544340385643_dp, dp)]

    reference(-1, :) = [cmplx(-0.29307595102112804_dp, 0.13885565528784205_dp, dp), &
                        cmplx(-7.48012773389838681E-002_dp, 2.47111487422850308E-002_dp, dp), &
                        cmplx(4.18132920773606545E-003_dp, 2.88245581086062734E-002_dp, dp)]

    reference(0, :) = [cmplx(2.3093859571793360_dp, -0.27778280813603101_dp, dp), &
                       cmplx(-0.47888510405820983_dp, -0.10687710560454694_dp, dp), &
                       cmplx(-0.14629353030178333_dp, 0.11757330280150616_dp, dp)]

    reference(1, :) = [cmplx(2.7740531868548701_dp, 0.14962894148476527_dp, dp), &
                       cmplx(-1.1137552636514965_dp, -0.62142636839955601_dp, dp), &
                       cmplx(-0.19978988247804791_dp, -0.15721755511764982_dp, dp)]

    reference(2, :) = [cmplx(1.9766151258990288_dp, -1.5339332987718250_dp, dp), &
                       cmplx(-1.2821663112964357_dp, -0.83921932331889137_dp, dp), &
                       cmplx(-0.22490376519340990_dp, -0.35729957620212965_dp, dp)]

    reference(3, :) = [cmplx(1.8018998814423435_dp, -1.5377626483618942_dp, dp), &
                       cmplx(-1.2332076419111702_dp, -0.83590823248728774_dp, dp), &
                       cmplx(-0.21308534763173001_dp, -0.34578979543778865_dp, dp)]

    reference(4, :) = [cmplx(1.8078323518182549_dp, -1.5055042908757781_dp, dp), &
                       cmplx(-1.2371411101622742_dp, -0.83643901794426512_dp, dp), &
                       cmplx(-0.21382739414861010_dp, -0.34454664746726749_dp, dp)]

    kpart = [10, 10, 1]

    call BC2N%construct(name="BC2N", &
                        from_file="./material_data/BC2N_tb.dat", &
                        fermi_energy=7.8461_dp)

    do ic_way = -2, 4

      call tsk%build_floquet_task(crys=BC2N, &
                                  Nharm=1, &
                                  axstart=[5.2E8_dp], axend=[5.2E8_dp], axsteps=[1], &
                                  pxstart=[0.0_dp], pxend=[0.0_dp], pxsteps=[1], &
                                  aystart=[0.0_dp], ayend=[0.0_dp], aysteps=[1], &
                                  pystart=[0.0_dp], pyend=[0.0_dp], pysteps=[1], &
                                  azstart=[0.0_dp], azend=[0.0_dp], azsteps=[1], &
                                  pzstart=[0.0_dp], pzend=[0.0_dp], pzsteps=[1], &
                                  omegastart=0.5_dp, omegaend=0.5_dp, omegasteps=1, &
                                  lambdastart=2.5_dp, lambdaend=2.5_dp, lambdasteps=1, &
                                  t0start=0.0_dp, t0end=0.0_dp, t0steps=1, &
                                  delta_smr=0.2_dp, &
                                  Nt=513, Ns=12, htk_calc_method=ic_way)

      call tsk%sample(BC2N, kpart, store_at, parallelization="none")

      test = store_at(:, 1)

      if (sum(abs(test - reference(ic_way, :))) > tol*epsilon(1.0_dp)) then
        allocate (error)
        return
      endif

    enddo

  end subroutine get_current_of_BC2N

end module C_Functionality_Suite
