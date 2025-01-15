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

  public :: get_current_FT_of_BC2N
  public :: get_current_FS_of_BC2N

contains

  subroutine get_current_FT_of_BC2N(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: BC2N
    type(currents_calc_tsk) :: tsk

    complex(dp) :: test(3), reference(-2:4, 3)
    integer :: kpart(3)
    complex(dp), allocatable :: store_at(:, :)
    integer :: ic_way

    reference(-2, :) = [cmplx(-6.4572648575015128_dp, 11.256121315897941_dp, dp), &
                        cmplx(-9.4484542780277199_dp, 109.12028314163328_dp, dp), &
                        cmplx(-6.4590805685544863_dp, 13.018492711390484_dp, dp)]

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

  end subroutine get_current_FT_of_BC2N

  subroutine get_current_FS_of_BC2N(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: BC2N
    type(currents_calc_tsk) :: tsk

    complex(dp) :: test(3), reference(-2:4, 3)
    integer :: kpart(3)
    complex(dp), allocatable :: store_at(:, :)
    integer :: ic_way

    reference(-2, :) = [cmplx(-0.30489629607544844_dp, 0.20059940732210838_dp, dp), &
                        cmplx(-0.12756888019115609_dp, 0.11614522288943543_dp, dp), &
                        cmplx(6.14314150058167807E-002_dp, -2.19264376328466640E-002_dp, dp)]

    reference(-1, :) = [cmplx(-1.85911565633049315E-002_dp, -5.98241262704674471E-003_dp, dp), &
                        cmplx(-6.69139363543799943E-003_dp, 3.37767376254179477E-003_dp, dp), &
                        cmplx(-1.76825455997346639E-003_dp, 4.16662852627827304E-003_dp, dp)]

    reference(0, :) = [cmplx(0.22298209132570049_dp, 1.04136291613247511E-002_dp, dp), &
                       cmplx(-3.84272319743509810E-002_dp, 8.36117256310328627E-003_dp, dp), &
                       cmplx(-1.45700112104078065E-002_dp, 2.41482394277682325E-002_dp, dp)]

    reference(1, :) = [cmplx(0.33386128486093930_dp, 3.22554507077790695E-002_dp, dp), &
                       cmplx(-0.18209534193236807_dp, 7.98984272500733866E-002_dp, dp), &
                       cmplx(2.33605198826365105E-003_dp, 1.78241284182735563E-002_dp, dp)]

    reference(2, :) = [cmplx(0.22796013336393295_dp, -0.15204354703498885_dp, dp), &
                       cmplx(-0.18556535390910361_dp, 4.92344305953587719E-002_dp, dp), &
                       cmplx(2.24619159388056433E-003_dp, 7.10552745165181701E-003_dp, dp)]

    reference(3, :) = [cmplx(0.19743384744533568_dp, -0.10205789159356374_dp, dp), &
                       cmplx(-0.17607488847049163_dp, 5.11761828120606693E-002_dp, dp), &
                       cmplx(3.39241549714349193E-003_dp, 4.28837338995650963E-003_dp, dp)]

    reference(4, :) = [cmplx(0.20243948477128046_dp, -0.11388164792638485_dp, dp), &
                       cmplx(-0.17882287453082590_dp, 4.92406354171959493E-002_dp, dp), &
                       cmplx(4.39119811134413077E-003_dp, 3.79271617321666867E-003_dp, dp)]

    kpart = [5, 5, 1]

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
                                  FS_component_calc=.true., FS_kpt_tolerance=0.001_dp, &
                                  Nt=513, Ns=12, htk_calc_method=ic_way)

      call tsk%sample(BC2N, kpart, store_at, parallelization="none")

      test = store_at(:, 1)

      if (sum(abs(test - reference(ic_way, :))) > tol*epsilon(1.0_dp)) then
        allocate (error)
        return
      endif

    enddo

  end subroutine get_current_FS_of_BC2N

end module C_Functionality_Suite
