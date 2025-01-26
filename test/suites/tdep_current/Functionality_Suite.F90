module TC_Functionality_Suite
  use, intrinsic :: iso_fortran_env, only: error_unit
  use MPI
  use OMP_LIB

  use Floquet_kinds, only: dp
  use Tdep_Currents, only: tdep_currents_calc_tsk
  use WannInt, only: crystal

  use testdrive, only: error_type

  implicit none
  private

  real(dp) :: tol = 1.0E5_dp

  public :: get_tdep_current_of_BC2N

contains

  subroutine get_tdep_current_of_BC2N(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: BC2N
    type(tdep_currents_calc_tsk) :: tsk

    complex(dp) :: test(3), reference(-2:4, 3)
    integer :: kpart(3)
    complex(dp), allocatable :: store_at(:, :)
    integer :: ic_way

    reference(-2, :) = [cmplx(-8.3871662933876525_dp, -6.08067481253897479E-002_dp, dp), &
                        cmplx(-21.560472264835447_dp, 3.44491467984179012E-002_dp, dp), &
                        cmplx(13.597224953364135_dp, 3.09384700768779945E-003_dp, dp)]

    reference(-1, :) = [cmplx(4.8506364039405714_dp, -0.20628607453554884_dp, dp), &
                        cmplx(-4.55462322809279674E-002_dp, -7.45985858677241254E-003_dp, dp), &
                        cmplx(6.80952351440239356E-002_dp, 2.12133864980772278E-003_dp, dp)]

    reference(0, :) = [cmplx(6.3328201937558282_dp, 0.19904163955110224_dp, dp), &
                       cmplx(-3.1030378411273158_dp, 4.82596225558391823E-003_dp, dp), &
                       cmplx(1.85890376490947772E-002_dp, -1.40550821061841808E-004_dp, dp)]

    reference(1, :) = [cmplx(-10.318714541016668_dp, 0.20478076088298910_dp, dp), &
                       cmplx(-1.7892034580162257_dp, 5.91207248422148321E-003_dp, dp), &
                       cmplx(1.5004008384783700_dp, 8.41596261399119592E-004_dp, dp)]

    reference(2, :) = [cmplx(-8.7929151022388314_dp, 0.19855891012879126_dp, dp), &
                       cmplx(-1.9680407654801257_dp, 4.45606966445461373E-003_dp, dp), &
                       cmplx(1.3889466891907676_dp, 9.95293531382268875E-004_dp, dp)]

    reference(3, :) = [cmplx(-7.5909189066773441_dp, 0.19912285038754682_dp, dp), &
                       cmplx(-1.9232484094358568_dp, 4.60234960180202265E-003_dp, dp), &
                       cmplx(1.3291186846448024_dp, 1.18361158908286950E-003_dp, dp)]

    reference(4, :) = [cmplx(-7.6121729722526457_dp, 0.19923179533027691_dp, dp), &
                       cmplx(-1.9391895407048714_dp, 4.59157279180482519E-003_dp, dp), &
                       cmplx(1.3363283259475964_dp, 1.18964443956968567E-003_dp, dp)]

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
                                  t0start=0.0_dp, t0end=0.0_dp, t0steps=1, &
                                  tstart=1.5_dp, tend=1.5_dp, tsteps=1, &
                                  Nt=513, Ns=12, htk_calc_method=ic_way)

      call tsk%sample(BC2N, kpart, store_at, parallelization="none")

      test = store_at(:, 1)

      if (sum(abs(test - reference(ic_way, :))) > tol*epsilon(1.0_dp)) then
        allocate (error)
        return
      endif

    enddo

  end subroutine get_tdep_current_of_BC2N

end module TC_Functionality_Suite
