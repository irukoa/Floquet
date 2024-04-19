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

  real(dp) :: tol = 1.0E3_dp

  public :: get_current_of_BC2N

contains

  subroutine get_current_of_BC2N(error)
    type(error_type), allocatable, intent(out) :: error

    type(crystal) :: BC2N
    type(currents_calc_tsk) :: tsk

    complex(dp) :: test(3), reference(-1:4, 3)
    integer :: kpart(3)
    complex(dp), allocatable :: store_at(:, :)
    integer :: ic_way

    reference(-1, :) = [cmplx(-1.11799602898074358E-006_dp, 5.29692288543098641E-007_dp, dp), &
                        cmplx(-2.85344228130279038E-007_dp, 9.42655515376473648E-008_dp, dp), &
                        cmplx(1.59505050954287165E-008_dp, 1.09956962999749273E-007_dp, dp)]

    reference(0, :) = [cmplx(8.80960829612478651E-006_dp, -1.05965731863415150E-006_dp, dp), &
                       cmplx(-1.82680169699939664E-006_dp, -4.07703802507579576E-007_dp, dp), &
                       cmplx(-5.58065530020840933E-007_dp, 4.48506556707405692E-007_dp, dp)]

    reference(1, :) = [cmplx(1.05821731065935903E-005_dp, 5.70789113940297207E-007_dp, dp), &
                       cmplx(-4.24863915882681451E-006_dp, -2.37055346832106020E-006_dp, dp), &
                       cmplx(-7.62137918388549457E-007_dp, -5.99737377615546476E-007_dp, dp)]

    reference(2, :) = [cmplx(7.54018831596004016E-006_dp, -5.85149116047601719E-006_dp, dp), &
                       cmplx(-4.89107632177900565E-006_dp, -3.20136765792423773E-006_dp, dp), &
                       cmplx(-8.57939778112067813E-007_dp, -1.36298971634723531E-006_dp, dp)]

    reference(3, :) = [cmplx(6.87370255066811934E-006_dp, -5.86609896988637637E-006_dp, dp), &
                       cmplx(-4.70431381954639496E-006_dp, -3.18873684878268333E-006_dp, dp), &
                       cmplx(-8.12856092955513038E-007_dp, -1.31908338713756046E-006_dp, dp)]

    reference(4, :) = [cmplx(6.89633312918950980E-006_dp, -5.74304310179053529E-006_dp, dp), &
                       cmplx(-4.71931881012830444E-006_dp, -3.19076163461404845E-006_dp, dp), &
                       cmplx(-8.15686775774422075E-007_dp, -1.31434115397364611E-006_dp, dp)]

    kpart = [10, 10, 1]

    call BC2N%construct(name="BC2N", &
                        from_file="./material_data/BC2N_tb.dat", &
                        fermi_energy=7.8461_dp)

    do ic_way = -1, 4

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
                                  delta_smr=0.1_dp, &
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
