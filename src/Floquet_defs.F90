module Floquet_defs

  use Floquet_kinds, only: dp

  implicit none
  private

  public :: cmplx_0, cmplx_1, cmplx_i
  public :: pi

  real(dp), parameter :: pi = acos(-1.0_dp)

  complex(dp), parameter :: cmplx_0 = cmplx(0.0_dp, 0.0_dp, dp), &
                            cmplx_1 = cmplx(1.0_dp, 0.0_dp, dp), &
                            cmplx_i = cmplx(0.0_dp, 1.0_dp, dp)

end module Floquet_defs
