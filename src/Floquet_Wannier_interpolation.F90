module Floquet_Wannier_interpolation

  use Floquet_kinds, only: dp
  use Floquet_defs, only: cmplx_0, cmplx_i
  use SsTC_driver, only: container, crystal

  implicit none
  private

  public :: wannier_momentum
  public :: wannier_2nd_momentum
  public :: NAD

contains

  pure function wannier_momentum(HW, DHW, AW)

    complex(dp), intent(in) :: HW(:, :), AW(:, :, :)
    type(container), intent(in) :: DHW

    complex(dp), dimension(size(HW(:, 1)), size(HW(1, :)), 3) :: rshpDHW, &
                                                                 wannier_momentum

    integer :: i

    rshpDHW = reshape(DHW%cdp_storage, shape(rshpDHW))

    do i = 1, 3
      wannier_momentum(:, :, i) = rshpDHW(:, :, i) - &
                                  cmplx_i*(matmul(AW(:, :, i), HW) - matmul(HW, AW(:, :, i)))
    enddo

  end function wannier_momentum

  pure function wannier_2nd_momentum(HW, DHW, DDHW, AW, DAW, PW)
    complex(dp), intent(in) :: HW(:, :), AW(:, :, :), PW(:, :, :)
    type(container), intent(in) :: DHW, DDHW, DAW

    complex(dp) :: rshpDHW(size(HW(:, 1)), size(HW(1, :)), 3), &
                   rshpDAW(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   rshpDDHW(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   wannier_2nd_momentum(size(HW(:, 1)), size(HW(1, :)), 3, 3)

    integer :: i, j

    rshpDHW = reshape(DHW%cdp_storage, shape(rshpDHW))
    rshpDDHW = reshape(DDHW%cdp_storage, shape(rshpDDHW))
    rshpDAW = reshape(DAW%cdp_storage, shape(rshpDAW))

    do i = 1, 3
      do j = 1, 3
        wannier_2nd_momentum(:, :, i, j) = rshpDDHW(:, :, i, j) - &
                                           cmplx_i*(matmul(rshpDAW(:, :, i, j), HW) + matmul(AW(:, :, i), rshpDHW(:, :, j)) - &
                                                    matmul(HW, rshpDAW(:, :, i, j)) - matmul(rshpDHW(:, :, j), AW(:, :, i))) - &
                                           cmplx_i*(matmul(AW(:, :, j), PW(:, :, i)) - matmul(PW(:, :, i), AW(:, :, j)))
      enddo
    enddo

  end function wannier_2nd_momentum

  pure function NAD(eig, rot, DHW, deg_thr, deg_offset)
    real(dp), intent(in) :: eig(:), deg_thr, deg_offset
    complex(dp), intent(in) :: rot(:, :)
    type(container), intent(in) :: DHW

    complex(dp) :: rshpDHW(size(eig), size(eig), 3), &
                   NAD(size(eig), size(eig), 3)

    integer :: i, n, m

    rshpDHW = reshape(DHW%cdp_storage, shape(rshpDHW))

    do i = 1, 3
      NAD(:, :, i) = matmul(matmul(transpose(conjg(rot)), rshpDHW(:, :, i)), rot)
      do n = 1, size(eig)
        do m = 1, size(eig)
          if (abs(eig(n) - eig(m)) < deg_thr) then
            NAD(n, m, i) = cmplx_0
          else
            NAD(n, m, i) = NAD(n, m, i)*((eig(m) - eig(n))/((eig(m) - eig(n))**2 + (deg_offset)**2))
          endif
        enddo
      enddo
    enddo

  end function NAD

end module Floquet_Wannier_interpolation
