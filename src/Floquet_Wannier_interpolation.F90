module Floquet_Wannier_interpolation

  use Floquet_kinds, only: dp
  use Floquet_defs, only: cmplx_0, cmplx_i
  use SsTC_driver, only: container

  implicit none
  private

  public :: wannier_momentum
  public :: wannier_2nd_momentum
  public :: wannier_3rd_momentum
  public :: wannier_4th_momentum
  public :: NAD

contains

  pure function wannier_momentum(HW, DHW, &
                                 AW)

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

  pure function Dwannier_momentum(HW, DHW, DDHW, &
                                  AW, DAW)

    complex(dp), intent(in) :: HW(:, :), AW(:, :, :)
    type(container), intent(in) :: DHW, DDHW, DAW

    complex(dp) :: rshpDHW(size(HW(:, 1)), size(HW(1, :)), 3), &
                   rshpDDHW(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   rshpDAW(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   Dwannier_momentum(size(HW(:, 1)), size(HW(1, :)), 3, 3)

    integer :: i, j

    rshpDHW = reshape(DHW%cdp_storage, shape(rshpDHW))
    rshpDDHW = reshape(DDHW%cdp_storage, shape(rshpDDHW))
    rshpDAW = reshape(DAW%cdp_storage, shape(rshpDAW))

    do j = 1, 3
      do i = 1, 3
        Dwannier_momentum(:, :, i, j) = rshpDDHW(:, :, i, j) - &
                                        cmplx_i*(matmul(rshpDAW(:, :, i, j), HW) + matmul(AW(:, :, i), rshpDHW(:, :, j)) - &
                                                 matmul(rshpDHW(:, :, j), AW(:, :, i)) - matmul(HW, rshpDAW(:, :, i, j)))
      enddo
    enddo

  end function Dwannier_momentum

  pure function DDwannier_momentum(HW, DHW, DDHW, DDDHW, &
                                   AW, DAW, DDAW)

    complex(dp), intent(in) :: HW(:, :), AW(:, :, :)
    type(container), intent(in) :: DHW, DDHW, DDDHW, DAW, DDAW

    complex(dp) :: rshpDHW(size(HW(:, 1)), size(HW(1, :)), 3), &
                   rshpDDHW(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   rshpDDDHW(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3), &
                   rshpDAW(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   rshpDDAW(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3), &
                   DDwannier_momentum(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3)

    integer :: i, j, l

    rshpDHW = reshape(DHW%cdp_storage, shape(rshpDHW))
    rshpDDHW = reshape(DDHW%cdp_storage, shape(rshpDDHW))
    rshpDDDHW = reshape(DDDHW%cdp_storage, shape(rshpDDDHW))
    rshpDAW = reshape(DAW%cdp_storage, shape(rshpDAW))
    rshpDDAW = reshape(DDAW%cdp_storage, shape(rshpDDAW))

    do l = 1, 3
      do j = 1, 3
        do i = 1, 3
          DDwannier_momentum(:, :, i, j, l) = rshpDDDHW(:, :, i, j, l) - &
                                              cmplx_i*(matmul(rshpDDAW(:, :, i, j, l), HW) + &
                                                       matmul(rshpDAW(:, :, i, j), rshpDHW(:, :, l)) + &
                                                       matmul(rshpDAW(:, :, i, l), rshpDHW(:, :, j)) + &
                                                       matmul(AW(:, :, i), rshpDDHW(:, :, j, l)) - &
                                                       matmul(rshpDDHW(:, :, j, l), AW(:, :, i)) - &
                                                       matmul(rshpDHW(:, :, j), rshpDAW(:, :, i, l)) - &
                                                       matmul(rshpDHW(:, :, l), rshpDAW(:, :, i, j)) - &
                                                       matmul(HW, rshpDDAW(:, :, i, j, l)))
        enddo
      enddo
    enddo
  end function DDwannier_momentum

  pure function wannier_2nd_momentum(HW, DHW, DDHW, &
                                     AW, DAW, &
                                     P1W)

    complex(dp), intent(in) :: HW(:, :), AW(:, :, :), &
                               P1W(:, :, :)
    type(container), intent(in) :: DHW, DDHW, DAW

    complex(dp) :: DP1W(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   wannier_2nd_momentum(size(HW(:, 1)), size(HW(1, :)), 3, 3)

    integer :: i, j

    DP1W = Dwannier_momentum(HW, DHW, DDHW, &
                             AW, DAW)

    do j = 1, 3
      do i = 1, 3
        wannier_2nd_momentum(:, :, i, j) = DP1W(:, :, i, j) - &
                                           cmplx_i*(matmul(AW(:, :, j), P1W(:, :, i)) - matmul(P1W(:, :, i), AW(:, :, j)))
      enddo
    enddo

  end function wannier_2nd_momentum

  pure function Dwannier_2nd_momentum(HW, DHW, DDHW, DDDHW, &
                                      AW, DAW, DDAW, &
                                      P1W, DP1W)

    complex(dp), intent(in) :: HW(:, :), AW(:, :, :), &
                               P1W(:, :, :), DP1W(:, :, :, :)
    type(container), intent(in) :: DHW, DDHW, DDDHW, DAW, DDAW

    complex(dp) :: rshpDHW(size(HW(:, 1)), size(HW(1, :)), 3), &
                   rshpDDHW(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   rshpDDDHW(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3), &
                   rshpDAW(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   rshpDDAW(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3), &
                   Dwannier_2nd_momentum(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3)

    integer :: i, j, l

    rshpDHW = reshape(DHW%cdp_storage, shape(rshpDHW))
    rshpDDHW = reshape(DDHW%cdp_storage, shape(rshpDDHW))
    rshpDDDHW = reshape(DDDHW%cdp_storage, shape(rshpDDDHW))
    rshpDAW = reshape(DAW%cdp_storage, shape(rshpDAW))
    rshpDDAW = reshape(DDAW%cdp_storage, shape(rshpDDAW))

    do l = 1, 3
      do j = 1, 3
        do i = 1, 3
          Dwannier_2nd_momentum(:, :, i, j, l) = rshpDDDHW(:, :, i, j, l) - &
                                                 cmplx_i*(matmul(rshpDDAW(:, :, i, j, l), HW) + &
                                                          matmul(rshpDAW(:, :, i, j), rshpDHW(:, :, l)) + &
                                                          matmul(rshpDAW(:, :, i, l), rshpDHW(:, :, j)) + &
                                                          matmul(AW(:, :, i), rshpDDHW(:, :, j, l)) - &
                                                          matmul(rshpDDHW(:, :, j, l), AW(:, :, i)) - &
                                                          matmul(rshpDHW(:, :, j), rshpDAW(:, :, i, l)) - &
                                                          matmul(rshpDHW(:, :, l), rshpDAW(:, :, i, j)) - &
                                                          matmul(HW, rshpDDAW(:, :, i, j, l))) - &
                                                 cmplx_i*(matmul(rshpDAW(:, :, j, l), P1W(:, :, i)) + &
                                                          matmul(AW(:, :, j), DP1W(:, :, i, l)) - &
                                                          matmul(DP1W(:, :, i, l), AW(:, :, j)) - &
                                                          matmul(P1W(:, :, i), rshpDAW(:, :, j, l)))
        enddo
      enddo
    enddo

  end function Dwannier_2nd_momentum

  pure function wannier_3rd_momentum(HW, DHW, DDHW, DDDHW, &
                                     AW, DAW, DDAW, &
                                     P1W, P2W)

    complex(dp), intent(in) :: HW(:, :), AW(:, :, :), &
                               P1W(:, :, :), P2W(:, :, :, :)
    type(container), intent(in) :: DHW, DDHW, DDDHW, DAW, DDAW

    complex(dp) :: DP1W(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   DP2W(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3), &
                   wannier_3rd_momentum(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3)

    integer :: i, j, l

    DP1W = Dwannier_momentum(HW, DHW, DDHW, &
                             AW, DAW)
    DP2W = Dwannier_2nd_momentum(HW, DHW, DDHW, DDDHW, &
                                 AW, DAW, DDAW, &
                                 P1W, DP1W)

    do l = 1, 3
      do j = 1, 3
        do i = 1, 3
          wannier_3rd_momentum(:, :, i, j, l) = DP2W(:, :, i, j, l) - &
                                                cmplx_i*(matmul(AW(:, :, l), P2W(:, :, i, j)) - &
                                                         matmul(P2W(:, :, i, j), AW(:, :, l)))
        enddo
      enddo
    enddo

  end function wannier_3rd_momentum

  pure function Dwannier_3rd_momentum(HW, DHW, DDHW, DDDHW, DDDDHW, &
                                      AW, DAW, DDAW, DDDAW, &
                                      P1W, DP1W, DDP1W, P2W, DP2W)

    complex(dp), intent(in) :: HW(:, :), AW(:, :, :), &
                               P1W(:, :, :), DP1W(:, :, :, :), DDP1W(:, :, :, :, :), &
                               P2W(:, :, :, :), DP2W(:, :, :, :, :)
    type(container), intent(in) :: DHW, DDHW, DDDHW, DDDDHW, DAW, DDAW, DDDAW

    complex(dp) :: rshpDHW(size(HW(:, 1)), size(HW(1, :)), 3), &
                   rshpDDHW(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   rshpDDDHW(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3), &
                   rshpDDDDHW(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3, 3), &
                   rshpDAW(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   rshpDDAW(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3), &
                   rshpDDDAW(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3, 3), &
                   Dwannier_3rd_momentum(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3, 3)

    integer :: i, j, l, o

    rshpDHW = reshape(DHW%cdp_storage, shape(rshpDHW))
    rshpDDHW = reshape(DDHW%cdp_storage, shape(rshpDDHW))
    rshpDDDHW = reshape(DDDHW%cdp_storage, shape(rshpDDDHW))
    rshpDDDDHW = reshape(DDDDHW%cdp_storage, shape(rshpDDDDHW))
    rshpDAW = reshape(DAW%cdp_storage, shape(rshpDAW))
    rshpDDAW = reshape(DDAW%cdp_storage, shape(rshpDDAW))
    rshpDDDAW = reshape(DDDAW%cdp_storage, shape(rshpDDDAW))

    do o = 1, 3
      do l = 1, 3
        do j = 1, 3
          do i = 1, 3
            Dwannier_3rd_momentum(:, :, i, j, l, o) = rshpDDDDHW(:, :, i, j, l, o) - &
                                                      cmplx_i*(matmul(rshpDDDAW(:, :, i, j, l, o), HW) + &
                                                               matmul(rshpDDAW(:, :, i, j, l), rshpDHW(:, :, o)) + &
                                                               matmul(rshpDDAW(:, :, i, j, o), rshpDHW(:, :, l)) + &
                                                               matmul(rshpDAW(:, :, i, j), rshpDDHW(:, :, l, o)) + &
                                                               matmul(rshpDDAW(:, :, i, l, o), rshpDHW(:, :, j)) + &
                                                               matmul(rshpDAW(:, :, i, l), rshpDDHW(:, :, j, o)) + &
                                                               matmul(rshpDAW(:, :, i, o), rshpDDHW(:, :, j, l)) + &
                                                               matmul(AW(:, :, i), rshpDDDHW(:, :, j, l, o)) - &
                                                               matmul(rshpDDDHW(:, :, j, l, o), AW(:, :, i)) - &
                                                               matmul(rshpDDHW(:, :, j, l), rshpDAW(:, :, i, o)) - &
                                                               matmul(rshpDDHW(:, :, j, o), rshpDAW(:, :, i, l)) - &
                                                               matmul(rshpDHW(:, :, j), rshpDDAW(:, :, i, l, o)) - &
                                                               matmul(rshpDDHW(:, :, l, o), rshpDAW(:, :, i, j)) - &
                                                               matmul(rshpDHW(:, :, l), rshpDDAW(:, :, i, j, o)) - &
                                                               matmul(rshpDHW(:, :, o), rshpDDAW(:, :, i, j, l)) - &
                                                               matmul(HW, rshpDDDAW(:, :, i, j, l, o))) - &
                                                      cmplx_i*(matmul(rshpDDAW(:, :, j, l, o), P1W(:, :, i)) + &
                                                               matmul(rshpDAW(:, :, j, l), DP1W(:, :, i, o)) + &
                                                               matmul(rshpDAW(:, :, j, o), DP1W(:, :, i, l)) + &
                                                               matmul(AW(:, :, j), DDP1W(:, :, i, l, o)) - &
                                                               matmul(DDP1W(:, :, i, l, o), AW(:, :, j)) - &
                                                               matmul(DP1W(:, :, i, l), rshpDAW(:, :, j, o)) - &
                                                               matmul(DP1W(:, :, i, o), rshpDAW(:, :, j, l)) - &
                                                               matmul(P1W(:, :, i), rshpDDAW(:, :, j, l, o))) - &
                                                      cmplx_i*(matmul(rshpDAW(:, :, l, o), P2W(:, :, i, j)) + &
                                                               matmul(AW(:, :, l), DP2W(:, :, i, j, o)) - &
                                                               matmul(DP2W(:, :, i, j, o), AW(:, :, l)) - &
                                                               matmul(P2W(:, :, i, j), rshpDAW(:, :, l, o)))
          enddo
        enddo
      enddo
    enddo

  end function Dwannier_3rd_momentum

  pure function wannier_4th_momentum(HW, DHW, DDHW, DDDHW, DDDDHW, &
                                     AW, DAW, DDAW, DDDAW, &
                                     P1W, P2W, P3W)

    complex(dp), intent(in) :: HW(:, :), AW(:, :, :), &
                               P1W(:, :, :), P2W(:, :, :, :), P3W(:, :, :, :, :)
    type(container), intent(in) :: DHW, DDHW, DDDHW, DDDDHW, DAW, DDAW, DDDAW

    complex(dp) :: DP1W(size(HW(:, 1)), size(HW(1, :)), 3, 3), &
                   DDP1W(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3), &
                   DP2W(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3), &
                   DP3W(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3, 3), &
                   wannier_4th_momentum(size(HW(:, 1)), size(HW(1, :)), 3, 3, 3, 3)

    integer :: i, j, l, o

    DP1W = Dwannier_momentum(HW, DHW, DDHW, &
                             AW, DAW)
    DDP1W = DDwannier_momentum(HW, DHW, DDHW, DDDHW, &
                               AW, DAW, DDAW)
    DP2W = Dwannier_2nd_momentum(HW, DHW, DDHW, DDDHW, &
                                 AW, DAW, DDAW, &
                                 P1W, DP1W)
    DP3W = Dwannier_3rd_momentum(HW, DHW, DDHW, DDDHW, DDDDHW, &
                                 AW, DAW, DDAW, DDDAW, &
                                 P1W, DP1W, DDP1W, P2W, DP2W)

    do o = 1, 3
      do l = 1, 3
        do j = 1, 3
          do i = 1, 3
            wannier_4th_momentum(:, :, i, j, l, o) = DP3W(:, :, i, j, l, o) - &
                                                     cmplx_i*(matmul(AW(:, :, o), P3W(:, :, i, j, l)) - &
                                                              matmul(P3W(:, :, i, j, l), AW(:, :, o)))
          enddo
        enddo
      enddo
    enddo

  end function wannier_4th_momentum

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
