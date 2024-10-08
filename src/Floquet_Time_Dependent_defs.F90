module Floquet_Time_Dependent_defs

  use SsTC_driver, only: expsh
  use Floquet_kinds, only: dp
  use Floquet_defs, only: cmplx_0, cmplx_i

  implicit none
  private

  public :: E_field
  public :: A_field
  public :: htk_l_no_intra_appr
  public :: htk_v_no_curv_appr
  public :: htk_v_2_terms
  public :: htk_v_3_terms
  public :: htk_v_4_terms
  public :: htk_v_5_terms

contains

  pure function E_field(amplitudes, phases, omega, t) result(q)

    real(dp), intent(in) :: amplitudes(:, :), phases(:, :), &
                            omega, t

    real(dp) :: q(3)

    integer :: icoord, iharm

    q = 0.0_dp

    do iharm = 1, size(amplitudes(:, 1))
      do icoord = 1, 3
        q(icoord) = q(icoord) + &
                    amplitudes(iharm, icoord)*cos(iharm*omega*t - phases(iharm, icoord))
      enddo
    enddo

    !q is now the driving electric field E(t), in V/m.
    !We want q to represent e*E(t) in eV/A, so we multiply it by
    !|e| to pass it to J/m and then divide it by
    !|e| to pass it to eV/m.
    !Lastly pass to eV/A by multiplying by 1^-10.
    q = q*1.0E-10_dp

  end function E_field

  pure function A_field(amplitudes, phases, omega, t, t0) result(q)

    real(dp), intent(in) :: amplitudes(:, :), phases(:, :), &
                            omega, t, t0

    real(dp) :: q(3)

    integer :: icoord, iharm

    q = 0.0_dp

    do iharm = 1, size(amplitudes(:, 1))
      do icoord = 1, 3
        q(icoord) = q(icoord) + &
                    amplitudes(iharm, icoord)* &
                    (sin(iharm*omega*t - phases(iharm, icoord)) - sin(iharm*omega*t0 - phases(iharm, icoord)))/ &
                    (iharm*omega)
      enddo
    enddo

    q = -q

    !q is now the integral of the driving electric field, in V/(m*eV).
    !We want q to represent e*Q(t)/m in A^-1. (E(t) = dQ(t)/dt), so
    !we have to multiply by e to obatin the integral of the driving force field in J/(m*eV)
    !then divide by |e| to obain it in 1/m. Lastly pass to 1/A by multiplying by 1^-10.
    q = q*1.0E-10_dp
    !That is, in reality q holds e*A/\hbar.

  end function A_field

  pure function htk_l_no_intra_appr(eig, AH, q) result(H_TK)
    real(dp), intent(in)    :: q(3), eig(:)
    complex(dp), intent(in) :: AH(:, :, :)
    complex(dp)             :: H_TK(size(eig), size(eig))

    integer :: n, m

    H_TK = cmplx_0
    do n = 1, size(eig)
      H_TK(n, n) = cmplx(eig(n), 0.0_dp, dp)
      do m = 1, size(eig)
        if (n == m) cycle
        H_TK(n, m) = sum(q*AH(n, m, :))
      enddo
    enddo
  end function htk_l_no_intra_appr

  function htk_v_no_curv_appr(HH, AH, q) result(H_TK)
    complex(dp), intent(in) :: HH(:, :), AH(:, :, :)
    real(dp), intent(in)    :: q(3)
    complex(dp)             :: H_TK(size(HH(:, 1)), size(HH(:, 1)))

    complex(dp) :: rtimesA(size(HH(:, 1)), size(HH(:, 1)))
    integer :: n, i

    rtimesA = cmplx_0
    do i = 1, 3
      rtimesA = rtimesA + AH(:, :, i)*q(i)
    enddo
    do n = 1, size(HH(:, 1))
      rtimesA(n, n) = cmplx_0
    enddo

    H_TK = matmul(expsh(-cmplx_i*rtimesA), &
                  matmul(HH, expsh(cmplx_i*rtimesA)))

  end function htk_v_no_curv_appr

  pure function htk_v_2_terms(HW, P1W, q) result(H_TK)
    real(dp), intent(in)    :: q(3)
    complex(dp), intent(in) :: HW(:, :), &
                               P1W(:, :, :)
    complex(dp)             :: H_TK(size(HW(:, 1)), size(HW(:, 1)))

    integer :: n, m

    H_TK = HW

    do n = 1, size(HW(:, 1))
      do m = 1, size(HW(:, 1))
        H_TK(n, m) = H_TK(n, m) + &
                     sum(q*P1W(n, m, :))
      enddo
    enddo
  end function htk_v_2_terms

  pure function htk_v_3_terms(HW, P1W, P2W, q) result(H_TK)
    real(dp), intent(in)    :: q(3)
    complex(dp), intent(in) :: HW(:, :), &
                               P1W(:, :, :), &
                               P2W(:, :, :, :)
    complex(dp)             :: H_TK(size(HW(:, 1)), size(HW(:, 1)))

    integer :: n, m, i, j

    H_TK = HW

    do m = 1, size(HW(:, 1))
      do n = 1, size(HW(:, 1))
        H_TK(n, m) = H_TK(n, m) + &
                     sum(q*P1W(n, m, :))
      enddo
    enddo

    do j = 1, 3
      do i = 1, 3
        do m = 1, size(HW(:, 1))
          do n = 1, size(HW(:, 1))
            H_TK(n, m) = H_TK(n, m) + &
                         q(i)*q(j)*P2W(n, m, i, j)/(2.0_dp)
          enddo
        enddo
      enddo
    enddo
  end function htk_v_3_terms

  pure function htk_v_4_terms(HW, P1W, P2W, P3W, q) result(H_TK)
    real(dp), intent(in)    :: q(3)
    complex(dp), intent(in) :: HW(:, :), &
                               P1W(:, :, :), &
                               P2W(:, :, :, :), &
                               P3W(:, :, :, :, :)
    complex(dp)             :: H_TK(size(HW(:, 1)), size(HW(:, 1)))

    integer :: n, m, i, j, l

    H_TK = HW

    do m = 1, size(HW(:, 1))
      do n = 1, size(HW(:, 1))
        H_TK(n, m) = H_TK(n, m) + &
                     sum(q*P1W(n, m, :))
      enddo
    enddo

    do j = 1, 3
      do i = 1, 3
        do m = 1, size(HW(:, 1))
          do n = 1, size(HW(:, 1))
            H_TK(n, m) = H_TK(n, m) + &
                         q(i)*q(j)*P2W(n, m, i, j)/(2.0_dp)
          enddo
        enddo
      enddo
    enddo

    do l = 1, 3
      do j = 1, 3
        do i = 1, 3
          do m = 1, size(HW(:, 1))
            do n = 1, size(HW(:, 1))
              H_TK(n, m) = H_TK(n, m) + &
                           q(i)*q(j)*q(l)*P3W(n, m, i, j, l)/(6.0_dp)
            enddo
          enddo
        enddo
      enddo
    enddo
  end function htk_v_4_terms

  pure function htk_v_5_terms(HW, P1W, P2W, P3W, P4W, q) result(H_TK)
    real(dp), intent(in)    :: q(3)
    complex(dp), intent(in) :: HW(:, :), &
                               P1W(:, :, :), &
                               P2W(:, :, :, :), &
                               P3W(:, :, :, :, :), &
                               P4W(:, :, :, :, :, :)
    complex(dp)             :: H_TK(size(HW(:, 1)), size(HW(:, 1)))

    integer :: n, m, i, j, l, o

    H_TK = HW

    do m = 1, size(HW(:, 1))
      do n = 1, size(HW(:, 1))
        H_TK(n, m) = H_TK(n, m) + &
                     sum(q*P1W(n, m, :))
      enddo
    enddo

    do j = 1, 3
      do i = 1, 3
        do m = 1, size(HW(:, 1))
          do n = 1, size(HW(:, 1))
            H_TK(n, m) = H_TK(n, m) + &
                         q(i)*q(j)*P2W(n, m, i, j)/(2.0_dp)
          enddo
        enddo
      enddo
    enddo

    do l = 1, 3
      do j = 1, 3
        do i = 1, 3
          do m = 1, size(HW(:, 1))
            do n = 1, size(HW(:, 1))
              H_TK(n, m) = H_TK(n, m) + &
                           q(i)*q(j)*q(l)*P3W(n, m, i, j, l)/(6.0_dp)
            enddo
          enddo
        enddo
      enddo
    enddo

    do o = 1, 3
      do l = 1, 3
        do j = 1, 3
          do i = 1, 3
            do m = 1, size(HW(:, 1))
              do n = 1, size(HW(:, 1))
                H_TK(n, m) = H_TK(n, m) + &
                             q(i)*q(j)*q(l)*q(o)*P4W(n, m, i, j, l, o)/(24.0_dp)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  end function htk_v_5_terms

end module Floquet_Time_Dependent_defs
