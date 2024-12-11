module Quasienergies

  use SsTC_driver, only: container, task_specifier, crystal, &
    diagonalize, expsh, logu
  use Floquet_kinds, only: dp
  use Floquet_defs, only: cmplx_0, cmplx_i, pi
  use Floquet_Time_Dependent_defs, only: E_field, A_field, &
    htk_l_no_intra_appr, htk_v_no_curv_appr, &
    htk_v_2_terms, htk_v_3_terms, &
    htk_v_4_terms, htk_v_5_terms, &
    htk_diagonal_TB
  use Floquet_Wannier_interpolation, only: wannier_momentum, &
    wannier_2nd_momentum, &
    wannier_3rd_momentum, &
    wannier_4th_momentum, &
    NADHW => NAD

  implicit none
  private

  public :: quasienergies_calc_tsk

  type, extends(task_specifier) :: quasienergies_calc_tsk
    private
    integer :: Nt = 513
    integer :: c_way = 1
    logical :: f_initialized = .false.
  contains
    private
    procedure, pass(self), public :: build_floquet_task => constructor
    procedure, pass(self), public :: ntsteps => get_Nt
    procedure, pass(self), public :: htk_calc_method => get_c_way
    procedure, pass(self), public :: is_floquet_initialized => f_initialized
  end type

contains

  subroutine constructor(self, crys, &
                         Nharm, &
                         axstart, axend, axsteps, &
                         pxstart, pxend, pxsteps, &
                         aystart, ayend, aysteps, &
                         pystart, pyend, pysteps, &
                         azstart, azend, azsteps, &
                         pzstart, pzend, pzsteps, &
                         omegastart, omegaend, omegasteps, &
                         t0start, t0end, t0steps, &
                         Nt, htk_calc_method)

    class(quasienergies_calc_tsk), intent(out) :: self

    type(crystal), intent(in)  :: crys

    integer, intent(in) :: Nharm
    real(dp), dimension(Nharm), intent(in) :: axstart, axend, &
                                              pxstart, pxend, &
                                              aystart, ayend, &
                                              pystart, pyend, &
                                              azstart, azend, &
                                              pzstart, pzend

    integer, dimension(Nharm), intent(in) :: axsteps, pxsteps, &
                                             aysteps, pysteps, &
                                             azsteps, pzsteps

    real(dp), intent(in) :: omegastart, omegaend, &
                            t0start, t0end

    integer, intent(in) :: omegasteps, t0steps

    integer, optional, intent(in) :: Nt, htk_calc_method

    real(dp), dimension(6*Nharm + 2) :: start, end
    integer, dimension(6*Nharm + 2) :: steps
    integer :: iharm

    if (.not. crys%initialized()) error stop "Floquet: Error #1: crys must be initialized."

    if (Nharm < 1) error stop "Floquet: Error #1: Nharm must be a positive integer."

    if (any(axsteps < 1)) error stop "Floquet: Error #1: axsteps must be a positive integer array."
    if (any(pxsteps < 1)) error stop "Floquet: Error #1: pxsteps must be a positive integer array."
    if (any(aysteps < 1)) error stop "Floquet: Error #1: aysteps must be a positive integer array."
    if (any(pysteps < 1)) error stop "Floquet: Error #1: pysteps must be a positive integer array."
    if (any(azsteps < 1)) error stop "Floquet: Error #1: azsteps must be a positive integer array."
    if (any(pzsteps < 1)) error stop "Floquet: Error #1: pzsteps must be a positive integer array."

    if (omegasteps < 1) error stop "Floquet: Error #1: omegasteps must be a positive integer."
    if (t0steps < 1) error stop "Floquet: Error #1: t0steps must be a positive integer."

    do iharm = 1, Nharm
      start(6*(iharm - 1) + 1) = axstart(iharm)
      end(6*(iharm - 1) + 1) = axend(iharm)
      steps(6*(iharm - 1) + 1) = axsteps(iharm)

      start(6*(iharm - 1) + 2) = pxstart(iharm)
      end(6*(iharm - 1) + 2) = pxend(iharm)
      steps(6*(iharm - 1) + 2) = pxsteps(iharm)

      start(6*(iharm - 1) + 3) = aystart(iharm)
      end(6*(iharm - 1) + 3) = ayend(iharm)
      steps(6*(iharm - 1) + 3) = aysteps(iharm)

      start(6*(iharm - 1) + 4) = pystart(iharm)
      end(6*(iharm - 1) + 4) = pyend(iharm)
      steps(6*(iharm - 1) + 4) = pysteps(iharm)

      start(6*(iharm - 1) + 5) = azstart(iharm)
      end(6*(iharm - 1) + 5) = azend(iharm)
      steps(6*(iharm - 1) + 5) = azsteps(iharm)

      start(6*(iharm - 1) + 6) = pzstart(iharm)
      end(6*(iharm - 1) + 6) = pzend(iharm)
      steps(6*(iharm - 1) + 6) = pzsteps(iharm)
    enddo

    start(6*Nharm + 1) = t0start
    end(6*Nharm + 1) = t0end
    steps(6*Nharm + 1) = t0steps

    start(6*Nharm + 2) = omegastart
    end(6*Nharm + 2) = omegaend
    steps(6*Nharm + 2) = omegasteps

    call self%construct(name="quasienergies", &
                        int_ind=[crys%num_bands()], &
                        cont_data_start=start, &
                        cont_data_end=end, &
                        cont_data_steps=steps, &
                        calculator=quasienergy)

    if (present(Nt)) then
      if (Nt < 1) error stop "Floquet: Error #1: Nt must be a positive integer."
      self%Nt = Nt
    endif

    if (present(htk_calc_method)) then
      select case (htk_calc_method)
      case (-2:4)
        self%c_way = htk_calc_method
      case default
        error stop "Floquet: Error #1: htk_calc_method must be an integer in the [-2, 4] range."
      end select
    endif

    self%f_initialized = .true.

  end subroutine

  pure elemental integer function get_Nt(self)
    class(quasienergies_calc_tsk), intent(in) :: self
    if (.not. self%f_initialized) error stop "Floquet: Error #2: Floquet task is not initialized."
    get_Nt = self%Nt
  end function get_Nt

  pure elemental integer function get_c_way(self)
    class(quasienergies_calc_tsk), intent(in) :: self
    if (.not. self%f_initialized) error stop "Floquet: Error #2: Floquet task is not initialized."
    get_c_way = self%c_way
  end function get_c_way

  pure elemental logical function f_initialized(self)
    class(quasienergies_calc_tsk), intent(in) :: self
    f_initialized = self%f_initialized
  end function f_initialized

  function quasienergy(self, crys, k, other) result(qsn)
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    real(dp), intent(in) :: k(3)
    class(*), optional, intent(in) :: other

    complex(dp) :: qsn(self%idims%size(), self%cdims%size())

    integer :: r_mem, r_arr(self%cdims%rank()), &
               i_mem

    real(dp) :: omega, t0, dt, &
                tper, q(3), &
                quasi(crys%num_bands()), &
                eig(crys%num_bands())

    real(dp), allocatable :: amplitudes(:, :), &
                             phases(:, :)

    integer :: Nharm, iharm, it, i

    complex(dp) :: H_TK(crys%num_bands(), crys%num_bands()), &
                   expH_TK(crys%num_bands(), crys%num_bands()), &
                   tev(crys%num_bands(), crys%num_bands()), &
                   hf(crys%num_bands(), crys%num_bands()), &
                   rot(crys%num_bands(), crys%num_bands()), &
                   AW(crys%num_bands(), crys%num_bands(), 3), &
                   NAD(crys%num_bands(), crys%num_bands(), 3), &
                   AH(crys%num_bands(), crys%num_bands(), 3), &
                   HW(crys%num_bands(), crys%num_bands()), &
                   HH(crys%num_bands(), crys%num_bands()), &
                   P1W(crys%num_bands(), crys%num_bands(), 3), &
                   P2W(crys%num_bands(), crys%num_bands(), 3, 3), &
                   P3W(crys%num_bands(), crys%num_bands(), 3, 3, 3), &
                   P4W(crys%num_bands(), crys%num_bands(), 3, 3, 3, 3)

    type(container) :: DHW, DDHW, DDDHW, DDDDHW, &
                       DAW, DDAW, DDDAW

    qsn = cmplx_0

    HW = crys%hamiltonian(kpt=k)
    AW = crys%berry_connection(kpt=k)

    select type (self)
    type is (quasienergies_calc_tsk)
      if (.not. self%f_initialized) error stop "Floquet: Error #2: Floquet task is not initialized."
      !Gather required quantities depending on the method to calculate H(k, t).
      select case (self%c_way)
      case (-2) !Velocity gauge: tight-binding approximation.
      case (-1, 0)
        !-1: Length gauge: no-intraband appr.
        ! 0: Velocity gauge: no k-curvature appr.
        call diagonalize(matrix=HW, P=rot, D=HH, eig=eig)
        NAD = NADHW(eig=eig, &
                    rot=rot, &
                    DHW=crys%hamiltonian(kpt=k, derivative=1), &
                    deg_thr=0.01_dp, deg_offset=0.04_dp)
        do i = 1, 3
          AH(:, :, i) = matmul(matmul(transpose(conjg(rot)), AW(:, :, i)), rot) + &
                        cmplx_i*NAD(:, :, i)
        enddo
      case (1)  !Velocity gauge: 2 terms in expansion.
        DHW = crys%hamiltonian(kpt=k, derivative=1)
        P1W = wannier_momentum(HW=HW, &
                               DHW=DHW, &
                               AW=AW)
      case (2)  !Velocity gauge: 3 terms in expansion.
        DHW = crys%hamiltonian(kpt=k, derivative=1)
        P1W = wannier_momentum(HW=HW, &
                               DHW=DHW, &
                               AW=AW)
        DDHW = crys%hamiltonian(kpt=k, derivative=2)
        DAW = crys%berry_connection(kpt=k, derivative=1)
        P2W = wannier_2nd_momentum(HW=HW, &
                                   DHW=DHW, &
                                   DDHW=DDHW, &
                                   AW=AW, &
                                   DAW=DAW, &
                                   P1W=P1W)
      case (3)  !Velocity gauge: 4 terms in expansion.
        DHW = crys%hamiltonian(kpt=k, derivative=1)
        P1W = wannier_momentum(HW=HW, &
                               DHW=DHW, &
                               AW=AW)
        DDHW = crys%hamiltonian(kpt=k, derivative=2)
        DAW = crys%berry_connection(kpt=k, derivative=1)
        P2W = wannier_2nd_momentum(HW=HW, &
                                   DHW=DHW, &
                                   DDHW=DDHW, &
                                   AW=AW, &
                                   DAW=DAW, &
                                   P1W=P1W)
        DDDHW = crys%hamiltonian(kpt=k, derivative=3)
        DDAW = crys%berry_connection(kpt=k, derivative=2)
        P3W = wannier_3rd_momentum(HW=HW, &
                                   DHW=DHW, &
                                   DDHW=DDHW, &
                                   DDDHW=DDDHW, &
                                   AW=AW, &
                                   DAW=DAW, &
                                   DDAW=DDAW, &
                                   P1W=P1W, &
                                   P2W=P2W)
      case (4)  !Velocity gauge: 5 terms in expansion.
        DHW = crys%hamiltonian(kpt=k, derivative=1)
        P1W = wannier_momentum(HW=HW, &
                               DHW=DHW, &
                               AW=AW)
        DDHW = crys%hamiltonian(kpt=k, derivative=2)
        DAW = crys%berry_connection(kpt=k, derivative=1)
        P2W = wannier_2nd_momentum(HW=HW, &
                                   DHW=DHW, &
                                   DDHW=DDHW, &
                                   AW=AW, &
                                   DAW=DAW, &
                                   P1W=P1W)
        DDDHW = crys%hamiltonian(kpt=k, derivative=3)
        DDAW = crys%berry_connection(kpt=k, derivative=2)
        P3W = wannier_3rd_momentum(HW=HW, &
                                   DHW=DHW, &
                                   DDHW=DDHW, &
                                   DDDHW=DDDHW, &
                                   AW=AW, &
                                   DAW=DAW, &
                                   DDAW=DDAW, &
                                   P1W=P1W, &
                                   P2W=P2W)
        DDDDHW = crys%hamiltonian(kpt=k, derivative=4)
        DDDAW = crys%berry_connection(kpt=k, derivative=3)
        P4W = wannier_4th_momentum(HW=HW, &
                                   DHW=DHW, &
                                   DDHW=DDHW, &
                                   DDDHW=DDDHW, &
                                   DDDDHW=DDDDHW, &
                                   AW=AW, &
                                   DAW=DAW, &
                                   DDAW=DDAW, &
                                   DDDAW=DDDAW, &
                                   P1W=P1W, &
                                   P2W=P2W, &
                                   P3W=P3W)
      end select

      Nharm = (self%cdims%rank() - 2)/6
      allocate (amplitudes(Nharm, 3), phases(Nharm, 3))

      !Loop for each external variable.
      do r_mem = 1, self%cdims%size()

        r_arr = self%cdims%ind(r_mem)

        omega = self%cdt(var=self%cdims%rank(), step=r_arr(self%cdims%rank()))

        t0 = self%cdt(var=self%cdims%rank() - 1, step=r_arr(self%cdims%rank() - 1))

        do iharm = 1, Nharm
          amplitudes(iharm, 1) = self%cdt(var=6*(iharm - 1) + 1, step=r_arr(6*(iharm - 1) + 1))
          phases(iharm, 1) = self%cdt(var=6*(iharm - 1) + 2, step=r_arr(6*(iharm - 1) + 2))
          amplitudes(iharm, 2) = self%cdt(var=6*(iharm - 1) + 3, step=r_arr(6*(iharm - 1) + 3))
          phases(iharm, 2) = self%cdt(var=6*(iharm - 1) + 4, step=r_arr(6*(iharm - 1) + 4))
          amplitudes(iharm, 3) = self%cdt(var=6*(iharm - 1) + 5, step=r_arr(6*(iharm - 1) + 5))
          phases(iharm, 3) = self%cdt(var=6*(iharm - 1) + 6, step=r_arr(6*(iharm - 1) + 6))
        enddo

        dt = (2*pi/omega)/real(self%Nt - 1, dp)
        tev = cmplx_0
        do i = 1, crys%num_bands()
          tev(i, i) = cmplx(1.0_dp, 0.0_dp, dp)
        enddo

        do it = 2, self%Nt

          tper = t0 + dt*real(it - 1, dp) !In eV^-1.

          !Calculate H(k, t) depending on the chosen method.
          select case (self%c_way)

          case (-2) !Velocity gauge: tight-binding approximation.
            !The Hamiltonian is defined in the Wannier state basis
            !as H_nm(R, t) = H0_nm(R)e^(i*q*A*rtb_nm(R)/hbar),
            !with rtb_nm(R) = A_{nn}(R=0) - A_{mm}(R=0) - R.

            q = A_field(amplitudes, phases, omega, tper, t0) !In A^-1

            H_TK = htk_diagonal_TB(crys, k, q)

          case (-1) !Length gauge: no-intraband appr.
            !No intraband approximation.
            !The interaction term is V(t) = e*E(t)*r.
            !We calculate H_nm(k) and r_nm(k) = A_nm(1-delta_nm)
            !in the Hamiltonian basis.
            !The diagonal Berry connection does not appear never,
            !so this approach is inexact but U(1)-gauge covariant.

            q = E_field(amplitudes, phases, omega, tper) !In eV * A^-1

            H_TK = htk_l_no_intra_appr(eig, AH, q) !In eV.

          case (0) !Velocity gauge: no k-curvature appr.
            !The interaction Hamiltonian is
            !exp(-i*q*A*r/hbar)*H_0*exp(i*q*A*r/hbar).
            !We calculate H_nm(k) and r_nm(k) = A_nm(1-delta_nm)
            !in the Hamiltonian basis.
            !The diagonal Berry connection does not appear never,
            !so this approach is inexact but U(1)-gauge covariant.

            q = A_field(amplitudes, phases, omega, tper, t0) !In A^-1

            H_TK = htk_v_no_curv_appr(HH, AH, q)

          case (1) !Velocity gauge: 2 terms in expansion.
            !The interaction term is V(t) = e*Q(t)*p/m
            !Q(t) = -\int dt E(t).
            !We calculate the momentum:
            !p_{nm}^a(k) = dH_{0\;nm}/dk^a - i[A^a, H_0]_{nm}
            !and then add the interaction term to H_0 in the Wannier basis.

            q = A_field(amplitudes, phases, omega, tper, t0) !In A^-1

            H_TK = htk_v_2_terms(HW, P1W, q) !In eV.

          case (2) !Velocity gauge: 3 terms in expansion.
            !The interaction term is V(t) = e*Q(t)*p/m + \sum_{ij} e^2*Q_j(t)*Q_i(t)*[D^j, [D^i, H]]/2\hbar^2
            !Q(t) = -\int dt E(t).

            q = A_field(amplitudes, phases, omega, tper, t0) !In A^-1

            H_TK = htk_v_3_terms(HW, P1W, P2W, q) !In eV.

          case (3) !Velocity gauge: 4 terms in expansion.
            !The interaction term is V(t) = e*Q(t)*p/m + \sum_{ij} e^2*Q_j(t)*Q_i(t)*[D^j, [D^i, H]]/2\hbar^2 +
            !\sum_{ijl} e^3*Q_l(t)*Q_j(t)*Q_i(t)*[D^l, [D^j, [D^i, H]]]/6\hbar^3
            !Q(t) = -\int dt E(t).

            q = A_field(amplitudes, phases, omega, tper, t0) !In A^-1

            H_TK = htk_v_4_terms(HW, P1W, P2W, P3W, q) !In eV.

          case (4) !Velocity gauge: 4 terms in expansion.
            !The interaction term is V(t) = e*Q(t)*p/m + \sum_{ij} e^2*Q_j(t)*Q_i(t)*[D^j, [D^i, H]]/2\hbar^2 +
            !\sum_{ijl} e^3*Q_l(t)*Q_j(t)*Q_i(t)*[D^l, [D^j, [D^i, H]]]/6\hbar^3 +
            !\sum_{ijlo} e^4*Q_o(t)*Q_l(t)*Q_j(t)*Q_i(t)*[D^o, [D^l, [D^j, [D^i, H]]]]/24\hbar^4
            !Q(t) = -\int dt E(t).

            q = A_field(amplitudes, phases, omega, tper, t0) !In A^-1

            H_TK = htk_v_5_terms(HW, P1W, P2W, P3W, P4W, q) !In eV.

          end select

          !Get exponential for each time instant.
          expH_TK = expsh(-cmplx_i*dt*H_TK)

          !Accumulate.
          tev = matmul(tev, expH_TK)

        enddo

        !Get effective Floquet Hamiltonian in:
        ! - Hamiltonian basis if c_way = -1, 0.
        ! - Wannier basis if c_way = -2, 1, 2, 3, 4.
        hf = cmplx_i*omega*logu(tev)/(2*pi)
        call diagonalize(matrix=hf, P=rot, eig=quasi)

        do i_mem = 1, self%idims%size()
          qsn(i_mem, r_mem) = quasi(i_mem)
        enddo

      enddo

      deallocate (amplitudes, phases)

    end select

  end function quasienergy

end module Quasienergies
