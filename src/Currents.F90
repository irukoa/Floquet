module Currents

  use SsTC_driver, only: container, task_specifier, crystal, &
    diagonalize, expsh, logu, dirac_delta
  use Extrapolation_Integration, only: extrapolation
  use Floquet_kinds, only: dp
  use Floquet_defs, only: cmplx_0, cmplx_1, cmplx_i, pi
  use Floquet_Time_Dependent_defs, only: E_field, A_field, &
    htk_l_no_intra_appr, htk_v_no_curv_appr, &
    htk_v_2_terms, htk_v_3_terms, &
    htk_v_4_terms, htk_v_5_terms
  use Floquet_Wannier_interpolation, only: wannier_momentum, &
    wannier_2nd_momentum, &
    wannier_3rd_momentum, &
    wannier_4th_momentum, &
    NADHW => NAD

  implicit none
  private

  public :: currents_calc_tsk

  type, extends(task_specifier) :: currents_calc_tsk
    private
    integer :: Nt = 513
    integer :: Ns = 10
    integer :: c_way = 1
    real(dp) :: delta_smr = 0.04_dp
    logical :: f_initialized = .false.
  contains
    private
    procedure, pass(self), public :: build_floquet_task => constructor
    procedure, pass(self), public :: ntsteps => get_Nt
    procedure, pass(self), public :: nsharms => get_Ns
    procedure, pass(self), public :: htk_calc_method => get_c_way
    procedure, pass(self), public :: smr => get_delta_smr
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
                         lambdastart, lambdaend, lambdasteps, &
                         delta_smr, &
                         Nt, Ns, htk_calc_method)

    class(currents_calc_tsk), intent(out) :: self

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
                            t0start, t0end, &
                            lambdastart, lambdaend

    integer, intent(in) :: omegasteps, t0steps, lambdasteps

    real(dp), optional, intent(in) :: delta_smr

    integer, optional, intent(in) :: Nt, Ns, htk_calc_method

    real(dp), dimension(6*Nharm + 3) :: start, end
    integer, dimension(6*Nharm + 3) :: steps
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

    start(6*Nharm + 3) = lambdastart
    end(6*Nharm + 3) = lambdaend
    steps(6*Nharm + 3) = lambdasteps

    call self%construct(name="currents", &
                        int_ind=[3], &
                        cont_data_start=start, &
                        cont_data_end=end, &
                        cont_data_steps=steps, &
                        calculator=current)

    if (present(delta_smr)) then
      if (delta_smr < 0.0_dp) error stop "Floquet: Error #1: delta_smr must be positive and real."
      self%delta_smr = delta_smr
    endif

    if (present(Nt)) then
      if (Nt < 1) error stop "Floquet: Error #1: Nt must be a positive integer."
      self%Nt = Nt
    endif

    if (present(Ns)) then
      if (Nt < 1) error stop "Floquet: Error #1: Ns must be a positive integer."
      self%Ns = Ns
    endif

    if (present(htk_calc_method)) then
      select case (htk_calc_method)
      case (-1:4)
        self%c_way = htk_calc_method
      case default
        error stop "Floquet: Error #1: htk_calc_method must be an integer in the [-1, 4] range."
      end select
    endif

    self%f_initialized = .true.

  end subroutine

  pure elemental integer function get_Nt(self)
    class(currents_calc_tsk), intent(in) :: self
    if (.not. self%f_initialized) error stop "Floquet: Error #2: Floquet task is not initialized."
    get_Nt = self%Nt
  end function get_Nt

  pure elemental integer function get_Ns(self)
    class(currents_calc_tsk), intent(in) :: self
    if (.not. self%f_initialized) error stop "Floquet: Error #2: Floquet task is not initialized."
    get_Ns = self%Ns
  end function get_Ns

  pure elemental integer function get_c_way(self)
    class(currents_calc_tsk), intent(in) :: self
    if (.not. self%f_initialized) error stop "Floquet: Error #2: Floquet task is not initialized."
    get_c_way = self%c_way
  end function get_c_way

  pure elemental function get_delta_smr(self)
    class(currents_calc_tsk), intent(in) :: self
    real(dp) :: get_delta_smr
    if (.not. self%f_initialized) error stop "Floquet: Error #2: Floquet task is not initialized."
    get_delta_smr = self%delta_smr
  end function get_delta_smr

  pure elemental logical function f_initialized(self)
    class(currents_calc_tsk), intent(in) :: self
    f_initialized = self%f_initialized
  end function f_initialized

  function current(self, crys, k, other) result(crr)
    class(task_specifier), intent(in) :: self
    class(crystal), intent(in) :: crys
    real(dp), intent(in) :: k(3)
    class(*), optional, intent(in) :: other

    complex(dp) :: crr(self%idims%size(), self%cdims%size())

    integer :: r_mem, rp_mem, r_arr(self%cdims%rank()), &
               i_mem

    real(dp) :: omega, t0, lambda, dt, &
                tper, q(3), &
                quasi(crys%num_bands()), &
                eig(crys%num_bands()), &
                smr

    real(dp), allocatable :: amplitudes(:, :), &
                             phases(:, :)

    integer :: Nharm, iharm, it, is, ir, i, n, m, l, p, ilambda, &
               cdims_shp(self%cdims%rank())

    complex(dp) :: H_TK(crys%num_bands(), crys%num_bands()), &
                   expH_TK(crys%num_bands(), crys%num_bands()), &
                   hf(crys%num_bands(), crys%num_bands()), &
                   rotH(crys%num_bands(), crys%num_bands()), &
                   rotF(crys%num_bands(), crys%num_bands()), &
                   rho(crys%num_bands(), crys%num_bands()), &
                   AW(crys%num_bands(), crys%num_bands(), 3), &
                   NAD(crys%num_bands(), crys%num_bands(), 3), &
                   AH(crys%num_bands(), crys%num_bands(), 3), &
                   HW(crys%num_bands(), crys%num_bands()), &
                   HH(crys%num_bands(), crys%num_bands()), &
                   P1W(crys%num_bands(), crys%num_bands(), 3), &
                   P2W(crys%num_bands(), crys%num_bands(), 3, 3), &
                   P3W(crys%num_bands(), crys%num_bands(), 3, 3, 3), &
                   P4W(crys%num_bands(), crys%num_bands(), 3, 3, 3, 3), &
                   tmp

    complex(dp), allocatable :: tev(:, :, :), pt(:, :, :, :), qs(:, :, :)

    integer, allocatable :: dims(:), pperm(:, :)

    type(container) :: DHW, DDHW, DDDHW, DDDDHW, &
                       DAW, DDAW, DDDAW

    cdims_shp = self%cdims%shape()

    crr = cmplx_0

    HW = crys%hamiltonian(kpt=k)
    AW = crys%berry_connection(kpt=k)

    DHW = crys%hamiltonian(kpt=k, derivative=1)
    P1W = wannier_momentum(HW=HW, &
                           DHW=DHW, &
                           AW=AW)

    call diagonalize(matrix=HW, P=rotH, D=HH, eig=eig)

    rho = cmplx_0
    !Get occupation matrix in the H gauge.
    do i = 1, crys%num_bands()
      if (eig(i) <= crys%fermi_energy()) rho(i, i) = cmplx_1
    enddo

    select type (self)
    type is (currents_calc_tsk)
      if (.not. self%f_initialized) error stop "Floquet: Error #2: Floquet task is not initialized."
      !Gather required quantities depending on the method to calculate H(k, t).
      select case (self%c_way)
      case (-1, 0)
        !-1: Length gauge: no-intraband appr.
        ! 0: Velocity gauge: no k-curvature appr.
        NAD = NADHW(eig=eig, &
                    rot=rotH, &
                    DHW=crys%hamiltonian(kpt=k, derivative=1), &
                    deg_thr=0.01_dp, deg_offset=0.04_dp)
        do i = 1, 3
          AH(:, :, i) = matmul(matmul(transpose(conjg(rotH)), AW(:, :, i)), rotH) + &
                        cmplx_i*NAD(:, :, i)
        enddo
      case (1)  !Velocity gauge: 2 terms in expansion.
        !All quantities already calculated.
        !Rotate rho back to W gauge if c_way >= 1.
        rho = matmul(rotH, matmul(rho, transpose(conjg(rotH))))
      case (2)  !Velocity gauge: 3 terms in expansion.
        !Rotate rho back to W gauge if c_way >= 1.
        rho = matmul(rotH, matmul(rho, transpose(conjg(rotH))))
        DDHW = crys%hamiltonian(kpt=k, derivative=2)
        DAW = crys%berry_connection(kpt=k, derivative=1)
        P2W = wannier_2nd_momentum(HW=HW, &
                                   DHW=DHW, &
                                   DDHW=DDHW, &
                                   AW=AW, &
                                   DAW=DAW, &
                                   P1W=P1W)
      case (3)  !Velocity gauge: 4 terms in expansion.
        !Rotate rho back to W gauge if c_way >= 1.
        rho = matmul(rotH, matmul(rho, transpose(conjg(rotH))))
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
        !Rotate rho back to W gauge if c_way >= 1.
        rho = matmul(rotH, matmul(rho, transpose(conjg(rotH))))
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

      Nharm = (self%cdims%rank() - 3)/6
      allocate (amplitudes(Nharm, 3), phases(Nharm, 3))
      allocate (tev(crys%num_bands(), crys%num_bands(), self%Nt))
      allocate (pt(self%Nt, -self%Ns:self%Ns, crys%num_bands(), crys%num_bands()))
      allocate (qs(-self%Ns:self%Ns, crys%num_bands(), crys%num_bands()))

      !Here we create a partial permutation dictionary,
      !iterating for all variables except lambda.
      allocate (dims(self%cdims%rank() - 1))
      do i = 1, size(dims)
        dims(i) = i
      enddo
      pperm = self%cdims%partial_permutation(dims)

      !Loop for each external variable except lambda.
      do rp_mem = 1, size(pperm(:, 1))

        r_arr = pperm(rp_mem, :)

        omega = self%cdt(var=self%cdims%rank() - 1, step=r_arr(self%cdims%rank() - 1))

        smr = self%delta_smr*omega !Adaptive delta smearing, proportional to omega.

        t0 = self%cdt(var=self%cdims%rank() - 2, step=r_arr(self%cdims%rank() - 2))

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
          tev(i, i, 1) = cmplx(1.0_dp, 0.0_dp, dp)
        enddo

        do it = 2, self%Nt

          tper = t0 + dt*real(it - 1, dp) !In eV^-1.

          !Calculate H(k, t) depending on the chosen method.
          select case (self%c_way)

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
          tev(:, :, it) = matmul(tev(:, :, it - 1), expH_TK)

        enddo

        !get the effective Floquet Hamiltonian.
        hf = cmplx_i*omega*logu(tev(:, :, self%Nt))/(2*pi)

        !Get P(t) = U(t)e^{-i*H_F*(tper - t0)}, and multiply by
        !e^{-i*s*w*(tper - t0)}, so we also get the integrand of Q_s.
        do is = -self%Ns, self%Ns
          do it = 1, self%Nt
            tper = t0 + dt*real(it - 1, dp) !In eV^-1.
            pt(it, is, :, :) = matmul(tev(:, :, it), expsh(cmplx_i*(tper - t0)*hf))
            pt(it, is, :, :) = pt(it, is, :, :)*exp(-cmplx_i*real(is, dp)*omega*(tper - t0))
          enddo
        enddo

        !Get Q_s.
        do m = 1, crys%num_bands()
          do n = 1, crys%num_bands()
            do is = -self%Ns, self%Ns
              qs(is, n, m) = extrapolation(pt(:, is, n, m))
            enddo
          enddo
        enddo

        !Right now the effective Floquet Hamiltonian, rho and Q_s are in:
        ! - Hamiltonian basis if c_way = -1, 0.
        ! - Wannier basis if c_way = 1, 2, 3, 4.
        !And, the momentum is in the Wannier basis.
        !If c_way <= 0 then, we rotate the momentum to the Hamiltonian basis.
        if (self%c_way <= 0) then
          do i = 1, 3
            P1W(:, :, i) = matmul(matmul(transpose(conjg(rotH)), P1W(:, :, i)), rotH)
          enddo
        endif

        call diagonalize(matrix=hf, P=rotF, eig=quasi)
        !Now, rotF holds the unitary rotation from the Wannier or Hamiltonian basis to the
        !Floquet basis.

        !Rotate Q_s, rho and the momentum to the Floquet basis
        do is = -self%Ns, self%Ns
          qs(is, :, :) = matmul(matmul(transpose(conjg(rotF)), qs(is, :, :)), rotF)
        enddo
        rho = matmul(matmul(transpose(conjg(rotF)), rho), rotF)
        do i = 1, 3
          P1W(:, :, i) = matmul(matmul(transpose(conjg(rotF)), P1W(:, :, i)), rotF)
        enddo

        !Compute the Fourier transform of the current exactly.
        do i = 1, 3
          do ilambda = 1, cdims_shp(self%cdims%rank())
            !Iterate now over lambda.
            r_arr(self%cdims%rank()) = ilambda
            !Get memory layout index and lambda.
            r_mem = self%cdims%ind(r_arr)
            lambda = self%cdt(var=self%cdims%rank(), step=r_arr(self%cdims%rank()))
            !Get current.
            tmp = cmplx_0
!$OMP             SIMD COLLAPSE(6) REDUCTION (+: tmp)
            do n = 1, crys%num_bands()
            do m = 1, crys%num_bands()
            do l = 1, crys%num_bands()
            do p = 1, crys%num_bands()
            do ir = -self%Ns, self%Ns
            do is = -self%Ns, self%Ns
              tmp = tmp + rho(n, m)*P1W(l, p, i)*conjg(qs(is, l, m))*qs(ir, p, n)* &
                    dirac_delta(quasi(n) - quasi(m) + real(ir - is, dp)*omega - lambda, smr)
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
!$OMP             END SIMD
            crr(i, r_mem) = tmp
          enddo
        enddo

      enddo

      deallocate (tev, pt, qs, amplitudes, phases, dims, pperm)

    end select

  end function current

end module Currents
