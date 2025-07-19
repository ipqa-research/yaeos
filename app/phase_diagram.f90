program phase_diagram
   !! Program for calculation of phase diagrams.
   use forsus, only: Substance, forsus_dir, forsus_default_dir
   use yaeos, only: pr, &
      SoaveRedlichKwong, PengRobinson76, PengRobinson78, RKPR, &
      EquilibriumState, ArModel, PTEnvel2, &
      pt_envelope_2ph, saturation_pressure, saturation_temperature, k_wilson
   implicit none

   ! ===========================================================================
   ! Variables definition
   ! ---------------------------------------------------------------------------
   integer, parameter :: nc=4
   class(ArModel), allocatable :: model ! Thermodynamic model to be used
   type(EquilibriumState) :: sat_point  ! Init
   type(PTEnvel2) :: envelope           ! PT Phase envelope
   real(pr) :: tc(nc), pc(nc), w(nc)    ! Component's critical constants
   real(pr) :: kij(nc, nc)              ! Binary interaction parameters
   real(pr) :: n(nc)                    ! Termodynamic variables
   type(Substance) :: sus(nc)           ! Substances to use
   ! ===========================================================================

   ! forsus database directory
   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir

   ! Find the selected substances on the database and extract their
   ! critical constants
   sus(1) = Substance("methane")
   sus(2) = Substance("n-hexane")
   sus(3) = Substance("carbon dioxide")
   sus(4) = Substance("n-tricosane")
   call get_critical_constants(sus, tc, pc, w)

   ! Model definition
   kij(3, :) = 0.08_pr
   kij(:, 3) = 0.08_pr

   model = PengRobinson76(tc, pc, w, kij)

   ! Composition vector
   n = [3, 1, 120, 1]
   n = n/sum(n)

   ! Calculate a dew point at low pressure to later
   ! initialize the phase envelope
   sat_point = saturation_temperature(model, n, P=0.01_pr, kind="dew", t0=350._pr)

   ! Calculate phase envelope
   envelope = pt_envelope_2ph(model, n, sat_point, delta_0=0.1_pr)

   ! Write the phase envelope to screen
   write(1, *) envelope

   sat_point = saturation_pressure(model, n, T=200._pr, kind="bubble", p0=10._pr)
   envelope = pt_envelope_2ph(model, n, sat_point)
   write(2, *) envelope

   find_hpl:block
      integer :: i
      real(pr) :: T, P
      real(pr) :: z(nc), y(nc)
      real(pr) :: lnphi_y(nc), lnphi_z(nc)
      type(EquilibriumState) :: fr
      real(pr) :: diffs(nc), Ts(nc)
      integer :: ncomp

      z = n/sum(n)

      P = 1000

      do ncomp=1,nc
         T = maxval(envelope%points%T)
         y = 0
         y(ncomp) = 1
         do i=500, 100, -10
            T = real(i, pr)
            call model%lnphi_pt(n, P, T, root_type="liquid", lnPhi=lnphi_z)
            call model%lnphi_pt(y, P, T, root_type="liquid", lnPhi=lnphi_y)
            diffs(ncomp) = log(z(ncomp)) + lnphi_z(ncomp) - log(y(ncomp)) - lnphi_y(ncomp)
            if (diffs(ncomp) > 0) exit
         end do
         Ts(ncomp) = T
      end do

      T = maxval(Ts, mask=diffs>0)
      ncomp = findloc(Ts, T, dim=1)

      y=0
      y(ncomp) = 1

      fr%x = z
      fr%y = y + 1e-5
      fr%y = fr%y/sum(fr%y)
      fr%T = T
      fr%P = P
      fr%kind = "liquid-liquid"
      envelope = pt_envelope_2ph(model, n, fr, specified_variable_0=nc+2, delta_0=-5.0_pr, iterations=1000)
      write(3, *) envelope
   end block find_hpl

contains

   subroutine get_critical_constants(subs, tc_in, pc_in, w_in)
      type(Substance) :: subs(:)
      real(pr), intent(out) :: tc_in(:), pc_in(:), w_in(:)

      tc_in = subs%critical%critical_temperature%value
      pc_in = subs%critical%critical_pressure%value/1e5
      w_in = subs%critical%acentric_factor%value
   end subroutine get_critical_constants
end program phase_diagram
