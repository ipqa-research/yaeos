program phase_diagram
   use yaeos

   implicit none
   type(CubicEoS) :: model ! Thermodynamic model to be used
   type(HV_NRTL) :: mr
   type(EquilibriumState) :: sat
   integer, parameter :: nc = 3
   real(pr) :: Tc(nc), Pc(nc), w(nc)    ! Component's critical constants

   real(pr) :: a, z(nc), zi(nc), z0(nc)

   character(len=*), parameter :: path = "inputs/pt_diagram.nml"

   call read_composition(path, z0, zi)
   call read_model(path)

   ! call calculate_px_envelope(390._pr)
   ! call calculate_px_envelope(323.15_pr)

   a = 0.805
   z = zi*a + z0*(1 - a)
   sat = saturation_temperature(model, z, P=1e-1_pr, kind="dew", t0=250._pr)
   print *, sat
   ! call calculate_dew_envelope

contains
   subroutine read_model(path)
      character(len=*), intent(in) :: path
      call critical_constants(nc, path, Tc, Pc, w)
      model = PengRobinson76(Tc, Pc, w)
      mr = mixrule_hvrtl(nc, path)
      call model%set_mixrule(mr)
   end subroutine read_model

   subroutine read_composition(path, z0, zi)
      character(len=*), intent(in) :: path
      real(pr), intent(out) :: z0(:)
      real(pr), intent(out) :: zi(:)
      namelist /compositions/ z0, zi

      integer :: unit

      open(newunit=unit, file=path)
      read(unit, nml=compositions)
      close(unit)
   end subroutine read_composition

   type(HV_NRTL) function mixrule_hvrtl(nc, filepath) result(mr)
      integer, intent(in) :: nc
      character(len=*), intent(in) :: filepath
      real(pr) :: gji0(nc, nc), gjiT(nc, nc), alpha(nc, nc), kij(nc, nc)
      logical :: use_kij(nc, nc)
      integer :: unit
      integer :: i

      namelist /parameters_hvnrtl/ gji0, gjiT, alpha, kij, use_kij
      open(newunit=unit, file=filepath)
      read(unit, nml=parameters_hvnrtl)
      close(unit)
      mr = init_hvnrtl(b=model%b, del1=model%del1, gji0=gji0, gjiT=gjiT, alpha=alpha, kij=kij, use_kij=use_kij)
   end function mixrule_hvrtl

   subroutine critical_constants(nc, filepath, Tc, Pc, w)
      integer, intent(in) :: nc
      character(len=*), intent(in) :: filepath
      real(pr) :: Pc(nc), Tc(nc), w(nc)
      integer :: unit

      namelist /parameters_critical/ tc, pc, w
      open(newunit=unit, file=filepath, status='old', action='read')
      read(unit, nml=parameters_critical)
      close(unit)
   end subroutine critical_constants

   subroutine calculate_dew_envelope
      type(PTEnvelMP) :: envelope
      integer, parameter :: np=1

      real(pr) :: x_l0(np, nc), w0(nc), betas0(np), P0, T0
      character(len=14) :: kinds_x(np), kind_w
      integer :: ns0
      real(pr) :: ds0
      real(pr) :: beta_w=0

      integer :: i

      kinds_x = "vapor"
      kind_w = "liquid"

      x_l0(1, :) = sat%y
      w0 = sat%x
      P0 = sat%P
      T0 = sat%T
      betas0 = 1

      ns0 = nc*np+np+1
      ds0 = 0.1_pr

      envelope = pt_envelope(&
         model, z, np, kinds_x, kind_w, &
         x_l0, w0, betas0, p0, t0, &
         ns0, ds0, beta_w)

      do i=1,size(envelope%points)
         write(*, *) envelope%points(i)%T, envelope%points(i)%P, envelope%points(i)%iters
      end do

   end subroutine calculate_dew_envelope

   subroutine calculate_px_envelope(T)
      real(pr), intent(in) :: T
      type(PXEnvelMP) :: envelope
      type(EquilibriumState) :: sat
      integer, parameter :: np=1
      real(pr) :: x_l0(np, nc), w0(nc), betas0(np), a0, P0
      character(len=14) :: kinds_x(np), kind_w
      integer :: ns0
      real(pr) :: ds0
      real(pr) :: beta_w=0
      integer :: i

      type(PurePsat) :: pure_psat
      kinds_x = "liquid"
      kind_w = "vapor"


      ! ==============================================================
      ! Calculation of the saturation line of the pure component 2 to
      ! get the initial guess of pressure p0
      ! --------------------------------------------------------------
      pure_psat = pure_saturation_line(model, 2, 1e-10_pr, T*0.9_pr)
      p0 = pure_psat%get_P(T)

      a0 = 1e-5_pr
      z = a0 * zi + (1 - a0) * z0

      ! Calculate a saturation point at T to get intial compositions
      ! and pressure of the line.
      sat = saturation_pressure(model, z, T=T, kind="bubble", p0=p0)


      x_l0(1, :) = z
      w0 = sat%y/sum(sat%y)
      betas0 = [1._pr]
      ns0 = nc*np + np + 2
      ds0 = 1e-5_pr

      p0 = sat%P

      print *, sum(z)
      print *, sum(w0)
      print *, sum(x_l0(1, :))


      envelope = px_envelope(&
         model=model, z0=z0, zi=zi, np=np, T=T, &
         kinds_x=kinds_x, kind_w=kind_w, &
         x_l0=x_l0, w0=w0, betas0=betas0, p0=p0, alpha0=a0, &
         ns0=ns0, ds0=ds0, beta_w=beta_w &
         )
      write(1, *) "ENVEL"
      do i=1, size(envelope%points)
         write(1, *) envelope%alpha(i), envelope%points(i)%P, envelope%points(i)%kind_w
         write(2, *) envelope%alpha(i), log(envelope%points(i)%w/envelope%points(i)%x_l(1, :))
      end do
      write(1, *) "ENVEL"

   end subroutine calculate_px_envelope

end program phase_diagram
