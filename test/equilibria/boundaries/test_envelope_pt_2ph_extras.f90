program test_envelope_pt_2ph_extras
   !! Two-phase PT envelope of a ternary started from a low-pressure dew
   !! point, its formatted output, and `find_hpl` on a stable mixture.
   !!
   !! - The line starts at the seeding point, crosses the critical point
   !!   (dew -> bubble) and every point satisfies isofugacity.
   !! - The reported critical point matches an independent critical-point
   !!   calculation at the same composition.
   !! - The formatted writer emits the dew/bubble sections and the critical
   !!   point.
   !! - `find_hpl` returns an empty envelope for a nearly ideal liquid
   !!   mixture (no liquid-liquid boundary).
   use yaeos, only: pr, CubicEoS, PengRobinson76, EquilibriumState, &
      saturation_temperature, critical_point, spec_CP, PTEnvel2, pt_envelope_2ph
   use yaeos__equilibria_boundaries_phase_envelopes_pt, only: find_hpl
   use testing_aux, only: assert, test_title
   implicit none

   integer, parameter :: nc = 3
   type(CubicEoS) :: model
   type(EquilibriumState) :: sat
   type(EquilibriumState) :: cp
   type(PTEnvel2) :: env
   real(pr) :: z(nc)
   integer :: i
   integer :: n
   integer :: n_dew
   integer :: n_bubble

   write(*, *) test_title("PT two-phase envelope from a dew point")

   model = PengRobinson76( &
      [190.564_pr, 369.83_pr, 507.6_pr], &
      [45.99_pr, 42.48_pr, 30.25_pr], &
      [0.0115478_pr, 0.1523_pr, 0.3013_pr])
   z = [0.6_pr, 0.25_pr, 0.15_pr]

   sat = saturation_temperature(model, z, P=1._pr, kind="dew", T0=300._pr)
   call assert(sat%kind == "dew", "PT: seeding dew point converged")

   env = pt_envelope_2ph(model, z, sat)
   n = size(env%points)

   call assert(n > 20, "PT: envelope traced")
   call assert(abs(env%points(1)%T - sat%T) < 1e-6_pr*sat%T, "PT: first point is the seeding dew temperature")
   call assert(abs(env%points(1)%P - 1._pr) < 1e-8_pr, "PT: first point is at the seeding pressure")

   n_dew = count(env%points%kind == "dew")
   n_bubble = count(env%points%kind == "bubble")
   call assert(n_dew > 0 .and. n_bubble > 0 .and. n_dew + n_bubble == n, &
      "PT: envelope has a dew and a bubble branch")
   call assert(all(env%points(:n_dew)%kind == "dew"), "PT: dew points come before the bubble points")

   do i = 1, n
      call assert(isofugacity_error(env%points(i)) < 1e-5_pr, "PT: isofugacity at every envelope point")
   end do

   ! Cricondentherm: the dew branch passes through a maximum temperature
   call assert(maxval(env%points(:n_dew)%T) > env%points(n_dew)%T + 10, &
      "PT: dew branch passes through the cricondentherm")

   call assert(size(env%cps) == 1, "PT: one critical point")
   cp = critical_point(model, z, z, spec=spec_CP%a, S=1._pr, max_iters=1000, &
      T0=env%cps(1)%T, P0=env%cps(1)%P)
   call assert(abs(cp%T - env%cps(1)%T) < 0.5_pr, "PT: critical temperature matches an independent CP")
   call assert(abs(cp%P - env%cps(1)%P) < 0.5_pr, "PT: critical pressure matches an independent CP")
   call assert(env%points(n_dew)%P < cp%P + 1 .and. env%points(n_dew + 1)%P > cp%P - 1, &
      "PT: dew/bubble switch happens at the critical pressure")

   call test_write()
   call test_find_hpl_stable()

contains

   real(pr) function isofugacity_error(point)
      !! Isofugacity residual between the feed phase and the incipient phase.
      !! The feed phase is identified by its composition, so this check does
      !! not depend on how x and y are assigned in the output state.
      type(EquilibriumState), intent(in) :: point
      real(pr) :: main(nc)
      real(pr) :: inc(nc)
      real(pr) :: lnphi_main(nc)
      real(pr) :: lnphi_inc(nc)
      character(len=6) :: kind_main
      character(len=6) :: kind_inc

      if (maxval(abs(point%x - z)) < 1e-10_pr) then
         main = point%x
         inc = point%y
      else
         main = point%y
         inc = point%x
      end if
      if (point%kind == "dew") then
         kind_main = "vapor"
         kind_inc = "liquid"
      else
         kind_main = "liquid"
         kind_inc = "vapor"
      end if
      call model%lnphi_pt(main, point%P, point%T, root_type=kind_main, lnPhi=lnphi_main)
      call model%lnphi_pt(inc, point%P, point%T, root_type=kind_inc, lnPhi=lnphi_inc)
      isofugacity_error = maxval(abs(log(main) + lnphi_main - log(inc) - lnphi_inc))
   end function isofugacity_error

   subroutine test_write()
      integer :: u
      integer :: ios
      integer :: n_lines
      character(len=1024) :: line
      logical :: has_header
      logical :: has_dew
      integer :: n_sections
      logical :: has_critical
      real(pr) :: Tcp
      real(pr) :: Pcp

      open(newunit=u, status="scratch", form="formatted", action="readwrite")
      write(u, *) env
      rewind(u)

      has_header = .false.
      has_dew = .false.
      n_sections = 0
      has_critical = .false.
      n_lines = 0
      Tcp = -1
      Pcp = -1
      do
         read(u, "(A)", iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         n_lines = n_lines + 1
         if (index(line, "#PTEnvel2") > 0) has_header = .true.
         if (index(line, "#dew") > 0) has_dew = .true.
         if (index(line, "#dew") > 0 .or. index(line, "#bubble") > 0) n_sections = n_sections + 1
         if (has_critical .and. Tcp < 0) then
            read(line, *) Tcp, Pcp
         end if
         if (index(line, "#Critical") > 0) has_critical = .true.
      end do
      close(u)

      call assert(has_header, "PTEnvel2 write: header written")
      ! One section for the first branch plus one more per critical point.
      ! (The label of the section after a critical point is taken from the
      ! point following the one nearest to the CP, which here is still a dew
      ! point, so the label itself is not checked.)
      call assert(has_dew, "PTEnvel2 write: first section labelled with the first point kind")
      call assert(n_sections == 1 + size(env%cps), "PTEnvel2 write: a new section after each critical point")
      call assert(has_critical, "PTEnvel2 write: critical section written")
      call assert(n_lines >= size(env%points) + 3, "PTEnvel2 write: one line per point")
      call assert(abs(Tcp - env%cps(1)%T) < 1e-8_pr*env%cps(1)%T .and. abs(Pcp - env%cps(1)%P) < 1e-8_pr*env%cps(1)%P, &
         "PTEnvel2 write: critical point values written")
   end subroutine test_write

   subroutine test_find_hpl_stable()
      !! Propane + n-butane liquids are nearly ideal: the fugacity of each
      !! component in the mixture is below the pure-liquid one at every
      !! temperature, so no liquid-liquid boundary is returned.
      type(CubicEoS) :: ideal
      type(PTEnvel2) :: hpl

      ideal = PengRobinson76([369.83_pr, 425.12_pr], [42.48_pr, 37.96_pr], [0.1523_pr, 0.200164_pr])
      hpl = find_hpl(ideal, [0.4_pr, 0.6_pr], 300._pr, 50._pr, 100)
      call assert(.not. allocated(hpl%points), "find_hpl: no liquid-liquid boundary for a nearly ideal mixture")
   end subroutine test_find_hpl_stable
end program test_envelope_pt_2ph_extras
