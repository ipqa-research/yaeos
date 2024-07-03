program phase_diagram
   use forsus, only: Substance, forsus_dir
   use yaeos
   use yaeos__phase_equilibria_stability, only: tm, min_tpd
   use yaeos, only: flash
   implicit none

   integer, parameter :: nc=2

   class(ArModel), allocatable :: model
   type(Substance) :: sus(nc)
   real(pr) :: tc(nc), pc(nc), ac(nc)
   real(pr) :: z(nc), T, P
   real(pr) :: w(nc), mintpd, mins(nc, nc)

   forsus_dir = "build/dependencies/forsus/data/json"
   sus(1) = Substance("methane")
   sus(2) = Substance("hydrogen sulfide")

   z = [0.13, 1-0.13]
   z = z/sum(z)
   
   P = 20.0_pr
   T = 190._pr

   tc = sus%critical%critical_temperature%value
   pc = sus%critical%critical_pressure%value/1e5_pr
   ac = sus%critical%acentric_factor%value
   
   model = SoaveRedlichKwong(tc, pc, ac)

   call min_tpd(model, z, P, T, mintpd, w)
   print *, mintpd, w/sum(w)

   P = 15
   call min_tpd(model, z, P, T, mintpd, w)
   print *, mintpd, w/sum(w)

end program phase_diagram
