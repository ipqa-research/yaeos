program phase_diagram
   use forsus, only: Substance, forsus_dir
   use yaeos
   use yaeos__phase_equilibria_stability, only: tpd, min_tpd
   use yaeos, only: flash
   implicit none

   integer, parameter :: nc=3

   class(ArModel), allocatable :: model
   type(EquilibriaState) :: bub, fr
   type(PTEnvel2) :: env
   type(Substance) :: sus(nc)
   real(pr) :: tc(nc), pc(nc), ac(nc), kij(nc, nc), lij(nc, nc), T, P
   real(pr) :: z(nc), w(nc), mintpd, lnphiw(nc), d(nc), wold(nc)
   integer :: i, j, k

   forsus_dir = "build/dependencies/forsus/data/json"
   sus(1) = Substance("methane")
   sus(2) = Substance("hydrogen sulfide")
   sus(3) = Substance("ethane")

   kij = 0
   lij = 0

   tc = sus%critical%critical_temperature%value
   pc = sus%critical%critical_pressure%value/1e5_pr
   ac = sus%critical%acentric_factor%value

   model = SoaveRedlichKwong(tc, pc, ac, kij, lij)

   z = [1., 0.1, 1.]
   P = 15.6_pr
   T = 200._pr
   z = z/sum(z)

   ! do j=1,3
   !    w = 0.01
   !    w(j) = 0.98
   !    do i=1,50 
   !       wold = w
   !       mintpd = tpd(model, z, w, p, t, lnphiw=lnphiw, outd=d)
   !       w = exp(d - lnphiw)
   !       if (maxval(abs(w - wold)) < 1e-5) exit
   !    end do
   !    print *, j, mintpd, w
   ! end do

   call min_tpd(model, z, P, T, mintpd, w)

   print *, z, P, T
   print *, mintpd, w
   ! print *, w
   ! print *, w/sum(w)
end program phase_diagram
