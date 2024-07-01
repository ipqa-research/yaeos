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
   real(pr) :: z(nc), w(nc), mintpd
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

   do i=1,200, 10
      w(1) = real(i, pr)/100
      do j=i,200, 10
         w(2) = real(200-j,pr)/100
         w(3) = 2 - w(1) - w(2)
         mintpd = tpd(model, z, w, p, t)
         write(4, *) w(1), w(2), mintpd
      end do
      write(4, *) ""
   end do

   ! write(3, *) env

   call min_tpd(model, z, P, T, mintpd, w)
   print *, z, P, T
   print *, mintpd
   print *, w
   print *, w/sum(w)
   
   fr = flash(model, z, t, p_spec=p, iters=i)
   print *, "FLASH", i
   write (*, *) fr%x
   write (*, *) fr%y

   bub = saturation_pressure(model, z, t-50, kind="bubble")
   env = pt_envelope_2ph(model, z, first_point=bub)

   write (1, *) env

end program phase_diagram
