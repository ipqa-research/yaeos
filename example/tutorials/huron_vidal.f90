program main
   !! Example of using CubicEoS with Huron-Vidal mixing rules with an
   !! NRTL model as the Ge model

   use numerical_differentiation_module
   use forsus, only: Substance, forsus_dir
   use yaeos, only: pr, SoaveRedlichKwong, CubicEoS, NRTL
   use yaeos__models_cubic_mixing_rules_huron_vidal, only: HV

   implicit none
   integer, parameter :: nc = 2

   type(numdiff_type) :: numdiff

   real(pr) :: T = 150, n(nc)

   real(pr) :: a(nc, nc), b(nc, nc), c(nc, nc) ! NRTL parameters
   real(pr) :: tc(nc), pc(nc), w(nc) ! Cubic EoS parameters

   real(pr) :: alpha(nc), Tr(nc)
   type(NRTL) :: ge_model ! Excess Gibbs model that will be used
   type(CubicEoS) :: model ! Main model
   type(HV) :: mixrule

   real(pr), dimension(nc) :: ai, daidt, daidt2
   real(pr) ::  D, dDdT, dDdT2, dDi(nc), dDidT(nc), dDij(nc, nc)
   real(pr) ::dx=0.001, Tdx, vardx, vardx2, dfn(nc), dfn2(nc, nc)
   real(pr) :: Ge, Gen(nc), Gen2(nc,nc)
   type(Substance)  :: sus(nc)
   real(pr) :: aux, dT, dT2

   integer :: i, j

   forsus_dir = "./build/dependencies/forsus/data/json"

   sus(1) = Substance("water")
   sus(2) = Substance("ethanol")

   tc = sus%critical%critical_temperature%value
   w = sus%critical%acentric_factor%value
   pc = sus%critical%critical_pressure%value/1e5

   a=0; b=0; c=0


   a(1, 2) = 3.458
   a(2, 1) = -0.801

   b(1, 2) = -586.1
   b(2, 1) = 246.2

   c(1, 2) = 0.3
   c(2, 1) = 0.3

   ge_model = NRTL(a, b, c)

   n = [0.2, 0.8]
   T = 150
   Tr = T/Tc
   
   model = SoaveRedlichKwong(tc, pc, w)
   mixrule = HV(ge_model, model%b)
   mixrule%q = 0.53

   call model%alpha%alpha(Tr, ai, daidt, daidt2)
   ai = ai*model%ac
   daidt = daidt*model%ac/Tc
   daidt2 = daidt2*model%ac/Tc**2

   dx = 0.01_pr

   n = [0.2, 0.8]

   do i = 1, nc
      dfn(i) = df(i, dx=0.00001_pr)
      do j = 1, nc
         dfn2(i, j) = df2(i, j, dx=1e-2_pr)
      end do
   end do

   call mixrule%Dmix(n, T, ai, daidt, daidt2, D, dDdT, dDdT2, dDi, dDidT, dDij)
   print *, "D", D
  
   print *, "dDdT"
   aux = df(nc + 1, dx)
   print *, "anal", dDdT
   print *, "numd", aux
   
   print *, "dDdT2"
   print *,"anal", dDdT2
   aux = df2(nc + 1, nc + 1, dx=1._pr)
   print *,"numd", aux
   
   print *, "dDi"
   print *, "anal", dDi
   print *, "numd", dfn
   
   print *, "dDij"
   print *, "anal", dDij
   print *, "numd", dfn2

   j = 2
   do i=1,99
      n(j) = real(i, pr)/100
      call ge_model%excess_gibbs(n, t, ge=ge, Gen=Gen, Gen2=Gen2)
      call mixrule%Dmix(n, T, ai, daidt, daidt2, D, dDdT, dDdT2, dDi, dDidT, dDij)
      print *, "D", n(j), D, dDi, dDij
      print *, "Ge", n(j), Ge, Gen, Gen2
   end do
contains
   real(pr) function df(i, dx)
      integer :: i
      real(pr) :: ndx(nc), tdx, dx
      real(pr) :: f0, f1, f2

      real(pr) ::  D, dDdT, dDdT2, dDi(nc), dDidT(nc), dDij(nc, nc)
      if (i == nc + 1) then
         tdx = T + dx
         f0 = DfromT(T)
         f1 = DfromT(Tdx)
      else

         ndx = n
         call mixrule%Dmix(ndx, T, ai, daidt, daidt2, f0, dDdT, dDdT2, dDi, dDidT, dDij)

         ndx(i) = n(i) + dx
         call mixrule%Dmix(ndx, T, ai, daidt, daidt2, f1, dDdT, dDdT2, dDi, dDidT, dDij)
      end if

      df = (f1 - f0)/dx
   end function

   real(pr) function df2(i, j, dx)
      integer :: i, j
      real(pr) :: ndx(nc), dx

      real(pr) :: f0, f1, f2, f3, f4
      real(pr) :: Tdx
      real(pr) ::  D, dDdT, dDdT2, dDi(nc), dDidT(nc), dDij(nc, nc)

      df2 = 0

      call mixrule%Dmix(ndx, T, ai, daidt, daidt2, f0, dDdT, dDdT2, dDi, dDidT, dDij)

      if (i == j .and. i == nc + 1) then
         Tdx = T + dx
         f1 = DfromT(Tdx)
         Tdx = T - dx
         f2 = DfromT(Tdx)
         df2 = (f1 - 2*f0 + f2)/(dx**2)
      
      else if (i == j) then
         ndx(i) = n(i) + dx
         call mixrule%Dmix(ndx, T, ai, daidt, daidt2, f1, dDdT, dDdT2, dDi, dDidT, dDij)
         ndx(i) = n(i) - dx
         call mixrule%Dmix(ndx, T, ai, daidt, daidt2, f2, dDdT, dDdT2, dDi, dDidT, dDij)
         df2 = (f1 - 2*f0 + f2)/(dx**2)
      
      else
         ndx(i) = n(i) + dx
         ndx(j) = n(j) + dx
         call mixrule%Dmix(ndx, T, ai, daidt, daidt2, f1, dDdT, dDdT2, dDi, dDidT, dDij)

         ndx(i) = n(i) + dx
         ndx(j) = n(j) - dx
         call mixrule%Dmix(ndx, T, ai, daidt, daidt2, f2, dDdT, dDdT2, dDi, dDidT, dDij)

         ndx(i) = n(i) - dx
         ndx(j) = n(j) + dx
         call mixrule%Dmix(ndx, T, ai, daidt, daidt2, f3, dDdT, dDdT2, dDi, dDidT, dDij)

         ndx(i) = n(i) - dx
         ndx(j) = n(j) - dx
         call mixrule%Dmix(ndx, T, ai, daidt, daidt2, f4, dDdT, dDdT2, dDi, dDidT, dDij)
      end if

      df2 = (f1 - f2 - f3 + f4)/(4*dx**2)
   end function

   real(pr) function DfromT(T, dT, dT2)
      real(pr) :: T
      real(pr), optional :: dT, dT2
      real(pr) :: ai(nc), daidt(nc), daidt2(nc), Tr(nc), Tdx

      real(pr) ::  D, dDdT, dDdT2, dDi(nc), dDidT(nc), dDij(nc, nc)
      Tr = T/Tc
      call model%alpha%alpha(Tr, ai, daidt, daidt2)
      ai = ai*model%ac
      daidt = daidt*model%ac/Tc
      daidt2 = daidt2*model%ac/Tc**2
      call mixrule%Dmix(n, T, ai, daidt, daidt2, DfromT, dDdT, dDdT2, dDi, dDidT, dDij)

      if (present(dT)) dT = dDdT
      if (present(dT2)) dT2 = dDdT2
   end function
end program