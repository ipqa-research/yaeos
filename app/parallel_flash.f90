program flasher
   !! Example program to make flash calculations with a PengRobinson76
   !! EoS model.
   !! In this case the program is parallelized using OpenMP.
   !! To run it in parallel you must use:
   !!
   !! ```
   !! fpm run --profile release --flag "-fopenmp"
   !! ```
   !!

   ! Import the relevant parts of yaeos
   use yaeos, only: pr, EquilibriumState, flash, PengRobinson76, ArModel, PTEnvel2, pt_envelope_2ph
   use yaeos, only: saturation_pressure, saturation_temperature, k_wilson, p_wilson
   use forsus, only: Substance, forsus_dir, forsus_default_dir
   use fortime, only: timer
   implicit none

   ! Variables definition:
   class(ArModel), allocatable :: model !! Model to use
   type(EquilibriumState) :: flash_result !! Result of Flash calculation
   type(EquilibriumState) :: sat

   integer, parameter :: nc=7

   type(Substance) :: sus(nc)
   type(PTEnvel2) :: env

   real(pr) :: Tc(nc), Pc(nc), w(nc)
   real(pr) :: n(nc), T, P, k0(nc)
   integer :: iter

   real(pr) :: dP, dT
   integer, parameter :: nP=100, nT=100
   real(pr) :: T0=200, Tf=600
   real(pr) :: P0=1, Pf=130

   real(pr) :: ts(nt, np), ps(nt, np), betas(nt, np)

   type(Timer) :: tim

   integer :: i, j, envpoint

   forsus_dir = "build/dependencies/forsus/" // forsus_default_dir


   ! ==========================================================================
   ! Set up the model and the substances
   ! --------------------------------------------------------------------------
   n = [(i, i=1, nc)]
   n = n/sum(n)
   sus(1) = Substance("methane")
   sus(2) = Substance("carbon dioxide")
   sus(3) = Substance("n-decane")
   sus(4) = Substance("n-eicosane")
   sus(5) = Substance("water")
   sus(6) = Substance("n-dodecane")
   sus(7) = Substance("n-heptane")
   Tc = sus%critical%critical_temperature%value
   Pc = sus%critical%critical_pressure%value/1e5
   w = sus%critical%acentric_factor%value
   model = PengRobinson76(Tc, Pc, w)


   ! ==========================================================================
   ! Calculate the phase-envelope to determine the upper limits
   ! of pressure and temperature
   ! --------------------------------------------------------------------------
   sat = saturation_temperature(model, n, P=0.01_pr, kind="dew", T0=500._pr)
   print *, sat

   env = pt_envelope_2ph(model, n, sat)
   write(2, *) env

   ! ==========================================================================
   ! Perform the flash calculations
   ! --------------------------------------------------------------------------
   T0 = 200
   P0 = 10
   Tf = maxval(env%points%T)
   Pf = maxval(env%points%P)

   dp = (Pf - P0)/(np-1)
   dt = (Tf - T0)/(nt-1)
   call tim%timer_start()
   !$OMP PARALLEL DO PRIVATE(i, j, t, p, k0, flash_result) shared(model, ts, ps, betas)
   do i=1,np
      T = T0 + (i-1)*dt
      P = P0
      do j=1,nt
         P = P0 + (j-1)*dp
         flash_result = flash(model, n, T=T, P_spec=P, iters=iter)
         ts(i, j) = flash_result%T
         ps(i, j) = flash_result%P
         betas(i, j) = flash_result%beta
      end do
   end do
   !$OMP END PARALLEL DO
   call tim%timer_stop()

   ! ==========================================================================
   ! Print the results
   ! Since parallelization is used, the output may not be in order, so we
   ! print it after it is calculated.
   ! --------------------------------------------------------------------------
   do i=1,nt
      do j=1,np
         if (Ps(i, j) > 0) then
            write(3, "(3(F10.3, 1x))") ts(i, j), ps(i, j), betas(i, j)
         end if
      end do
      write(3, *) ""
   end do
end program flasher
