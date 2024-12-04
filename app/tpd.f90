program phase_diagram
   use forsus, only: Substance, forsus_dir
   use yaeos
   implicit none

   integer, parameter :: nc=12


   class(ArModel), allocatable :: model
   type(Substance) :: sus(nc)
   real(pr) :: tc(nc), pc(nc), ac(nc)
   real(pr) :: z(nc), T, P
   real(pr) :: w(nc), mintpd, allmin(nc, nc), di(nc), dw(nc)
   real(pr) :: val, lnphi_w(nc)
   type(PTEnvel2) :: env
   type(EquilibriumState) :: sat
   integer :: i, j

   P = 3
   T = 140

   model = get_model()

   call min_tpd(model, z, P, T, mintpd, w, allmin)
   ! print *, mintpd, w/sum(w)

   do i=1,nc
      print *, allmin(i,:)/sum(allmin(i,:))
   end do

   print *, "SS"

   call model%lnphi_pt(z, P, T, root_type="stable", lnPhi=di)
   di = log(z) + di

   sat = saturation_temperature(model, z, P=0.001_pr, kind="dew", T0=500._pr)
   env = pt_envelope_2ph(model, z, sat)
   write(1, *) env
   env = find_hpl(model, z, T0=300._pr, p0=maxval(env%points%P))
   write(1, *) env

   ! w =[0.001, 0.999]

contains

   type(CubicEoS) function get_model()
      real(pr) :: tc(nc), pc(nc), w(nc), kij(nc,nc)
      z=[0.0656,0.3711,0.0538,0.0373,0.0261,0.0187,0.0218,0.1791,0.091,0.0605,0.0447,0.0302]

      tc=[304.088888888889,190.6,305.4,369.8,425.2,469.6,507.4,616.2,698.9,770.4,853.1,1001.2]
      pc=[73.7343491450634,45.9196083838941,48.7516547159404,42.3795504688362, &
          37.9291919470491,33.6811224489796,29.6353419746277,28.8261858797573,&
          19.3186017650303,16.5876999448428,15.2728212906784,14.6659542195256]
      w= [0.228,0.008,0.098,0.152,0.193,0.251,0.296,0.454,0.787,1.048,1.276,1.299]
      
      kij = 0
      kij(1, 2) = 0.12
      kij(1, 3:) = 0.15
      kij(:, 1) = kij(1, :)

      get_model = PengRobinson78(tc, pc, w, kij=kij)
   end function get_model


end program phase_diagram
