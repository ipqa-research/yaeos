module fixtures_models
   use yaeos, only: pr, R, CubicEoS, NRTL
   use yaeos__tapenade_ar_api, only: ArModelTapenade

contains

   type(CubicEoS) function binary_PR76() result(eos)
      use yaeos, only: PengRobinson76
      integer, parameter :: n = 2
      real(pr) :: tc(n), pc(n), w(n)
      real(pr) :: kij(n, n), lij(n, n)
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]

      kij = reshape([0., 0.1, 0.1, 0.], [n, n])
      lij = kij/2
      eos = PengRobinson76(tc, pc, w, kij, lij)
   end function binary_PR76

   type(CubicEoS) function binary_PR78() result(eos)
      use yaeos, only: PengRobinson78
      integer, parameter :: n = 2
      real(pr) :: tc(n), pc(n), w(n)
      real(pr) :: kij(n, n), lij(n, n)
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]

      kij = reshape([0., 0.1, 0.1, 0.], [n, n])
      lij = kij/2

      eos = PengRobinson78(tc, pc, w, kij, lij)
   end function binary_PR78

   type(CubicEoS) function binary_SRK() result(eos)
      use yaeos, only: SoaveRedlichKwong
      integer, parameter :: n = 2
      real(pr) :: tc(n), pc(n), w(n)
      real(pr) :: kij(n, n), lij(n, n)
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]

      kij = reshape([0., 0.1, 0.1, 0.], [n, n])
      lij = kij/2

      eos = SoaveRedlichKwong(tc, pc, w, kij, lij)
   end function binary_SRK

   type(CubicEoS) function binary_RKPR() result(eos)
      use yaeos, only: RKPR
      integer, parameter :: n = 2
      real(pr) :: tc(n), pc(n), w(n), zc(n)
      real(pr) :: kij(n, n), lij(n, n)
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]
      zc = [0.23, 0.26]

      kij = reshape([0., 0.1, 0.1, 0.], [n, n])
      lij = kij/2

      eos = RKPR(tc, pc, w, zc, kij, lij)
   end function binary_RKPR

   type(hdPR76) function binary_PR76_hd() result(eos)
      use autodiff_hyperdual_pr76, only: setup, hdPR76
      integer, parameter :: n = 2
      real(pr) :: tc(n), pc(n), w(n)
      real(pr) :: kij(n, n), lij(n, n)
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]

      kij = reshape([0., 0.1, 0.1, 0.], [n, n])
      lij = kij/2

      eos = setup(tc, pc, w, kij, lij)
   end function binary_PR76_hd

   type(TPR76) function binary_PR76_tape() result(eos)
      use autodiff_tapenade_pr76, only: setup_model, TPR76
      integer, parameter :: n = 2
      real(pr) :: tc(n), pc(n), w(n)
      real(pr) :: kij(n, n), lij(n, n)
      tc = [190._pr, 310._pr]
      pc = [14._pr, 30._pr]
      w = [0.001_pr, 0.03_pr]

      kij = reshape([0., 0.1, 0.1, 0.], [n, n])
      lij = kij/2

      eos = setup_model(tc, pc, w, kij, lij)
   end function binary_PR76_tape

   type(NRTL) function binary_NRTL_tape() result(model)
      integer, parameter :: n = 2
      real(pr) :: a(n, n), b(n, n), c(n, n)

      real(pr) :: ge, GeT, GeT2, Gen(n), Gen2(n, n), GeTn(n)
      real(pr) :: lngamma(n)

      a = 0; b = 0; c = 0

      a(1, 2) = 3.458
      a(2, 1) = -0.801

      b(1, 2) = -586.1
      b(2, 1) = 246.2

      c(1, 2) = 0.3
      c(2, 1) = 0.3

      model = NRTL(a, b, c)
   end function binary_NRTL_tape

   type(CubicEoS) function binary_NRTL_SRK() result(model)
      use yaeos, only: NRTL, SoaveRedlichKwong, pr, CubicEoS, MHV

      integer, parameter :: nc = 2
      real(pr), dimension(nc,nc) :: a, b, c
      real(pr), dimension(nc) :: tc, pc, w
      type(MHV) :: mixrule
      type(NRTL) :: ge_model

      tc = [647.13999999999999, 513.91999999999996]
      pc = [220.63999999999999, 61.479999999999997]
      w =  [0.34399999999999997, 0.64900000000000002]
      a = 0; b = 0; c = 0

      ! NRTL model parameters
      a(1, 2) = 3.458
      a(2, 1) = -0.801

      b(1, 2) = -586.1
      b(2, 1) = 246.2

      c(1, 2) = 0.3
      c(2, 1) = 0.3

      ge_model = NRTL(a, b, c)

      ! Define the model to be SRK
      model = SoaveRedlichKwong(tc, pc, w)
      mixrule = MHV(ge=ge_model, q=0.593_pr, b=model%b)
      deallocate (model%mixrule)
      model%mixrule = mixrule
   end function binary_NRTL_SRK

   type(CubicEoS) function binary_NRTL_SRK_HV() result(model)
      use yaeos, only: NRTL, SoaveRedlichKwong, pr, CubicEoS, HV

      integer, parameter :: nc = 2
      real(pr), dimension(nc,nc) :: a, b, c
      real(pr), dimension(nc) :: tc, pc, w
      type(HV) :: mixrule
      type(NRTL) :: ge_model

      tc = [647.13999999999999, 513.91999999999996]
      pc = [220.63999999999999, 61.479999999999997]
      w =  [0.34399999999999997, 0.64900000000000002]
      a = 0; b = 0; c = 0

      ! NRTL model parameters
      a(1, 2) = 3.458
      a(2, 1) = -0.801

      b(1, 2) = -586.1
      b(2, 1) = 246.2

      c(1, 2) = 0.3
      c(2, 1) = 0.3

      ge_model = NRTL(a, b, c)

      ! Define the model to be SRK
      model = SoaveRedlichKwong(tc, pc, w)
      mixrule%ge = ge_model
      mixrule%bi = model%b
      mixrule%del1 = model%del1
      call model%set_mixrule(mixrule)
   end function binary_NRTL_SRK_HV

   type(CubicEoS) function multicomponent_PR(z0, zi) result(model)
      use yaeos, only: PengRobinson78
      integer, parameter :: nc = 12
      real(pr), intent(out) :: z0(nc), zi(nc)

      real(pr) :: tc(nc), pc(nc), w(nc), kij(nc,nc)
      z0=[0.0656,0.3711,0.0538,0.0373,0.0261,0.0187,&
         0.0218,0.1791,0.091,0.0605,0.0447,0.0302]
      zi = 0*z0
      zi(1) = 1

      tc=[304.088888888889,190.6,305.4,369.8,425.2,469.6,507.4,616.2,&
         698.9,770.4,853.1,1001.2]
      pc=[73.7343491450634,45.9196083838941,48.7516547159404,42.3795504688362, &
         37.9291919470491,33.6811224489796,29.6353419746277,28.8261858797573,&
         19.3186017650303,16.5876999448428,15.2728212906784,14.6659542195256]
      w= [0.228,0.008,0.098,0.152,0.193,0.251,0.296,&
         0.454,0.787,1.048,1.276,1.299]
      kij = 0
      kij(1, 2) = 0.12
      kij(1, 3:) = 0.15
      kij(:, 1) = kij(1, :)

      model = PengRobinson78(tc, pc, w, kij=kij)
   end function multicomponent_PR

   type(CubicEoS) function asphaltenes_srk(z0, zi) result(model)
      use yaeos, only: SoaveRedlichKwong
      integer, parameter :: nc = 16
      real(pr), intent(out) :: z0(nc), zi(nc)

      real(pr) :: tc(nc), pc(nc), w(nc), kij(nc,nc)

      !names="N2" "CO2" "H2S" "C1" "C2" "C3" "iC4" "nC4" "iC5" "nC5" "C6" "C7-C25" "C26-C49" "C50-C64-PN" "C65-C80-PN" "C50-C80-A"
      z0=[0.0049, 0.11369, 0.0322, 0.27357, 0.09409, 0.06699, 0.0081, 0.0317, 0.0122, 0.0198, 0.0249, 0.25278, 0.05524, 0.00545, 0.00183, 0.00256]
      zi = 0*z0
      zi(2) = 1

      tc=[126.2, 304.2, 373.25, 190.6, 305.4, 369.8, 408.1, 425.2, 460.4, 469.6, 507.4, 666.25, 893.95, 956.95, 1118.45, 1286.75]
      pc=[33.94, 73.76, 89.37, 46., 48.84, 42.46, 36.48, 38., 33.84, 33.74, 29.69, 20.39, 14.19, 13.72, 14.1, 18.11]
      w=[0.04, 0.225, 0.1, 0.008, 0.098, 0.152, 0.176, 0.193, 0.227, 0.251, 0.296, 0.756, 1.262, 1.313, 0.876, 1.274]

      kij(1, :) =  [0., -0.032, 0.17, 0.028, 0.041, 0.076, 0.094, 0.07, 0.087, 0.088, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08]
      kij(2, :) =  [-0.032, 0., 0.099, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.12, 0.1, 0.1, 0.1, 0.1, 0.1]
      kij(3, :) =  [0.17, 0.099, 0., 0.08, 0.085, 0.089, 0.051, 0.06, 0.06, 0.069, 0.05, 0., 0., 0., 0., -0.01]
      kij(4, :) =  [0.028, 0.12, 0.08, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(5, :) =  [0.041, 0.12, 0.085, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(6, :) =  [0.076, 0.12, 0.089, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(7, :) =  [0.094, 0.12, 0.051, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(8, :) =  [0.07, 0.12, 0.06, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(9, :) =  [0.087, 0.12, 0.06, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(10, :) = [0.088, 0.12, 0.069, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(11, :) = [0.08, 0.12, 0.05, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.017]
      kij(12, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(13, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(14, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(15, :) = [0.08, 0.1, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
      kij(16, :) = [0.08, 0.1, 0.0, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.0, 0.0, 0.0, 0.0, 0.0]

      model = SoaveRedlichKwong(tc, pc, w, kij=kij)
   end function asphaltenes_srk

end module fixtures_models
