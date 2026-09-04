module yaeos__extra_fluids
   !! # Predefined Fluids
   !! Collection of predefined fluids to be used as examples and in tests.
   !!
   !! # Description
   !! This module collects a set of different fluids intended to be used for tests
   !! or in practice when training on using the library.
   !!
   !! The fluids are stored as functions that return a derived type `FixtureFluid`
   !! which contains the instantiated ArModel of the fluid, and its relevant
   !! compositions `z0` and `zi`
   !!
   !!
   !! # Examples
   !!
   !! ```fortran
   !!  A basic code example
   !! ```
   use yaeos__constants, only: pr, R
   use yaeos__models

   type :: FixtureFluid
      integer :: nc
      !! Number of components
      character(len=:), allocatable :: reference
      !! Reference of where the fluid was taken from
      real(pr), allocatable :: z0(:)
      !! Original composition
      real(pr), allocatable :: zi(:)
      !! Injection composition (if exists)
      class(ArModel), allocatable :: ar_model
   end type FixtureFluid

contains

   type(FixtureFluid) function oil_gao()
      !! J. Gao, R. Okuno, and H. A. Li,
      !! “An Experimental Study of Multiphase Behavior for n-Butane/Bitumen/Water Mixtures,”
      !! SPE Journal, 2016, doi: 10.2118/180736-pa.
      integer, parameter :: nc=6

      type(CubicEoS) :: ar_model

      real(pr) :: tc(nc), pc(nc), w(nc)

      type(QMR) :: mixrule
      real(pr) :: kij(nc, nc), lij(nc, nc)


      tc = [647.05_pr, 425.15_pr, 708.15_pr, 768.25_pr, 998.15_pr, 1346.05_pr]
      pc = [220.64_pr, 37.96_pr, 21.46_pr, 15.07_pr, 13.64_pr, 10.45_pr]
      w = [0.3433_pr, 0.2014_pr, 0.8423_pr, 0.9429_pr, 1.0225_pr, 1.1486_pr]

      ar_model = PengRobinson78(tc, pc, w)

      kij(1, :) = [0.0_pr, 0.636_pr, 0.2006_pr, 0.1694_pr, 0.1694_pr, 0.1694_pr]
      kij(2, :) = [0.636_pr, 0.0_pr, -0.0005_pr, -0.0011_pr, -0.0018_pr, -0.0031_pr]
      kij(3, :) = [0.2006_pr, -0.0005_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(4, :) = [0.1694_pr, -0.0011_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(5, :) = [0.1694_pr, -0.0018_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(6, :) = [0.1694_pr, -0.0031_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]

      lij(1, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(2, :) = [0.0_pr, 0.0_pr, -0.0_pr, -0.0_pr, -0.0_pr, -0.0_pr]
      lij(3, :) = [0.0_pr, -0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(4, :) = [0.0_pr, -0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(5, :) = [0.0_pr, -0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(6, :) = [0.0_pr, -0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]

      mixrule = QMR(k=kij, l=lij)

      call ar_model%set_mixrule(mixrule)

      oil_gao%nc = nc
      oil_gao%reference = "doi: 10.2118/180736-pa"
      oil_gao%ar_model = ar_model
      oil_gao%z0 = [0.6202_pr, 0.3702_pr, 0.0051_pr, 0.0023_pr, 0.0014_pr, 0.0008_pr]
   end function oil_gao

   type(FixtureFluid) function co2_h2o_isop() result(fluid)
      use yaeos__models, only: CubicEoS, SoaveRedlichKwong
      integer, parameter :: nc = 3
      real(pr) :: Tc(nc), Pc(nc), w(nc), kij(nc, nc)
      Tc = [304.21, 647.13, 508.3]
      Pc = [73.3, 220.55, 47.64]
      w = [0.2236, 0.3449, 0.6669]

      kij(1, :) = [0._pr, 0.19_pr, 0.1215_pr]
      kij(2, :) = [0.19_pr, 0._pr, -0.1727_pr]
      kij(3, :) = [0.1215_pr, -0.1727_pr, 0._pr]

      fluid%nc = nc
      fluid%ar_model = SoaveRedlichKwong(Tc, Pc, w, kij=kij)
      fluid%z0 = [0.997493267163008, 0.000567120831996, 0.00193961200499]
      fluid%zi = [0.122878702720943, 0.376116603168507, 0.50100469411055]
   end function co2_h2o_isop
   
   type(FixtureFluid) function oil_b71() result(fluid)
      integer, parameter :: nc=16

      type(CubicEoS) :: ar_model

      real(pr) :: tc(nc), pc(nc), w(nc)

      type(QMR) :: mixrule
      real(pr) :: kij(nc, nc), lij(nc, nc)

      tc = [304.2_pr, 126.2_pr, 190.6_pr, 305.4_pr, 369.8_pr, 408.1_pr, 425.2_pr, 460.4_pr, 469.6_pr, 506.35_pr, 566.55_pr, 647.06_pr, 719.44_pr, 784.93_pr, 846.33_pr, 919.39_pr]
      pc = [72.8_pr, 33.5_pr, 45.4_pr, 48.2_pr, 41.9_pr, 36.0_pr, 37.5_pr, 33.4_pr, 33.3_pr, 33.9_pr, 25.3_pr, 19.1_pr, 14.2_pr, 10.5_pr, 7.5_pr, 4.76_pr]
      w = [0.225_pr, 0.04_pr, 0.008_pr, 0.098_pr, 0.152_pr, 0.176_pr, 0.193_pr, 0.227_pr, 0.251_pr, 0.299_pr, 0.3884_pr, 0.5289_pr, 0.6911_pr, 0.8782_pr, 1.1009_pr, 1.4478_pr]

      ar_model = PengRobinson78(tc, pc, w)

      kij(1, :) = [0.0_pr, -0.02_pr, 0.075_pr, 0.08_pr, 0.08_pr, 0.085_pr, 0.085_pr, 0.085_pr, 0.085_pr, 0.095_pr, 0.095_pr, 0.095_pr, 0.095_pr, 0.095_pr, 0.095_pr, 0.095_pr]
      kij(2, :) = [-0.02_pr, 0.0_pr, 0.08_pr, 0.07_pr, 0.07_pr, 0.06_pr, 0.06_pr, 0.06_pr, 0.06_pr, 0.05_pr, 0.1_pr, 0.12_pr, 0.12_pr, 0.12_pr, 0.12_pr, 0.12_pr]
      kij(3, :) = [0.075_pr, 0.08_pr, 0.0_pr, 0.003_pr, 0.01_pr, 0.018_pr, 0.018_pr, 0.025_pr, 0.026_pr, 0.036_pr, 0.049_pr, 0.073_pr, 0.098_pr, 0.124_pr, 0.149_pr, 0.181_pr]
      kij(4, :) = [0.08_pr, 0.07_pr, 0.003_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(5, :) = [0.08_pr, 0.07_pr, 0.01_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(6, :) = [0.085_pr, 0.06_pr, 0.018_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(7, :) = [0.085_pr, 0.06_pr, 0.018_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(8, :) = [0.085_pr, 0.06_pr, 0.025_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(9, :) = [0.085_pr, 0.06_pr, 0.026_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(10, :) = [0.095_pr, 0.05_pr, 0.036_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(11, :) = [0.095_pr, 0.1_pr, 0.049_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(12, :) = [0.095_pr, 0.12_pr, 0.073_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(13, :) = [0.095_pr, 0.12_pr, 0.098_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(14, :) = [0.095_pr, 0.12_pr, 0.124_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(15, :) = [0.095_pr, 0.12_pr, 0.149_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      kij(16, :) = [0.095_pr, 0.12_pr, 0.181_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]

      lij(1, :) = [0.0_pr, -0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(2, :) = [-0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(3, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(4, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(5, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(6, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(7, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(8, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(9, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(10, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(11, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(12, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(13, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(14, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(15, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]
      lij(16, :) = [0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr]

      mixrule = QMR(k=kij, l=lij)
      call ar_model%set_mixrule(mixrule)
      
      fluid%z0=[0.710319_pr,   0.001392_pr,   0.04727_pr,    0.011687_pr,   0.008613_pr,   0.001044_pr, 0.009541_pr,   0.004582_pr,   0.006235_pr,   0.009628_pr,   0.05340408_pr, 0.04752839_pr, 0.03690424_pr, 0.02809752_pr, 0.01705432_pr, 0.00670045_pr]
      fluid%ar_model = ar_model
      fluid%nc = nc
   end function oil_b71

end module yaeos__extra_fluids
