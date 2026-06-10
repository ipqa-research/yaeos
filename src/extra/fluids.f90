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
end module yaeos__extra_fluids
