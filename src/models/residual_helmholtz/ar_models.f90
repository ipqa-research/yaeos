module yaeos__models_ar
   !! Module that defines the basics of a residual Helmholtz energy.
   !!
   !! All the residual properties that are calculated in this library are
   !! based on residual Helmholtz Equations of State. Following the book by
   !! Michelsen and Mollerup.
   !!
   !! In this library up to second derivatives of residual Helmholtz energy
   !! are used. Because they're the fundamentals for phase equilibria
   !! calculation.
   !!
   !! @note Later on, third derivative with respect to volume will be included
   !! since it's importance on calculation of critical points.
   use yaeos__constants, only: pr
   use yaeos__models_base, only: BaseModel
   implicit none

   type, abstract, extends(BaseModel) :: ArModel
      !! Abstract residual Helmholtz model.
      !!
      !! This derived type defines the basics needed for the calculation
      !! of residual properties.
      !! The basics of a residual Helmholtz model is a routine that calculates
      !! all the needed derivatives of \(Ar\) `residual_helmholtz` and
      !! a volume initializer function, that is used to initialize a Newton
      !! solver of volume when specifying pressure.
      character(len=:), allocatable :: name !! Name of the model
   contains
      procedure(abs_residual_helmholtz), deferred :: residual_helmholtz
      !! Method to calculate residual helmholtz energy and derivatives.
      procedure(abs_volume_initializer), deferred :: get_v0
      !! Volume initializer
   end type ArModel

   abstract interface
      subroutine abs_residual_helmholtz(&
         self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
         )
         !! Residual Helmholtz model generic interface.
         !!
         !! This interface represents how an Ar model should be implemented.
         !! By our standard, a Resiudal Helmholtz model takes as input:
         !!
         !! - The mixture's number of moles vector.
         !! - Volume, by default in liters.
         !! - Temperature, by default in Kelvin.
         !!
         !! All the output arguments are optional. While this keeps a long
         !! signature for the implementation, this is done this way to take
         !! advantage of any inner optimizations to calculate derivatives
         !! inside the procedure.
         !!
         !! Once the model is implemented, the signature can be short like
         !! `model%residual_helmholtz(n, v, t, ArT2=dArdT2)`
         import ArModel, pr
         class(ArModel), intent(in) :: self !! ArModel
         real(pr), intent(in) :: n(:) !! Moles vector
         real(pr), intent(in) :: v !! Volume [L]
         real(pr), intent(in) :: t !! Temperature [K]
         real(pr), optional, intent(out) :: Ar !! Residual Helmoltz energy
         real(pr), optional, intent(out) :: ArV !! \(\frac{dAr}{dV}\)
         real(pr), optional, intent(out) :: ArT !! \(\frac{dAr}{dT}\)
         real(pr), optional, intent(out) :: ArT2 !! \(\frac{d^2Ar}{dT^2}\)
         real(pr), optional, intent(out) :: ArTV !! \(\frac{d^2Ar}{dTV}\)
         real(pr), optional, intent(out) :: ArV2 !! \(\frac{d^2Ar}{dV^2}\)
         real(pr), optional, intent(out) :: Arn(size(n)) !! \(\frac{dAr}{dn_i}\)
         real(pr), optional, intent(out) :: ArVn(size(n)) !! \(\frac{d^2Ar}{dVn_i}\)
         real(pr), optional, intent(out) :: ArTn(size(n)) !! \(\frac{d^2Ar}{dTn_i}\)
         real(pr), optional, intent(out) :: Arn2(size(n), size(n))!! \(\frac{d^2Ar}{dn_{ij}}\)
      end subroutine abs_residual_helmholtz

      function abs_volume_initializer(self, n, p, t)
         !! Initialization of volume.
         import ArModel, pr
         class(ArModel), intent(in) :: self !! Ar Model
         real(pr), intent(in) :: n(:) !! Moles vector
         real(pr), intent(in) :: p !! Pressure [bar]
         real(pr), intent(in) :: t !! Temperature [K]
         real(pr) :: abs_volume_initializer !! Initial volume [L]
      end function abs_volume_initializer
   end interface

contains

   subroutine volume(eos, n, T, P, V, root_type)
      !! Generic volume solver
      use yaeos__constants, only: pr, R
      use yaeos__math, only: newton
      class(ArModel), intent(in) :: eos
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: P, T
      real(pr), intent(out) :: V
      character(len=*), intent(in) :: root_type

      integer :: max_iters=30
      real(pr) :: tol=1e-7

      real(pr) :: totnRT, GrL, GrV, Gr
      real(pr) :: Vliq, Vvap

      GrL = HUGE(GrL)
      GrV = HUGE(GrV)

      totnRT = sum(n) * R * T
      select case(root_type)
       case("liquid")
         Vliq = eos%get_v0(n, P, T)! *1.001_pr
         call newton(foo, Vliq, tol=tol, max_iters=max_iters)
         GrL = Gr
       case("vapor")
         Vvap = R * T / P
         call newton(foo, Vvap, tol=tol, max_iters=max_iters)
         GrV = Gr
       case("stable")
         Vliq = eos%get_v0(n, P, T)*1.00001_pr
         call newton(foo, Vliq, tol=tol, max_iters=max_iters)
         GrL = Gr
         
         Vvap = R * T / P
         call newton(foo, Vvap, tol=tol, max_iters=max_iters)
         GrV = Gr
      end select

      if (GrL < GrV) then
         V = Vliq
      else
         V = Vvap
      end if

   contains
      subroutine foo(x, f, df)
         real(pr), intent(in) :: x
         real(pr), intent(out) :: f, df
         real(pr) :: Ar, ArV, ArV2, Pcalc, dPcalcdV, Vin
         Vin = x
         call eos%residual_helmholtz(n, Vin, T, Ar=Ar, ArV=ArV, ArV2=ArV2)
         Pcalc = totnRT / Vin - ArV
         dPcalcdV = -totnRT / Vin**2 - ArV2
         f = Pcalc - p
         df = dPcalcdV
         Gr = Ar + P * Vin - totnRT - totnRT * log(P*Vin/(R*T))
      end subroutine foo
   end subroutine volume
end module yaeos__models_ar
