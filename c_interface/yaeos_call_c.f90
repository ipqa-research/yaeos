module yaeos__call_c
   use yaeos, only: ArModel, pr
   use iso_c_binding
   implicit none

   type, extends(ArModel) :: CArModel
      integer(c_int) :: id !! ID of the model
      procedure(abs_c_residual_helmholtz), nopass, pointer :: Ar 
         !! Residual Helmholtz energy function
         !! (C function pointer)
   contains
      procedure :: residual_helmholtz => residual_helmholtz
      procedure :: get_v0 => get_v0
   end type CArModel

   abstract interface
      subroutine abs_c_residual_helmholtz(&
         id, nc, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
         ) bind(C)
         import c_double, c_int
         integer(c_int), intent(in) :: id
         integer(c_int), intent(in) :: nc
         real(c_double), intent(in) :: n(nc) !! Moles vector
         real(c_double), intent(in) :: v !! Volume [L]
         real(c_double), intent(in) :: t !! Temperature [K]
         real(c_double), intent(out) :: Ar !! Residual Helmoltz energy
         real(c_double), intent(out) :: ArV !! \(\frac{dAr}{dV}\)
         real(c_double), intent(out) :: ArT !! \(\frac{dAr}{dT}\)
         real(c_double), intent(out) :: ArT2 !! \(\frac{d^2Ar}{dT^2}\)
         real(c_double), intent(out) :: ArTV !! \(\frac{d^2Ar}{dTV}\)
         real(c_double), intent(out) :: ArV2 !! \(\frac{d^2Ar}{dV^2}\)
         real(c_double), intent(out) :: Arn(size(n)) !! \(\frac{dAr}{dn_i}\)
         real(c_double), intent(out) :: ArVn(size(n)) !! \(\frac{d^2Ar}{dVn_i}\)
         real(c_double), intent(out) :: ArTn(size(n)) !! \(\frac{d^2Ar}{dTn_i}\)
         real(c_double), intent(out) :: Arn2(size(n), size(n))!! \(\frac{d^2Ar}{dn_{ij}}\)
      end subroutine abs_c_residual_helmholtz
   end interface
contains
   subroutine residual_helmholtz(&
      self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      class(CArModel), intent(in) :: self
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

      integer :: i

      call self%Ar(&
         self%id, size(n), n, V, T, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2&
         )
   end subroutine residual_helmholtz

   function get_v0(self, n, P, T)
      class(CArModel), intent(in) :: self
      real(pr), intent(in) :: n(:) !! Moles vector
      real(pr), intent(in) :: P !! Pressure [Pa]
      real(pr), intent(in) :: T !! Temperature [K]
      real(pr) :: get_v0 !! Volume [L]

      get_v0 = 0.01_pr
   end function

   subroutine set_model(id, nc, Tc, Pc, w, ar_foo_c) bind(C, name="set_yaeos_model")
      integer(c_int), intent(in), value :: id
      integer(c_int), intent(in), value :: nc
      real(c_double), intent(in) :: Tc(nc) !! Critical temperature [K]
      real(c_double), intent(in) :: Pc(nc) !! Critical pressure [Pa]
      real(c_double), intent(in) :: w(nc) !! Acentric factor
      type(c_funptr), value :: ar_foo_c
      procedure(abs_c_residual_helmholtz), pointer :: ar_foo

      type(CArModel) :: model

      real(c_double) :: n(nc), v, t, Ar
      real(c_double) :: ArV, ArT, ArTV, ArV2, ArT2 
      real(c_double) :: Arn(nc), ArVn(nc), ArTn(nc), Arn2(nc, nc)

      integer(c_int) :: c_id, c_nc

      c_id = id
      c_nc = nc
      n = [1.0_pr, 2.0_pr]
      v = 3.0_pr
      t = 300.0_pr

      call c_f_procpointer(ar_foo_c, ar_foo)
      model%components%Tc = Tc
      print *, "Tc", model%components%Tc
      model%components%Pc = Pc
      model%components%w = w
      model%id = id
      model%Ar => ar_foo
      call model%residual_helmholtz(n, v, t, Ar=Ar, &
         ArV=ArV, ArT=ArT, ArTV=ArTV, ArV2=ArV2, ArT2=ArT2, &
         Arn=Arn, ArVn=ArVn, ArTn=ArTn, Arn2=Arn2&
         )
   end subroutine
end module yaeos__call_c
