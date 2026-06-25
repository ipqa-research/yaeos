module yaeos__models_ar_c_interface
   !! Module to inject externally defined models inside yaeos.
   use iso_c_binding
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel
   implicit none

   abstract interface
      subroutine c_ar_func(nc, n, v, t, &
         calc_Ar,  Ar,  &
         calc_ArV, ArV, &
         calc_ArT, ArT, &
         calc_ArTV, ArTV, &
         calc_ArV2, ArV2, &
         calc_ArT2, ArT2, &
         calc_Arn,  Arn,  &
         calc_ArVn, ArVn, &
         calc_ArTn, ArTn, &
         calc_Arn2, Arn2  &
         ) bind(C)
         !! Definition of an ArModel on the C-side.
         import c_int, c_double, c_bool
         integer(c_int), value, intent(in) :: nc
         real(c_double), intent(in)  :: n(nc)
         real(c_double), value, intent(in)  :: v, t
         logical(c_bool), value, intent(in) :: &
            calc_Ar, calc_ArV, calc_ArT, calc_ArTV, &
            calc_ArV2, calc_ArT2, calc_Arn, calc_ArVn, &
            calc_ArTn, calc_Arn2
         real(c_double), intent(out) :: Ar, ArV, ArT, ArTV, ArV2, ArT2
         real(c_double), intent(out) :: Arn(nc), ArVn(nc), ArTn(nc)
         real(c_double), intent(out) :: Arn2(nc, nc)
      end subroutine c_ar_func
   end interface

   type, extends(ArModel) :: CInjectedModel
      !! ArModel injected from C.
      !! Holds a pointer to a C-defined function that calculates the residual
      !! Helmholtz energy and its derivatives.
      type(c_funptr) :: c_func = c_null_funptr
      integer :: nc = 0
   contains
      procedure :: residual_helmholtz => call_c_func
      procedure :: get_v0 => get_v0
   end type CInjectedModel

contains

   subroutine call_c_func(&
      self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      !! Call the C-side function
      class(CInjectedModel), intent(in) :: self
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

      procedure(c_ar_func), pointer :: f
      integer(c_int) :: nc
      real(c_double) :: Ar_, ArV_, ArT_, ArTV_, ArV2_, ArT2_
      real(c_double), allocatable :: Arn_(:), ArVn_(:), ArTn_(:), Arn2_(:,:)

      ! Convert the C pointer to a callable subroutine
      call c_f_procpointer(self%c_func, f)
      nc = size(n)
      allocate(Arn_(nc), ArVn_(nc), ArTn_(nc), Arn2_(nc,nc))

      ! Call the C procedure
      call f(nc, n, v, t, &
         logical(present(Ar),   c_bool), Ar_,   &
         logical(present(ArV),  c_bool), ArV_,  &
         logical(present(ArT),  c_bool), ArT_,  &
         logical(present(ArTV), c_bool), ArTV_, &
         logical(present(ArV2), c_bool), ArV2_, &
         logical(present(ArT2), c_bool), ArT2_, &
         logical(present(Arn),  c_bool), Arn_,  &
         logical(present(ArVn), c_bool), ArVn_, &
         logical(present(ArTn), c_bool), ArTn_, &
         logical(present(Arn2), c_bool), Arn2_  &
         )

      ! Copy requested results
      if (present(Ar))   Ar   = Ar_
      if (present(ArV))  ArV  = ArV_
      if (present(ArT))  ArT  = ArT_
      if (present(ArTV)) ArTV = ArTV_
      if (present(ArV2)) ArV2 = ArV2_
      if (present(ArT2)) ArT2 = ArT2_
      if (present(Arn))  Arn  = Arn_
      if (present(ArVn)) ArVn = ArVn_
      if (present(ArTn)) ArTn = ArTn_
      if (present(Arn2)) Arn2 = Arn2_
   end subroutine call_c_func

   real(pr) function get_v0(self, n, P, T)
      class(CInjectedModel), intent(in) :: self
      real(pr), intent(in) :: n(:), P, T
      get_v0 = 0.01
   end function get_v0

   ! Factory: crea un CInjectedModel desde C/Python
   subroutine make_c_model(f_ptr, nc, model_ptr) bind(C, name="make_c_model")
      !! Construct the thermodynamic model
      type(c_funptr), value, intent(in) :: f_ptr
      integer(c_int), value, intent(in) :: nc
      type(c_ptr), intent(out) :: model_ptr

      type(CInjectedModel), pointer :: model
      allocate(model)
      model%c_func = f_ptr
      model%nc = nc
      model_ptr = c_loc(model)
   end subroutine make_c_model

end module yaeos__models_ar_c_interface
