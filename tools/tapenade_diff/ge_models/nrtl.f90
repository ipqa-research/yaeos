module nrtl_mod
   use yaeos___tapenade_ge_api, only: GeModelTapenade
   ! Replace pr and R parameters with the following line on the generated
   ! code. Generated code sometimes misses those variables and defines them
   ! (with no value) inside.
   ! use yaeos___constants, only: pr, R
   implicit none

   integer, parameter :: pr = 8
   real(8), parameter :: R=0.08314472

   type :: NRTL
      real(8), allocatable :: a(:, :)
      real(8), allocatable :: b(:, :)
      real(8), allocatable :: c(:, :)
   contains
      ! Point to the generated procedures procedures 
      ! and replate `type :: NRTL` with `type, extends(GeModelTapenade) :: NRTL`
      procedure :: ge => excess_gibbs
      procedure :: ge_b => excess_gibbs_b
      procedure :: ge_d => excess_gibbs_d
      procedure :: ge_d_b => excess_gibbs_d_b
      procedure :: ge_d_d => excess_gibbs_d_d
   end type NRTL

contains

   type(NRTL) function setup(a_mat, b_mat, c_mat) result(model)
      ! A setup function for the NRTL model.
      ! Sets the attributes with provided values.
      real(8), intent(in) :: a_mat(:, :)
      real(8), intent(in) :: b_mat(:, :)
      real(8), intent(in) :: c_mat(:, :)

      model%a = a_mat
      model%b = b_mat
      model%c = c_mat
   end function setup

   subroutine excess_gibbs(model, n, T, ge)
      class(NRTL) :: model
      real(8), intent(in) :: n(:)
      real(8), intent(in) :: T
      real(8), intent(out) :: ge

      real(8) :: x(size(n)), G(size(n), size(n)), tau(size(n), size(n))

      ! Due to a tapenade bug, the model attributes are defined as variables,
      ! though they aren't used. This makes it possible to define sizes for
      ! generated auxiliar variables.
      real(8) :: a(size(n), size(n)), b(size(n), size(n)), c(size(n), size(n))
      real(8) :: down
      integer :: i, j

      x = n/sum(n)

      tau = model%a(:, :) + model%b(:, :)/T
      G = exp(-model%c * tau)

      ge = 0
      do i=1,size(n)
         ge = ge + x(i) * sum(x(:) * tau(:, i) * G(:, i)) / sum(x(:) *  g(:, i))
      end do
      ge = sum(n) * R * T * ge
   end subroutine excess_gibbs
end module nrtl_mod
