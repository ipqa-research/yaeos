module yaeos__adiff_hyperdual_ge_api
   !! # Hyperdual excess-Gibbs models
   !! Base type and derivative adapter for GE models using hyperdual AD.
   !!
   !! # Description
   !! Extend `GeModelAdiff` and implement the hyperdual excess-Gibbs function.
   !! The adapter extracts composition, temperature, and mixed derivatives to
   !! provide the standard `GeModel` activity-coefficient interface.
   !!
   !! # Examples
   !! See the group-contribution and NRTL model implementations for examples.
   !!
   !! # References
   !! Hyperdual arithmetic is provided by `hyperdual_mod`.
   use yaeos__constants, only: pr
   use yaeos__models_ge, only: GeModel
   use hyperdual_mod

   implicit none

   type, abstract, extends(GeModel) :: GeModelAdiff
      !! # GeModelAdiff
      !! Base type for hyperdual excess-Gibbs models.
      !!
      !! # Description
      !! Concrete extensions implement the hyperdual-valued `Ge` function.
      !! This adapter evaluates it with derivative seeds in composition and
      !! temperature, then provides the `GeModel` excess-Gibbs derivative API.
      !!
      !! # Examples
      !! See the UNIFAC implementation for a model using hyperdual AD.
   contains
      procedure(hyperdual_ge), deferred :: Ge
      procedure :: excess_gibbs => excess_gibbs
   end type GeModelAdiff

   abstract interface
      type(hyperdual) function hyperdual_Ge(self, n, t)
         !! Evaluate the model's primal excess-Gibbs energy in hyperdual form.
         import hyperdual, GeModelAdiff
         class(GeModelAdiff) :: self !! Concrete GE model; passed-object dummy.
         type(hyperdual), intent(in) :: n(:) !! Component mole amounts as hyperdual values.
         type(hyperdual), intent(in) :: t !! Temperature as a hyperdual value.
      end function hyperdual_Ge
   end interface
contains

   subroutine excess_gibbs(self, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      !! # excess_gibbs
      !! Evaluate excess Gibbs energy and requested derivatives by AD.
      !!
      !! # Description
      !! Seeds hyperdual composition and temperature variables, evaluates the
      !! concrete model's `Ge` function, and returns requested first and second
      !! derivatives through the standard `GeModel` interface.
      class(GeModelAdiff), intent(in) :: self !! Model providing the hyperdual `Ge` function.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
      real(pr), intent(in) :: t !! Temperature [K].
      real(pr), optional, intent(out) :: Ge !! Excess Gibbs energy.
      real(pr), optional, intent(out) :: GeT !! Temperature derivative of excess Gibbs energy.
      real(pr), optional, intent(out) :: GeT2 !! Second temperature derivative of excess Gibbs energy.
      real(pr), optional, dimension(size(n)), intent(out) :: Gen !! Mole-number gradient of excess Gibbs energy.
      real(pr), optional, dimension(size(n)), intent(out) :: GeTn !! Mixed temperature/mole-number derivative of excess Gibbs energy.
      real(pr), optional, intent(out) :: Gen2(size(n), size(n)) !! Hessian with respect to component mole amounts.

      type(hyperdual) :: d_t !! Hyperdual temperature input.
      type(hyperdual) :: d_n(size(n)) !! Hyperdual component mole amounts.
      type(hyperdual) :: d_Ge

      real(pr) :: dGe(size(n)+1, size(n)+1)
      integer :: nc !! Number of components.

      logical :: any_deriv

      any_deriv = .false.

      nc = size(n)

      if (present(GeT) .or. present(GeT2)) then
         any_deriv = .true.
         if (present(GeT2)) then
            call get_dgedt2
         end if

         if (.not. (present(GeT2) .and. .not. present(GeTn))) then
            call get_dgedt
         end if
      end if

      if (present(GeTn)) then
         any_deriv = .true.
         call get_dgedtn
      end if


      if (present(Gen) .or. present(Gen2)) then
         any_deriv = .true.
         if (present(Gen2)) then
            call get_dgedn2
         else
            call get_dgedn
         end if
      end if

      if (present(Ge)) then
         if (.not. any_deriv) then
            call reset_vars
            d_ge = self%Ge(d_n, d_t)
         end if
         Ge = d_ge%f0
      end if

   contains

      subroutine get_dgedn()
         !! # get_dgedn
         !! Extract first mole-number derivatives of the excess Gibbs energy.
         integer :: i !! Loop index `i`.
         integer :: j !! Loop index `j`.

         do i=2, size(n), 2
            call reset_vars
            d_n(i-1)%f1 = 1
            d_n(i )%f2 = 1

            d_Ge = self%Ge(d_n, d_t)

            Gen(i-1) = d_Ge%f1
            Gen(i) = d_Ge%f2
         end do

         if (mod(size(n), 2) /= 0) then
            call reset_vars
            d_n(size(n))%f1 = 1
            d_Ge = self%Ge(d_n, d_t)
            Gen(size(n)) = d_Ge%f1
         end if

      end subroutine get_dgedn

      subroutine get_dgedn2()
         !! # get_dgedn2
         !! Extract the symmetric mole-number Hessian of the excess Gibbs energy.
         integer :: i !! Loop index `i`.
         integer :: j !! Loop index `j`.

         do i=1,size(n)
            do j=i,size(n)
               call reset_vars
               d_n(i)%f1 = 1
               d_n(j)%f2 = 1

               d_Ge = self%Ge(d_n, d_t)

               if(present(Gen)) Gen(i) = d_Ge%f1
               Gen2(i, j) = d_Ge%f12
               Gen2(j, i) = d_Ge%f12
            end do
         end do
      end subroutine get_dgedn2

      subroutine get_dgedtn()
         !! # get_dgedtn
         !! Extract mixed temperature-mole derivatives of the excess Gibbs energy.
         integer :: i !! Loop/index variable.

         do i=1,size(n)
            call reset_vars
            d_n(i)%f1 = 1
            d_t%f2 = 1
            d_Ge = self%Ge(d_n, d_t)
            if (present(Gen)) Gen(i) = d_Ge%f1
            if (present(GeT)) GeT = d_Ge%f2
            GeTn(i) = d_Ge%f12
         end do
      end subroutine get_dgedtn

      subroutine get_dgedt()
         !! # get_dgedt
         !! Extract the first temperature derivative of the excess Gibbs energy.
         call reset_vars
         d_t%f1 = 1
         d_Ge = self%Ge(d_n, d_t)
         GeT = d_Ge%f1
      end subroutine get_dgedt

      subroutine get_dgedt2()
         !! # get_dgedt2
         !! Extract the second temperature derivative of the excess Gibbs energy.
         call reset_vars
         d_t%f1 = 1
         d_t%f2 = 1
         d_Ge = self%Ge(d_n, d_t)
         if (present(GeT)) GeT = d_Ge%f1
         if (present(GeT2)) GeT2 = d_Ge%f12
      end subroutine get_dgedt2

      subroutine reset_vars()
         !! # reset_vars
         !! Reset hyperdual mole and temperature values before seeding a derivative.
         d_n = n
         d_t = t
      end subroutine reset_vars
   end subroutine excess_gibbs

end module yaeos__adiff_hyperdual_ge_api
