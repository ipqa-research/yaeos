module yaeos__adiff_hyperdual_ar_api
   !! # Hyperdual residual Helmholtz models
   !! Base type and derivative adapter for models defined with hyperdual numbers.
   !!
   !! # Description
   !! Extend `ArModelAdiff`, implement its hyperdual `Ar` function and the
   !! inherited volume initializer, and the adapter extracts the first and
   !! second derivatives required by the standard `ArModel` interface.
   !! This is the intended user entry point for operator-overloaded AD; the
   !! adapter routines are its implementation details.
   !!
   !! # Examples
   !! All that is needed is a derived type extending `ArModelAdiff`, a
   !! hyperdual-valued `Ar` function, and the inherited volume initializer.
   !!
   !! ```fortran
   !! module newmodel
   !! use yaeos__adiff_hyperdual_ar_api, only: ArModelAdiff
   !! use yaeos__constants, only: pr
   !! use hyperdual_mod, only: hyperdual
   !! implicit none
   !!
   !! type, extends(ArModelAdiff) :: YourNewModel
   !!    real(pr) :: a = 1.0_pr
   !! contains
   !!    procedure :: Ar => arfun
   !!    procedure :: get_v0 => v0
   !! end type
   !! contains
   !! type(hyperdual) function arfun(self, n, v, t) result(arval)
   !!    class(YourNewModel) :: self
   !!    type(hyperdual), intent(in) :: n(:), v, t
   !!
   !!    arval = sum(n) * t + self%a * v
   !! end function arfun
   !!
   !! function v0(self, n, p, t) result(vol)
   !!    class(YourNewModel), intent(in) :: self
   !!    real(pr), intent(in) :: n(:), p, t
   !!    real(pr) :: vol
   !!
   !!    vol = self%a * sum(n)
   !! end function v0
   !! end module newmodel
   !! ```
   !!
   !! # References
   !! See `example/extra/adiff_pr76.f90` and the user guide's automatic-
   !! differentiation page for a complete implementation.
   !!
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel
   use hyperdual_mod

   implicit none

   type, abstract, extends(ArModel) :: ArModelAdiff
      !! # ArModelAdiff
      !! Base type for residual Helmholtz models defined with hyperdual arithmetic.
      !!
      !! # Description
      !! A concrete extension implements `Ar` as a function returning the
      !! residual Helmholtz energy in hyperdual form. The adapter seeds the
      !! mole amounts, volume, and temperature to extract requested first and
      !! second derivatives, then exposes them through the `ArModel` API.
      !! A concrete type must also implement the inherited `get_v0` initializer.
      !!
      !! # Examples
      !! See the module-level example and `example/extra/adiff_pr76.f90`.
   contains
      procedure(hyperdual_ar), deferred :: Ar
      procedure :: residual_helmholtz => residual_helmholtz
   end type ArModelAdiff

   abstract interface
      type(hyperdual) function hyperdual_Ar(self, n, v, t)
         !! Evaluate the model's residual Helmholtz energy \(A^r(n,V,T)\) [bar L] in hyperdual arithmetic.
         import hyperdual, ArModelAdiff
         class(ArModelAdiff) :: self !! Concrete EOS model; passed-object dummy.
         type(hyperdual), intent(in) :: n(:) !! Component mole amounts as hyperdual values.
         type(hyperdual), intent(in) :: v !! Total volume as a hyperdual value.
         type(hyperdual), intent(in) :: t !! Temperature as a hyperdual value.
      end function hyperdual_Ar
   end interface
contains

   subroutine residual_helmholtz(&
      self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      !! # residual_helmholtz
      !! Evaluate residual Helmholtz energy and requested derivatives by AD.
      !!
      !! # Description
      !! Converts real inputs into hyperdual variables, seeds their derivative
      !! slots, evaluates the concrete model's `Ar` function, and copies the
      !! requested first- and second-order derivatives to optional outputs.
      !!
      !! # Examples
      !! Usually called through inherited methods such as `pressure` or
      !! `lnphi_pt` rather than invoked directly.
      class(ArModelAdiff), intent(in) :: self !! Model instance.
      real(pr), intent(in) :: n(:) !! Component mole amounts [mol].
      real(pr), intent(in) :: v !! Total volume [L].
      real(pr), intent(in) :: t !! Temperature [K].
      real(pr), optional, intent(out) :: Ar !! Residual Helmholtz energy [bar L].
      real(pr), optional, intent(out) :: ArV !! \(\partial A^r/\partial V\) [bar].
      real(pr), optional, intent(out) :: ArT !! \(\partial A^r/\partial T\) [bar L/K].
      real(pr), optional, intent(out) :: ArT2 !! \(\partial^2 A^r/\partial T^2\) [bar L/K^2].
      real(pr), optional, intent(out) :: ArTV !! \(\partial^2 A^r/\partial T\partial V\) [bar/K].
      real(pr), optional, intent(out) :: ArV2 !! \(\partial^2 A^r/\partial V^2\) [bar/L].
      real(pr), optional, dimension(size(n)), intent(out) :: Arn !! \(\partial A^r/\partial n_i\) [bar L/mol].
      real(pr), optional, dimension(size(n)), intent(out) :: ArVn !! \(\partial^2 A^r/\partial V\partial n_i\) [bar/mol].
      real(pr), optional, dimension(size(n)), intent(out) :: ArTn !! \(\partial^2 A^r/\partial T\partial n_i\) [bar L/(K mol)].
      real(pr), optional, intent(out) :: Arn2(size(n), size(n)) !! \(\partial^2 A^r/\partial n_i\partial n_j\) [bar L/mol^2].

      type(hyperdual) :: d_v !! Hyperdual total-volume input.
      type(hyperdual) :: d_t !! Hyperdual temperature input.
      type(hyperdual) :: d_n(size(n)) !! Hyperdual component mole amounts.
      type(hyperdual) :: d_Ar

      logical :: any_deriv

      any_deriv = .false.

      ! 1. Evaluate Second and Cross Derivatives Independently
      if (present(ArV2)) then
         any_deriv = .true.
         call get_dardv2()
      end if

      if (present(ArT2)) then
         any_deriv = .true.
         call get_dardt2()
      end if

      if (present(ArTV)) then
         any_deriv = .true.
         call get_dardvt()
      end if

      if (present(ArVn)) then
         any_deriv = .true.
         call get_dardvn()
      end if

      if (present(ArTn)) then
         any_deriv = .true.
         call get_dardtn()
      end if

      if (present(Arn2)) then
         any_deriv = .true.
         call get_dardn2()
      else if (present(Arn)) then
         ! Only calculate Arn if Arn2 is not present, because get_dardn2 also yields Arn
         any_deriv = .true.
         call get_dardn()
      end if

      ! 2. Evaluate Pure First Derivatives (Only if not already computed!)
      if (present(ArV)) then
         any_deriv = .true.
         if (.not. (present(ArV2) .or. present(ArVn) .or. present(ArTV))) then
            call get_dardv()
         end if
      end if

      if (present(ArT)) then
         any_deriv = .true.
         if (.not. (present(ArT2) .or. present(ArTn) .or. present(ArTV))) then
            call get_dardt()
         end if
      end if

      ! 3. Evaluate Base Function
      if (present(Ar)) then
         if (.not. any_deriv) then
            call reset_vars()
            d_ar = self%Ar(d_n, d_v, d_t)
         end if
         Ar = d_Ar%f0
      end if

   contains

      subroutine get_dardn()
         !! # get_dardn
         !! Extract the first composition derivatives from hyperdual lanes.
         integer :: i !! Loop/index variable.
         do i=2, size(n), 2
            call reset_vars()
            d_n(i-1)%f1 = 1
            d_n(i)%f2 = 1

            d_Ar = self%Ar(d_n, d_v, d_t)

            if(present(Arn)) then
               Arn(i-1) = d_Ar%f1
               Arn(i)   = d_Ar%f2
            end if
         end do

         if (mod(size(n), 2) /= 0) then
            call reset_vars()
            d_n(size(n))%f1 = 1
            d_Ar = self%Ar(d_n, d_v, d_t)
            if(present(Arn)) Arn(size(n)) = d_Ar%f1
         end if
      end subroutine get_dardn

      subroutine get_dardn2()
         !! # get_dardn2
         !! Extract the symmetric mole-number Hessian and requested gradient.
         integer :: i !! Loop index `i`.
         integer :: j !! Loop index `j`.
         do i=1,size(n)
            do j=i,size(n)
               call reset_vars()
               d_n(i)%f1 = 1
               d_n(j)%f2 = 1

               d_Ar = self%Ar(d_n, d_v, d_t)

               if (present(Arn)) then
                  Arn(i) = d_Ar%f1
                  Arn(j) = d_Ar%f2
               end if
               Arn2(i, j) = d_Ar%f12
               Arn2(j, i) = d_Ar%f12
            end do
         end do
      end subroutine get_dardn2

      subroutine get_dardvn()
         !! # get_dardvn
         !! Extract mixed volume-mole derivatives and any requested first derivatives.
         integer :: i !! Loop/index variable.
         do i=1,size(n)
            call reset_vars()
            d_n(i)%f1 = 1
            d_v%f2 = 1
            d_Ar = self%Ar(d_n, d_v, d_t)
            if (present(Arn)) Arn(i) = d_Ar%f1
            if (present(ArV)) ArV = d_Ar%f2
            ArVn(i) = d_Ar%f12
         end do
      end subroutine get_dardvn

      subroutine get_dardtn()
         !! # get_dardtn
         !! Extract mixed temperature-mole derivatives and any requested first derivatives.
         integer :: i !! Loop/index variable.
         do i=1,size(n)
            call reset_vars()
            d_n(i)%f1 = 1
            d_t%f2 = 1
            d_Ar = self%Ar(d_n, d_v, d_t)
            if (present(Arn)) Arn(i) = d_Ar%f1
            if (present(ArT)) ArT = d_Ar%f2
            ArTn(i) = d_Ar%f12
         end do
      end subroutine get_dardtn

      subroutine get_dardv()
         !! # get_dardv
         !! Extract the first derivative with respect to volume.
         call reset_vars()
         d_v%f1 = 1
         d_Ar = self%Ar(d_n, d_v, d_t)
         ArV = d_Ar%f1
      end subroutine get_dardv

      subroutine get_dardt()
         !! # get_dardt
         !! Extract the first derivative with respect to temperature.
         call reset_vars()
         d_t%f1 = 1
         d_Ar = self%Ar(d_n, d_v, d_t)
         ArT = d_Ar%f1
      end subroutine get_dardt

      subroutine get_dardt2()
         !! # get_dardt2
         !! Extract the second temperature derivative and optional first derivative.
         call reset_vars()
         d_t%f1 = 1
         d_t%f2 = 1
         d_Ar = self%Ar(d_n, d_v, d_t)
         if (present(ArT)) ArT = d_Ar%f1
         ArT2 = d_Ar%f12
      end subroutine get_dardt2

      subroutine get_dardv2()
         !! # get_dardv2
         !! Extract the second volume derivative and optional first derivative.
         call reset_vars()
         d_v%f1 = 1
         d_v%f2 = 1
         d_Ar = self%Ar(d_n, d_v, d_t)
         if (present(ArV)) ArV = d_Ar%f1
         ArV2 = d_Ar%f12
      end subroutine get_dardv2

      subroutine get_dardvt()
         !! # get_dardvt
         !! Extract the mixed volume-temperature derivative and first derivatives.
         call reset_vars()
         d_v%f1 = 1
         d_t%f2 = 1
         d_Ar = self%Ar(d_n, d_v, d_t)
         if (present(ArV)) ArV = d_Ar%f1
         if (present(ArT)) ArT = d_Ar%f2
         ArTV = d_Ar%f12
      end subroutine get_dardvt

      subroutine reset_vars()
         !! # reset_vars
         !! Restore hyperdual inputs to the unperturbed state before a new seed.
         d_n = n
         d_v = v
         d_t = t
      end subroutine reset_vars
   end subroutine residual_helmholtz

end module yaeos__adiff_hyperdual_ar_api
