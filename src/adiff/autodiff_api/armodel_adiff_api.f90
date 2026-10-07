module yaeos__adiff_hyperdual_ar_api
   !! Module that contains the automatic differentiation logic for an Ar model.
   !!
   !! All that is needed to define an Ar model that uses automatic
   !! differentiation with hyperdual numbers is to define a new derived type
   !! that overloads the method to the Ar function that you want to use.
   !! A minimal example follows:
   !!
   !! ```fortran
   !! module newmodel
   !! use yaeos__adiff_hyperdual_ar_api, only: ArModelAdiff
   !!
   !! type, extends(ArModelAdiff) :: YourNewModel
   !!       type(Substances) :: composition
   !!       real(8) :: parameters(:)
   !!     contains
   !!       procedure :: Ar => arfun
   !!       procedure :: get_v0 => v0
   !! end type
   !! contains
   !! subroutine arfun(self, n, v, t, Ar)
   !!    class(YourNewModel), intent(in) :: self
   !!    type(hyperdual), intent(in) :: n(:) ! Number of moles
   !!    type(hyperdual), intent(in) :: v ! Volume [L]
   !!    type(hyperdual), intent(in) :: t ! Temperature [K]
   !!    type(hyperdual), intent(out) :: ar_value ! Residual Helmholtz Energy
   !!
   !!    ! A very complicated residual helmholtz function of a mixture
   !!    Ar = sum(n) * v * t
   !! end subroutine
   !!
   !! function v0(self, n, p, t)
   !!    class(YourNewModel), intent(in) :: self
   !!    real(pr), intent(in) :: n(:) ! Number of moles
   !!    real(pr), intent(in) :: p ! Pressure [bar]
   !!    real(pr), intent(in) :: t ! Temperature [K]
   !!    real(pr) :: v0
   !!
   !!    v0 = self%parameters(3)
   !! end function
   !! ```
   !!
   !! A complete implementation of the PR76 Equation of State can me found in
   !! `example/adiff/adiff_pr76.f90`
   !!
   use yaeos__constants, only: pr
   use yaeos__models_ar, only: ArModel
   use hyperdual_mod

   implicit none

   type, abstract, extends(ArModel) :: ArModelAdiff
   contains
      procedure(hyperdual_ar), deferred :: Ar
      procedure :: residual_helmholtz => residual_helmholtz
   end type ArModelAdiff

   abstract interface
      type(hyperdual) function hyperdual_Ar(self, n, v, t)
         import hyperdual, ArModelAdiff
         class(ArModelAdiff) :: self
         type(hyperdual), intent(in) :: n(:)
         type(hyperdual), intent(in) :: v
         type(hyperdual), intent(in) :: t
      end function hyperdual_Ar
   end interface
contains

   subroutine residual_helmholtz(&
      self, n, v, t, Ar, ArV, ArT, ArTV, ArV2, ArT2, Arn, ArVn, ArTn, Arn2 &
      )
      class(ArModelAdiff), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: v, t
      real(pr), optional, intent(out) :: Ar, ArV, ArT, ArT2, ArTV, ArV2
      real(pr), optional, dimension(size(n)), intent(out) :: Arn, ArVn, ArTn
      real(pr), optional, intent(out) :: Arn2(size(n), size(n))

      type(hyperdual) :: d_v, d_t, d_n(size(n))
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
         integer :: i
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
         integer :: i, j
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
         integer :: i
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
         integer :: i
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
         call reset_vars()
         d_v%f1 = 1
         d_Ar = self%Ar(d_n, d_v, d_t)
         ArV = d_Ar%f1
      end subroutine get_dardv

      subroutine get_dardt()
         call reset_vars()
         d_t%f1 = 1
         d_Ar = self%Ar(d_n, d_v, d_t)
         ArT = d_Ar%f1
      end subroutine get_dardt

      subroutine get_dardt2()
         call reset_vars()
         d_t%f1 = 1
         d_t%f2 = 1
         d_Ar = self%Ar(d_n, d_v, d_t)
         if (present(ArT)) ArT = d_Ar%f1
         ArT2 = d_Ar%f12
      end subroutine get_dardt2

      subroutine get_dardv2()
         call reset_vars()
         d_v%f1 = 1
         d_v%f2 = 1
         d_Ar = self%Ar(d_n, d_v, d_t)
         if (present(ArV)) ArV = d_Ar%f1
         ArV2 = d_Ar%f12
      end subroutine get_dardv2

      subroutine get_dardvt()
         call reset_vars()
         d_v%f1 = 1
         d_t%f2 = 1
         d_Ar = self%Ar(d_n, d_v, d_t)
         if (present(ArV)) ArV = d_Ar%f1
         if (present(ArT)) ArT = d_Ar%f2
         ArTV = d_Ar%f12
      end subroutine get_dardvt

      subroutine reset_vars()
         d_n = n
         d_v = v
         d_t = t
      end subroutine reset_vars
   end subroutine residual_helmholtz

end module yaeos__adiff_hyperdual_ar_api
