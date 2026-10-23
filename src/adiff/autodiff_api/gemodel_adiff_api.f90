module yaeos__adiff_hyperdual_ge_api
   !! Module that contains the automatic differentiation logic for an Ge model.
   !!
   !! All that is needed to define an Ge model that uses automatic
   !! differentiation with hyperdual numbers is to define a new derived type
   !! that overloads the method to the Ge function that you want to use.
   use yaeos__constants, only: pr
   use yaeos__models_ge, only: GeModel
   use hyperdual_mod

   implicit none

   type, abstract, extends(GeModel) :: GeModelAdiff
   contains
      procedure(hyperdual_ge), deferred :: Ge
      procedure :: excess_gibbs => excess_gibbs
   end type GeModelAdiff

   abstract interface
      type(hyperdual) function hyperdual_Ge(self, n, t)
         import hyperdual, GeModelAdiff
         class(GeModelAdiff) :: self
         type(hyperdual), intent(in) :: n(:)
         type(hyperdual), intent(in) :: t
      end function hyperdual_Ge
   end interface
contains

   subroutine excess_gibbs(self, n, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      class(GeModelAdiff), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: t
      real(pr), optional, intent(out) :: Ge, GeT, GeT2
      real(pr), optional, dimension(size(n)), intent(out) :: Gen, GeTn
      real(pr), optional, intent(out) :: Gen2(size(n), size(n))

      type(hyperdual) :: d_t, d_n(size(n))
      type(hyperdual) :: d_Ge

      real(pr) :: dGe(size(n)+1, size(n)+1)
      integer :: nc

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
         integer :: i, j

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
         integer :: i, j

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
         integer :: i

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
         call reset_vars
         d_t%f1 = 1
         d_Ge = self%Ge(d_n, d_t)
         GeT = d_Ge%f1
      end subroutine get_dgedt

      subroutine get_dgedt2()
         call reset_vars
         d_t%f1 = 1
         d_t%f2 = 1
         d_Ge = self%Ge(d_n, d_t)
         if (present(GeT)) GeT = d_Ge%f1
         if (present(GeT2)) GeT2 = d_Ge%f12
      end subroutine get_dgedt2

      subroutine reset_vars()
         d_n = n
         d_t = t
      end subroutine reset_vars
   end subroutine excess_gibbs

end module yaeos__adiff_hyperdual_ge_api
