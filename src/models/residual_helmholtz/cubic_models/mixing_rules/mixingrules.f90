module mixrule_classicvdw
   !!  Mixing rules module
   
   !! Cubic EoS Mixing Rule
   !! Contains a single subroutine as an attribute that receives the
   !! set of parameters of pure compounds and returns the mixture's
   !! parameters and their derivatives
   use constants, only: pr
   use hyperdual_mod
   
   implicit none

   real(pr), allocatable :: kij(:, :)
   real(pr), allocatable :: lij(:, :)

contains
   ! ===========================================================================
   ! Constructors
   ! ---------------------------------------------------------------------------
   subroutine setup_ClassicVdW(kij_in, lij_in)
      real(pr) :: kij_in(:, :)
      real(pr) :: lij_in(:, :)

      kij = kij
      lij = lij
      
   end subroutine
   ! ===========================================================================

   ! ===========================================================================
   ! Mixing subroutines
   ! ---------------------------------------------------------------------------
   pure subroutine mix_ClassicVdW(z, v, t, a, b, c, amix, bmix, cmix)
      type(hyperdual), intent(in) :: z(:), v, t, a(size(z)), b(size(z)), c(size(z))
      type(hyperdual), intent(out) :: amix, bmix, cmix
      
      call a_mix_classic(z, v, t, a, amix)
      call b_mix_classic(z, v, t, b, bmix)
      call c_mix_classic(z, v, t, c, cmix)
   end subroutine

   pure subroutine a_mix_classic(z, v, t, a, amix)
      type(hyperdual), intent(in) :: z(:), v, t, a(size(z))
      type(hyperdual), intent(out) :: amix

      type(hyperdual) :: aij(size(z), size(z))

      integer :: i, j

      aij = 0.0_pr
      amix = 0.0_pr
      do j = 1, size(z)
         aij(:, j) = sqrt(a(:) * a(j)) * (1.0_pr - kij(:, j))
         amix = amix + sum(z(:) * z(j) * aij(:, j))
      end do
   end subroutine

   pure subroutine b_mix_classic(z, v, t, b, bmix)
      type(hyperdual), intent(in) :: z(:), v, t, b(size(z))
      type(hyperdual), intent(out) :: bmix

      type(hyperdual) :: bij(size(z), size(z))

      integer :: i

      bmix = 0.0_pr
      do i = 1, size(z)
         bij(:, i) = (b(:) + b(i))/2.0_pr*(1.0_pr - lij(:, i))
         bmix = bmix + sum(z(:)*z(i)*bij(:, i))
      end do

      bmix = bmix/sum(z)
   end subroutine

   pure subroutine c_mix_classic(z, v, t, c, cmix)
      type(hyperdual), intent(in) :: z(:), v, t, c(size(z))
      type(hyperdual), intent(out) :: cmix

      cmix = sum(z*c)
   end subroutine
   
   pure subroutine del_mix_classic(z, v, t, del, delmix)
      type(hyperdual), intent(in) :: z(:), v, t, del(size(z))
      type(hyperdual), intent(out) :: delmix

      delmix = del(1)
   end subroutine
   ! ===========================================================================
end module