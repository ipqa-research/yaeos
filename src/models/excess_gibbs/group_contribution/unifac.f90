module yaeos_models_ge_group_contribution_unifac
   use yaeos_constants, only: pr, R
   use yaeos_models_ge, only: GeModel
   implicit none

   type :: Groups
      integer, allocatable :: groups_ids(:) !! Indexes (ids) of each group
      integer, allocatable :: number_of_groups(:) !! \(\nu\)
      real(pr) :: surface_area !! q
      real(pr) :: volume !! r
   end type Groups

   type, extends(GeModel) :: UNIFAC
      !! UNIFAC parameters
      integer :: ngroups !! Total number of individual groups
      real(pr) :: z = 10 !! Model constant
      
      real(pr), allocatable :: group_area(:) !! Q_k
      real(pr), allocatable :: group_volume(:) !! R_k

      class(PsiFunction), allocatable :: psi_function

      type(Groups), allocatable :: molecules(:) !! Substances present in the system
      type(Groups) :: groups_stew !! All the groups present in the system
   contains
      procedure :: excess_gibbs
   end type UNIFAC

   type, abstract :: PsiFunction
   contains
      procedure(temperature_dependence), deferred :: psi
   end type PsiFunction

   type, extends(PsiFunction) :: UNIFACPsi
      real(pr), allocatable :: Eij(:, :)
   contains
      procedure :: psi => UNIFAC_temperature_dependence
   end type

   abstract interface
      subroutine temperature_dependence(self, systems_groups, T, psi, dpsidt, dpsidt2)
         import pr, PsiFunction, Groups
         class(PsiFunction) :: self
         class(Groups) :: systems_groups
         real(pr), intent(in) :: T
         real(pr), intent(out) :: psi(:, :)
         real(pr), optional, intent(out) :: dpsidt(:, :)
         real(pr), optional, intent(out) :: dpsidt2(:, :)
      end subroutine temperature_dependence
   end interface

contains

   subroutine excess_gibbs(self, n, t, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      !! Excess Gibbs and derivs procedure
      class(UNIFAC), intent(in) :: self !! Model
      real(pr), intent(in) ::n(:) !! Moles vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), optional, intent(out) :: Ge !! Excess Gibbs
      real(pr), optional, intent(out) :: GeT !! \(\frac{dG^E}{dT}\)
      real(pr), optional, intent(out) :: GeT2 !! \(\frac{d^2G^E}{dT^2}\)
      real(pr), optional, intent(out) :: Gen(size(n))
      real(pr), optional, intent(out) :: GeTn(size(n))
      real(pr), optional, intent(out) :: Gen2(size(n), size(n))

      real(pr) :: x(size(n))
      real(pr) :: ln_gamma_c(size(n)), ln_gamma_r(size(n)), ln_activity(size(n))


      integer :: i, j, nc

      x = n/sum(n)

      nc = size(self%molecules)

      ln_activity = 0

      call combinatorial_activity(self, x, ln_gamma_c)
      call residual_activity(self, x, T, ln_gamma_r)
      ln_activity = ln_gamma_c + ln_gamma_r

      if (present(Ge)) Ge = sum(x * ln_activity)
   end subroutine excess_gibbs

   subroutine residual_activity(self, x, T, ln_gamma_r)
      class(UNIFAC) :: self
      real(pr), intent(in) :: x(:)
      real(pr), intent(in) :: T
      real(pr), intent(out) :: ln_gamma_r(:)

      real(pr) :: ln_Gamma(size(self%groups_stew%groups_ids))
      real(pr) :: ln_Gamma_i(size(self%groups_stew%groups_ids))

      real(pr) :: xpure(size(x))
      integer :: i, k, nc

      ln_gamma = 0
      nc = size(self%molecules)

      call group_big_gamma(self, x, T, ln_Gamma)
      do i=1,nc
         xpure = 0
         xpure(i) = 1

         ln_gamma_i = 0
         call group_big_gamma(self, xpure, T, ln_Gamma_i)
         do k=1,size(self%molecules(i)%groups_ids)
            ln_gamma_r(i) = ln_gamma_r(i) + self%molecules(i)%number_of_groups(k) * (ln_gamma(k) - ln_gamma_i(k))
         end do
      end do
   end subroutine residual_activity

   subroutine combinatorial_activity(self, x, ln_gamma_c)
      class(UNIFAC) :: self
      real(pr), intent(in) :: x(:)
      real(pr), intent(out) :: ln_gamma_c(:)

      real(pr) :: theta(size(x)), phi(size(x)), L(size(x))

      associate (&
         q => self%molecules%surface_area,&
         r => self%molecules%volume,&
         z => self%z &
         )

         theta = x * q / sum(x * q)
         phi = x * r / sum(x * r)
         L = 0.5_pr * z * (r - q) - (r - 1)
         ln_gamma_c = log(phi/x) + z/2*q * log(theta/phi) + L - phi/x * sum(x*L)
      end associate
   end subroutine combinatorial_activity

   subroutine ln_activity_coefficient(&
      self, x, T, ln_gamma, &
      dln_gammadt, dln_gammadt2, dln_gammadx, dln_gammadtx, dln_gammadx2 &
      )
      class(UNIFAC), intent(in) :: self
      real(pr), intent(in) :: x(:)
      real(pr), intent(in) :: T
      real(pr), intent(out) :: ln_gamma(:)
      real(pr), optional, intent(out) :: dln_gammadt(:)
      real(pr), optional, intent(out) :: dln_gammadt2(:)
      real(pr), optional, intent(out) :: dln_gammadx(:, :)
      real(pr), optional, intent(out) :: dln_gammadtx(:, :)
      real(pr), optional, intent(out) :: dln_gammadx2(:, :, :)

      real(pr) :: ln_gamma_c(size(x))
      real(pr) :: dln_gamma_c_dt(size(x))
      real(pr) :: dln_gamma_c_dt2(size(x))
      real(pr) :: dln_gamma_c_dx (size(x), size(x))
      real(pr) :: dln_gamma_c_dtx(size(x), size(x))
      real(pr) :: dln_gamma_c_dx2(size(x), size(x), size(x))

      real(pr) :: ln_gamma_r(size(x))
      real(pr) :: dln_gamma_r_dt(size(x))
      real(pr) :: dln_gamma_r_dt2(size(x))
      real(pr) :: dln_gamma_r_dx (size(x), size(x))
      real(pr) :: dln_gamma_r_dtx(size(x), size(x))
      real(pr) :: dln_gamma_r_dx2(size(x), size(x), size(x))

      call combinatorial_activity(self, x, ln_gamma_c)
      call residual_activity(self, x, T, ln_gamma_r)

      ln_gamma = ln_gamma_c + ln_gamma_r
   end subroutine ln_activity_coefficient

   subroutine UNIFAC_temperature_dependence(self, systems_groups, T, psi, dpsidt, dpsidt2)
      class(UNIFACPsi) :: self !! \(\psi\) function
      class(Groups) :: systems_groups !! Groups in the system
      real(pr), intent(in) :: T !! Temperature
      real(pr), intent(out) :: psi(:, :) !! \(\psi\)
      real(pr), optional, intent(out) :: dpsidt(:, :) 
      real(pr), optional, intent(out) :: dpsidt2(:, :)

      integer :: i, j
      integer :: ig, jg
      integer :: ngroups

      ngroups = size(systems_groups%groups_ids)

      do concurrent(i=1:ngroups, j=1:ngroups)
         ig = systems_groups%groups_ids(i)
         jg = systems_groups%groups_ids(j)
         psi(i, j) = exp(-self%Eij(ig, jg) / T)
      end do

   end subroutine UNIFAC_temperature_dependence

   subroutine group_area_fraction(&
      self, x, &
      theta, dthetadx, dthetadx2 &
      )
      type(UNIFAC) :: self
      real(pr), intent(in) :: x(:)


      real(pr), intent(out) :: theta(:) !! Group area fraction
      real(pr), optional, intent(out) :: dthetadx(:, :)
      real(pr), optional, intent(out) :: dthetadx2(:, :, :)

      real(pr) :: gf(size(theta)) !! Group fraction
      real(pr) :: dgfdx(size(theta), size(x))
      real(pr) :: dgfdx2(size(theta), size(x), size(x))

      real(pr) :: ga(size(theta)) !! Area of group j in the system
      real(pr) :: dgadx(size(theta), size(x))

      real(pr) :: total_groups_area, dtga_dx(size(x))
      integer :: i, l, j, m, gi

      associate(&
         group_area => self%group_area, stew => self%groups_stew, &
         molecs => self%molecules &
         )

         theta = 0
         dgadx = 0

         total_groups_area = 0
         dtga_dx = 0
         ga = 0

         do l=1,size(molecs)
            do m=1,size(molecs(l)%number_of_groups)
               gi = molecs(l)%groups_ids(m)
               j = findloc(stew%groups_ids, gi, dim=1)

               ga(j) = ga(j) + x(l) * group_area(gi) * molecs(l)%number_of_groups(m)

               if (present(dthetadx)) then
                  dgadx(j, l) = group_area(gi) * molecs(l)%number_of_groups(m)
                  dtga_dx(l) = dtga_dx(l) + molecs(l)%number_of_groups(m) * group_area(gi)
               end if

               total_groups_area = total_groups_area + x(l) * molecs(l)%number_of_groups(m) * group_area(gi)
            end do
         end do

         theta = ga/total_groups_area
         if (present(dthetadx)) then
            do i=1,size(x)
               dthetadx(:, i) = &
                  dgadx(:, i)/total_groups_area &
                  - ga(:)/total_groups_area**2 * dtga_dx(i)

               if (present(dthetadx2)) then
                  do j=i,size(x)
                     dthetadx2(:, i, j) = &
                        dtga_dx(i) * dtga_dx(j) * ga(:)/total_groups_area**3 &
                        - ( dgadx(:, i) * dgadx(:, j)**2 )/total_groups_area**2
                  end do
               end if
            end do
         end if
      end associate
   end subroutine group_area_fraction

   subroutine group_big_gamma(self, x, T, ln_Gamma, dln_Gammadx)
      class(UNIFAC) :: self
      real(pr), intent(in) :: x(:) !! Molar fractions
      real(pr), intent(in) :: T !! Temperature
      real(pr), intent(out) :: ln_Gamma(:) !! \(\ln \Gamma_i\)
      real(pr), optional, intent(out) :: dln_Gammadx(:, :) !!

      real(pr) :: gf(size(self%groups_stew%groups_ids))
      real(pr) :: theta(size(self%groups_stew%groups_ids))
      real(pr) :: dthetadx(size(self%groups_stew%groups_ids), size(x))
      real(pr) :: psi(&
         size(self%groups_stew%groups_ids),&
         size(self%groups_stew%groups_ids) &
         )

      integer :: ng
      integer :: n, m, k, gi
      integer :: i

      logical :: dx

      real(pr) :: up, updx(size(x))
      real(pr) :: updown, updowndx(size(x))
      real(pr) :: down(size(theta)), downdx(size(theta), size(x))

      ng = size(self%groups_stew%groups_ids)

      dx = present(dln_Gammadx)


      if (dx) then
         dln_Gammadx = 0
         call group_area_fraction(self, x, theta=theta, dthetadx=dthetadx)
      else
         call group_area_fraction(self, x, theta=theta)
      end if

      
      call self%psi_function%psi(self%groups_stew, T, psi)

      down = 0
      do concurrent(m=1:ng, n=1:ng)
         down(m) = down(m) + theta(n) * psi(n, m)
         if(dx) then
            downdx(m, :) = downdx(m, :) + dthetadx(n, :) * psi(n, m)
         end if
      end do

      do k=1,ng
         gi = self%groups_stew%groups_ids(k)
         updown = sum(theta(:) * psi(k, :)/down(:))

         ln_Gamma(k) = self%group_area(gi) * (&
            1 - log(sum(theta(:) * psi(:, k))) - updown &
         )

         if (dx) then
            do i=1,size(x)
               dln_Gammadx(k, i) = self%group_area(gi) * (&
                  - sum(psi(:, k) * dthetadx(:, i))/sum(theta(:) * psi(:, k)) &
                  - sum(psi(k, :) * dthetadx(:, i)/down(:)) &
                  + sum(downdx(:, i) * theta(:) * psi(k, :)/down**2) &
               )
            end do
         end if
      end do
   end subroutine group_big_gamma

   type(UNIFAC) function setup_unifac(molecules, Eij, Qk, Rk)
      !! UNIFAC model initialization.
      type(Groups), intent(in) :: molecules(:) !! Molecules
      real(pr), intent(in) :: Eij(:, :) !! Interaction Matrix
      real(pr), intent(in) :: Qk(:) !! Group areas
      real(pr), intent(in) :: Rk(:) !! Group volumes

      type(Groups) :: soup
      type(UNIFACPsi) :: psi_function
      integer :: gi, i, j

      setup_unifac%molecules = molecules

      allocate(soup%groups_ids(0))
      allocate(soup%number_of_groups(0))

      ! ========================================================================
      ! Count all the individual groups and each molecule volume and area
      ! ------------------------------------------------------------------------
      associate(&
         r => setup_unifac%molecules%volume, &
         q => setup_unifac%molecules%surface_area &
         )
         ! Get all the groups indexes and counts into a single stew of groups.
         do i=1,size(molecules)
            r(i) = 0
            q(i) = 0

            do j=1,size(molecules(i)%groups_ids)
               gi = molecules(i)%groups_ids(j)

               ! Calculate molecule i volume and area
               r(i) = r(i) + molecules(i)%number_of_groups(j) * Rk(gi)
               q(i) = q(i) + molecules(i)%number_of_groups(j) * Qk(gi)

               if (all(soup%groups_ids - gi  /= 0)) then
                  ! Add group if it wasn't included yet
                  soup%groups_ids = [soup%groups_ids, gi]
                  soup%number_of_groups = [soup%number_of_groups, 0]
               end if

               ! Find where is the group located in the main soup of
               ! groups.
               gi = findloc(soup%groups_ids - gi, 0, dim=1)

               soup%number_of_groups(gi) = soup%number_of_groups(gi) &
                  + molecules(i)%number_of_groups(j)
            end do
         end do
      end associate
      ! ========================================================================

      psi_function%Eij = Eij
      setup_unifac%groups_stew = soup
      setup_unifac%ngroups = size(soup%number_of_groups)
      setup_unifac%psi_function = psi_function
      setup_unifac%group_area = Qk
      setup_unifac%group_volume = Rk
   end function setup_unifac
end module yaeos_models_ge_group_contribution_unifac