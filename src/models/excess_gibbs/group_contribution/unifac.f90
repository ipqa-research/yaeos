module yaeos__models_ge_group_contribution_unifac
   !! UNIFAC module
   use yaeos__constants, only: pr, R
   use yaeos__models_ge, only: GeModel
   implicit none

   type :: Groups
      !! Group derived type.
      !!
      !! This type represent a molecule and it's groups
      !!
      integer, allocatable :: groups_ids(:)
      !! Indexes (ids) of each group in the main group matrix
      integer, allocatable :: number_of_groups(:)
      !! Count of each group in the molecule
      real(pr) :: surface_area
      !! Molecule surface area \(q\)
      real(pr) :: volume
      !! Molecule volume \(r\)
   end type Groups

   type, extends(GeModel) :: UNIFAC
      !! UNIFAC model derived type
      !!
      !! This type holds the needed parameters for using a UNIFAC \(G^E\) model
      !! mainly group areas, volumes and what temperature dependence function
      !! \(\psi(T)\) to use.
      !! It also holds the individual molecules of a particular system and
      !! the set of all groups in the system as a "stew" of groups instead of
      !! being them included in particular molecules.
      integer :: ngroups
      !! Total number of individual groups
      real(pr) :: z = 10
      !! Model constant
      real(pr), allocatable :: group_area(:)
      !! Group areas \(Q_k\)
      real(pr), allocatable :: group_volume(:)
      !! Group volumes \(R_k\)
      class(PsiFunction), allocatable :: psi_function
      !! Temperature dependance function of the model
      type(Groups), allocatable :: molecules(:)
      !! Substances present in the system
      type(Groups) :: groups_stew
      !! All the groups present in the system
   contains
      procedure :: excess_gibbs
   end type UNIFAC

   type, abstract :: PsiFunction
   contains
      procedure(temperature_dependence), deferred :: psi
   end type PsiFunction

   type, extends(PsiFunction) :: UNIFACPsi
      !! Original UNIFAC \(psi\) function
      !! \[
      !!    \psi_{ij}(T) = \exp(-\frac{A_{ij}}{T})
      !! \]
      real(pr), allocatable :: Eij(:, :)
   contains
      procedure :: psi => UNIFAC_temperature_dependence
   end type UNIFACPsi

   abstract interface
      subroutine temperature_dependence(&
         self, systems_groups, T, psi, dpsidt, dpsidt2&
         )
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
      use iso_fortran_env, only: error_unit
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

      write(error_unit, *) "WARN: UNIFAC not fully implemented yet"

      x = n/sum(n)

      nc = size(self%molecules)

      ln_activity = 0

      call combinatorial_activity(self, x, ln_gamma_c)
      call residual_activity(self, x, T, ln_gamma_r)
      ln_activity = ln_gamma_c + ln_gamma_r

      if (present(Ge)) Ge = sum(x * ln_activity)
      if (present(GeN)) Gen = ln_activity * (R*T)
   end subroutine excess_gibbs

   subroutine residual_activity(&
      self, x, T, ln_gamma_r, dln_gamma_dn, dln_gamma_dt, dln_gamma_dtn&
      )
      class(UNIFAC) :: self
      real(pr), intent(in) :: x(:)
      real(pr), intent(in) :: T
      real(pr), intent(out) :: ln_gamma_r(:)
      real(pr), optional, intent(out) :: dln_gamma_dn(:, :)
      real(pr), optional, intent(out) :: dln_gamma_dT(:)
      real(pr), optional, intent(out) :: dln_gamma_dTn(:, :)

      real(pr) :: ln_G(self%ngroups)
      real(pr) :: dln_Gdn(self%ngroups, size(x))

      real(pr) :: ln_G_i(size(x), self%ngroups)
      real(pr) :: dln_G_idn(size(x), self%ngroups, size(x))

      real(pr) :: xpure(size(x))
      integer :: i, j, k, nc, gi

      ln_G = 0
      ln_G_i = 0
      ln_gamma_r = 0

      if (present(dln_gamma_dn)) dln_gamma_dn = 0
      nc = size(self%molecules)

      call group_big_gamma(self, x, T, ln_G, dln_gammadx=dln_Gdn)

      do i=1,nc
         xpure = 0
         xpure(i) = 1
         call group_big_gamma(&
            self, xpure, T, ln_G_i(i, :), dln_Gammadx=dln_G_idn(i, :, :) &
            )
      end do

      do i=1,nc
         do k=1,size(self%molecules(i)%groups_ids)
            ! Get the index of the group k of molecule i in the whole stew of
            ! groups.
            ! TODO: These kind of `findloc` maybe can be optimized.
            gi = findloc(&
               self%groups_stew%groups_ids - self%molecules(i)%groups_ids(k), &
               0, dim=1 &
               )

            ln_gamma_r(i) = ln_gamma_r(i) &
               + self%molecules(i)%number_of_groups(k) &
               * (ln_G(gi) - ln_G_i(i, gi))

            if (present(dln_gamma_dn)) then
               ! Compositional derivative
               do j=1,i
                  dln_gamma_dn(i, j) = dln_gamma_dn(i, j) &
                     + self%molecules(i)%number_of_groups(k) &
                     * (dln_Gdn(gi, j) - dln_G_idn(i, gi, j))
                  dln_gamma_dn(j, i) = dln_gamma_dn(i, j)
               end do
            end if
         end do
      end do
   end subroutine residual_activity

   subroutine combinatorial_activity(self, x, ln_gamma_c, dln_gamma_dx)
      class(UNIFAC) :: self
      real(pr), intent(in) :: x(:)
      real(pr), intent(out) :: ln_gamma_c(:)
      real(pr), optional, intent(out) :: dln_gamma_dx(:, :)

      real(pr) :: theta(size(x)), phi(size(x)), L(size(x))
      real(pr) :: V(size(x)), dVdx(size(x), size(x)), Vsum
      real(pr) :: F(size(x)), dFdx(size(x), size(x)), Fsum
      real(pr) :: dthetadx(size(x), size(x))
      real(pr) :: dphidx(size(x), size(x))

      real(pr) :: xq, xr
      integer :: i, j


      associate (&
         q => self%molecules%surface_area,&
         r => self%molecules%volume,&
         z => self%z &
         )
         xq = dot_product(x, q)
         xr = dot_product(x, r)

         V = r / xr
         Vsum = 1/xr
         F = q / xq
         Fsum = 1/xq

         theta = x * q / xq
         phi = x * r / xr
         L = 0.5_pr * z * (r - q) - (r - 1)
         ln_gamma_c = log(phi/x) + z/2*q * log(theta/phi) + L - phi/x * sum(x*L)

         if (present(dln_gamma_dx)) then
            dln_gamma_dx = 0
            do concurrent(i=1:size(x), j=1:size(x))
               dVdx(i, j) = -r(i)*r(j)*Vsum**2
               dFdx(i, j) = -q(i)*q(j)*Fsum**2

               dln_gamma_dx(i, j) = -0.5_pr * z * q(i) * ( &
                  (dVdx(i,j)/F(i) - V(i)*dFdx(i,j)/F(i)**2) * F(i)/V(i) &
                  - dVdx(i, j)/F(i) + V(i) * dFdx(i,j)/F(i)**2 &
                  ) &
                  - dVdx(i, j) + dVdx(i, j)/V(i)
            end do
         end if
      end associate
   end subroutine combinatorial_activity

   subroutine ln_activity_coefficient(&
      self, n, T, lngamma, &
      dln_gammadt, dln_gammadt2, dln_gammadn, dln_gammadtn, dln_gammadn2 &
      )
      class(UNIFAC), intent(in) :: self
      real(pr), intent(in) :: n(:)
      real(pr), intent(in) :: T
      real(pr), intent(out) :: lngamma(:)
      real(pr), optional, intent(out) :: dln_gammadt(:)
      real(pr), optional, intent(out) :: dln_gammadt2(:)
      real(pr), optional, intent(out) :: dln_gammadn(:, :)
      real(pr), optional, intent(out) :: dln_gammadtn(:, :)
      real(pr), optional, intent(out) :: dln_gammadn2(:, :, :)

      real(pr) :: ln_gamma_c(size(n))
      real(pr) :: dln_gamma_c_dt(size(n))
      real(pr) :: dln_gamma_c_dt2(size(n))
      real(pr) :: dln_gamma_c_dn (size(n), size(n))
      real(pr) :: dln_gamma_c_dtn(size(n), size(n))
      real(pr) :: dln_gamma_c_dn2(size(n), size(n), size(n))

      real(pr) :: ln_gamma_r(size(n))
      real(pr) :: dln_gamma_r_dt(size(n))
      real(pr) :: dln_gamma_r_dt2(size(n))
      real(pr) :: dln_gamma_r_dn (size(n), size(n))
      real(pr) :: dln_gamma_r_dtn(size(n), size(n))
      real(pr) :: dln_gamma_r_dn2(size(n), size(n), size(n))

      call combinatorial_activity(self, n, ln_gamma_c)
      call residual_activity(self, n, T, ln_gamma_r)

      lngamma = ln_gamma_c + ln_gamma_r
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
         if (present(dpsidt)) &
            dpsidt(i, j) = self%Eij(ig, jg) * psi(i, j) / T**2
         if (present(dpsidt2)) &
            dpsidt2(i, j) = &
            self%Eij(ig, jg) * (self%Eij(ig, jg) - 2*T) * psi(i, j) / T**4
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

               total_groups_area = total_groups_area &
                  + x(l) * molecs(l)%number_of_groups(m) * group_area(gi)
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

   subroutine group_big_gamma(&
      self, x, T, ln_Gamma, dln_Gammadx, dln_Gammadt, dln_Gammadt2, dln_Gammadx2&
      )
      class(UNIFAC) :: self
      real(pr), intent(in) :: x(:) !< Molar fractions
      real(pr), intent(in) :: T !! Temperature
      real(pr), intent(out) :: ln_Gamma(:) !! \(\ln \Gamma_i\)
      real(pr), optional, intent(out) :: dln_Gammadt(:) !! \(\ln \Gamma_i\)
      real(pr), optional, intent(out) :: dln_Gammadt2(:) !! \(\ln \Gamma_i\)
      real(pr), optional, intent(out) :: dln_Gammadx(:, :) !!
      real(pr), optional, intent(out) :: dln_Gammadx2(:, :, :) !!

      real(pr) :: theta(self%ngroups)
      !! Group area fractions
      real(pr) :: dthetadx(self%ngroups, size(x))
      !! First derivative of group area fractions with composition
      real(pr) :: dthetadx2(self%ngroups, size(x), size(x))
      !! Second derivative of group area fractions with composition
      real(pr) :: psi(self%ngroups,self%ngroups)
      !! \(\psi(T)\) parameter
      real(pr) :: dpsidt(self%ngroups, self%ngroups)
      real(pr) :: dpsidt2(self%ngroups,self%ngroups)

      ! Number of groups and components
      integer :: ng, nc

      ! Indexes used for groups
      integer :: n, m, k, gi

      ! Indexes for components
      integer :: i, j

      logical :: dx, dx2, dt, dt2

      ! Auxiliar variables to ease calculations
      real(pr) :: updown, updowndt, updowndx(size(x))
      real(pr) :: down(size(theta)), downdt(size(theta)), &
         downdt2(size(theta)), downdx(size(theta), size(x))

      ng = self%ngroups
      nc = size(x)
      dt = present(dln_Gammadt)
      dt2 = present(dln_Gammadt2)
      dx = present(dln_Gammadx)
      dx2 = present(dln_gammadx2)

      ! Initializate as zero
      ln_gamma = 0

      ! ========================================================================
      ! Calculate only the needed derivatives of the area fractions
      ! ------------------------------------------------------------------------
      if (dx2) then
         dln_Gammadx = 0
         dln_Gammadx2 = 0
         call group_area_fraction(&
            self, x, theta=theta, dthetadx=dthetadx, dthetadx2=dthetadx2 &
            )
      else if (dx) then
         dln_Gammadx = 0
         call group_area_fraction(self, x, theta=theta, dthetadx=dthetadx)
      else
         call group_area_fraction(self, x, theta=theta)
      end if

      ! ========================================================================
      ! Temperature dependance
      ! ------------------------------------------------------------------------
      if (dt2) then
         dln_Gammadt = 0
         dln_Gammadt2 = 0
         call self%psi_function%psi(self%groups_stew, T, psi, dpsidt, dpsidt2)
      else if(dt) then
         dln_Gammadt = 0
         call self%psi_function%psi(self%groups_stew, T, psi, dpsidt)
      else
         call self%psi_function%psi(self%groups_stew, T, psi)
      end if
      ! ========================================================================
      ! This parameters are used to avoid repeated calculations
      ! ------------------------------------------------------------------------
      down = 0
      downdx = 0
      do concurrent(m=1:ng, n=1:ng)
         down(m) = down(m) + theta(n) * psi(n, m)
         if (dt)  downdt(m) = downdt(m) + theta(n) * dpsidt(n, m)
         if (dt2) downdt2(m) = downdt(m) + theta(n) * dpsidt2(n, m)
         if(dx)   downdx(m, :) = downdx(m, :) + dthetadx(n, :) * psi(n, m)
      end do

      do k=1,ng
         if (theta(k) == 0) cycle

         ! Get the group index on the main matrix
         gi = self%groups_stew%groups_ids(k)
         updown = sum(theta(:) * psi(k, :)/down(:))

         ln_Gamma(k) = self%group_area(gi) * (&
            1 - log(sum(theta(:) * psi(:, k))) - updown &
            )

         if (dt) then
            temp_deriv: block
               real(pr) :: F(ng), Z(ng)
               do i=1,ng
                  F(i) = sum(theta(:) * dpsidt(:, i))
                  Z(i) = 1/sum(theta(:) * psi(:, i))
               end do

               dln_Gammadt(k) = self%group_area(gi) * (&
                  sum(&
                  Z(:) * (&
                  theta(:) * dpsidt(:, k) &
                  + theta(:) * psi(:, k) * F(:)* Z(:) &
                  ))  - F(k) * Z(k) &
                  )
            end block temp_deriv
         end if

         if (dx) then
            do concurrent(i=1:nc)
               dln_Gammadx(k, i) = self%group_area(gi) * ( &
                  - sum(dthetadx(:, i) * psi(:, k))/sum(theta(:) * psi(:, k)) &
                  - sum(dthetadx(:, i) * psi(k, :) / down(:))                 &
                  + sum(downdx(:, i) * theta(:) * psi(k, :) / down(:)**2)     &
                  )
            end do
            if (dx2) then
               do concurrent(j=1:nc)
               end do
            end if
         end if
      end do
   end subroutine group_big_gamma

   type(UNIFAC) function setup_unifac(molecules, Eij, Qk, Rk)
      !! UNIFAC model initialization.
      type(Groups), intent(in) :: molecules(:) !! Molecules
      real(pr), intent(in) :: Eij(:, :) !! Interaction Matrix
      real(pr), intent(in) :: Qk(:) !! Group k areas
      real(pr), intent(in) :: Rk(:) !! Group k volumes

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
end module yaeos__models_ge_group_contribution_unifac
