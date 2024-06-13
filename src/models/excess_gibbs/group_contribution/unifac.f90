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
      integer :: nmolecules
      !! Total number of molecules
      real(pr) :: z = 10
      !! Model constant
      real(pr), allocatable :: group_area(:)
      !! Group areas \(Q_k\)
      real(pr), allocatable :: group_volume(:)
      !! Group volumes \(R_k\)
      real(pr), allocatable :: thetas_ij(:, :)
      !! Area fractions of groups j (row) on molecules i (column)
      real(pr), allocatable :: vij(:,:)
      !! Ocurrences of each group j on each molecule i
      real(pr), allocatable :: qk(:)
      !! Area of each group k
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
         real(pr), optional, intent(out) :: psi(:, :)
         real(pr), optional, intent(out) :: dpsidt(:, :)
         real(pr), optional, intent(out) :: dpsidt2(:, :)
      end subroutine temperature_dependence
   end interface

contains

   subroutine excess_gibbs(self, n, t, Ge, GeT, GeT2, Gen, GeTn, Gen2)
      !! Excess Gibbs and derivs procedure
      class(UNIFAC), intent(in) :: self !! Model
      real(pr), intent(in) :: n(:) !! Moles vector
      real(pr), intent(in) :: t !! Temperature [K]
      real(pr), optional, intent(out) :: Ge !! Excess Gibbs
      real(pr), optional, intent(out) :: GeT !! \(\frac{dG^E}{dT}\)
      real(pr), optional, intent(out) :: GeT2 !! \(\frac{d^2G^E}{dT^2}\)
      real(pr), optional, intent(out) :: Gen(size(n))
      real(pr), optional, intent(out) :: GeTn(size(n))
      real(pr), optional, intent(out) :: Gen2(size(n), size(n))

      real(pr) :: x(size(n))
      real(pr) :: ln_gamma_c(size(n)), ln_gamma_r(size(n)), ln_activity(size(n))
      real(pr) :: nt, dxidni(size(n), size(n))
      real(pr) :: Ge_c, Ge_r
      real(pr) :: dGe_c_dn(size(n)), dGe_r_dn(size(n))

      integer :: i, j, nc

      nt = sum(n)

      x = n/nt

      nc = self%nmolecules

      ln_activity = 0

      !call combinatorial_activity(self, x, ln_gamma_c)
      !call residual_activity(self, x, T, ln_gamma_r)
      !ln_activity = ln_gamma_c + ln_gamma_r

      if (present(Ge)) then
         call Ge_combinatorial(self, n, Ge_c)
         call Ge_residual(self, n, T, Ge_r)
         Ge =  (Ge_c + Ge_r) * (R*T)
      end if

      if (present(GeN)) then
         call Ge_combinatorial(self, n, Ge_c, dGe_c_dn)
         call Ge_residual(self, n, T, Ge_r, dGe_dn=dGe_r_dn)
         Gen = (dGe_c_dn + dGe_r_dn)
      end if
   end subroutine excess_gibbs

   subroutine Ge_combinatorial(self, n, Ge, dGe_dn, dGe_dn2)
      class(UNIFAC) :: self

      real(pr), intent(in) :: n(self%nmolecules)
      real(pr), optional, intent(out) :: Ge
      real(pr), optional, intent(out) :: dGe_dn(self%nmolecules)
      real(pr), optional, intent(out) :: dGe_dn2(self%nmolecules,self%nmolecules)

      ! Flory-Huggins variables
      real(pr) :: Ge_fh
      real(pr) :: dGe_fh_dn(self%nmolecules)
      real(pr) :: dGe_fh_dn2(self%nmolecules,self%nmolecules)

      ! Staverman-Guggenheim variables
      real(pr) :: Ge_sg
      real(pr) :: dGe_sg_dn(self%nmolecules)
      real(pr) :: dGe_sg_dn2(self%nmolecules,self%nmolecules)

      ! utility
      real(pr) :: nq, nr, n_t
      integer :: i, j

      associate(&
         q => self%molecules%surface_area,&
         r => self%molecules%volume,&
         z => self%z &
         )

         nr = dot_product(n, r)
         nq = dot_product(n, q)
         n_t = sum(n)

         if (present(Ge)) then
            Ge_fh = dot_product(n, log(r)) - n_t * log(nr) + n_t * log(n_t)
            Ge_sg = z/2 * sum(n * q * (log(q/r) - log(nq) + log(nr)))
         end if

         if (present(dGe_dn)) then
            dGe_fh_dn = log(r) - log(nr) + log(n_t) + 1 - n_t * r / nr
            dGe_sg_dn = z/2*q*(-log((r*nq)/(q*nr)) - 1 + (r*nq)/(q*nr))
         end if

         if (present(dGe_dn2)) then
            dGe_fh_dn2 = 0.0_pr
            dGe_sg_dn2 = 0.0_pr
            do concurrent(i=1:size(n), j=1:size(n))
               dGe_fh_dn2(i,j) = -(r(i) + r(j))/nr + 1/n_t + n_t*r(i)*r(j)/ nr**2
               dGe_sg_dn2(i,j) = z/2*(-q(i)*q(j)/nq + (q(i)*r(j) + q(j)*r(i))/nr - r(i)*r(j)*nq/nr**2)
            end do
         end if
      end associate

      if (present(Ge)) Ge = Ge_fh + Ge_sg
      if (present(dGe_dn)) dGe_dn = dGe_fh_dn + dGe_sg_dn
      if (present(dGe_dn2)) dGe_dn2 = dGe_fh_dn2 + dGe_sg_dn2
   end subroutine Ge_combinatorial

   subroutine Ge_residual(self, n, T, Ge, dGe_dn, dGe_dn2, dGe_dT, dGe_dT2, dGe_dTn)
      class(UNIFAC) :: self

      real(pr), intent(in) :: n(self%nmolecules)
      real(pr), intent(in) :: T
      real(pr), optional, intent(out) :: Ge
      real(pr), optional, intent(out) :: dGe_dn(self%nmolecules)
      real(pr), optional, intent(out) :: dGe_dn2(self%nmolecules, self%nmolecules)
      real(pr), optional, intent(out) :: dGe_dT
      real(pr), optional, intent(out) :: dGe_dT2
      real(pr), optional, intent(out) :: dGe_dTn(self%nmolecules, self%nmolecules)

      ! Thetas variables
      real(pr) :: theta_j(self%ngroups)

      ! Ejk variables
      real(pr) :: Ejk(self%ngroups, self%ngroups)
      real(pr) :: dEjk_dt(self%ngroups, self%ngroups)
      real(pr) :: dEjk_dt2(self%ngroups, self%ngroups)

      ! Lambdas variables
      real(pr) :: lambda_k(self%ngroups)
      real(pr) :: dlambda_k_dT(self%ngroups)
      real(pr) :: dlambda_k_dT2(self%ngroups)
      real(pr) :: dlambda_k_dn(self%nmolecules, self%ngroups)
      real(pr) :: dlambda_k_dn2(self%nmolecules, self%nmolecules, self%ngroups)
      real(pr) :: dlambda_k_dndT(self%nmolecules, self%ngroups)

      real(pr) :: lambda_ik(self%nmolecules, self%ngroups)
      real(pr) :: dlambda_ik_dT(self%nmolecules, self%ngroups)
      real(pr) :: dlambda_ik_dT2(self%nmolecules, self%ngroups)

      ! Auxiliars
      real(pr) :: sum_vij_Qj_Ejk(self%nmolecules, self%ngroups)
      real(pr) :: sum_ni_vij_Qj_Ejk(self%ngroups)
      real(pr) :: sum_vik_Qk(self%nmolecules)
      real(pr) :: sum_vQ_Lambda(self%nmolecules)
      real(pr) :: sum_nl_vlj(self%ngroups)
      real(pr) :: sum_ni_vik_Qk
      real(pr) :: aux_sum(self%nmolecules)
      real(pr) :: sum_Q_v_dlambda_k_dn(self%nmolecules, self%nmolecules)
      real(pr) :: aux_sum2
      real(pr) :: sum_vij_Qj_dEjk_dT(self%nmolecules, self%ngroups)
      real(pr) :: sum_vij_Qj_dEjk_dT2(self%nmolecules, self%ngroups)
      real(pr) :: sum_ni_vij_Qj_dEjk_dT(self%ngroups)
      real(pr) :: sum_vij_Qj_dlambdas_dT(self%nmolecules)

      ! Indexes used for groups
      integer :: j, k

      ! Indexes used for components
      integer :: i, l

      ! logicals
      logical :: pge, dn, dn2, dt, dt2, dtn

      pge = present(Ge)
      dn = present(dGe_dn)
      dn2 = present(dGe_dn2)
      dt = present(dGe_dT)
      dt2 = present(dGe_dT2)
      dtn = present(dGe_dTn)

      ! ========================================================================
      ! Ejk
      ! ------------------------------------------------------------------------
      if (pge .and. .not. (dt .or. dt2 .or. dtn)) then
         call self%psi_function%psi(self%groups_stew, T, psi=Ejk)
      elseif ((dt .or. dtn) .and. .not. (pge .or. dt2)) then
         call self%psi_function%psi(self%groups_stew, T, dpsidt=dEjk_dt)
      elseif (dt2 .and. .not. (pge .or. dt .or. dtn)) then
         call self%psi_function%psi(self%groups_stew, T, dpsidt2=dEjk_dt2)
      elseif ((pge .and. (dt .or. dtn)) .and. .not. dt2) then
         call self%psi_function%psi(self%groups_stew, T, psi=Ejk, dpsidt=dEjk_dt)
      elseif ((pge .and. dt2) .and. .not. (dt .or. dtn)) then
         call self%psi_function%psi(self%groups_stew, T, psi=Ejk, dpsidt2=dEjk_dt2)
      elseif (((dt .or. dtn) .and. dt2) .and. .not. pge) then
         call self%psi_function%psi(self%groups_stew, T, dpsidt=dEjk_dt, dpsidt2=dEjk_dt2)
      else
         call self%psi_function%psi(self%groups_stew, T, psi=Ejk, dpsidt=dEjk_dt, dpsidt2=dEjk_dt2)
      end if

      ! ========================================================================
      ! Auxiliars
      ! ------------------------------------------------------------------------
      do i=1,self%nmolecules
         sum_vik_Qk(i) = sum(self%vij(i,:) * self%qk)
      end do
      sum_ni_vik_Qk = sum(n * sum_vik_Qk)

      if (dtn .or. dt2 .or. dt) then
         do concurrent(i=1:self%nmolecules, k=1:self%ngroups)
            sum_vij_Qj_dEjk_dT(i,k) = sum(self%vij(i,:) * self%qk * dEjk_dT(:,k))
            sum_vij_Qj_dEjk_dT2(i,k) = sum(self%vij(i,:) * self%qk * dEjk_dT2(:,k))
         end do
      end if

      ! ========================================================================
      ! Thetas
      ! ------------------------------------------------------------------------
      do j=1,self%ngroups
         sum_nl_vlj(j) = sum(n * self%vij(:,j))
         theta_j(j) = sum_nl_vlj(j) * self%qk(j) / sum_ni_vik_Qk
      end do

      ! ========================================================================
      ! Lambda_k
      ! ------------------------------------------------------------------------
      ! Lambda_k
      if (pge .or. dn) then
         do k=1,self%ngroups
            lambda_k(k) = log(sum(theta_j * Ejk(:,k)))
         end do
      end if

      ! Lambda_k first compositional derivatives
      if (dn .or. dtn .or. dn2) then
         do concurrent (i=1:self%nmolecules, k=1:self%ngroups)
            sum_vij_Qj_Ejk(i,k) = sum(self%vij(i,:) * self%qk * Ejk(:,k))
         end do

         do k=1,self%ngroups
            sum_ni_vij_Qj_Ejk(k) = sum(n * sum_vij_Qj_Ejk(:,k))
         end do

         do i=1,self%nmolecules
            dlambda_k_dn(i,:) = sum_vij_Qj_Ejk(i,:) / sum_ni_vij_Qj_Ejk - sum_vik_Qk(i) / sum_ni_vik_Qk
         end do
      end if

      ! Lambda_k second compositional derivatives
      if (dn2) then
         do concurrent (i=1:self%nmolecules,l=1:self%nmolecules)
            sum_Q_v_dlambda_k_dn(i,l) = sum(self%qk * self%vij(l,:) * dlambda_k_dn(i,:))
            dlambda_k_dn2(i,l,:) = (&
               - sum_vij_Qj_Ejk(i,:) * sum_vij_Qj_Ejk(l,:) / sum_ni_vij_Qj_Ejk**2 &
               + sum_vik_Qk(i) * sum_vik_Qk(l) / sum_ni_vik_Qk**2 &
               )
         end do
      end if

      ! Temperature derivatives
      if (dt .or. dtn .or. dt2) then
         do k=1,self%nmolecules
            sum_ni_vij_Qj_dEjk_dT(k) = sum(n * sum_vij_Qj_dEjk_dT(:,k))
            dlambda_k_dT(k) = sum_ni_vij_Qj_dEjk_dT(k) / sum_ni_vij_Qj_Ejk(k)
            dlambda_k_dT2(k) = sum(n * sum_vij_Qj_dEjk_dT2(:,k)) / sum_ni_vij_Qj_Ejk(k) - dlambda_k_dT(k)**2
         end do
      end if

      ! ========================================================================
      ! Lambda_ik
      ! ------------------------------------------------------------------------
      do concurrent (k=1:self%ngroups, i=1:self%nmolecules)
         lambda_ik(i, k) = sum(self%thetas_ij(i, :) * Ejk(:, k))
      end do
      lambda_ik = log(lambda_ik)

      ! Temperature derivatives
      if (dt .or. dt2) then
         dlambda_ik_dT = sum_vij_Qj_dEjk_dT / sum_vij_Qj_Ejk
         if (dt2) dlambda_ik_dT2 = sum_vij_Qj_dEjk_dT2 / sum_vij_Qj_Ejk - dlambda_ik_dT * dlambda_ik_dT
      end if

      if (dtn) then
         do i=1,self%nmolecules
            dlambda_k_dndT(i,:) = (&
               sum_vij_Qj_dEjk_dT(i,:) / sum_ni_vij_Qj_Ejk &
               - sum_vij_Qj_Ejk(i,:) * sum_ni_vij_Qj_dEjk_dT / sum_ni_vij_Qj_Ejk**2 &
               )
         end do
      end if

      ! ========================================================================
      ! Ge
      ! ------------------------------------------------------------------------
      if (pge .or. dn) then
         do i=1,self%nmolecules
            sum_vQ_Lambda(i) = sum(self%vij(i,:) * self%qk * (lambda_k - lambda_ik(i,:)))
         end do
      end if

      if (pge) Ge = - sum(n * sum_vQ_Lambda)

      ! ========================================================================
      ! dGe_dn
      ! ------------------------------------------------------------------------
      if (dn) then
         do i=1,self%nmolecules
            aux_sum(i) = sum(sum_nl_vlj * self%qk * dlambda_k_dn(i,:))
         end do
         dGe_dn = -sum_vQ_Lambda - aux_sum
      end if

      ! ========================================================================
      ! dGe_dn2
      ! ------------------------------------------------------------------------
      if (dn2) then
         do concurrent (i=1:self%nmolecules,l=1:self%nmolecules)
            aux_sum2 = sum(sum_nl_vlj * dlambda_k_dn2(i,l,:) * self%qk)
            dGe_dn2(i,l) = -(sum_Q_v_dlambda_k_dn(i,l) + sum_Q_v_dlambda_k_dn(l,i)) - aux_sum2
         end do
      end if

      ! ========================================================================
      ! dGe_dT, dGe_dT2, dGE_dnT
      ! ------------------------------------------------------------------------
      do i=1,self%nmolecules
         sum_vij_Qj_dlambdas_dT(i) = sum(self%vij(i,:) * self%qk * (dlambda_k_dT - dlambda_ik_dT(i,:)))
      end do

      if (dt) then
         dGe_dT = -sum(n * sum_vij_Qj_dlambdas_dT)
      end if

   end subroutine Ge_residual

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
   end subroutine ln_activity_coefficient

   subroutine UNIFAC_temperature_dependence(self, systems_groups, T, psi, dpsidt, dpsidt2)
      class(UNIFACPsi) :: self !! \(\psi\) function
      class(Groups) :: systems_groups !! Groups in the system
      real(pr), intent(in) :: T !! Temperature
      real(pr), optional, intent(out) :: psi(:, :) !! \(\psi\)
      real(pr), optional, intent(out) :: dpsidt(:, :)
      real(pr), optional, intent(out) :: dpsidt2(:, :)

      integer :: i, j
      integer :: ig, jg
      integer :: ngroups

      ngroups = size(systems_groups%groups_ids)

      do concurrent(i=1:ngroups, j=1:ngroups)
         ig = systems_groups%groups_ids(i)
         jg = systems_groups%groups_ids(j)
         if (present(psi)) &
            psi(i, j) = exp(-self%Eij(ig, jg) / T)
         if (present(dpsidt)) &
            dpsidt(i, j) = self%Eij(ig, jg) * psi(i, j) / T**2
         if (present(dpsidt2)) &
            dpsidt2(i, j) = &
            self%Eij(ig, jg) * (self%Eij(ig, jg) - 2*T) * psi(i, j) / T**4
      end do

   end subroutine UNIFAC_temperature_dependence

   function thetas_i(nm, ng, group_area, stew, molecules) result(thetas_ij)
      integer, intent(in) :: nm !! Number of molecules
      integer, intent(in) :: ng !! Number of groups
      real(pr), intent(in) :: group_area(:) !! Group k areas
      type(Groups), intent(in) :: stew !! All the groups present in the system
      type(Groups), intent(in) :: molecules(:) !! Molecules
      real(pr) :: thetas_ij(nm, ng) !! Group j area fraction on molecule i

      real(pr) :: total_area_i(nm)
      real(pr) :: qki_contribution

      integer :: gi !! group k id
      integer :: i, j, k

      thetas_ij = 0.0_pr
      total_area_i = 0.0_pr

      ! Obtain the total area of each molecule
      do i=1,size(molecules)
         do k=1,size(molecules(i)%number_of_groups)
            gi = molecules(i)%groups_ids(k)

            ! Locate group k in the stew ordering (position j of group k).
            j = findloc(stew%groups_ids, gi, dim=1)

            ! Contribution of the group k to the molecule i area.
            qki_contribution = (&
               group_area(gi) * molecules(i)%number_of_groups(k) &
               )

            ! Adding to the total area of each molecule
            total_area_i(i) = total_area_i(i) + qki_contribution
         end do
      end do

      ! Calculate the fraction of each group on each molecule
      thetas_ij = 0.0_pr

      do i=1,size(molecules)
         do k=1,size(molecules(i)%number_of_groups)
            gi = molecules(i)%groups_ids(k)

            j = findloc(stew%groups_ids, gi, dim=1)

            thetas_ij(i, j) = group_area(gi) * molecules(i)%number_of_groups(k) / total_area_i(i)
         end do
      end do
   end function thetas_i

   type(UNIFAC) function setup_unifac(molecules, Eij, Qk, Rk)
      !! UNIFAC model initialization.
      type(Groups), intent(in) :: molecules(:) !! Molecules
      real(pr), intent(in) :: Eij(:, :) !! Interaction Matrix
      real(pr), intent(in) :: Qk(:) !! Group k areas
      real(pr), intent(in) :: Rk(:) !! Group k volumes

      type(Groups) :: soup
      type(UNIFACPsi) :: psi_function

      integer, allocatable :: vij(:, :)
      real(pr), allocatable :: qij(:,:)
      real(pr), allocatable :: qks(:)

      integer :: gi, i, j, k

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
      ! Build a matrix vij and vector qk
      ! ------------------------------------------------------------------------
      allocate(vij(size(molecules), size(soup%number_of_groups)))
      allocate(qij(size(molecules), size(soup%number_of_groups)))
      allocate(qks(size(soup%number_of_groups)))

      vij = 0
      qij = 0.0_pr
      do i=1,size(molecules)
         do k=1,size(molecules(i)%number_of_groups)
            gi = molecules(i)%groups_ids(k)

            ! Index of group for Area
            j = findloc(soup%groups_ids, gi, dim=1)

            vij(i,j) = molecules(i)%number_of_groups(k)
            qij(i,j) = Qk(gi)

            if (Qk(gi) /= 0) then
               qks(j) = Qk(gi)
            end if
         end do
      end do
      ! ========================================================================

      psi_function%Eij = Eij
      setup_unifac%groups_stew = soup
      setup_unifac%ngroups = size(soup%number_of_groups)
      setup_unifac%nmolecules = size(molecules)
      setup_unifac%psi_function = psi_function
      setup_unifac%group_area = Qk
      setup_unifac%group_volume = Rk
      setup_unifac%thetas_ij = thetas_i(&
         size(molecules), size(soup%number_of_groups), Qk, soup, molecules)
      setup_unifac%vij = vij
      setup_unifac%qk = qks
   end function setup_unifac
end module yaeos__models_ge_group_contribution_unifac
