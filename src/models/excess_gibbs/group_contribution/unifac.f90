module yaeos__models_ge_group_contribution_unifac
   !! UNIFAC module
   use yaeos__constants, only: pr, R
   use yaeos__models_ge, only: GeModel
   use stdlib_io_npy, only: load_npy
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
      real(pr), allocatable :: Aij(:, :)
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

      real(pr) :: n_t

      ! Combinatorial
      real(pr) :: Ge_c
      real(pr) :: dGe_c_dn(self%nmolecules)
      real(pr) :: dGe_c_dn2(self%nmolecules, self%nmolecules)

      ! Residual
      real(pr) :: Ge_r
      real(pr) :: dGe_r_dn(self%nmolecules)
      real(pr) :: dGe_r_dn2(self%nmolecules, self%nmolecules)
      real(pr) :: dGe_r_dT, dGe_r_dT2, dGe_r_dTn(self%nmolecules)

      ! logical
      logical :: pge, dt, dt2, dtn, dn, dn2

      n_t = sum(n)

      !if (present(Ge) .and. .not.(GeT GeT2 Gen GeTn Gen2))
      call Ge_combinatorial(self, n, Ge_c, dGe_c_dn, dGe_c_dn2)
      call Ge_residual(self, n, T, Ge_r, dGe_r_dn, dGe_r_dn2, dGe_r_dT, dGe_r_dT2, dGe_r_dTn)

      if (present(Ge)) Ge = Ge_c * R * T + Ge_r
      if (present(Gen)) Gen = dGe_c_dn * R * T + dGe_r_dn
      if (present(Gen2)) Gen2 = (dGe_c_dn2 * R * T + dGe_r_dn2) * n_t
      if (present(GeT)) GeT = (dGe_r_dT + Ge_c * R)
      if (present(GeT2)) GeT2 = dGe_r_dT2
      if (present(GeTn)) GeTn = R * dGe_c_dn + dGe_r_dTn
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
      real(pr), optional, intent(out) :: dGe_dTn(self%nmolecules)

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
      real(pr) :: sum_vij_Qj_dlambdas_dT2(self%nmolecules)

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
      if ((dt .or. dtn) .and. .not. dt2) then
         call self%psi_function%psi(self%groups_stew, T, psi=Ejk, dpsidt=dEjk_dt)
      elseif (dt2 .and. .not. (dt .or. dtn)) then
         call self%psi_function%psi(self%groups_stew, T, psi=Ejk, dpsidt2=dEjk_dt2)
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
      if (dn .or. dt .or. dtn .or. dn2) then
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

      if (dt .or. dtn .or. dt2) then
         do k=1,self%ngroups
            sum_ni_vij_Qj_dEjk_dT(k) = sum(n * sum_vij_Qj_dEjk_dT(:,k))
            dlambda_k_dT(k) = sum(theta_j * dEjk_dt(:, k)) / sum(theta_j * Ejk(:, k))
            dlambda_k_dT2(k) = sum(n * sum_vij_Qj_dEjk_dT2(:,k)) / sum_ni_vij_Qj_Ejk(k) - dlambda_k_dT(k)**2
         end do
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
      ! Lambda_ik
      ! ------------------------------------------------------------------------
      if (pge .or. dn) then
         lambda_ik = 0.0_pr
         do concurrent (i=1:self%nmolecules, k=1:self%ngroups)
            if (self%vij(i,k) /= 0) then
               lambda_ik(i,k) = log(sum(self%thetas_ij(i, :) * Ejk(:, k)))
            end if
         end do
      end if

      ! Temperature derivatives
      if (dt .or. dt2) then
         dlambda_ik_dT = 0.0_pr
         do concurrent (i=1:self%nmolecules, k=1:self%ngroups)
            if (self%vij(i,k) /= 0) then
               dlambda_ik_dT(i,k) = sum(self%thetas_ij(i,:) * dEjk_dt(:, k)) / sum(self%thetas_ij(i,:) * Ejk(:, k))
            end if
         end do

         if (dt2) dlambda_ik_dT2 = sum_vij_Qj_dEjk_dT2 / sum_vij_Qj_Ejk - dlambda_ik_dT * dlambda_ik_dT
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
      if (dt .or. dt2) then
         do i=1,self%nmolecules
            sum_vij_Qj_dlambdas_dT(i) = sum(self%vij(i,:) * self%qk * (dlambda_k_dT - dlambda_ik_dT(i,:)))
         end do

         dGe_dT = -sum(n * sum_vij_Qj_dlambdas_dT)
      end if

      if (dt2) then
         do i=1,self%nmolecules
            sum_vij_Qj_dlambdas_dT2(i) = sum(self%vij(i,:) * self%qk * (dlambda_k_dT2 - dlambda_ik_dT2(i,:)))
         end do

         dGe_dT2 = -sum(n * sum_vij_Qj_dlambdas_dT2)
      end if

      if (dtn) then
         do i=1,self%nmolecules
            aux_sum(i) = sum(sum_nl_vlj * self%qk * dlambda_k_dndT(i,:))
         end do

         dGe_dTn = - sum_vij_Qj_dLambdas_dT - aux_sum
      end if

      ! ========================================================================
      ! From reduced Ge to Ge
      ! ------------------------------------------------------------------------
      if (present(dGe_dT2)) then
         dGe_dT2 = R * (2.0 * dGe_dT + T * dGe_dT2)
      end if

      if (present(dGe_dT)) then
         dGe_dT = R * (Ge + dGe_dT * T)
      end if

      if (present(dGe_dTn)) then
         dGe_dTn = R * (dGe_dn + dGe_dTn * T)
      end if

      Ge = Ge * R * T

      if (present(dGe_dn)) then
         dGe_dn = dGe_dn * R * T
      end if

      if (present(dGe_dn2)) then
         dGe_dn2 = dGe_dn2 * R * T
      end if
   end subroutine Ge_residual

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

      real(pr) :: Aij
      real(pr) :: Eij

      ngroups = size(systems_groups%groups_ids)

      do concurrent(i=1:ngroups, j=1:ngroups)
         ig = systems_groups%groups_ids(i)
         jg = systems_groups%groups_ids(j)
         
         Aij = self%Aij(ig, jg)
         Eij = exp(- Aij / T)

         if (present(psi)) &
            psi(i, j) = Eij
         if (present(dpsidt)) &
            dpsidt(i, j) = Aij * Eij / T**2
         if (present(dpsidt2)) &
            dpsidt2(i, j) = Aij * (Aij - 2_pr*T) * Eij / T**4
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

   subroutine get_unifac_default_parameters(Aij, Qk, Rk)
      real(pr), allocatable, optional, intent(out) :: Aij(:,:)
      real(pr), allocatable, optional, intent(out) :: Qk(:)
      real(pr), allocatable, optional, intent(out) :: Rk(:)

      if (present(Aij)) then
         call load_npy("data/unifac_aij.npy", Aij)
      end if

      if (present(Qk)) then
         call load_npy("data/unifac_Qk.npy", Qk)
      end if

      if (present(Rk)) then
         call load_npy("data/unifac_Rk.npy", Rk)
      end if
   end subroutine get_unifac_default_parameters

   type(UNIFAC) function setup_unifac(molecules, Aij, Qk, Rk)
      !! UNIFAC model initialization.
      type(Groups), intent(in) :: molecules(:) !! Molecules
      real(pr), optional, intent(in) :: Aij(:, :) !! Interaction Matrix
      real(pr), optional, intent(in) :: Qk(:) !! Group k areas
      real(pr), optional, intent(in) :: Rk(:) !! Group k volumes

      type(Groups) :: soup
      type(UNIFACPsi) :: psi_function

      ! Default UNIFAC parameters
      real(pr), allocatable:: Aij_default(:, :)
      real(pr), allocatable:: Qk_default(:)
      real(pr), allocatable:: Rk_default(:)

      ! Parameters to perform the setup
      real(pr), allocatable:: Aij_final(:, :)
      real(pr), allocatable:: Qk_final(:)
      real(pr), allocatable:: Rk_final(:)

      !
      integer, allocatable :: vij(:, :)
      real(pr), allocatable :: qks(:)

      integer :: gi, i, j, k

      setup_unifac%molecules = molecules

      allocate(soup%groups_ids(0))
      allocate(soup%number_of_groups(0))

      ! ========================================================================
      ! Load default UNIFAC parameters if not provided
      ! ------------------------------------------------------------------------
      ! Aij
      if (.not. present(Aij)) then
         call get_unifac_default_parameters(Aij=Aij_default)
         
         allocate(Aij_final(size(Aij_default, 1), size(Aij_default, 2)))
         Aij_final(:,:) = Aij_default(:,:)
      else
         allocate(Aij_final(size(Aij, 1), size(Aij, 2)))
         Aij_final(:,:) = Aij(:,:) 
      end if

      ! Qk
      if (.not. present(Qk)) then
         call get_unifac_default_parameters(Qk=Qk_default)
         
         allocate(Qk_final(size(Qk_default)))
         Qk_final(:) = Qk_default(:)
      else
         allocate(Qk_final(size(Qk)))
         Qk_final(:) = Qk(:) 
      end if

      ! Rk
      if (.not. present(Rk)) then
         call get_unifac_default_parameters(Rk=Rk_default)
         
         allocate(Rk_final(size(Rk_default)))
         Rk_final(:) = Rk_default(:)
      else
         allocate(Rk_final(size(Rk)))
         Rk_final(:) = Rk(:) 
      end if

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
               r(i) = r(i) + molecules(i)%number_of_groups(j) * Rk_final(gi)
               q(i) = q(i) + molecules(i)%number_of_groups(j) * Qk_final(gi)

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
      allocate(qks(size(soup%number_of_groups)))

      vij = 0
      qks = 0.0_pr
      do i=1,size(molecules)
         do k=1,size(molecules(i)%number_of_groups)
            gi = molecules(i)%groups_ids(k)

            ! Index of group for Area
            j = findloc(soup%groups_ids, gi, dim=1)

            vij(i,j) = molecules(i)%number_of_groups(k)

            qks(j) = Qk_final(gi)
         end do
      end do
      ! ========================================================================

      psi_function%Aij = Aij_final
      setup_unifac%groups_stew = soup
      setup_unifac%ngroups = size(soup%number_of_groups)
      setup_unifac%nmolecules = size(molecules)
      setup_unifac%psi_function = psi_function
      setup_unifac%group_area = Qk_final
      setup_unifac%group_volume = Rk_final
      setup_unifac%thetas_ij = thetas_i(&
         size(molecules), size(soup%number_of_groups), Qk_final, soup, molecules)
      setup_unifac%vij = vij
      setup_unifac%qk = qks
   end function setup_unifac
end module yaeos__models_ge_group_contribution_unifac
