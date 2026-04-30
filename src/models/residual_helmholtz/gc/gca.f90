! module yaeos__models_ar_gc_gca
!    !! Group Contribution with Association Equation of State (GCA-EoS)
!    !! Refactored to be compatible with the yaeos framework
!    use yaeos__constants, only: pr, R
!    use yaeos__models_ar, only: ArModel
!    implicit none
! 
!    private
!    public :: GCAEOS
! 
!    ! Parameters
!    integer, parameter :: NCM = 30, NGM = 30, NSM = 24, NGAM = 14
!    real(pr), parameter :: RGC = 82.05_pr  ! Gas constant in atm·cm³/(mol·K)
!    real(pr), parameter :: PI = 3.14159265358979323846_pr
! 
!    type, extends(ArModel) :: GCAEOS
!       !! GCA-EoS model compatible with yaeos framework
!       integer :: nc  !! Number of components
!       integer :: ng  !! Number of dispersive groups
!       integer :: nga !! Number of associating groups
!       integer :: nst !! Total number of association sites
! 
!       ! Component properties
!       real(pr), allocatable :: tc(:)     !! Critical temperature [K]
!       real(pr), allocatable :: pc(:)     !! Critical pressure [atm]
!       real(pr), allocatable :: omega(:)  !! Acentric factor
!       real(pr), allocatable :: dc(:)     !! Critical diameter [cm/mol^(1/3)]
!       real(pr), allocatable :: d(:)      !! Temperature-dependent diameter
!       real(pr), allocatable :: dt(:)     !! Temperature derivative of diameter
! 
!       ! Group properties - Dispersive
!       integer, allocatable :: ny(:,:)           !! Group occurrence matrix (nc, ng)
!       real(pr), allocatable :: q(:)             !! Surface area parameter
!       real(pr), allocatable :: gstr(:)          !! Energy parameter g*
!       real(pr), allocatable :: g1(:)            !! Temperature coefficient g'
!       real(pr), allocatable :: g2(:)            !! Temperature coefficient g''
!       real(pr), allocatable :: tstr(:)          !! Reference temperature T*
!       real(pr), allocatable :: tspl(:)          !! Spline temperature
!       real(pr), allocatable :: epx(:)           !! Exponent for high-T region
!       real(pr), allocatable :: g(:,:)           !! Group interaction energies
!       real(pr), allocatable :: dgdt(:,:)        !! Temperature derivative of g
!       real(pr), allocatable :: kstr(:,:)        !! Binary interaction parameter k*
!       real(pr), allocatable :: kp(:,:)          !! Temperature dependence k'
!       real(pr), allocatable :: alpha(:,:)       !! Damping factor alpha
! 
!       ! Association properties
!       integer, allocatable :: sigma(:,:)        !! Site occurrence matrix (nst, nc)
!       real(pr), allocatable :: delta(:,:)       !! Association strength
!       real(pr), allocatable :: ddelt(:,:)       !! Temperature derivative of Delta
!       real(pr), allocatable :: eps_r(:,:)       !! Association energy/R [K]
!       real(pr), allocatable :: kappa(:,:)       !! Association volume [cm³/mol]
! 
!    contains
!       procedure :: residual_helmholtz => gca_residual_helmholtz
!       procedure :: get_v0 => gca_get_b
!       procedure :: update_diameters
!       procedure :: update_group_energies
!       procedure :: update_association_parameters
!    end type GCAEOS
! 
!    interface GCAEOS
!       module procedure :: init_gcaeos
!    end interface GCAEOS
! 
! contains
! 
!    type(GCAEOS) function init_gcaeos(nc, ng, nga, nst) result(eos)
!       !! Initialize GCA-EoS model
!       integer, intent(in) :: nc   !! Number of components
!       integer, intent(in) :: ng   !! Number of dispersive groups
!       integer, intent(in) :: nga  !! Number of associating groups (0 if no association)
!       integer, intent(in) :: nst  !! Number of association sites (0 if no association)
! 
!       eos%nc = nc
!       eos%ng = ng
!       eos%nga = nga
!       eos%nst = nst
! 
!       ! Allocate component properties
!       allocate(eos%tc(nc), eos%pc(nc), eos%omega(nc))
!       allocate(eos%dc(nc), eos%d(nc), eos%dt(nc))
! 
!       ! Allocate group properties
!       allocate(eos%ny(nc, ng))
!       allocate(eos%q(ng), eos%gstr(ng), eos%g1(ng), eos%g2(ng))
!       allocate(eos%tstr(ng), eos%tspl(ng), eos%epx(ng))
!       allocate(eos%g(ng, ng), eos%dgdt(ng, ng))
!       allocate(eos%kstr(ng, ng), eos%kp(ng, ng), eos%alpha(ng, ng))
! 
!       ! Allocate association properties if needed
!       if (nst > 0) then
!          allocate(eos%sigma(nst, nc))
!          allocate(eos%delta(nst, nst), eos%ddelt(nst, nst))
!          allocate(eos%eps_r(nst, nst), eos%kappa(nst, nst))
! 
!          eos%sigma = 0
!          eos%delta = 0.0_pr
!          eos%ddelt = 0.0_pr
!          eos%eps_r = 0.0_pr
!          eos%kappa = 0.0_pr
!       end if
! 
!       ! Initialize matrices
!       eos%ny = 0
!       eos%g = 0.0_pr
!       eos%dgdt = 0.0_pr
!       eos%kstr = 1.0_pr  ! Default to 1 (no interaction)
!       eos%kp = 0.0_pr
!       eos%alpha = 0.0_pr
! 
!    end function init_gcaeos
! 
!    subroutine update_diameters(self, T, d, dt)
!       !! Update temperature-dependent soft-sphere diameters
!       !! d(T) = dc * (1 - 0.12 * exp(-2*Tc/(3*T))) * 1.065655
!       class(GCAEOS), intent(in) :: self
!       real(pr), intent(in) :: T  !! Temperature [K]
!       real(pr), intent(out) :: d(:)   !! Diameters [cm/mol^(1/3)]
!       real(pr), intent(out) :: dt(:)  !! Temperature derivatives [cm/mol^(1/3)/K]
! 
!       integer :: i
!       real(pr) :: arg, darg, darg_dt
! 
!       do i = 1, self%nc
!          arg = -2.0_pr * self%tc(i) / (3.0_pr * T)
!          darg = exp(arg)
! 
!          ! Derivada del argumento: d(-2Tc/3T)/dT = 2Tc/(3T^2)
!          ! Que es lo mismo que: -arg / T
!          darg_dt = darg * (2.0_pr * self%tc(i) / (3.0_pr * T**2))
! 
!          ! dt = dc * (-0.12 * darg_dt) * 1.065655
!          dt(i) = -self%dc(i) * 0.12_pr * darg_dt * 1.065655_pr
!          d(i) = self%dc(i) * (1.0_pr - 0.12_pr * darg) * 1.065655_pr
!       end do
! 
!    end subroutine update_diameters
! 
!    subroutine update_group_energies(self, T, g, dgdt)
!       !! Update temperature-dependent group interaction energies
!       class(GCAEOS), intent(in) :: self
!       real(pr), intent(in) :: T  !! Temperature [K]
!       real(pr), intent(out) :: g(:,:)    !! Group energies [atm·cm⁶/mol²]
!       real(pr), intent(out) :: dgdt(:,:) !! Temperature derivatives
! 
!       integer :: i, j
!       real(pr) :: tr, help, k_ij, dkij
! 
!       g = 0.0_pr
!       dgdt = 0.0_pr
! 
!       ! Update pure group energies
!       do i = 1, self%ng
!          if (T < self%tspl(i)) then
!             ! Low temperature region
!             tr = T / self%tstr(i)
!             dgdt(i, i) = self%gstr(i) * (self%g1(i) / self%tstr(i) + &
!                self%g2(i) / T)
!             g(i, i) = self%gstr(i) * (1.0_pr + self%g1(i) * (tr - 1.0_pr) + &
!                self%g2(i) * log(tr))
!          else
!             ! High temperature region
!             g(i, i) = self%gstr(i) / 4.0_pr * (self%tspl(i) / T)**self%epx(i)
!             dgdt(i, i) = -self%epx(i) * g(i, i) / T
!          end if
!       end do
! 
!       ! Update cross-interaction energies
!       if (self%ng > 1) then
!          do i = 1, self%ng - 1
!             do j = i + 1, self%ng
!                tr = 2.0_pr * T / (self%tstr(i) + self%tstr(j))
!                k_ij = self%kstr(i, j) * (1.0_pr + self%kp(i, j) * log(tr))
!                help = sqrt(g(i, i) * g(j, j))
!                g(i, j) = help * k_ij
!                g(j, i) = g(i, j)
! 
!                dkij = self%kstr(i, j) * self%kp(i, j) / T
!                dgdt(i, j) = dkij * help + k_ij / (2.0_pr * help) * &
!                   (g(i, i) * dgdt(j, j) + g(j, j) * dgdt(i, i))
!                dgdt(j, i) = dgdt(i, j)
!             end do
!          end do
!       end if
! 
!    end subroutine update_group_energies
! 
!    subroutine update_association_parameters(self, T, delta, ddelt)
!       !! Update temperature-dependent association parameters
!       !! Δ(T) = κ * (exp(ε/T) - 1)
!       class(GCAEOS), intent(in) :: self
!       real(pr), intent(in) :: T  !! Temperature [K]
!       real(pr), intent(out) :: delta(:,:)  !! Association strength [cm³/mol]
!       real(pr), intent(out) :: ddelt(:,:)  !! Temperature derivatives
! 
!       integer :: i, j
! 
!       delta = 0.0_pr
!       ddelt = 0.0_pr
! 
!       if (self%nst > 0) then
!          do i = 1, self%nst
!             do j = i, self%nst
!                delta(i, j) = self%kappa(i, j) * &
!                   (exp(self%eps_r(i, j) / T) - 1.0_pr)
!                delta(j, i) = delta(i, j)
! 
!                ddelt(i, j) = -(delta(i, j) + self%kappa(i, j)) * &
!                   self%eps_r(i, j) / T**2
!                ddelt(j, i) = ddelt(i, j)
!             end do
!          end do
!       end if
! 
!    end subroutine update_association_parameters
! 
!    subroutine update_temperature_parameters(self, T)
!       !! Update all temperature-dependent parameters
!       !! This is a convenience wrapper - updates internal arrays
!       !! Note: Cannot be used inside residual_helmholtz due to intent(in)
!       class(GCAEOS), intent(inout) :: self
!       real(pr), intent(in) :: T  !! Temperature [K]
! 
!       call self%update_diameters(T, self%d, self%dt)
!       call self%update_group_energies(T, self%g, self%dgdt)
! 
!       if (self%nst > 0) then
!          call self%update_association_parameters(T, self%delta, self%ddelt)
!       end if
! 
!    end subroutine update_temperature_parameters
! 
!    subroutine gca_residual_helmholtz(self, n, v, t, Ar, ArV, ArT, ArTV, &
!       ArV2, ArT2, Arn, ArVn, ArTn, Arn2)
!       !! Calculate residual Helmholtz energy and derivatives
!       !! Following yaeos ArModel interface
!       class(GCAEOS), intent(in) :: self
!       real(pr), intent(in) :: n(:)   !! Moles vector [mol]
!       real(pr), intent(in) :: v      !! Volume [L]
!       real(pr), intent(in) :: t      !! Temperature [K]
!       real(pr), optional, intent(out) :: Ar      !! Residual Helmholtz [bar·L]
!       real(pr), optional, intent(out) :: ArV     !! dAr/dV [bar]
!       real(pr), optional, intent(out) :: ArT     !! dAr/dT [bar·L/K]
!       real(pr), optional, intent(out) :: ArTV    !! d²Ar/dVdT [bar/K]
!       real(pr), optional, intent(out) :: ArV2    !! d²Ar/dV² [bar/L]
!       real(pr), optional, intent(out) :: ArT2    !! d²Ar/dT² [bar·L/K²]
!       real(pr), optional, intent(out) :: Arn(size(n))  !! dAr/dn [bar·L/mol]
!       real(pr), optional, intent(out) :: ArVn(size(n)) !! d²Ar/dVdn [bar/mol]
!       real(pr), optional, intent(out) :: ArTn(size(n)) !! d²Ar/dTdn [bar·L/(mol·K)]
!       real(pr), optional, intent(out) :: Arn2(size(n),size(n)) !! d²Ar/dndn [bar·L/mol²]
! 
!       ! Local variables
!       real(pr) :: vgc, rho, rt_local, rg_ratio, ntot
!       real(pr) :: ar_local, arv_local, art_local, artv_local
!       real(pr) :: arv2_local, art2_local
!       real(pr), allocatable :: arn_local(:), arvn_local(:), artn_local(:)
!       real(pr), allocatable :: arn2_local(:,:)
! 
!       integer :: ider, itemp
! 
!       ! Conversion factors
!       rg_ratio = R / RGC  ! Convert atm·cm³ to bar·L
!       vgc = 1000.0_pr * v ! Convert L to cm³
!       ntot = sum(n)
!       rho = ntot / vgc
!       rt_local = t * RGC
! 
!       ! Determine which derivatives to calculate
!       ider = 0
!       itemp = 0
!       if (present(Arn2)) ider = 2
!       if (present(Arn) .or. present(ArVn) .or. present(ArTn)) ider = max(ider, 1)
!       if (present(ArT) .or. present(ArTV) .or. present(ArT2) .or. present(ArTn)) itemp = 1
! 
!       ! Allocate local arrays
!       allocate(arn_local(self%nc), arvn_local(self%nc), artn_local(self%nc))
!       allocate(arn2_local(self%nc, self%nc))
! 
!       ! Call the core GCA-EoS calculation
!       ! This is a simplified wrapper - the actual implementation would call
!       ! the detailed calculation similar to GCAEOS_Ar
!       call gca_core_calculation(self, self%nc, ider, itemp, t, vgc, n, &
!          ar_local, arv_local, art_local, artv_local, &
!          arv2_local, art2_local, arn_local, &
!          arvn_local, artn_local, arn2_local)
! 
!       ! Assign outputs if requested
!       if (present(Ar)) Ar = ar_local
!       if (present(ArV)) ArV = arv_local
!       if (present(ArT)) ArT = art_local
!       if (present(ArTV)) ArTV = artv_local
!       if (present(ArV2)) ArV2 = arv2_local
!       if (present(ArT2)) ArT2 = art2_local
!       if (present(Arn)) Arn = arn_local
!       if (present(ArVn)) ArVn = arvn_local
!       if (present(ArTn)) ArTn = artn_local
!       if (present(Arn2)) Arn2 = arn2_local
! 
!    end subroutine gca_residual_helmholtz
! 
!    subroutine gca_core_calculation(self, nc, ider, itemp, t, vgc, n, &
!       ar, arv, art, artv, arv2, art2, &
!       arn, arvn, artn, arn2)
!       !! Core GCA-EoS calculation
!       !! Combines attractive, free-volume, and association contributions
!       type(GCAEOS), intent(in) :: self
!       integer, intent(in) :: nc, ider, itemp
!       real(pr), intent(in) :: t, vgc, n(:)
!       real(pr), intent(out) :: ar, arv, art, artv, arv2, art2
!       real(pr), intent(out) :: arn(:), arvn(:), artn(:), arn2(:,:)
! 
!       ! Contribution variables
!       real(pr) :: f_att, dfv_att, dfvdv_att, dadt_att, dartv_att
!       real(pr) :: f_fv, dfv_fv, dfvdv_fv, dadt_fv, dartv_fv
!       real(pr) :: f_assoc, dfv_assoc, dfvdv_assoc, dadt_assoc, dartv_assoc
!       real(pr) :: f_total, dfv_total, dfvdv_total
!       real(pr) :: rg_ratio, rt_local, ntot, v2
! 
!       real(pr), allocatable :: dfdn_att(:), dfvdn_att(:), dfdnt_att(:), dfdn2_att(:,:)
!       real(pr), allocatable :: dfdn_fv(:), dfvdn_fv(:), dfdnt_fv(:), dfdn2_fv(:,:)
!       real(pr), allocatable :: dfdn_assoc(:), dfvdn_assoc(:), dfdnt_assoc(:), dfdn2_assoc(:,:)
!       real(pr), allocatable :: dpdn(:)
! 
!       ! Temperature-dependent parameters (local copies)
!       real(pr), allocatable :: d(:), dt(:)
!       real(pr), allocatable :: g(:,:), dgdt(:,:)
!       real(pr), allocatable :: delta(:,:), ddelt(:,:)
! 
!       real(pr) :: dadt_att_proper
!       ! Allocate derivative arrays
!       allocate(dfdn_att(nc), dfvdn_att(nc), dfdnt_att(nc), dfdn2_att(nc,nc))
!       allocate(dfdn_fv(nc), dfvdn_fv(nc), dfdnt_fv(nc), dfdn2_fv(nc,nc))
!       allocate(dfdn_assoc(nc), dfvdn_assoc(nc), dfdnt_assoc(nc), dfdn2_assoc(nc,nc))
!       allocate(dpdn(nc))
! 
!       ! Allocate and calculate temperature-dependent parameters
!       allocate(d(nc), dt(nc))
!       allocate(g(self%ng, self%ng), dgdt(self%ng, self%ng))
! 
!       call self%update_diameters(t, d, dt)
!       call self%update_group_energies(t, g, dgdt)
! 
!       if (self%nst > 0) then
!          allocate(delta(self%nst, self%nst), ddelt(self%nst, self%nst))
!          call self%update_association_parameters(t, delta, ddelt)
!       end if
! 
!       ! Initialize
!       rg_ratio = R / RGC
!       ntot = sum(n)
!       v2 = vgc * vgc
!       rt_local = t * RGC
! 
!       ! Initialize all outputs
!       f_att = 0.0_pr; dfv_att = 0.0_pr; dfvdv_att = 0.0_pr
!       dadt_att = 0.0_pr; dartv_att = 0.0_pr
!       dfdn_att = 0.0_pr; dfvdn_att = 0.0_pr; dfdnt_att = 0.0_pr
!       dfdn2_att = 0.0_pr
! 
!       f_fv = 0.0_pr; dfv_fv = 0.0_pr; dfvdv_fv = 0.0_pr
!       dadt_fv = 0.0_pr; dartv_fv = 0.0_pr
!       dfdn_fv = 0.0_pr; dfvdn_fv = 0.0_pr; dfdnt_fv = 0.0_pr
!       dfdn2_fv = 0.0_pr
! 
!       f_assoc = 0.0_pr; dfv_assoc = 0.0_pr; dfvdv_assoc = 0.0_pr
!       dadt_assoc = 0.0_pr; dartv_assoc = 0.0_pr
!       dfdn_assoc = 0.0_pr; dfvdn_assoc = 0.0_pr; dfdnt_assoc = 0.0_pr
!       dfdn2_assoc = 0.0_pr
! 
!       ! Calculate attractive contribution
!       call calc_attractive_contribution(self, nc, ider, itemp, t, vgc, n, &
!          g, dgdt, &
!          f_att, dfv_att, dfvdv_att, dadt_att, &
!          dartv_att, dfdn_att, dfvdn_att, &
!          dfdnt_att, dfdn2_att)
! 
!       ! Calculate free-volume contribution
!       call calc_freevolume_contribution(self, nc, ider, itemp, t, vgc, n, &
!          d, dt, &
!          f_fv, dfv_fv, dfvdv_fv, dadt_fv, &
!          dartv_fv, dfdn_fv, dfvdn_fv, &
!          dfdnt_fv, dfdn2_fv)
! 
!       ! Calculate association contribution
!       if (self%nst > 0) then
!          call calc_association_contribution(self, nc, ider, itemp, t, vgc, n, &
!             delta, ddelt, &
!             f_assoc, dfv_assoc, dfvdv_assoc, &
!             dadt_assoc, dartv_assoc, dfdn_assoc, &
!             dfvdn_assoc, dfdnt_assoc, dfdn2_assoc)
!       end if
! 
!       ! Combine contributions
!       f_total = f_att + f_fv + f_assoc
!       dfv_total = dfv_att + dfv_fv + dfv_assoc
!       dfvdv_total = dfvdv_att + dfvdv_fv + dfvdv_assoc
! 
!       ! Convert to Helmholtz energy and derivatives
!       ar = rt_local * f_total * rg_ratio
!       arv = rt_local * dfv_total * 1000.0_pr * rg_ratio
!       arv2 = rt_local * dfvdv_total * 1.0e6_pr * rg_ratio
! 
!       if (itemp > 0) then
! 
!          ! Free-volume: already returns dA^fv/dT [atm·cm³/K]
!          ! Attractive: returns dF^att/dT (dimensionless), needs conversion
!          ! Association: already returns dA^assoc/dT [atm·cm³/K]
! 
!          ! Convert attractive contribution
!          dadt_att_proper = RGC * (f_att + t * dadt_att)
! 
!          art = (dadt_fv + dadt_att_proper + dadt_assoc) * rg_ratio
!          print *, "dart calculation:", [dadt_fv, dadt_att_proper, dadt_assoc] * rg_ratio
!          artv = -(dartv_fv + dartv_att + dartv_assoc) * rg_ratio * 1000.0_pr
!       end if
! 
!       if (ider >= 1) then
!          arn = rt_local * (dfdn_att + dfdn_fv + dfdn_assoc) * rg_ratio
!          arvn = rt_local * (dfvdn_att + dfvdn_fv + dfvdn_assoc) * rg_ratio * 1000.0_pr
! 
!          ! Convert pressure derivative
!          dpdn = rt_local / vgc + rt_local * (dfvdn_att + dfvdn_fv + dfvdn_assoc)
!       end if
! 
!       if (itemp > 0 .and. ider >= 1) then
!          artn = rt_local * (dfdnt_att + dfdnt_fv + dfdnt_assoc + &
!             (dfdn_att + dfdn_fv + dfdn_assoc) / t) * rg_ratio
!       end if
! 
!       if (ider >= 2) then
!          arn2 = rt_local * (dfdn2_att + dfdn2_fv + dfdn2_assoc) * rg_ratio / ntot
!       end if
! 
!    end subroutine gca_core_calculation
! 
!    subroutine calc_attractive_contribution(self, nc, ider, itemp, t, vgc, n, &
!       g, dgdt, &
!       f_att, dfv_att, dfvdv_att, dadt_att, &
!       dartv_att, dfdn_att, dfvdn_att, &
!       dfdnt_att, dfdn2_att)
!       !! Calculate attractive contribution to Helmholtz energy
!       type(GCAEOS), intent(in) :: self
!       integer, intent(in) :: nc, ider, itemp
!       real(pr), intent(in) :: t, vgc, n(:)
!       real(pr), intent(in) :: g(:,:), dgdt(:,:)  !! Temperature-dependent energies
!       real(pr), intent(out) :: f_att, dfv_att, dfvdv_att, dadt_att, dartv_att
!       real(pr), intent(out) :: dfdn_att(:), dfvdn_att(:), dfdnt_att(:), dfdn2_att(:,:)
! 
!       ! Local variables
!       integer :: i, j, k
!       real(pr) :: dnom, qt, rt_local, ntot
!       real(pr) :: aah, argu, tet, thetae, thetaa, thetau
!       real(pr) :: h2, h3, h4, h5, h6, h7, h11, h12
!       real(pr) :: help_val, dhel, dat, ddel, dcrt, darg, alf
!       real(pr) :: akj, dakj, aa, arg
!       real(pr) :: xmsi, xmsk
!       real(pr) :: d2, d3, d4, d5, d6, d7, delik
! 
!       real(pr), allocatable :: nyt(:), theta(:)
!       real(pr), allocatable :: help2(:), help4(:), help5(:), help6(:)
!       real(pr), allocatable :: help11(:), help12(:)
!       real(pr), allocatable :: ms(:), ps(:,:)
!       real(pr), allocatable :: help3(:,:), help7(:,:), help8(:,:)
!       real(pr), allocatable :: help9(:,:), help10(:,:)
!       real(pr), allocatable :: aat_jk(:,:), deltaaat_jk(:,:), e(:,:)
!       real(pr), allocatable :: dh2dt(:), dh4dt(:), dh5dt(:), dh6dt(:)
!       real(pr), allocatable :: dh3dt(:,:), dh7dt(:,:)
! 
!       real(pr), parameter :: zz = 10.0_pr
!       real(pr), parameter :: macheps = epsilon(1.0_pr)
! 
!       ! Allocate arrays
!       allocate(nyt(self%ng), theta(self%ng))
!       allocate(help2(self%ng), help4(self%ng), help5(self%ng), help6(self%ng))
!       allocate(help11(self%ng), help12(self%ng))
!       allocate(ms(nc), ps(nc, self%ng))
!       allocate(help3(nc, self%ng), help7(nc, self%ng), help8(nc, self%ng))
!       allocate(help9(nc, self%ng), help10(nc, self%ng))
!       allocate(aat_jk(self%ng, self%ng), deltaaat_jk(self%ng, self%ng))
!       allocate(e(self%ng, self%ng))
!       allocate(dh2dt(self%ng), dh4dt(self%ng), dh5dt(self%ng), dh6dt(self%ng))
!       allocate(dh3dt(nc, self%ng), dh7dt(nc, self%ng))
! 
!       ! Initialize
!       ntot = sum(n)
!       rt_local = t * RGC
! 
!       ! Calculate group moles: ny_j = sum_i(n_i * nu_ji)
!       nyt = matmul(n(:nc), real(self%ny(:nc, :self%ng), pr))
! 
!       ! Total surface area
!       qt = dot_product(nyt(:self%ng), self%q(:self%ng))
! 
!       ! Surface fraction of each group
!       do j = 1, self%ng
!          if (qt > macheps) then
!             theta(j) = self%q(j) * nyt(j) / qt
!          else
!             theta(j) = 0.0_pr
!          end if
!       end do
! 
!       dnom = qt / (rt_local * vgc)
! 
!       ! Calculate auxiliary variables for each group k
!       dfv_att = 0.0_pr
!       dfvdv_att = 0.0_pr
! 
!       do k = 1, self%ng
!          help2(k) = 0.0_pr  ! H2_k/H4_k
!          help4(k) = 0.0_pr  ! H4_k
!          help5(k) = 0.0_pr  ! H5_k/H4_k
!          help6(k) = 0.0_pr  ! H6_k/H4_k
!          help11(k) = 0.0_pr
!          help12(k) = 0.0_pr
! 
!          do j = 1, self%ng
!             aah = g(j, k) * dnom
!             argu = self%alpha(j, k) * (g(j, k) - g(k, k)) * dnom
! 
!             aat_jk(j, k) = aah
!             deltaaat_jk(j, k) = argu
!             e(j, k) = exp(argu)
! 
!             if (theta(j) >= macheps) then
!                thetae = theta(j) * e(j, k)
!                help4(k) = help4(k) + thetae
! 
!                thetaa = thetae * aah
!                help2(k) = help2(k) + thetaa
! 
!                thetaa = thetaa * argu
!                help5(k) = help5(k) + thetaa
! 
!                thetaa = thetaa * argu
!                help11(k) = help11(k) + thetaa
! 
!                thetae = thetae * argu
!                help6(k) = help6(k) + thetae
! 
!                thetae = thetae * argu
!                help12(k) = help12(k) + thetae
!             end if
!          end do
! 
!          if (help4(k) > macheps) then
!             help2(k) = help2(k) / help4(k)
!             help5(k) = help5(k) / help4(k)
!             help6(k) = help6(k) / help4(k)
!             help11(k) = help11(k) / help4(k)
!             help12(k) = help12(k) / help4(k)
!          end if
! 
!          if (nyt(k) >= macheps) then
!             dfv_att = dfv_att + nyt(k) * self%q(k) * &
!                (help5(k) + help2(k) - help2(k) * help6(k))
!             dfvdv_att = dfvdv_att - (zz / 2.0_pr) * nyt(k) * self%q(k) * &
!                ((2.0_pr * help2(k) * help6(k)**2 - &
!                help2(k) * (help12(k) + 2.0_pr * help6(k)) - &
!                2.0_pr * help5(k) * help6(k)) + &
!                help11(k) + 2.0_pr * help5(k)) / (vgc * vgc)
!          end if
!       end do
! 
!       dfv_att = (zz / 2.0_pr) * dfv_att / vgc
!       dfvdv_att = dfvdv_att - 2.0_pr * dfv_att / vgc
! 
!       ! Helmholtz energy
!       f_att = -(zz / 2.0_pr) * ntot * &
!          sum(nyt(:self%ng) * self%q(:self%ng) * help2(:self%ng))
! 
!       ! Calculate MS and PS arrays for composition derivatives
!       do i = 1, nc
!          ms(i) = 0.0_pr
!          do k = 1, self%ng
!             ps(i, k) = 0.0_pr
!             if (self%ny(i, k) /= 0) then
!                ps(i, k) = real(self%ny(i, k), pr) * self%q(k)
!                ms(i) = ms(i) + ps(i, k)
!             end if
!          end do
!       end do
! 
!       ! Calculate HELP3, HELP7, HELP8, HELP9, HELP10 for composition derivatives
!       do i = 1, nc
!          xmsi = ms(i)
!          do k = 1, self%ng
!             help3(i, k) = 0.0_pr
!             help7(i, k) = 0.0_pr
!             help9(i, k) = 0.0_pr
!             help10(i, k) = 0.0_pr
! 
!             do j = 1, self%ng
!                if (self%ny(i, j) > 0) then
!                   tet = e(j, k)
!                   akj = aat_jk(j, k)
! 
!                   help_val = ps(i, j) * tet
!                   help7(i, k) = help7(i, k) + help_val
! 
!                   help_val = help_val * deltaaat_jk(j, k)
!                   help10(i, k) = help10(i, k) + help_val
! 
!                   help_val = ps(i, j) * tet * akj
!                   help3(i, k) = help3(i, k) + help_val
! 
!                   help9(i, k) = help9(i, k) + help_val * deltaaat_jk(j, k)
!                end if
!             end do
! 
!             if (help4(k) > macheps) then
!                help3(i, k) = help3(i, k) / help4(k)
!                help7(i, k) = help7(i, k) / help4(k)
!                help9(i, k) = help9(i, k) / help4(k)
!                help10(i, k) = help10(i, k) / help4(k)
!             end if
! 
!             help8(i, k) = help7(i, k) - xmsi * (1.0_pr - help6(k))
!          end do
!       end do
! 
!       ! Temperature derivatives
!       if (itemp /= 0) then
!          call calc_attractive_temperature_derivatives(self, nc, t, vgc, n, &
!             theta, nyt, help2, help4, &
!             help5, help6, aat_jk, &
!             deltaaat_jk, e, ps, ms, &
!             g, dgdt, &
!             dadt_att, dartv_att, &
!             dh2dt, dh4dt, dh5dt, dh6dt, &
!             dh3dt, dh7dt)
!       end if
! 
!       ! Composition derivatives
!       if (ider >= 1) then
!          call calc_attractive_composition_derivatives(self, nc, ider, itemp, &
!             t, vgc, qt, ntot, n, &
!             theta, nyt, help2, help3, &
!             help4, help5, help6, help7, &
!             help8, help9, help10, help11, &
!             help12, ms, ps, dh2dt, dh3dt, &
!             dh4dt, dh5dt, dh6dt, dh7dt, &
!             dfdn_att, dfvdn_att, &
!             dfdnt_att, dfdn2_att)
!       end if
! 
!    end subroutine calc_attractive_contribution
! 
!    subroutine calc_attractive_temperature_derivatives(self, nc, t, vgc, n, &
!       theta, nyt, help2, help4, &
!       help5, help6, aat_jk, &
!       deltaaat_jk, e, ps, ms, &
!       g, dgdt, &
!       dadt_att, dartv_att, &
!       dh2dt, dh4dt, dh5dt, dh6dt, &
!       dh3dt, dh7dt)
!       !! Calculate temperature derivatives for attractive contribution
!       type(GCAEOS), intent(in) :: self
!       integer, intent(in) :: nc
!       real(pr), intent(in) :: t, vgc, n(:)
!       real(pr), intent(in) :: theta(:), nyt(:), help2(:), help4(:)
!       real(pr), intent(in) :: help5(:), help6(:)
!       real(pr), intent(in) :: aat_jk(:,:), deltaaat_jk(:,:), e(:,:)
!       real(pr), intent(in) :: ps(:,:), ms(:)
!       real(pr), intent(in) :: g(:,:), dgdt(:,:)  !! Temperature-dependent energies
!       real(pr), intent(out) :: dadt_att, dartv_att
!       real(pr), intent(out) :: dh2dt(:), dh4dt(:), dh5dt(:), dh6dt(:)
!       real(pr), intent(out) :: dh3dt(:,:), dh7dt(:,:)
! 
!       integer :: i, j, k
!       real(pr) :: aa, dat, dakj, akj, ddel, darg, dcrt, alf, arg
!       real(pr) :: tet, thetau, help_val, helpi, helpi1, helpi2, help1
!       real(pr) :: h30, h31, h32, h33, h34, h35, h36
!       real(pr) :: h2, h4, h5, h6, dh2t, dh4t, dh5t, dh6t
!       real(pr) :: dhel, hel, rt_local, dnom, psik
! 
!       real(pr), parameter :: zz = 10.0_pr
! 
!       rt_local = t * RGC
!       dnom = sum(nyt * self%q) / (rt_local * vgc)
! 
!       dadt_att = 0.0_pr
!       dhel = 0.0_pr
!       dh3dt = 0.0_pr
!       dh7dt = 0.0_pr
! 
!       do j = 1, self%ng
!          h30 = 0.0_pr; h31 = 0.0_pr; h32 = 0.0_pr
!          h33 = 0.0_pr; h34 = 0.0_pr; h35 = 0.0_pr; h36 = 0.0_pr
! 
!          aa = g(j, j)
!          dat = dgdt(j, j)
! 
!          do k = 1, self%ng
!             tet = e(k, j)
!             thetau = theta(k) * tet
!             dakj = dgdt(k, j)
!             arg = deltaaat_jk(k, j)
!             akj = g(k, j)
!             ddel = dakj - dat
! 
!             darg = 0.0_pr
!             dcrt = abs(akj - aa)
!             if (dcrt > 1.0e-2_pr) then
!                darg = arg * (ddel / (akj - aa) - 1.0_pr / t)
!             end if
! 
!             alf = self%alpha(k, j)
!             helpi = tet * darg
!             help_val = thetau * darg
! 
!             h30 = h30 + help_val
! 
!             helpi1 = helpi * akj * dnom
!             help1 = help_val * akj * dnom
!             h31 = h31 + help1
!             h32 = h32 + help_val * arg
!             h33 = h33 + help1 * arg
! 
!             helpi2 = tet * dakj * dnom
!             help_val = thetau * dakj * dnom
!             h34 = h34 + help_val
!             h35 = h35 + help_val * arg
!             h36 = h36 + thetau * ddel * alf * dnom
! 
!             do i = 1, nc
!                if (self%ny(i, k) /= 0) then
!                   psik = ps(i, k)
!                   dh3dt(i, j) = dh3dt(i, j) + helpi1 * psik + helpi2 * psik
!                   dh7dt(i, j) = dh7dt(i, j) + helpi * psik
!                end if
!             end do
!          end do
! 
!          h2 = help2(j)
!          h4 = help4(j)
!          h5 = help5(j)
!          h6 = help6(j)
! 
!          dh2t = (h31 + h34) / h4
!          dh4t = h30 / h4
!          dh5t = (h33 + h35 + h31) / h4
!          dh6t = -h6 / t + (h32 + h36) / h4
! 
!          dh2dt(j) = dh2t
!          dh4dt(j) = dh4t
!          dh5dt(j) = dh5t
!          dh6dt(j) = dh6t
! 
!          dadt_att = dadt_att + nyt(j) * self%q(j) * (dh2t - dh4t * h2)
! 
!          hel = dh5t + dh2t - dh2t * h6 - h2 * dh6t - (h5 + h2 - 2.0_pr * h2 * h6) * dh4t
!          dhel = dhel + hel * nyt(j) * self%q(j)
!       end do
! 
!       dadt_att = -(zz / 2.0_pr) * dadt_att
!       dartv_att = -(zz / 2.0_pr) * rt_local / vgc * dhel
! 
!    end subroutine calc_attractive_temperature_derivatives
! 
!    subroutine calc_attractive_composition_derivatives(self, nc, ider, itemp, &
!       t, vgc, qt, ntot, n, &
!       theta, nyt, help2, help3, &
!       help4, help5, help6, help7, &
!       help8, help9, help10, help11, &
!       help12, ms, ps, dh2dt, dh3dt, &
!       dh4dt, dh5dt, dh6dt, dh7dt, &
!       dfdn_att, dfvdn_att, &
!       dfdnt_att, dfdn2_att)
!       !! Calculate composition derivatives for attractive contribution
!       type(GCAEOS), intent(in) :: self
!       integer, intent(in) :: nc, ider, itemp
!       real(pr), intent(in) :: t, vgc, qt, ntot, n(:)
!       real(pr), intent(in) :: theta(:), nyt(:), help2(:), help3(:,:)
!       real(pr), intent(in) :: help4(:), help5(:), help6(:), help7(:,:)
!       real(pr), intent(in) :: help8(:,:), help9(:,:), help10(:,:)
!       real(pr), intent(in) :: help11(:), help12(:), ms(:), ps(:,:)
!       real(pr), intent(in) :: dh2dt(:), dh3dt(:,:), dh4dt(:)
!       real(pr), intent(in) :: dh5dt(:), dh6dt(:), dh7dt(:,:)
!       real(pr), intent(out) :: dfdn_att(:), dfvdn_att(:)
!       real(pr), intent(out) :: dfdnt_att(:), dfdn2_att(:,:)
! 
!       integer :: i, j, k, nc1
!       real(pr) :: xmsi, xmsk, tet, psij, pskj
!       real(pr) :: h12, h11
!       real(pr) :: h2, h3ij, h4, h5, h6, h7ij, h8ij
!       real(pr) :: h9ij, h9kj, h10ij, h10kj, h20, h21
!       real(pr) :: deli, delik, hder, dp, helpt
!       real(pr) :: dh2t, dh4t, dh5t, dh6t
!       real(pr) :: h2h4, h2h8, h3h5, h7h6, h8, d2, d3, d4, d5, d6, d7
!       real(pr) :: rt_local
! 
!       real(pr), parameter :: zz = 10.0_pr
! 
!       rt_local = t * RGC
! 
!       do i = 1, nc
!          xmsi = ms(i)
!          nc1 = i
!          dfdn_att(i) = 0.0_pr
!          dfvdn_att(i) = 0.0_pr
!          dfdnt_att(i) = 0.0_pr
! 
!          if (ider >= 2) then
!             dfdn2_att(i, nc1:nc) = 0.0_pr
!          end if
! 
!          do j = 1, self%ng
!             tet = theta(j)
!             psij = ps(i, j)
!             h3ij = help3(i, j)
!             h7ij = help7(i, j)
!             h2 = help2(j)
!             h4 = help4(j)
!             h5 = help5(j)
!             h6 = help6(j)
!             h8ij = help8(i, j)
! 
!             ! Calculate dF/dni attractive contribution
!             hder = -h2 * (h7ij - xmsi * (1.0_pr - h6)) + h3ij + xmsi * h5
!             deli = -(hder * tet + psij * h2)
!             dfdn_att(i) = dfdn_att(i) + deli
! 
!             ! Temperature-composition cross derivatives
!             if (itemp > 0) then
!                dh2t = dh2dt(j)
!                dh4t = dh4dt(j)
!                dh5t = dh5dt(j)
!                dh6t = dh6dt(j)
! 
!                helpt = psij * (dh2t - h2 * dh4t) + &
!                   tet * (dh3dt(i, j) / h4 + xmsi * dh5t - &
!                   (h3ij + xmsi * h5) * dh4t) - &
!                   tet * h8ij * (dh2t - 2.0_pr * h2 * dh4t) - &
!                   tet * h2 * (dh7dt(i, j) / h4 - dh4t * xmsi + dh6t * xmsi)
! 
!                dfdnt_att(i) = dfdnt_att(i) - deli / t - helpt
!             end if
! 
!             ! Volume-composition derivatives
!             if (ider >= 1) then
!                h9ij = help9(i, j)
!                h10ij = help10(i, j)
!                h11 = help11(j)
!                h12 = help12(j)
!                h20 = h5 + h2 - h2 * h6
!                h21 = nyt(j) * self%q(j)
! 
!                dp = -psij * h20 - h21 / qt * (xmsi * (3.0_pr * h5 + h2 + h11) + &
!                   h9ij + h3ij) + h21 / qt * (h6 * (h3ij + 2.0_pr * xmsi * h5) + &
!                   h2 * (h10ij + xmsi * h12) + h7ij * (h5 + h2) + &
!                   3.0_pr * xmsi * h2 * h6) - h21 / qt * 2.0_pr * &
!                   (h7ij * h2 * h6 + xmsi * h2 * h6 * h6)
! 
!                dfvdn_att(i) = dfvdn_att(i) + (zz / 2.0_pr) * dp * rt_local / vgc
!             end if
! 
!             ! Second composition derivatives
!             if (ider >= 2) then
!                h2h8 = h2 * h8ij / qt
!                h3h5 = (h3ij + xmsi * h5) / qt
!                h7h6 = tet * h2 * (h7ij + h6 * xmsi)
!                h8 = tet * h8ij
!                h2h4 = tet * h2
! 
!                do k = nc1, nc
!                   xmsk = ms(k)
!                   pskj = ps(k, j)
!                   h9kj = help9(k, j)
!                   h10kj = help10(k, j)
! 
!                   d2 = (help3(k, j) + xmsk * h5) / qt
!                   d3 = xmsk / qt * (h9ij + h3ij)
!                   d4 = (help7(k, j) - xmsk * (1.0_pr - h6)) / qt
!                   d5 = (h9kj + xmsk * (h11 + h5)) / qt
!                   d6 = (h10kj + xmsk * h12) / qt
!                   d7 = xmsk * h10ij / qt
! 
!                   delik = -(psij * d2 + (pskj - tet * xmsk) * h3h5) + &
!                      (pskj - tet * xmsk) * h2h8 - d4 * h7h6 + d2 * h8 - &
!                      deli * d4 + h2h4 * (d7 + xmsi * d6) - &
!                      tet * (d3 + xmsi * d5)
! 
!                   dfdn2_att(i, k) = dfdn2_att(i, k) + delik
!                end do
!             end if
!          end do
!       end do
! 
!       ! Apply zz factor
!       dfdn_att = dfdn_att * (zz / 2.0_pr)
!       dfdnt_att = dfdnt_att * (zz / 2.0_pr)
! 
!       if (ider >= 2) then
!          dfdn2_att = dfdn2_att * (zz / 2.0_pr)
!       end if
! 
!    end subroutine calc_attractive_composition_derivatives
! 
!    subroutine calc_freevolume_contribution(self, nc, ider, itemp, t, vgc, n, &
!       d, dt, &
!       f_fv, dfv_fv, dfvdv_fv, dadt_fv, &
!       dartv_fv, dfdn_fv, dfvdn_fv, &
!       dfdnt_fv, dfdn2_fv)
!       !! Calculate free-volume (repulsive) contribution
!       type(GCAEOS), intent(in) :: self
!       integer, intent(in) :: nc, ider, itemp
!       real(pr), intent(in) :: t, vgc, n(:)
!       real(pr), intent(in) :: d(:), dt(:)  !! Temperature-dependent diameters
!       real(pr), intent(out) :: f_fv, dfv_fv, dfvdv_fv, dadt_fv, dartv_fv
!       real(pr), intent(out) :: dfdn_fv(:), dfvdn_fv(:), dfdnt_fv(:), dfdn2_fv(:,:)
! 
!       integer :: i, j, k
!       real(pr) :: ntot, b_gc, rt_local
!       real(pr) :: lambda1, lambda2, lambda3
!       real(pr) :: y, y2, ln_y, dydt, dydv, dydvv, dydvt
!       real(pr) :: r0, r1, r2, r3, r4, r6, r33, r34
!       real(pr) :: prep, dfvr, dfvdvr, dffvdt, dpdtr
!       real(pr) :: dia, dla, diav, diat, deldt, xdt
!       real(pr) :: tlam1, tlam2, tlam3
!       real(pr) :: r20, r21, r25, r30, r31
!       real(pr) :: dl1, dl2, dl3, dl1tn, dl2tn, dl3tn
!       real(pr) :: dfdnr, dydni, dydtni, dyvdni
!       real(pr) :: dprep, dyik, dydnk, d2fdn2_fv
!       real(pr) :: dl1ik, dl2ik, dl3ik
!       real(pr) :: dlam1k, dlam2k, dlam3k
! 
!       real(pr), allocatable :: dlam1(:), dlam2(:), dlam3(:)
!       real(pr), allocatable :: dlamt1(:), dlamt2(:), dlamt3(:)
!       real(pr), allocatable :: dydn(:), dydtn(:), dyvdn(:)
! 
!       real(pr), parameter :: pi = acos(-1.0_pr)
!       real(pr), parameter :: macheps = epsilon(1.0_pr)
! 
!       ! Allocate arrays
!       allocate(dlam1(nc), dlam2(nc), dlam3(nc))
!       allocate(dlamt1(nc), dlamt2(nc), dlamt3(nc))
!       allocate(dydn(nc), dydtn(nc), dyvdn(nc))
! 
!       ! Initialize
!       ntot = sum(n)
!       rt_local = t * RGC
! 
!       ! Calculate lambda parameters
!       lambda1 = 0.0_pr
!       lambda2 = 0.0_pr
!       lambda3 = 0.0_pr
!       tlam1 = 0.0_pr
!       tlam2 = 0.0_pr
!       tlam3 = 0.0_pr
! 
!       do i = 1, nc
!          dia = d(i)
!          dla = n(i) * dia
!          lambda1 = lambda1 + dla
! 
!          diav = dia * dia
!          dla = dla * dia
!          lambda2 = lambda2 + dla
! 
!          lambda3 = lambda3 + dla * dia
! 
!          if (itemp /= 0) then
!             diat = dt(i)
!             deldt = diat
!             xdt = n(i) * diat
! 
!             dlamt1(i) = deldt
!             tlam1 = tlam1 + xdt
! 
!             deldt = 2.0_pr * deldt * dia
!             dlamt2(i) = deldt
!             xdt = xdt * dia
!             tlam2 = tlam2 + xdt
! 
!             dlamt3(i) = 3.0_pr * deldt * dia / 2.0_pr
!             tlam3 = tlam3 + xdt * dia
!          end if
! 
!          dlam1(i) = dia / lambda1
!          dlam2(i) = diav / lambda2
!          dlam3(i) = diav * dia
!       end do
! 
!       ! Net co-volume
!       b_gc = pi * lambda3 / 6.0_pr
! 
!       ! Y function and derivatives
!       y = 1.0_pr / (1.0_pr - b_gc / vgc)
!       y2 = y * y
!       dydv = -y2 * b_gc / (vgc * vgc)
!       dydvv = 2.0_pr * dydv**2 / y - 2.0_pr * dydv / vgc
! 
!       ! Temperature derivatives of Y
!       if (itemp /= 0) then
!          tlam2 = tlam2 * 2.0_pr
!          tlam3 = tlam3 * 3.0_pr
!          dydt = y2 * pi * tlam3 / vgc / 6.0_pr
!          dydvt = dydv * tlam3 / lambda3 + 2.0_pr * dydv * dydt / y
! 
!          tlam3 = tlam3 / lambda3
!          tlam2 = tlam2 / lambda2
!          tlam1 = tlam1 / lambda1
!       end if
! 
!       ! Composition derivatives of lambda
!       do i = 1, nc
!          dlam3(i) = dlam3(i) / lambda3
!       end do
! 
!       ! R auxiliary variables
!       r0 = lambda2 / lambda3
!       r1 = r0 * lambda1
!       r2 = r0 * r0 * lambda2
!       r3 = r2
!       ln_y = log(y)
!       r4 = y2 - y - ln_y
!       r34 = y - 1.0_pr
! 
!       if (r34 <= macheps) then
!          r34 = b_gc * y / vgc
!       end if
! 
!       r6 = 2.0_pr * y - 1.0_pr - 1.0_pr / y
!       r33 = 2.0_pr + 1.0_pr / y2
! 
!       ! Free-volume Helmholtz energy
!       f_fv = r1 * r34 + r2 * r4 + ntot * ln_y
! 
!       ! First volume derivative
!       dfvr = (r1 + r2 * r6 + 1.0_pr / y) * dydv
!       dfv_fv = dfvr
! 
!       ! Second volume derivative
!       dfvdvr = dfvr * dydvv / dydv + (r2 * r33 - 1.0_pr / y2) * dydv**2
!       dfvdv_fv = dfvdvr
! 
!       ! Pressure (repulsive)
!       prep = -rt_local * dfvr
! 
!       ! Temperature derivatives
!       if (itemp /= 0) then
!          r30 = tlam1 + tlam2 - tlam3
!          r31 = 3.0_pr * tlam2 - 2.0_pr * tlam3
! 
!          dffvdt = 3.0_pr * r1 * r30 * r34 + 3.0_pr * r1 * dydt + &
!             r3 * r6 * dydt + r3 * r31 * r4 + ntot / y * dydt
! 
!          dadt_fv = -RGC * (t * dffvdt + f_fv)
!          ! FROM (incorrect):
!          dadt_fv = -RGC * (t * dffvdt + f_fv)
! 
!          ! TO (correct):
!          dadt_fv = RGC * (f_fv + t * dffvdt)
! 
!          dpdtr = prep / t + dydvt / dydv * prep - &
!             rt_local * dydv * (3.0_pr * r1 * r30 + r3 * r31 * r6 + &
!             r3 * r33 * dydt - ntot / y2 * dydt)
! 
!          dartv_fv = -dpdtr
!       end if
! 
!       ! Composition derivatives
!       if (ider >= 1) then
!          do i = 1, nc
!             dydn(i) = y2 * pi * dlam3(i) / vgc / 6.0_pr
!             dyvdn(i) = (2.0_pr * dydv / y - 1.0_pr / vgc) * dydn(i)
! 
!             if (itemp /= 0) then
!                dydtn(i) = 2.0_pr * dydn(i) / y * dydt + &
!                   dlamt3(i) / dlam3(i) * dydn(i)
!             end if
! 
!             dl1 = dlam1(i)
!             dl2 = dlam2(i)
!             dl3 = dlam3(i)
!             r2 = dl1 + dl2 - dl3
!             r20 = r2 + dydn(i) / r34
!             r20 = r20 * r1 * 3.0_pr
!             r21 = 3.0_pr * dl2 - 2.0_pr * dl3
! 
!             dfdnr = r34 * r20 + r3 * r4 * r21 + r3 * r6 * dydn(i) + &
!                ln_y + ntot / y * dydn(i)
! 
!             dfdn_fv(i) = dfdnr
! 
!             if (itemp /= 0) then
!                dl1tn = dlamt1(i) / lambda1 - tlam1 * dl1
!                dl2tn = dlamt2(i) / lambda2 - tlam2 * dl2
!                dl3tn = dlamt3(i) / lambda3 - tlam3 * dl3
! 
!                dfdnt_fv(i) = dydt * r20 + r34 * r30 * r20 + &
!                   3.0_pr * r34 * r1 * (dl1tn + dl2tn - dl3tn + &
!                   (dydtn(i) - dydt * dydn(i) / r34) / r34) + &
!                   r3 * (r31 * r4 * r21 + dydt * r6 * r21 + &
!                   r4 * (3.0_pr * dl2tn - 2.0_pr * dl3tn) + &
!                   r6 * r31 * dydn(i) + r33 * dydt * dydn(i) + &
!                   r6 * dydtn(i)) + dydt / y - &
!                   ntot / y * (dydn(i) * dydt / y - dydtn(i))
!             end if
! 
!             if (ider >= 1) then
!                dprep = prep * dyvdn(i) / dydv - dydv * rt_local * &
!                   (3.0_pr * r1 * r2 + r3 * r6 * r21 + &
!                   r3 * r33 * dydn(i) + 1.0_pr / y - ntot / y2 * dydn(i))
! 
!                dfvdn_fv(i) = dprep / rt_local
!             end if
!          end do
!       end if
! 
!       ! Second composition derivatives
!       if (ider >= 2) then
!          do i = 1, nc
!             dydni = dydn(i)
!             dl1 = dlam1(i)
!             dl2 = dlam2(i)
!             dl3 = dlam3(i)
!             r2 = dl1 + dl2 - dl3
!             r20 = r2 + dydni / r34
!             r20 = r20 * r1 * 3.0_pr
!             r21 = 3.0_pr * dl2 - 2.0_pr * dl3
! 
!             do k = i, nc
!                dydnk = dydn(k)
!                dlam1k = dlam1(k)
!                dlam2k = dlam2(k)
!                dlam3k = dlam3(k)
!                r25 = 3.0_pr * dlam2k - 2.0_pr * dlam3k
! 
!                d2fdn2_fv = r20 * (dydnk + r34 * (dlam1k + dlam2k - dlam3k))
! 
!                dl1ik = dlam1k * dl1
!                dl2ik = dlam2k * dl2
!                dl3ik = dlam3k * dl3
!                dyik = dydni * dydnk
! 
!                d2fdn2_fv = d2fdn2_fv + 3.0_pr * r34 * r1 * &
!                   (-dl1ik - dl2ik + dl3ik) - &
!                   3.0_pr * r1 * (dyik / r34 - 2.0_pr * dyik / y) + &
!                   r4 * r21 * r25 * r3 + r3 * r21 * r6 * dydnk + &
!                   r3 * r4 * (2.0_pr * dlam3k * dl3 - 3.0_pr * dlam2k * dl2) + &
!                   r3 * r25 * r6 * dydni + r3 * (r33 + 2.0_pr * r6 / y) * dyik + &
!                   (dydni + dydnk) / y + ntot * dyik / y2
! 
!                dfdn2_fv(i, k) = d2fdn2_fv
!                dfdn2_fv(k, i) = d2fdn2_fv
!             end do
!          end do
!       end if
! 
!    end subroutine calc_freevolume_contribution
! 
!    subroutine calc_association_contribution(self, nc, ider, itemp, t, vgc, n, &
!       delta, ddelt, &
!       f_assoc, dfv_assoc, dfvdv_assoc, &
!       dadt_assoc, dartv_assoc, dfdn_assoc, &
!       dfvdn_assoc, dfdnt_assoc, dfdn2_assoc)
!       !! Calculate association contribution
!       type(GCAEOS), intent(in) :: self
!       integer, intent(in) :: nc, ider, itemp
!       real(pr), intent(in) :: t, vgc, n(:)
!       real(pr), intent(in) :: delta(:,:), ddelt(:,:)  !! Temperature-dependent association
!       real(pr), intent(out) :: f_assoc, dfv_assoc, dfvdv_assoc
!       real(pr), intent(out) :: dadt_assoc, dartv_assoc
!       real(pr), intent(out) :: dfdn_assoc(:), dfvdn_assoc(:)
!       real(pr), intent(out) :: dfdnt_assoc(:), dfdn2_assoc(:,:)
! 
!       integer :: i, j, k, m, in
!       integer, parameter :: maxIt = 200
!       real(pr) :: ntot, rho, rt_local, sum_s
!       real(pr) :: aux_d2fn, dfast, dpdtas
! 
!       real(pr), allocatable :: sm(:), xs(:), sm_xs(:)
!       real(pr), allocatable :: dxs_dv(:), dxs_dt(:), aux(:)
!       real(pr), allocatable :: dxs_dn(:,:), dxs_ds(:,:)
!       real(pr), allocatable :: h(:,:), d2fassoc_dsdt(:)
!       integer, allocatable :: pivot_vector(:)
! 
!       ! Allocate arrays
!       allocate(sm(self%nst), xs(self%nst), sm_xs(self%nst))
!       allocate(dxs_dv(self%nst), dxs_dt(self%nst), aux(self%nst))
!       allocate(dxs_dn(self%nst, nc), dxs_ds(self%nst, self%nst))
!       allocate(h(self%nst, self%nst), d2fassoc_dsdt(self%nst))
!       allocate(pivot_vector(self%nst))
! 
!       ! Initialize
!       ntot = sum(n)
!       rho = 1.0_pr / vgc
!       rt_local = t * RGC
! 
!       xs = 1.0_pr
!       dxs_dv = 0.0_pr
!       dxs_ds = 0.0_pr
!       sm_xs = 0.0_pr
!       sm = 0.0_pr
! 
!       ! Calculate site moles
!       sm = matmul(real(self%sigma(:self%nst, :nc), pr), n(:nc))
! 
!       ! Initial guess with direct substitution
!       do i = 1, 3
!          do j = 1, self%nst
!             aux(j) = 1.0_pr / (1.0_pr + sum(sm(:self%nst) * xs(:self%nst) * &
!                delta(:self%nst, j) / vgc))
!          end do
!          xs(:self%nst) = aux(:self%nst)
!       end do
! 
!       ! Total number of association sites
!       sum_s = sum(sm(:self%nst))
! 
!       ! Solve for non-bonded fractions using Newton-Raphson
!       call optinewton(maxIt, self%nst, sm, delta, rho, xs, h, pivot_vector, in)
! 
!       ! Number of non-bonded sites
!       sm_xs(:self%nst) = sm(:self%nst) * xs(:self%nst)
! 
!       ! Association Helmholtz energy
!       f_assoc = dot_product(sm(:self%nst), log(xs(:self%nst))) + &
!          (sum_s - sum(sm_xs(:self%nst))) / 2.0_pr
! 
!       ! First volume derivative
!       dxs_dv(:self%nst) = -sm(:self%nst) * &
!          matmul(delta(:self%nst, :self%nst), sm_xs(:self%nst)) / &
!          (vgc * vgc)
! 
!       call lubksb(h(:self%nst, :self%nst), self%nst, self%nst, &
!          pivot_vector, dxs_dv(:self%nst))
! 
!       dfv_assoc = (sum_s - sum(sm_xs(:self%nst))) / vgc / 2.0_pr
!       dfvdv_assoc = -(dfv_assoc + dot_product(sm(:self%nst), dxs_dv(:self%nst)) / 2.0_pr) / vgc
! 
!       ! Composition derivatives
!       if (ider >= 2 .or. itemp >= 1) then
!          do i = 1, self%nst
!             if (sm(i) >= epsilon(1.0_pr)) then
!                dxs_ds(i, :self%nst) = -sm(i) * xs(:self%nst) * &
!                   delta(i, :self%nst) / vgc / ntot
!                dxs_ds(i, i) = dxs_ds(i, i) + &
!                   (1.0_pr / xs(i) - 1.0_pr - &
!                   dot_product(sm_xs(:self%nst), delta(i, :self%nst)) / vgc) / ntot
!             end if
!          end do
! 
!          do i = 1, self%nst
!             dxs_ds(:self%nst, i) = -dxs_ds(:self%nst, i)
!             call lubksb(h(:self%nst, :self%nst), self%nst, self%nst, &
!                pivot_vector(:self%nst), dxs_ds(:self%nst, i))
!          end do
! 
!          dxs_dn(:self%nst, :nc) = matmul(dxs_ds(:self%nst, :self%nst), &
!             real(self%sigma(:self%nst, :nc), pr))
! 
!          if (itemp >= 1) then
!             aux(:self%nst) = matmul(ddelt(:self%nst, :self%nst), sm_xs(:self%nst))
! 
!             dfast = -dot_product(aux(:self%nst), sm_xs(:self%nst)) / 2.0_pr / vgc
!             dadt_assoc = RGC * (t * dfast + f_assoc)
! 
!             dpdtas = -RGC * dfv_assoc
! 
!             do i = 1, self%nst
!                dpdtas = dpdtas + rt_local * sm_xs(i) * &
!                   sum(sm_xs(:self%nst) * ddelt(i, :self%nst) * &
!                   (dxs_dv(:self%nst) / xs(:self%nst) - 1.0_pr / vgc / 2.0_pr)) / vgc
!             end do
! 
!             dartv_assoc = -dpdtas
! 
!             do k = 1, self%nst
!                aux(:self%nst) = sm(:self%nst) * dxs_ds(:self%nst, k)
!                aux(:self%nst) = matmul(ddelt(:self%nst, :self%nst), aux(:self%nst))
!                d2fassoc_dsdt(k) = dot_product(sm_xs(:self%nst), aux(:self%nst))
!                d2fassoc_dsdt(k) = -(d2fassoc_dsdt(k) + xs(k) * &
!                   dot_product(sm_xs(:self%nst), &
!                   ddelt(:self%nst, k))) / vgc
!             end do
!          end if
!       end if
! 
!       ! Composition derivatives of F
!       do i = 1, nc
!          dfdn_assoc(i) = 0.0_pr
!          dfvdn_assoc(i) = 0.0_pr
!          dfdnt_assoc(i) = 0.0_pr
! 
!          if (maxval(self%sigma(:self%nst, i)) > 0) then
!             do j = 1, self%nst
!                dfdn_assoc(i) = dfdn_assoc(i) + &
!                   log(xs(j)) * real(self%sigma(j, i), pr)
!             end do
! 
!             if (itemp > 0 .or. ider == 2) then
!                dfvdn_assoc(i) = -rt_local * &
!                   sum(real(self%sigma(:self%nst, i), pr) * &
!                   dxs_dv(:self%nst) / xs(:self%nst)) / rt_local
!             end if
! 
!             if (itemp > 0) then
!                dfdnt_assoc(i) = dot_product(real(self%sigma(:self%nst, i), pr), &
!                   d2fassoc_dsdt(:self%nst))
!             end if
! 
!             if (ider == 2) then
!                do k = i, nc
!                   aux_d2fn = sum(real(self%sigma(:self%nst, i), pr) * &
!                      dxs_dn(:self%nst, k) / xs(:self%nst))
!                   dfdn2_assoc(i, k) = dfdn2_assoc(i, k) + aux_d2fn
!                   dfdn2_assoc(k, i) = dfdn2_assoc(i, k)
!                end do
!             end if
!          end if
!       end do
! 
!    end subroutine calc_association_contribution
! 
!    ! ========================================================================
!    ! Helper subroutines for linear algebra and association calculations
!    ! ========================================================================
! 
!    subroutine ludcmp(a, n, np, indx, d)
!       !! LU decomposition of matrix A
!       !! Finds U such that A = U · L where U is upper triangular
!       !! and L is lower triangular with unit diagonal
!       integer, intent(in) :: n, np
!       real(pr), intent(inout) :: a(np, np)
!       integer, intent(out) :: indx(n)
!       real(pr), intent(out) :: d
! 
!       integer, parameter :: nmax = 500
!       real(pr), parameter :: tiny = 1.0e-20_pr
! 
!       integer :: i, imax, j, k
!       real(pr) :: aamax, dum, sum_val, vv(nmax)
! 
!       d = 1.0_pr
! 
!       do i = 1, n
!          aamax = 0.0_pr
!          do j = 1, n
!             if (abs(a(i, j)) > aamax) aamax = abs(a(i, j))
!          end do
! 
!          if (aamax == 0.0_pr) then
!             ! Singular matrix
!             aamax = tiny
!          end if
! 
!          vv(i) = 1.0_pr / aamax
!       end do
! 
!       do j = 1, n
!          do i = 1, j - 1
!             sum_val = a(i, j)
!             do k = 1, i - 1
!                sum_val = sum_val - a(i, k) * a(k, j)
!             end do
!             a(i, j) = sum_val
!          end do
! 
!          aamax = 0.0_pr
!          imax = j
! 
!          do i = j, n
!             sum_val = a(i, j)
!             do k = 1, j - 1
!                sum_val = sum_val - a(i, k) * a(k, j)
!             end do
!             a(i, j) = sum_val
!             dum = vv(i) * abs(sum_val)
! 
!             if (dum >= aamax) then
!                imax = i
!                aamax = dum
!             end if
!          end do
! 
!          if (j /= imax) then
!             do k = 1, n
!                dum = a(imax, k)
!                a(imax, k) = a(j, k)
!                a(j, k) = dum
!             end do
!             d = -d
!             vv(imax) = vv(j)
!          end if
! 
!          indx(j) = imax
! 
!          if (a(j, j) == 0.0_pr) a(j, j) = tiny
! 
!          if (j /= n) then
!             dum = 1.0_pr / a(j, j)
!             do i = j + 1, n
!                a(i, j) = a(i, j) * dum
!             end do
!          end if
!       end do
! 
!    end subroutine ludcmp
! 
!    subroutine lubksb(a, n, np, indx, b)
!       !! Solve the system U · x = b by back substitution
!       !! Complements LUDcmp which decomposes A into U*L
!       integer, intent(in) :: n, np
!       real(pr), intent(in) :: a(np, np)
!       integer, intent(in) :: indx(n)
!       real(pr), intent(inout) :: b(n)
! 
!       integer :: i, ii, j, ll
!       real(pr) :: sum_val
! 
!       ii = 0
! 
!       do i = 1, n
!          ll = indx(i)
!          sum_val = b(ll)
!          b(ll) = b(i)
! 
!          if (ii /= 0) then
!             do j = ii, i - 1
!                sum_val = sum_val - a(i, j) * b(j)
!             end do
!          else if (sum_val /= 0.0_pr) then
!             ii = i
!          end if
! 
!          b(i) = sum_val
!       end do
! 
!       do i = n, 1, -1
!          sum_val = b(i)
!          do j = i + 1, n
!             sum_val = sum_val - a(i, j) * b(j)
!          end do
!          b(i) = sum_val / a(i, i)
!       end do
! 
!    end subroutine lubksb
! 
!    subroutine optinewton(maxit, nts, sm, delta, rho, x, h, indx, in)
!       !! Calculate non-bonded fraction using Newton-Raphson optimization
!       !! Maximizes Q = A/RT|assoc,eq
!       !! Based on Michelsen & Mollerup approach
!       integer, intent(in) :: maxit, nts
!       real(pr), intent(in out) :: sm(nts), rho
!       real(pr), intent(in) :: delta(nts, nts)
!       real(pr), intent(inout) :: x(nts)
!       real(pr), intent(out) :: h(nts, nts)
!       integer, intent(out) :: indx(nts), in
! 
!       integer :: i, im
!       real(pr) :: q, qnew, alpha, d
!       real(pr) :: xnew(nts), xold(nts), dx(nts)
!       real(pr) :: g(nts), gnew(nts), hnew(nts, nts)
! 
!       ! Initial calculation
!       call qfunction(nts, x, sm, rho, delta, q, g, h)
!       im = 0
! 
!       ! Convergence loop
!       do in = 1, maxit
!          im = im + 1
!          dx = -g
! 
!          ! Solve linear system
!          call ludcmp(h, nts, nts, indx, d)
!          call lubksb(h, nts, nts, indx, dx)
! 
!          ! Initialize step size
!          alpha = 1.0_pr
! 
!          ! Calculate new X
!          do
!             do i = 1, nts
!                xnew(i) = x(i) + alpha * dx(i)
! 
!                ! If fraction is negative, reduce but don't make zero
!                if (xnew(i) <= 0.0_pr) then
!                   xnew(i) = x(i) / 5.0_pr
!                end if
!             end do
! 
!             ! Exit if step is small enough
!             if (maxval(abs(alpha * dx)) <= 1.0e-15_pr) exit
! 
!             ! Calculate new Q
!             call qfunction(nts, xnew, sm, rho, delta, qnew, gnew, hnew)
! 
!             ! Did Q increase?
!             if (qnew > q * (1.0_pr + 1.0e-14_pr)) then
!                q = qnew
!                x = xnew
!                g = gnew
!                h = hnew
!                exit
!             else
!                ! Reduce step size
!                alpha = alpha / 3.0_pr
!                im = im + 1
!             end if
!          end do
! 
!          ! Check convergence
!          if (maxval(abs(alpha * dx)) <= 1.0e-15_pr) exit
!       end do
! 
!    end subroutine optinewton
! 
!    subroutine qfunction(nts, x, sm, rho, delta, q, g, h)
!       !! Calculate Q function, gradient, and Hessian
!       !! Q = Sum(sm(k)·(ln X(k) - X(k) + 1)) -
!       !!     1/2·Sum(Sum(sm(k)·sm(l)·Delta(k,l)·X(l)/V))
!       integer, intent(in) :: nts
!       real(pr), intent(in out) :: x(nts), sm(nts), rho
!       real(pr), intent(in) :: delta(nts, nts)
!       real(pr), intent(out) :: q, g(nts), h(nts, nts)
! 
!       integer :: i, j
!       real(pr) :: q1, q2, suma
!       real(pr), parameter :: macheps = epsilon(1.0_pr)
! 
!       ! Initialize
!       q1 = 0.0_pr
!       q2 = 0.0_pr
!       g = 0.0_pr
!       h = 0.0_pr
! 
!       ! Calculate Q1, gradient and Hessian
!       do i = 1, nts
!          if (sm(i) > 2.0_pr * macheps) then
!             suma = 0.0_pr
! 
!             do j = 1, nts
!                if (delta(j, i) > 0.0_pr .and. sm(j) > 2.0_pr * macheps) then
!                   h(j, i) = -sm(i) * sm(j) * delta(j, i) * rho
!                   suma = suma - h(j, i) * x(j)
!                   q2 = q2 + h(j, i) * x(i) * x(j) / 2.0_pr
!                end if
!             end do
! 
!             g(i) = suma
!             h(i, i) = h(i, i) - (sm(i) + suma) / x(i)
! 
!             ! Avoid null elements on diagonal
!             if (h(i, i) == 0.0_pr) h(i, i) = 1.0e-20_pr
! 
!             q1 = q1 + sm(i) * (log(x(i)) - x(i) + 1.0_pr)
!             g(i) = sm(i) / x(i) - sm(i) - g(i)
!          else
!             ! Trace site, deactivated but has its own X at infinite dilution
!             h(i, i) = 1.0_pr
!             x(i) = 1.0_pr / (1.0_pr + sum(rho * sm(:nts) * x(:nts) * delta(:nts, i)))
!          end if
!       end do
! 
!       ! Total Q
!       q = q1 + q2
! 
!    end subroutine qfunction
! 
!    real(pr) function gca_get_b(self, n, P, T) result(b_gc)
!       !! Calculate the covolume parameter b_gc
!       class(GCAEOS), intent(in) :: self
!       real(pr), intent(in) :: n(:)
!       real(pr), intent(in) :: P, T
! 
!       real(pr) :: dia, dla
! 
!       real(pr) :: lambda1, lambda2, lambda3
! 
!       integer :: i, nc
! 
!       nc = size(n)
! 
!       ! Calculate lambda parameters
!       lambda1 = 0.0_pr
!       lambda2 = 0.0_pr
!       lambda3 = 0.0_pr
! 
!       do i = 1, nc
!          dia = self%d(i)
!          dla = n(i) * dia
!          lambda1 = lambda1 + dla
!          dla = dla * dia
!          lambda2 = lambda2 + dla
!          lambda3 = lambda3 + dla * dia
!       end do
! 
!       ! Net co-volume
!       b_gc = PI * lambda3 / 6.0_pr
!    end function gca_get_b
! 
! end module yaeos__models_ar_gc_gca
! 
! 