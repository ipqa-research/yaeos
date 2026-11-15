module yaeos__models_ge_implementations
   !! # Excess Gibbs model constructors
   !!
   !! # Description
   !! Convenience module that re-exports the excess Gibbs models available in
   !! `yaeos`: [[NRTL]], [[NRTLHV]], [[UNIQUAC]], [[UNIFAC]], and the
   !! Dortmund and PSRK variants.
   !! This is a convenience re-export module; each concrete model's constructor
   !! and parameter requirements are documented in its implementation module.
   !!
   !! # Examples
   !! ```fortran
   !! use yaeos__models_ge_implementations, only: NRTL, setup_unifac, Groups
   !! ```
   use yaeos__models_ge_NRTL, only: NRTL
   use yaeos__models_ge_nrtlhv, only: NRTLHV
   use yaeos__models_ge_group_contribution_dortmund, only: setup_dortmund
   use yaeos__models_ge_group_contribution_unifac, only: &
      Groups, setup_unifac, UNIFAC, excess_gibbs
   use yaeos__models_ge_uniquac, only: setup_uniquac, UNIQUAC
   use yaeos__models_ge_group_contribution_psrk, only: setup_psrk

   implicit none
end module yaeos__models_ge_implementations
