module yaeos__models_ge_implementations
   use yaeos__models_ge_NRTL, only: NRTL
   use yaeos__models_ge_group_contribution_dortmund, only: setup_dortmund
   use yaeos__models_ge_group_contribution_unifac, only: &
      Groups, setup_unifac, UNIFAC, excess_gibbs
   use yaeos__models_ge_uniquac, only: setup_uniquac, UNIQUAC
   use yaeos__models_ge_group_contribution_psrk, only: setup_psrk

   implicit none
end module yaeos__models_ge_implementations
