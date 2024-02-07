module yaeos_equilibria_saturation_points
    use yaeos_constants, only: pr
    use yaeos_models, only: ArModel
    use yaeos_thermoprops, only: fugacity_vt
    use yaeos_equilibria_equilibria_state, only: EquilibriaState

    real(pr) :: tol = 1e-5_pr
    integer :: max_iterations = 100
    real(pr) :: step_tol = 0.1_pr
end module