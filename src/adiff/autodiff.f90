module yaeos__autodiff
   !! # Automatic-differentiation interfaces
   !! Re-exports the hyperdual number implementation and the `ArModelAdiff`
   !! abstract type used to build residual Helmholtz models with operator
   !! overloading.
   !!
   !! # Description
   !! Extend `ArModelAdiff`, implement the model's hyperdual residual Helmholtz
   !! function and volume initializer, then use the inherited `ArModel`
   !! properties and phase-equilibrium routines. Tapenade-based interfaces are
   !! provided by separate API modules.
   !!
   !! # Examples
   !! See `example/extra/adiff_pr76.f90` for a complete model implementation.
   use hyperdual_mod
   use yaeos__adiff_hyperdual_ar_api, only: ArModelAdiff
   implicit none
end module
