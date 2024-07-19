module yaeos__models
   !! `yaeos` thermodynamic models
   !!
   !! On `yaeos` there are implemented a series of both residual Helmholtz
   !! energy (\(A_r\)) and excess Gibbs energy (\(G^E\)) models.
   !!
   !! This module takes all the relevant procedures and derived types
   !! related to them.
   !!
   !! - Residual Helmholtz model base type `ArModel` base derived type
   !!   that provides the basic structure that a residual Helmholtz model
   !!   should provide.
   !! - **Cubic Equations of state**:
   !!    - `AlphaFunction` type
   !!    - `CubicMixRule` type
   !!    - `CubicEos` type that extends `ArModel` to use a generic
   !!      two-parameter EoS. Implemented models that use this type can be
   !!      seen at [[yaeos__models_ar_cubic_implementations(module)]]
   !!    - `QMR` (Quadratic Mixing Rule) type: extensible derived type that 
   !!       defaults to classic vdW mixing rules
   !!

   ! Base model structure
   use yaeos__models_base, only: BaseModel

   ! Residual Helmholtz Models
   use yaeos__models_ar, only: ArModel, size

   ! Cubic EoS models
   use yaeos__models_ar_genericcubic, only: &
      CubicEoS, GenericCubic_Ar, AlphaFunction, CubicMixRule

   ! Alpha functions
   use yaeos__models_ar_cubic_alphas

   ! Mixing Rules
   use yaeos__models_ar_cubic_quadratic_mixing
   use yaeos__models_cubic_mixing_rules_huron_vidal

   ! Implemented models
   use yaeos__models_ar_cubic_implementations

   ! Ge Models
   use yaeos__models_ge, only: GeModel

   ! Implemented models
   use yaeos__models_ge_implementations

end module
