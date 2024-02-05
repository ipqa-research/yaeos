module yaeos_models
    !! Ya_EoS thermodynamic models

    ! Residual Helmholtz Models
    use yaeos_models_ar, only: ArModel

    ! Cubic EoS models
    use yaeos_models_ar_genericcubic, only: CubicEoS, GenericCubic_Ar
    use yaeos_models_ar_cubic_alphas, only: AlphaSoave
    use yaeos_models_ar_genericcubic_quadratic_mixing, only: QMR

    use yaeos_models_ar_cubic_pengrobinson76, only: PengRobinson76
end module
