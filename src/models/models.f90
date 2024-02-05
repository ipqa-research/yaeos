module yaeos_models
    use yaeos_models_ar, only: ArModel
    use yaeos_models_ar_genericcubic, only: CubicEoS, GenericCubic_Ar
    use yaeos_models_ar_cubic_alphas, only: AlphaSoave
    use yaeos_models_ar_genericcubic_quadratic_mixing, only: QMR

    use yaeos_models_ar_cubic_pengrobinson76, only: PengRobinson76
    use yaeos_models_ar_cubic_pengrobinson78, only: PengRobinson78
    use yaeos_models_ar_cubic_srk, only: SoaveRedlichKwong
end module
