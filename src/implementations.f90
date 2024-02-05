module yaeos_models_ar_cubic_implementations
   !! Implemented Cubic Equations of State.
   !!
   !! - PengRobinson76
   !! - PengRobinson78
   !! - SoaveRedlichKwong
   use yaeos_models_ar_cubic_pengrobinson76, only: PengRobinson76
   use yaeos_models_ar_cubic_pengrobinson78, only: PengRobinson78
   use yaeos_models_ar_cubic_SRK, only: SoaveRedlichKwong
end module
