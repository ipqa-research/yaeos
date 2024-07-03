module yaeos__consistency
   !! # yaeos__consistency
   !! Subroutine to evaluate the consistency of thermodynamic models.
   !!
   !! # Description
   !! Tools to evaluate the consistency of \(A^r\) and \(G^E\) models. This
   !! module also provides subroutines for numerical evaluations of \(A^r\) and
   !! \(G^E\) derivatives using central finite differences. The purpose of the
   !! module is to assist in the development of new models and ensure the
   !! accuracy of the derivatives implementation.
   !!
   !! # Examples
   !! For detailed explanations and examples of each consistency test, please
   !! refer to the API documentation of each submodule.
   !!
   !! - \(A^r\) consistency tests: [[yaeos__consistency_armodel]]
   !! - \(G^E\) consistency tests: [[yaeos__consistency_gemodel]]
   !!
   !! # References
   !! 1. Michelsen, M. L., & Mollerup, J. M. (2007). Thermodynamic models:
   !! Fundamentals & computational aspects (2. ed). Tie-Line Publications.
   !!
   ! Consistency test for ArModels
   use yaeos__consistency_armodel
   use yaeos__consistency_gemodel
end module yaeos__consistency
