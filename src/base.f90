module yaeos__models_base
    !! Basic element of a thermodynamic model.
    use yaeos__substance, only: Substances
    implicit none

    type, abstract :: BaseModel
        !! Base model type.
        !!
        !! Contains the important parts of most models and other procedures.
        type(Substances) :: components 
            !! Substances contained in the module
    end type
end module