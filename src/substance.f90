module yaeos_substance
    use yaeos_constants, only: pr

    type :: Substances
        !! A set of substances
        real(pr), allocatable :: tc(:), pc(:), w(:)
    end type
end module