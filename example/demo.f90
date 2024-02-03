program examples
    use bench, only: benchmarks
    use hyperdual_pr76, only: adiff_pr76

    call benchmarks
    call adiff_pr76
end program