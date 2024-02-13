program examples
    use bench, only: benchmarks
    use hyperdual_pr76, only: adiff_pr76
    use flashing, only: run_flashes
    use TapeRobinson, only: run_tape => main

    call run_tape
    call benchmarks
    ! call adiff_pr76
    ! call run_flashes
    ! call run_tape
end program