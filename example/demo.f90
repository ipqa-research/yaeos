program examples
    use bench, only: benchmarks => main
    use hyperdual_pr76, only: run_hyperdual_pr76 => main
    use flashing, only: run_flashes => main
    use TapeRobinson, only: run_tape_pr76 => main

    print *, "Running Tapenade generated PR76"
    call run_tape_pr76
    print *, "Running Hyperdual generated PR76"
    call run_hyperdual_pr76
    print *, "Running bencharks O(f(N))"
    call benchmarks
    print *, "Flash example"
    call run_flashes
end program