program examples
    ! use bench, only: benchmarks => main
    use hyperdual_pr76, only: run_hyperdual_pr76 => main
    use flashing, only: run_flashes => main
    use tape_nrtl, only: run_tape_nrtl => main
    use bench, only: benchmarks => main

    print *, "========================================="
    print *, "YAEOS DEMO"
    ! print *, "Running Tapenade generated NRTL model"
    ! call run_tape_nrtl
    ! print *, "Running Hyperdual generated PR76"
    ! call run_hyperdual_pr76
     print *, "Running bencharks O(f(N))"
     call benchmarks
    ! print *, "Flash example"
    ! call run_flashes
    print *, "========================================="
end program