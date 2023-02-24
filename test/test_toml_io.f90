!program main
!    use toml_io, only: read_system
!    use models, only: ArModel, CubicEOS
!    use mixing_rules
!
!    type(CubicEOS) :: model
!
!    call read_system("input.toml", model)
!
!    print *, "MODEL", " ", "MIXING_RULE"
!    print *, model%thermo_model, "  ", model%mixing_rule
!    print *, "=========================================="
!    print *, "COMPONENT", " ", "z_i", " ", "Tc(i)"
!    do i=1,size(model%names)
!        print *, model%names(i), model%z(i), model%tc(i)
!    end do
!
!end program main