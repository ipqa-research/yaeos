module tapenade_interfaces
    implicit none

    interface pushreal8
        subroutine pushreal8(realnum)
            real(8) :: realnum
        end subroutine
        
        subroutine pushreal8array(realnum)
            real(8) :: realnum(:)
        end subroutine
    end interface

    interface popreal8
        subroutine popreal8(realnum)
            real(8) :: realnum
        end subroutine
    end interface

    interface pushinteger4
        subroutine pushinteger4(intnum)
            integer :: intnum
        end subroutine
    end interface

    interface popinteger4
        subroutine popinteger4(intnum)
            integer :: intnum
        end subroutine
    end interface

    interface pushcontrol1b
        subroutine pushcontrol1b(intnum)
            integer :: intnum
        end subroutine
    end interface
end module