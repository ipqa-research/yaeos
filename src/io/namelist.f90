module yaeos__io_namelist
   !! File IO to set up models from namelist files.
   use yaeos, only: &
      pr, &
      ArModel, &
      CubicEoS, PengRobinson76, PengRobinson78, SoaveRedlichKwong, RKPR, &
      QMR, MHV
   implicit none

contains

   type(CubicEoS) function read_cubic_eos(filepath, nc, model_name) result(model)
      character(len=*), intent(in) :: filepath
      integer, intent(in) :: nc
      character(len=*), intent(in) :: model_name

      real(pr) :: Tc(nc), Pc(nc), w(nc), Zc(nc)

      call read_critical(filepath, Tc, Pc, w, Zc)

      select case (model_name)
       case("SoaveRedlichKwong")
         model = SoaveRedlichKwong(Tc, Pc, w)
       case("PengRobinson76")
         model = PengRobinson76(Tc, Pc, w)
       case("PengRobinson78")
         model = PengRobinson78(Tc, Pc, w)
       case("RKPR")
         model = RKPR(Tc, Pc, w, Zc)
      end select

      call read_cubic_alpha_mixrule(filepath, nc, model)
   end function read_cubic_eos

   subroutine set_cubic_alpha_mixrule(filepath, nc, model)
      character(len=*), intent(in) :: filepath
      integer, intent(in) :: nc
      type(CubicEoS), intent(in out) :: model

      character(len=50) :: alpha, mixrule
      namelist /yaeos_cubic/ alpha, mixrule

      integer :: fu

      alpha = "AlphaSoave"
      mixrule = "QMR"

      open(newunit=fu, file=filepath)
      read(fu, nml=yaeos_cubic)
      close(fu)
   end subroutine set_cubic_alpha_mixrule

   subroutine read_setup(filepath, nc, model)
      character(len=*), intent(in) :: filepath
      integer, intent(out) :: nc
      character(len=50), intent(out) :: model

      namelist /yaeos_setup/ nc, model

      integer :: fu

      open(newunit=fu, file=filepath)
      read(fu, nml=yaeos_setup)
      close(fu)
   end subroutine read_setup

   subroutine read_critical(filepath, Tc, Pc, w, Zc)
      character(len=*), intent(in) :: filepath
      real(pr), intent(out) :: Tc(:), Pc(:), w(:), Zc(:)
      integer :: fu

      namelist /nml_critical/ Tc, Pc, w, Zc

      Tc = 0
      Pc = 0
      w = 0
      Zc = 0

      open(newunit=fu, file=filepath)
      read(fu, nml=nml_critical)
      close(fu)
   end subroutine read_critical
end module yaeos__io_namelist
