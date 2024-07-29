module yaeos__io_namelist
   !! File IO to set up models from namelist files.
   use yaeos__models_ar_genericcubic, only: CubicMixRule
   use yaeos, only: &
      pr, &
      ArModel, &
      CubicEoS , PengRobinson76, PengRobinson78, SoaveRedlichKwong, RKPR, &
      QMR, MHV
   implicit none

contains

   function model_from_namelist(filepath) result(model)
      character(len=*), intent(in) :: filepath
      class(ArModel), allocatable :: model

      integer :: nc
      character(len=50) :: model_name

      call read_setup(filepath, nc, model_name)

      select case(model_name)
       case("SoaveRedlichKwong", "PengRobinson76", "PengRobinson78", "RKPR")
         model = read_cubic_eos(filepath, nc, model_name)
      end select
   end function model_from_namelist

   type(CubicEoS) function read_cubic_eos(filepath, nc, model_name) result(model)
      character(len=*), intent(in) :: filepath
      integer, intent(in) :: nc
      character(len=*), intent(in) :: model_name
      character(len=50) :: mixrule
      class(CubicMixRule), allocatable :: mix

      real(pr) :: Tc(nc), Pc(nc), w(nc), Zc(nc)

      integer :: fu

      namelist /yaeos_cubic/ mixrule

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
       case default
         error stop 1
      end select

      open(newunit=fu, file=filepath)
      read(fu, nml=yaeos_cubic)
      close(fu)

      select case (mixrule)
       case("QMR")
         mix = QMR_from_namelist(filepath, nc)
         deallocate(model%mixrule)
         model%mixrule = mix
       case default
         error stop 1
      end select
   end function read_cubic_eos

   type(QMR) function QMR_from_namelist(filepath, nc)
      character(len=*), intent(in) :: filepath
      integer, intent(in) :: nc

      integer :: fu

      real(pr) :: kij(nc, nc), lij(nc, nc)

      namelist /yaeos_qmr/ kij, lij

      kij = 0
      lij = 0
      open(newunit=fu, file=filepath)
      read(fu, nml=yaeos_qmr)
      close(fu)

      QMR_from_namelist = QMR(kij, lij)
   end function QMR_from_namelist

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

      namelist /yaeos_critical/ Tc, Pc, w, Zc

      Tc = 0
      Pc = 0
      w = 0
      Zc = 0

      open(newunit=fu, file=filepath)
      read(fu, nml=yaeos_critical)
      close(fu)
   end subroutine read_critical
end module yaeos__io_namelist
