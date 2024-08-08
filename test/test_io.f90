program test_io
   use yaeos, only: ArModel, CubicEoS, QMR
   use yaeos__io, only: model_from_namelist
   implicit none

   class(ArModel), allocatable :: model
   model = model_from_namelist("test/data/model.nml")

   select type(model)
   type is (CubicEoS)
      print *, model%components%Tc
      print *, model%components%Pc
      print *, model%components%w

      associate (mr => model%mixrule)
         select type(mr)
         type is (QMR)
            print *, mr%k(1, :)
            print *, mr%k(2, :)
            print *, mr%k(3, :)
         end select
      end associate

   end select
end program