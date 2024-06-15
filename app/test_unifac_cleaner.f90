program main
   use iso_fortran_env, only: int64
   use yaeos, only: pr
   use yaeos__models_ge_group_contribution_unifac
   use stdlib_io_npy, only: load_npy
   implicit none

   integer :: i

   type(UNIFAC) :: model
   integer, parameter :: nc = 3, ng = 4
   real(pr) :: x(nc) = [0.2, 0.7, 0.1], T=150

   ! integer, parameter :: nc = 2, ng = 3
   ! real(pr) :: x(nc) = [0.3, 0.7], T=150

   real(pr), allocatable :: Aij(:, :)
   real(pr), allocatable :: Qk(:), Rk(:)
   real(pr) :: dx=0.01, dpsidt_num(ng, ng), dpsidt(ng, ng)

   type(Groups) :: molecules(nc)

   real(pr) :: Ge, Ge_dt, Gen(nc), Gen2(nc,nc), GeT, GeT2, GeTn(nc)
   real(pr) :: d2GedT2, dGedn(nc), d2GedndT(nc), d2Gedn2(nc,nc)

   integer(int64) :: rate, st, et
   call system_clock(count_rate=rate)

   call load_npy("data/unifac_aij.npy", Aij)
   call load_npy("data/unifac_Qk.npy", Qk)
   call load_npy("data/unifac_Rk.npy", Rk)

   ! ! Ethane [CH3]
   molecules(1)%groups_ids = [1]
   molecules(1)%number_of_groups = [2]

   ! ! Ethanol [CH3, CH2, OH]
   molecules(2)%groups_ids = [1, 2, 14]
   molecules(2)%number_of_groups = [1, 1, 1]

   ! ! Methylamine [H3C-NH2]
   molecules(3)%groups_ids = [28]
   molecules(3)%number_of_groups = [1]


   model = setup_unifac(molecules)

   call excess_gibbs(model, x, T, Ge, GeT, GeT2, Gen, GeTn, Gen2)
   call excess_gibbs(model, x, T + dx, Ge_dt)

   print *, "Ge: "
   print *, Ge
   print *, "Thermo: "
   print *, -3.223992676822129
   print *, " "

   print *, "Gen:"
   print *, Gen
   print *, "Thermo: "
   print *, [10.53032277,  -2.37758326, -36.65748951]
   print *, " "

   print *, "ln_gammas:"
   print *, Gen / R / T
   print *, "Thermo: "
   print *, [0.84433781, -0.19063836, -2.93925506]
   print *, " "

   print *, "Gen2"
   print *, Gen2(1,:)
   print *, Gen2(2,:)
   print *, Gen2(3,:)
   print *, "Thermo: "
   print *, [-0.75249927,  0.13440904,  0.56413529] * R * T
   print *, [ 0.13440904,  0.34708386, -2.69840507] * R * T
   print *, [ 0.56413529, -2.69840507, 17.76056492] * R * T
   print *, " "

   print *, "dln_gammas_dn"
   print *, Gen2(1,:) / R / T
   print *, Gen2(2,:) / R / T
   print *, Gen2(3,:) / R / T
   print *, "Thermo: "
   print *, [-0.752499273, 0.134409037, 0.564135287]
   print *, [0.13440903599999998, 0.34708385599999997, -2.698405064]
   print *, [0.5641352889999995, -2.6984050710000003, 17.760564919]
   print *, " "

   print *, "GeT"
   print *, GeT
   print *, "numeric: "
   print *, (Ge_dt - Ge) / dx
   print *, "Thermo: "
   print *, 0.03268447167877294
   print *, " "

   print *, "GeT2"
   print *, GeT2
   print *, "Thermo: "
   print *, -0.0003594405355829625
   print *, " "

   print *, "GeTn"
   print *, GeTn
   print *, "Thermo: "
   print *, [0.06015389, 0.02239722, 0.04975642]
end program main

