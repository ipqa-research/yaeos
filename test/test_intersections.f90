program main
   use yaeos__math, only: Point, intersect_one_line
   use yaeos__constants, only: pr
   use testing_aux, only: test_title, assert
   implicit none
   type(Point), allocatable :: intersections(:)
   real(pr), allocatable :: lx(:), ly(:)

   write(*, *) test_title("Testing Line Segment Intersections")
   lx = [1, 2, 3, 4, 5, 6, 7, 7, 6, 5, 4]
   ly = [2, 3, 4, 5, 6, 7, 8, 9, 8, 6, 4]
   intersections = intersect_one_line(lx, ly)

   call assert(size(intersections) == 1, "One intersection expected")
   call assert(intersections(1)%x == 5.0_pr, "Intersection x-coordinate should be 5")
   call assert(intersections(1)%y == 6.0_pr, "Intersection y-coordinate should be 6")
end program