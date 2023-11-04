
program main
  use, intrinsic :: iso_fortran_env
  use bem, only: calc_bem
  implicit none

  integer(int32) :: div_num
  integer(int32), parameter :: MESH_DIV_NUM = 40
  integer(int32), parameter :: fo = 11
  real(real64) :: inner_point(2)
  real(real64) :: result_value
  real(real64), allocatable :: result_array(:)
  real(real64), allocatable :: inner_points_xy(:)  ! 領域内の点を並べて保存しておくための一次元配列
  real(real64), allocatable :: inner_points(:, :)
  integer(int32) :: i, j, point_num
  ! CPU時間計測用
  integer(int32)   :: begin_time, end_time, time_count_per_sec, time_count_max

  inner_points_xy = [real(real64) ::]
  point_num = 0

  do i = 1, MESH_DIV_NUM + 1
    do j = 1, MESH_DIV_NUM + 1
      ! 順番に点を考え, 領域内に入っていればinner_points_xy配列に追加する. point_numはそのような点の数
      inner_point(:) = [2.0d0/(real(MESH_DIV_NUM, real64))*(j - 1), 2.0d0/(real(MESH_DIV_NUM, real64))*(i - 1)] - [1.0d0, 1.0d0]
      if (dot_product(inner_point, inner_point) < 1.0d0) then
        point_num = point_num + 1
        inner_points_xy = [inner_points_xy, inner_point]
      end if
    end do
  end do
  ! 列挙し終わったらうまくreshapeして平面の点の配列になるようにする（二次元配列）.
  inner_points = reshape(inner_points_xy, [point_num, 2], order=[2, 1])

  ! 円周の分割サイズを入力（計算精度）
  print *, "input div_num : "
  READ (*, *) div_num
  if (div_num <= 0) stop

  print *, "start to measure time"
  call system_clock(begin_time, time_count_per_sec, time_count_max)
  ! 結果を保存する配列を用意
  allocate (result_array(point_num))
  do i = 1, point_num
    call calc_bem(div_num, inner_points(i, :), result_value)
    ! 厳密値との差を見るか実際の値を見るかいずれか
    ! result_array(i) = result_value - exact_u(inner_points(i,:))
    result_array(i) = result_value
  end do

  print *, "sum of error : ", sum(result_array)
  call system_clock(end_time)
  print *, "time measurement end. elapsed time: ", real(end_time - begin_time)/time_count_per_sec, "sec"

  open (fo, file='plot.dat')
  do i = 1, point_num
    write (fo, '(3e20.7)') inner_points(i, 1), inner_points(i, 2), result_array(i)
  end do
  close (fo)
  ! なぜかVScodeのexecタスクで実行するとファイルが生成されないのでターミナルから実行する.
  ! plot.datができたら, gnuplotを起動し[sp "plot.dat" w p]と入力するとプロットされる.

end program main
