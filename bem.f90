module bem
  use, intrinsic :: iso_fortran_env
  use subprogram
  implicit none
contains
  subroutine calc_bem(div_num, inner_point, result_value)
    !! 境界要素法の計算を行う
    implicit none
    !> 円周の分割数（境界要素数, 点の数）
    integer(int32), intent(in) :: div_num
    !> 空間内点の座標
    real(real64), intent(in) :: inner_point(2)
    !> 結果の値
    real(real64), intent(inout) :: result_value

    integer(int32) :: i, m, n, info
    !> 小林本のUの値を並べる行列
    real(real64), allocatable :: capital_u_mat(:, :)
    !> 小林本のWの値を並べる行列
    real(real64), allocatable :: capital_w_mat(:, :)
    !> 境界上の法線微分の値を並べる行列
    real(real64), allocatable :: q(:, :)
    !> 境界上の値を並べる行列
    real(real64), allocatable :: u(:, :)
    !> 解析的な解と比べるための厳密な値
    real(real64), allocatable :: exact_q(:, :)
    !> 円周上の点列. 端点であり, 代表点（`points`）はそれぞれの要素の中点.
    real(real64), allocatable :: edge_points(:, :)
    !> 要素を表す点列
    real(real64), allocatable :: points(:, :)

    allocate (edge_points(div_num, 2))
    allocate (points(div_num, 2))
    allocate (capital_u_mat(div_num, div_num))
    allocate (capital_w_mat(div_num, div_num))
    allocate (q(div_num, 1))
    allocate (exact_q(div_num, 1))
    allocate (u(div_num, 1))

    ! 円周上の点を用意する. (1,0)から始めて一周する.
    do i = 1, div_num
      edge_points(i, 1) = cos(2.0d0*real(i - 1, real64)*PI/real(div_num, real64))
      edge_points(i, 2) = sin(2.0d0*real(i - 1, real64)*PI/real(div_num, real64))
    end do
    ! 代表点
    do i = 1, div_num
      points(i, :) = (edge_points(i, :) + edge_points(modulo(i, div_num) + 1, :))/2.0d0
    end do

    ! 行列の各要素を計算
    do m = 1, div_num
      do n = 1, div_num
        capital_u_mat(m, n) = U_component(m, n, edge_points)
        capital_w_mat(m, n) = W_component(m, n, edge_points)
      end do
      u(m, 1) = exact_u(points(m, :))
      ! 確認のため厳密なqも計算する.
      exact_q(m, 1) = exact_u_normal_drv(points(m, :))
    end do

    ! 求めた行列を用いて共役な法線微分：qを計算.
    call solve(capital_u_mat, matmul(capital_w_mat, u), q, info)
    ! if (info == 0) then
    !   print *,"sum of error about q :",dot_product(q(:,1)-exact_q(:,1), q(:,1)-exact_q(:,1))
    ! end if

    ! 台形公式で内点の値を計算
    result_value = integral_trapezoidal(inner_point, points, q, u)

    print '("at point [" e20.7 "," ,e20.7 " ],", " result : " e20.7)', inner_point(1), inner_point(2), result_value

  end subroutine calc_bem
end module bem
