module bem
  use, intrinsic :: iso_fortran_env
  use subprogram
  implicit none
contains
  subroutine calc_bem(DIV_NUM, inner_point, result_value)
    implicit none
    integer(int32), intent(in) :: DIV_NUM
    real(real64), intent(in) :: inner_point(2)
    real(real64), intent(out) :: result_value
    integer(int32) :: i, m, n, info
    real(real64), allocatable :: lU(:, :), lW(:, :), q(:, :), u(:, :), exact_q(:, :), edge_points(:, :), points(:, :)

    ! 割り付け
    allocate (edge_points(DIV_NUM, 2))
    allocate (points(DIV_NUM, 2))
    allocate (lU(DIV_NUM, DIV_NUM))
    allocate (lW(DIV_NUM, DIV_NUM))
    allocate (q(DIV_NUM, 1))
    allocate (exact_q(DIV_NUM, 1))
    allocate (u(DIV_NUM, 1))

    ! 円周上の点を用意する. (1,0)から始めて一周する. これは端点で, 代表点はそれぞれの中点とする.
    do i = 1, DIV_NUM
      edge_points(i, 1) = cos(2.0d0*real(i - 1, real64)*PI/real(DIV_NUM, real64))
      edge_points(i, 2) = sin(2.0d0*real(i - 1, real64)*PI/real(DIV_NUM, real64))
    end do
    ! 代表点
    do i = 1, DIV_NUM
      points(i, :) = (edge_points(i, :) + edge_points(modulo(i, DIV_NUM) + 1, :))/2.0d0
    end do

    ! 行列の各要素を計算
    do m = 1, DIV_NUM
      do n = 1, DIV_NUM
        lU(m, n) = U_component(m, n, edge_points)
        lW(m, n) = W_component(m, n, edge_points)
      end do
      u(m, 1) = exact_u(points(m, :))
      ! 確認のため厳密なqも計算する.
      exact_q(m, 1) = exact_u_normal_drv(points(m, :))
    end do

    ! 求めた行列を用いて共役な法線微分（q）を計算.
    call solve(lU, matmul(lW, u), q, info)
    if (info == 0) then
      ! print *,"sum of error about q :",dot_product(q(:,1)-exact_q(:,1), q(:,1)-exact_q(:,1))
    end if

    ! 台形公式で内点の値を計算
    result_value = integral_trapezoidal(inner_point, points, q, u)

    print '("at (" e20.7 "," ,e20.7 "),", "result : " e20.7)', inner_point(1), inner_point(2), result_value
    ! print *,"exact value is ",exact_u(inner_point(:))," ... error is ",exact_u(inner_point(:))-result_value

  end subroutine calc_bem
end module bem
