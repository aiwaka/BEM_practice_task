module subprogram
  use, intrinsic :: iso_fortran_env
  use constants, only: PI
  implicit none

contains
  subroutine solve(A, b, x, info)
    !! Ax=bの方程式を解き, xに結果を代入し, infoにステータスを保存する.

    !> 係数行列
    real(real64), intent(in) :: A(:, :)
    !> 定数ベクトル
    real(real64), intent(in) :: b(:, :)
    !> 解を保存するベクトル
    real(real64), allocatable, intent(inout) :: x(:, :)
    !> 計算結果ステータスを保存する
    integer(int32), intent(inout) :: info

    integer(int32), allocatable :: ipiv(:)
    integer(int32) :: an, am, bn, bm
    ! Aが正方であるかチェック
    an = size(A, 1)
    am = size(A, 2)
    if (an /= am) then
      print *, "Error : solve (an not equal to am) : an=", an, ", am=", am
      return
    end if
    ! 基本的にbmは1であることを想定する.
    bn = size(b, 1)
    bm = size(b, 2)

    if (allocated(x)) deallocate (x)
    allocate (x, source=b)
    ! ピボット行列（入れ替えた結果を保存する行列. Aと同サイズの一次元配列とする）
    allocate (ipiv(an))

    call dgesv(an, bm, A, an, ipiv, x, bn, info)  ! lapackのサブルーチンを呼び出して解く.
    if (info /= 0) then
      ! infoが0でなければbad statusなのでその旨を表示する
      print *, "Error : solve : info=", info
      return
    end if
  end subroutine

  function integral_trapezoidal(x, points, q_array, u_array) result(retval)
    !! 台形公式を用いて周回積分を行う.
    !! 添字1とN+1の点が一致していることを利用して,
    !! すべての点を一回ずつ足してh/2をかける代わりにhをかけることにする.

    implicit none

    !> 計算する内点座標
    real(real64) :: x(2)
    !> 点列. index: (axis, point)
    real(real64) :: points(:, :)
    !> 境界上の法線微分の値の列. index: (point, [1])
    real(real64) :: q_array(:, :)
    !> 境界上の値の列. index: (point, [1])
    real(real64) :: u_array(:, :)

    real(real64) :: retval

    real(real64) :: h
    integer(int32) :: points_num, i

    retval = 0.0d0
    points_num = size(points, 2)
    h = 2.0d0*PI/real(points_num, real64)
    do i = 1, points_num
      ! 一重層ポテンシャル
      ! q_arrayは添字1を書いているが配列の大きさ1なので特に意味はない.
      retval = retval + fund_gamma(x, points(:, i))*q_array(i, 1)
      ! 二重層ポテンシャル
      retval = retval - fund_gamma_derivative(x, points(:, i))*u_array(i, 1)
    end do
    retval = retval*h
  end function integral_trapezoidal

  function exact_u(x) result(retval)
    !! 与えられた点に対して関数x^3 - 3xy^2の厳密な値を計算する.
    implicit none
    real(real64) :: x(2)
    real(real64) :: retval

    retval = x(1)**3 - 3*x(1)*x(2)**2
  end function exact_u

  function exact_u_normal_drv(x) result(retval)
    !! 与えられた境界上の点に対して関数x^3-3xy^2の法線微分を計算する.
    !! 境界は半径1の単位円とする.
    implicit none
    real(real64) :: x(2)
    real(real64) :: retval
    real(real64) :: theta

    theta = atan2(x(2), x(1))
    retval = 3*(x(1)**2 - x(2)**2)*cos(theta) - 6*x(1)*x(2)*sin(theta)
  end function exact_u_normal_drv

  function fund_gamma(x, y) result(retval)
    !! 任意の二点に対する（二次元）基本解の値を返す.
    implicit none
    real(real64) :: x(2), y(2)
    real(real64) :: retval

    retval = -log(sqrt(dot_product(x - y, x - y)))/2.0d0/PI
  end function fund_gamma

  function fund_gamma_derivative(x, y) result(retval)
    !! 任意の二点（ただしyは境界上）に対する基本解のyについての法線微分を返す.
    implicit none
    real(real64) :: x(2), y(2), theta
    real(real64) :: retval

    theta = atan2(y(2), y(1))
    retval = dot_product(x - y, [cos(theta), sin(theta)])/2.0d0/PI/dot_product(x - y, x - y)
  end function fund_gamma_derivative

  subroutine component_values(m, n, edge_points, capital_x1, capital_x2, capital_y1, capital_y2, r1, r2, h, theta)
    !! 小林本の表記に従った諸量を計算する.
    implicit none

    integer(int32), intent(in) :: m, n
    real(real64), intent(in) :: edge_points(:, :)
    real(real64), intent(inout) :: capital_x1, capital_x2, capital_y1, capital_y2, r1, r2, h, theta
    integer(int32) :: points_num
    real(real64) :: x(2), x1(2), x2(2)
    real(real64) :: t_vec(2), n_vec(2)

    ! pointsのサイズを保存
    points_num = size(edge_points, 2)

    ! moduloを用いて周期的に添え字を扱って配列外参照が起きないようにする
    ! m番目の区間の中点
    x(:) = (edge_points(:, modulo(m - 1, points_num) + 1) + edge_points(:, modulo(m, points_num) + 1))/2.0d0
    ! n番目の区間の端点
    x1(:) = edge_points(:, modulo(n - 1, points_num) + 1)
    x2(:) = edge_points(:, modulo(n, points_num) + 1)

    ! 概ね小林本の表式に合わせ, XやYを導入している.
    h = sqrt(dot_product(x1 - x2, x1 - x2))  ! 区間の長さ

    t_vec(:) = [(x2(1) - x1(1))/h, (x2(2) - x1(2))/h]
    n_vec(:) = [(x1(2) - x2(2))/h, (x2(1) - x1(1))/h]

    capital_x1 = dot_product(x - x1, t_vec)
    capital_x2 = dot_product(x - x2, t_vec)
    r1 = sqrt(dot_product(x - x1, x - x1))
    r2 = sqrt(dot_product(x - x2, x - x2))
    capital_y1 = dot_product(x - x1, n_vec)
    capital_y2 = dot_product(x - x2, n_vec)
    theta = atan2(capital_y2, capital_x2) - atan2(capital_y1, capital_x1)
  end subroutine component_values

  function U_component(m, n, edge_points) result(retval)
    !! 離散化された一重層ポテンシャルの行列の要素の値を返す
    implicit none
    integer(int32) :: m, n
    real(real64) :: edge_points(:, :)
    real(real64) :: retval
    real(real64) :: capital_x1, capital_x2, capital_y1, capital_y2, r1, r2, h, theta

    call component_values(m, n, edge_points, capital_x1, capital_x2, capital_y1, capital_y2, r1, r2, h, theta)
    if (m == n) then
      ! m=nの場合は例外的な計算
      retval = (1 - log(h/2.0d0))*h/2.0d0/PI
      return
    end if

    retval = (capital_x2*log(r2) - capital_x1*log(r1) + h - capital_y1*theta)/2.0d0/PI
  end function U_component

  function W_component(m, n, edge_points) result(retval)
    !! 離散化された二重層ポテンシャルの行列の要素の値を返す
    implicit none
    integer(int32) :: m, n
    real(real64) :: edge_points(:, :)
    real(real64) :: retval
    real(real64) :: capital_x1, capital_x2, capital_y1, capital_y2, r1, r2, h, theta

    if (m == n) then
      retval = 0.5d0
      return
    end if
    call component_values(m, n, edge_points, capital_x1, capital_x2, capital_y1, capital_y2, r1, r2, h, theta)

    retval = theta/2.0d0/PI
  end function W_component

end module subprogram
