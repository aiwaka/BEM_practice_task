module subprogram
    use,intrinsic :: iso_fortran_env
    implicit none

    REAL(real64),PARAMETER :: PI = acos(-1.0_real64)

contains
    subroutine solve(A,b,x,info)
        ! Ax=bの方程式を解き, xに結果を代入し, infoにステータスを保存する.
        real(real64),INTENT(IN) :: A(:,:), b(:,:)
        real(real64),ALLOCATABLE,INTENT(INOUT) :: x(:,:)
        integer(int32),INTENT(OUT) :: info
        INTEGER(int32), ALLOCATABLE :: ipiv(:)
        INTEGER(int32) :: an, am, bn, bm
        ! Aが正方であるかチェック
        an = size(A,1)
        am = size(A,2)
        if (an /= am) then
            print *,'Error : solve (an not equal to am) : an=',an,", am=",am
            return
        endif
        ! 基本的にbmは1であることを想定する.
        bn = size(b,1)
        bm = size(b,2)

        if (allocated(x)) DEALLOCATE(x)
        ALLOCATE(x,source=b)
        ! ピボット行列（入れ替えた結果を保存する行列. Aと同サイズの一次元配列とする）
        ALLOCATE(ipiv(an))

        call dgesv(an, bm, A, an, ipiv, x, bn, info)  ! lapackのサブルーチンを呼び出して解く.
        if (info /= 0) then
            ! infoが0でなければbad statusなのでその旨を表示する
            WRITE(*,*) 'Error : solve : info=', info
            RETURN
        endif
    end subroutine

    function integral_trapezoidal(x, points, q_array, u_array) result(retval)
        ! 台形公式を用いて周回積分を行う.
        ! 添字1とN+1の点が一致していることを利用して,
        ! すべての点を一回ずつ足してh/2をかける代わりにhをかけることにする.
        implicit none
        REAL(real64) :: x(2), points(:,:), q_array(:,:), u_array(:,:)
        REAL(real64) :: h
        REAL(real64) :: retval
        INTEGER(int32) :: points_num, i

        retval = 0.0d0
        points_num = size(points,1)
        h = 2.0d0*PI/real(points_num,real64)
        do i = 1, points_num
            ! 一重層ポテンシャル
            ! q_arrayは添字1を書いているが配列の大きさ1なので特に意味はない.
            retval = retval + fund_gamma(x,points(i,:))*q_array(i,1)
            ! 二重層ポテンシャル
            retval = retval - fund_gamma_derivative(x,points(i,:))*u_array(i,1)
        end do
        retval = retval*h
    end function integral_trapezoidal

    function exact_u(x) result(retval)
        ! 与えられた点に対してx^3 - 3xy^2の厳密な値を計算する.
        implicit none
        REAL(real64) :: x(2)
        REAL(real64) :: retval

        retval = x(1)**3 - 3*x(1)*x(2)**2
    end function exact_u

    function exact_u_normal_drv(x) result(retval)
        ! 半径1の単位円の領域とする.
        ! 与えられた点に対してx^3-3xy^2の法線微分を計算する.
        implicit none
        REAL(real64) :: x(2)
        REAL(real64) :: retval
        REAL(real64) :: theta

        theta = atan2(x(2),x(1))
        retval = 3*(x(1)**2 - x(2)**2)*cos(theta) - 6*x(1)*x(2)*sin(theta)
    end function exact_u_normal_drv

    function fund_gamma(x,y) result(retval)
        ! 任意の二点に対する（二次元）基本解の値を返す.
        implicit none
        REAL(real64) :: x(2),y(2)
        REAL(real64) :: retval

        retval = -log(sqrt(dot_product(x-y,x-y)))/2.0d0/PI
    end function fund_gamma

    function fund_gamma_derivative(x,y) result(retval)
        ! 任意の二点（ただしyは境界上）に対する基本解のyについての法線微分を返す.
        implicit none
        REAL(real64) :: x(2),y(2),theta
        REAL(real64) :: retval

        theta = atan2(y(2),y(1))
        retval = dot_product(x-y,[cos(theta),sin(theta)])/2.0d0/PI/dot_product(x-y,x-y)
    end function fund_gamma_derivative

    subroutine component_values(m, n, edge_points, lX1, lX2, lY1, lY2, r1, r2, h, theta)
        implicit none
        INTEGER(int32),intent(in) :: m, n
        REAL(real64),intent(in) :: edge_points(:,:)
        REAL(real64),intent(out) :: lX1, lX2, lY1, lY2, r1, r2, h, theta
        INTEGER(int32) :: points_num
        REAL(real64) :: x(2),x1(2),x2(2)
        REAL(real64) :: t_vec(2),n_vec(2)

        ! pointsのサイズを保存
        points_num = size(edge_points,1)

        ! moduloを用いて周期的に添え字を扱って配列外参照が起きないようにする
        ! m番目の区間の中点
        x(:) = (edge_points(modulo(m-1,points_num)+1,:) + edge_points(modulo(m,points_num)+1,:))/2.0d0
        ! n番目の区間の端点
        x1(:) = edge_points(modulo(n-1,points_num)+1,:)
        x2(:) = edge_points(modulo(n,points_num)+1,:)

        ! 概ね小林本の表式に合わせ, XやYを導入している.
        h = sqrt(dot_product(x1-x2,x1-x2))  ! 区間の長さ

        t_vec(:) = [(x2(1)-x1(1))/h, (x2(2)-x1(2))/h]
        n_vec(:) = [(x1(2)-x2(2))/h, (x2(1)-x1(1))/h]

        lX1 = dot_product(x-x1,t_vec)
        lX2 = dot_product(x-x2,t_vec)
        r1 = sqrt(dot_product(x-x1,x-x1))
        r2 = sqrt(dot_product(x-x2,x-x2))
        lY1 = dot_product(x-x1,n_vec)
        lY2 = dot_product(x-x2,n_vec)
        theta = atan2(lY2,lX2) - atan2(lY1,lX1)
    end subroutine component_values

    function U_component(m,n,edge_points) result(retval)
        ! 任意の点xと区間端点x1,x2を入力すると基本解の値を返す
        implicit none
        INTEGER(int32) :: m,n
        REAL(real64) :: edge_points(:,:)
        REAL(real64) :: retval
        REAL(real64) :: lX1,lX2,lY1,lY2,r1,r2, h, theta

        call component_values(m,n,edge_points,lX1,lX2,lY1,lY2,r1,r2,h,theta)
        if ( m == n ) then
            ! m=nの場合は例外的な計算
            retval = (1-log(h/2.0d0))*h/2.0d0/PI
            return
        end if

        retval = (lX2*log(r2)-lX1*log(r1)+h-lY1*theta)/2.0d0/PI
    end function U_component

    function W_component(m,n,edge_points) result(retval)
        ! 基本解の法線微分の値を計算する.
        implicit none
        INTEGER(int32) :: m,n
        REAL(real64) :: edge_points(:,:)
        REAL(real64) :: retval
        REAL(real64) :: lX1,lX2,lY1,lY2,r1,r2, h, theta

        if ( m == n ) then
            retval = 0.5d0
            return
        end if
        call component_values(m,n,edge_points,lX1,lX2,lY1,lY2,r1,r2,h,theta)

        retval = theta/2.0d0/PI
    end function W_component

end module subprogram

module bem
    use,intrinsic :: iso_fortran_env
    use subprogram
    implicit none
contains
    subroutine bem_calc(DIV_NUM, inner_point, result_value)
        implicit none
        INTEGER(int32),intent(in) :: DIV_NUM
        REAL(real64),intent(in) :: inner_point(2)
        REAL(real64),intent(out) :: result_value
        INTEGER(int32) :: i,m,n,info
        REAL(real64),ALLOCATABLE :: lU(:,:), lW(:,:), q(:,:), u(:,:), exact_q(:,:), edge_points(:,:),points(:,:)

        ! 割り付け
        ALLOCATE(edge_points(DIV_NUM,2))
        ALLOCATE(points(DIV_NUM,2))
        ALLOCATE(lU(DIV_NUM,DIV_NUM))
        ALLOCATE(lW(DIV_NUM,DIV_NUM))
        ALLOCATE(q(DIV_NUM,1))
        ALLOCATE(exact_q(DIV_NUM,1))
        ALLOCATE(u(DIV_NUM,1))

        ! 円周上の点を用意する. (1,0)から始めて一周する. これは端点で, 代表点はそれぞれの中点とする.
        do i = 1, DIV_NUM
            edge_points(i,1) = cos(2.0d0*(i-1)*PI/DIV_NUM)
            edge_points(i,2) = sin(2.0d0*(i-1)*PI/DIV_NUM)
        end do
        ! 代表点
        do i = 1, DIV_NUM
            points(i,:) = (edge_points(i,:) + edge_points(modulo(i,DIV_NUM)+1,:))/2.0d0
        end do

        ! 行列の各要素を計算
        do m = 1, DIV_NUM
            do n = 1, DIV_NUM
                lU(m,n) = U_component(m,n,edge_points)
                lW(m,n) = W_component(m,n,edge_points)
            end do
            u(m,1) = exact_u(points(m,:))
            ! 確認のため厳密なqも計算する.
            exact_q(m,1) = exact_u_normal_drv(points(m,:))
        end do

        ! 求めた行列を用いて共役な法線微分（q）を計算.
        call solve(lU,matmul(lW,u),q,info)
        if (info == 0) then
            ! print *,"sum of error about q :",dot_product(q(:,1)-exact_q(:,1), q(:,1)-exact_q(:,1))
        endif

        ! 台形公式で内点の値を計算
        result_value = integral_trapezoidal(inner_point, points, q, u)

        print '("at (" e20.7 "," ,e20.7 "),", "result : " e20.7)', inner_point(1), inner_point(2), result_value
        ! print *,"exact value is ",exact_u(inner_point(:))," ... error is ",exact_u(inner_point(:))-result_value

    end subroutine bem_calc
end module bem

program main
    use,intrinsic :: iso_fortran_env
    use bem
    implicit none

    INTEGER(int32) :: DIV_NUM
    INTEGER(int32),PARAMETER :: MESH_DIV_NUM = 40
    INTEGER(int32),PARAMETER :: fo = 11
    REAL(real64) :: inner_point(2)
    REAL(real64) :: result_value
    REAL(real64),ALLOCATABLE :: result_array(:)
    REAL(real64),ALLOCATABLE :: inner_points_xy(:)  ! 領域内の点を並べて保存しておくための一次元配列
    REAL(real64),ALLOCATABLE :: inner_points(:,:)
    INTEGER(int32) :: i,j,point_num

    inner_points_xy = [REAL(real64) :: ]
    point_num = 0

    do i = 1, MESH_DIV_NUM+1
        do j = 1, MESH_DIV_NUM+1
            ! 順番に点を考え, 領域内に入っていればinner_points_xy配列に追加する. point_numはそのような点の数
            inner_point(:) = [2.0d0/(real(MESH_DIV_NUM,real64))*(j-1), 2.0d0/(real(MESH_DIV_NUM,real64))*(i-1)] - [1.0d0, 1.0d0]
            if (dot_product(inner_point,inner_point) < 1.0d0) then
                point_num = point_num + 1
                inner_points_xy = [inner_points_xy, inner_point]
            endif
        end do
    end do
    ! 列挙し終わったらうまくreshapeして平面の点の配列になるようにする（二次元配列）.
    inner_points = reshape(inner_points_xy, [point_num,2], order=[2,1])

    ! 円周の分割サイズを入力（計算精度）
    print *,"input DIV_NUM : "
    READ (*,*) DIV_NUM
    if ( DIV_NUM <= 0 ) return

    ! 結果を保存する配列を用意
    ALLOCATE(result_array(point_num))
    do i = 1, point_num
        call bem_calc(DIV_NUM, inner_points(i,:), result_value)
        ! 厳密値との差を見るか実際の値を見るかいずれか
        ! result_array(i) = result_value - exact_u(inner_points(i,:))
        result_array(i) = result_value
    end do

    print *,"sum of error : ", sum(result_array)

    open(fo, file='plot.csv')
    do i = 1, point_num
        write(fo,'(3e20.7)') inner_points(i,1),inner_points(i,2),result_array(i)
    end do
    close(fo)
    ! なぜかVScodeのexecタスクで実行するとファイルが生成されないのでターミナルから実行する.
    ! plot.csvができたら, gnuplotを起動し[sp "plot.csv" w p]と入力するとプロットされる.


end program main