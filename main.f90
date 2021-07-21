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
        INTEGER, ALLOCATABLE :: ipiv(:)
        INTEGER :: an, am, bn, bm
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

    function integral_trapezoidal(x, points, q_array) result(retval)
        implicit none
        REAL(real64) :: x(2), points(:,:), q_array(:,:)
        REAL(real64) :: h
        REAL(real64) :: retval
        INTEGER(int32) :: points_num, i, cyclic_num

        retval = 0.0d0
        points_num = size(points,1)
        h = 2.0d0*PI/real(points_num,real64)
        do i = 1, points_num+1
            cyclic_num = modulo(i-1, points_num) + 1
            retval = retval + fund_gamma(x,points(cyclic_num,:))*q_array(cyclic_num,1)
        end do
        retval = retval*h/2.0d0
        
    end function integral_trapezoidal

    function exact_u(x) result(retval)
        ! x^3 - 3xy^2
        implicit none
        REAL(real64) :: x(2)
        REAL(real64) :: retval

        retval = x(1)**3 - 3*x(1)*x(2)**2
    end function exact_u

    function exact_u_normal_drv(x) result(retval)
        ! 半径1の単位円の領域とする
        implicit none
        REAL(real64) :: x(2)
        REAL(real64) :: retval
        REAL(real64) :: theta

        theta = atan2(x(2),x(1))
        retval = 3*(x(1)**2 - x(2)**2)*cos(theta) - 6*x(1)*x(2)*sin(theta)
    end function exact_u_normal_drv

    function fund_gamma(x,y) result(retval)
        ! 任意の二点に対する基本解の値を返す
        implicit none
        REAL(real64) :: x(2),y(2)
        REAL(real64) :: retval

        retval = -log(sqrt(dot_product(x-y,x-y)))/2.0d0/PI
    end function fund_gamma

    function U_component(m,n,edge_points) result(retval)
        ! 任意の点xと区間端点x1,x2を入力すると基本解の値を返す
        implicit none
        INTEGER(int32) :: m,n,points_num
        REAL(real64) :: edge_points(:,:)
        REAL(real64) :: x(2),x1(2),x2(2)
        REAL(real64) :: retval
        REAL(real64) :: t_vec(2),n_vec(2)
        REAL(real64) :: lX1,lX2,lY1,lY2,r1,r2, h, theta

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
        if ( m == n ) then
            ! m=nの場合は例外的な計算
            retval = (1-log(h/2.0d0))*h/2.0d0/PI
            return
        end if

        t_vec(:) = [(x2(1)-x1(1))/h, (x2(2)-x1(2))/h]
        n_vec(:) = [(x1(2)-x2(2))/h, (x2(1)-x1(1))/h]

        lX1 = dot_product(x-x1,t_vec)
        lX2 = dot_product(x-x2,t_vec)
        r1 = sqrt(dot_product(x-x1,x-x1))
        r2 = sqrt(dot_product(x-x2,x-x2))
        lY1 = dot_product(x-x1,n_vec)
        lY2 = dot_product(x-x2,n_vec)
        theta = atan2(lY2,lX2) - atan2(lY1,lX1)

        retval = (lX2*log(r2)-lX1*log(r1)+h-lY1*theta)/2.0d0/PI
    end function U_component

    function W_component(m,n,edge_points) result(retval)
        ! 基本解の法線微分の値を計算する.
        implicit none
        INTEGER(int32) :: m,n,points_num
        REAL(real64) :: edge_points(:,:)
        REAL(real64) :: x(2),x1(2),x2(2)
        REAL(real64) :: retval
        REAL(real64) :: t_vec(2),n_vec(2)
        REAL(real64) :: lX1,lX2,lY1,lY2,r1,r2, h, theta

        points_num = size(edge_points,1)

        x(:) = (edge_points(modulo(m-1,points_num)+1,:) + edge_points(modulo(m,points_num)+1,:))/2.0d0
        x1(:) = edge_points(modulo(n-1,points_num)+1,:)
        x2(:) = edge_points(modulo(n,points_num)+1,:)

        if ( m == n ) then
            retval = 0.5d0
            return
        end if
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

        retval = theta/2.0d0/PI
    end function W_component

end module subprogram

program bem
    use,intrinsic :: iso_fortran_env
    use subprogram
    implicit none

    INTEGER(int32) :: i,m,n,info
    INTEGER(int32) :: DIV_NUM
    REAL(real64),ALLOCATABLE :: lU(:,:), lW(:,:), q(:,:), u(:,:), exact_q(:,:), edge_points(:,:),points(:,:)
    REAL(real64) :: inner_point(2)
    REAL(real64) :: result_value

    ! 分割サイズを入力
    print *,"input DIV_NUM : "
    READ (*,*) DIV_NUM
    if ( DIV_NUM <= 0 ) return
    ! 値を求める点を入力
    inner_point(:) = [0.5d0,0.5d0]
    if ( dot_product(inner_point(:),inner_point(:)) >= 1.0d0 ) then
        print *,"your point is out of domain."
        stop
    end if  

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
        print *,"sum of error about q :",sum(q(:,1)-exact_q(:,1))
    endif

    ! 台形公式で内点の値を計算
    result_value = integral_trapezoidal(inner_point, points, q)

    print *,"at (", inner_point(1),",",inner_point(2),"), result : ",result_value
    print *,"exact value is ",exact_u(inner_point(:))

end program bem