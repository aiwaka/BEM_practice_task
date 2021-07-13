module subprogram
    use,intrinsic :: iso_fortran_env
    implicit none

    REAL(real64),PARAMETER :: PI = acos(-1.0_real64)

contains
    subroutine solve(A,b,x,info)
        real(real64),INTENT(IN) :: A(:,:), b(:,:)
        real(real64),ALLOCATABLE,INTENT(INOUT) :: x(:,:)
        integer(int32),INTENT(OUT) :: info
        INTEGER, ALLOCATABLE :: ipiv(:)
        INTEGER :: an, am, bn, bm
        an = size(A,1)
        am = size(A,2)
        if (an /= am) then
            print *,'Error : solve (an not equal to am) : an=',an,", am=",am
            return
        endif
        bn = size(b,1)
        bm = size(b,2)
        if (allocated(x)) then
            DEALLOCATE(x)
        end if
        ALLOCATE(x,source=b)
        ALLOCATE(ipiv(an))
        
        call dgesv(an, bm, A, an, ipiv, x, bn, info)  ! bmは1を想定（bは普通の配列）
        if (info /= 0) then
            WRITE(*,*) 'Error : solve : info=', info
            RETURN
        endif
    end subroutine

    function exact_u(x) result(retval)
        ! x^3 - 3xy^2
        implicit none
        REAL(real64) :: x(2)
        REAL(real64) :: retval
        
        retval = x(1)**3 - 3*x(1)*x(2)*x(2)
    end function exact_u

    function exact_u_norm_drv(x) result(retval)
        ! 半径1の単位円領域とする
        implicit none
        REAL(real64) :: x(2)
        REAL(real64) :: retval
        REAL(real64) :: theta

        theta = atan2(x(2),x(1))
        retval = (3*x(1)*x(1) - 3*x(2)*x(2))*cos(theta) + 6*x(1)*x(2)*sin(theta)
    end function exact_u_norm_drv

    function fund_gamma(x,y) result(retval)
        ! 任意の二点に対する基本解の値を返す
        implicit none
        REAL(real64) :: x(2),y(2)
        REAL(real64) :: retval
        
        retval = -log(sqrt(dot_product(x-y,x-y)))/2/PI
    end function fund_gamma

    function U_component(m,n,points) result(retval)
        ! 任意の点xと区間端点x1,x2を入力すると基本解の値を返す
        implicit none
        INTEGER(int32) :: m,n,points_num
        REAL(real64) :: points(:,:)
        REAL(real64) :: x(2),x1(2),x2(2)
        REAL(real64) :: retval
        REAL(real64) :: tVec(2),nVec(2)
        REAL(real64) :: lX1,lX2,lY1,lY2,r1,r2, h, theta

        ! pointsのサイズを保存
        points_num = size(points,1)

        ! moduloを用いて周期的に添え字を扱って配列外参照が起きないようにする
        ! m番目の区間の中点
        x(:) = (points(modulo(m-1,points_num)+1,:) + points(modulo(m,points_num)+1,:))/2
        ! n番目の区間の端点
        x1(:) = points(modulo(n-1,points_num)+1,:)
        x2(:) = points(modulo(n,points_num)+1,:)

        ! 概ね小林本の表式に合わせ, XやYを導入している.
        h = sqrt(dot_product(x1-x2,x1-x2))  ! 区間の長さ
        if ( m == n ) then
            ! m=nの場合は例外的な計算
            retval = (1-log(h/2))*h/2/PI
            return
        end if

        tVec(:) = [(x2(1)-x1(1))/h, (x2(2)-x1(2))/h]
        nVec(:) = [(x1(2)-x2(2))/h, (x2(1)-x1(1))/h]

        lX1 = dot_product(x-x1,tVec)
        lX2 = dot_product(x-x2,tVec)
        r1 = sqrt(dot_product(x-x1,x-x1))
        r2 = sqrt(dot_product(x-x2,x-x2))
        lY1 = dot_product(x-x1,nVec)
        lY2 = dot_product(x-x2,nVec)
        theta = atan2(lY2,lX2) - atan2(lY1,lX1)

        retval = (lX2*log(r2)-lX1*log(r1)+h-lY1*theta)/2/PI
    end function U_component

    function W_component(m,n,points) result(retval)
        implicit none
        INTEGER(int32) :: m,n,points_num
        REAL(real64) :: points(:,:)
        REAL(real64) :: x(2),x1(2),x2(2)
        REAL(real64) :: retval
        REAL(real64) :: tVec(2),nVec(2)
        REAL(real64) :: lX1,lX2,lY1,lY2,r1,r2, h, theta

        points_num = size(points,1)

        x(:) = (points(modulo(m-1,points_num)+1,:) + points(modulo(m,points_num)+1,:))/2
        x1(:) = points(modulo(n-1,points_num)+1,:)
        x2(:) = points(modulo(n,points_num)+1,:)

        if ( m == n ) then
            retval = 0.5d0
            return
        end if
        h = sqrt(dot_product(x1-x2,x1-x2))  ! 区間の長さ

        tVec(:) = [(x2(1)-x1(1))/h, (x2(2)-x1(2))/h]
        nVec(:) = [(x1(2)-x2(2))/h, (x2(1)-x1(1))/h]

        lX1 = dot_product(x-x1,tVec)
        lX2 = dot_product(x-x2,tVec)
        r1 = sqrt(dot_product(x-x1,x-x1))
        r2 = sqrt(dot_product(x-x2,x-x2))
        lY1 = dot_product(x-x1,nVec)
        lY2 = dot_product(x-x2,nVec)
        theta = atan2(lY2,lX2) - atan2(lY1,lX1)

        retval = theta/2/PI

        
    end function W_component

end module subprogram

program bem
    use,intrinsic :: iso_fortran_env
    use subprogram
    implicit none

    INTEGER(int32) :: i,m,n,info
    INTEGER(int32) :: DIV_NUM
    REAL(real64),ALLOCATABLE :: lU(:,:), lW(:,:), q(:,:), u(:,:), exact_q(:,:), points(:,:)

    print *,"input DIV_NUM : "
    READ (*,*) DIV_NUM
    if ( DIV_NUM <= 0 ) return

    ! 割り付け
    ALLOCATE(points(DIV_NUM,2))
    ALLOCATE(lU(DIV_NUM,DIV_NUM))
    ALLOCATE(lW(DIV_NUM,DIV_NUM))
    ALLOCATE(q(DIV_NUM,1))
    ALLOCATE(exact_q(DIV_NUM,1))
    ALLOCATE(u(DIV_NUM,1))

    ! 円周上の点を用意する
    do i = 1, DIV_NUM
        points(i,1) = cos(2*(i-1)*PI/DIV_NUM)
        points(i,2) = sin(2*(i-1)*PI/DIV_NUM)
    end do

    ! 行列の各要素を計算
    do m = 1, DIV_NUM
        do n = 1, DIV_NUM
            lU(m,n) = U_component(m,n,points)
            lW(m,n) = W_component(m,n,points)
        end do
        u(m,1) = exact_u((points(m,:) + points(modulo(m,DIV_NUM)+1,:))/2)
        exact_q(m,1) = exact_u_norm_drv((points(m,:) + points(modulo(m,DIV_NUM)+1,:))/2)
    end do
    do m = 1, DIV_NUM
        print *,u(m,1)
        print *,exact_q(m,1)
    end do

    call solve(lU,matmul(lW,u),q,info)
    if (info == 0) then
        print *,"Err about q :",q(:,1)-exact_q(:,1)
    endif
    stop
end program bem