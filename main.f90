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
            print *,'Error : solve (an not equal to am)'
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

    function fund_gamma(x,y) result(retval)
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
        REAL(real64) :: t_vec(2),n_vec(2)
        REAL(real64) :: lX1,lX2,lY1,lY2,r1,r2
        REAL(real64) :: h,theta

        ! pointsのサイズを保存
        points_num = size(points,1)

        ! moduloを用いて周期的に添え字を扱って配列外参照が起きないようにする
        ! m番目の区間の中点
        x(:) = (points(modulo(m-1,points_num)+1,:) + points(modulo(m,points_num)+1,:))/2
        ! n番目の区間の端点
        x1(:) = points(modulo(n-1,points_num)+1,:)
        x2(:) = points(modulo(n,points_num)+1,:)
        print *,"m=",m,",n=",n,",x(:):",x(:)

        ! 概ね小林本の表式に合わせ, XやYを導入している.
        h = sqrt(dot_product(x1-x2,x1-x2))  ! 区間の長さ
        if ( m == n ) then
            ! m=nの場合は例外的な計算
            retval = (1-log(h/2))*h/2/PI
            return
        end if

        t_vec(1:2) = [(x2(1)-x1(1))/h, (x2(2)-x1(2))/h]
        n_vec(1:2) = [(x1(2)-x2(2))/h, (x2(1)-x1(1))/h]
        lX1 = dot_product(x-x1,t_vec)
        lX2 = dot_product(x-x2,t_vec)
        r1 = sqrt(dot_product(x-x1,x-x1))
        r2 = sqrt(dot_product(x-x2,x-x2))
        lY1 = dot_product(x-x1,n_vec)
        lY2 = dot_product(x-x2,n_vec)
        theta = atan2(lY2,lX2) - atan2(lY1,lX1)

        retval = (lX2*log(r2)-lX1*log(r1)+h-lY1*theta)/2/PI
    end function U_component

end module subprogram

program bem
    use,intrinsic :: iso_fortran_env
    use subprogram
    implicit none

    INTEGER(int32) :: i,m,n,info
    REAL(real64),ALLOCATABLE :: A(:,:), b(:,:), x(:,:)
    ! REAL(real64),ALLOCATABLE :: D(:,:)
    INTEGER(int32) :: DIV_NUM
    REAL(real64),ALLOCATABLE :: points(:,:)

    print *,"input DIV_NUM : "
    READ (*,*) DIV_NUM
    if ( DIV_NUM <= 0 ) then
        return
    end if

    ALLOCATE(points(DIV_NUM,2))
    do i = 1, DIV_NUM
        points(i,1) = cos(2*i*PI/DIV_NUM)
        points(i,2) = sin(2*i*PI/DIV_NUM)
    end do
    do i = 1, DIV_NUM
        print *,i,";",points(i,:)
    end do
    
    do m = 1, DIV_NUM
        do n = 1, DIV_NUM
            A(m,n) = U_component(m,n,points)
        end do
    end do

    ! Ax=bを解くのを試したときのもの
    ! ALLOCATE(A, source=reshape([ 2.0d0, 4.0d0, 6.0d0, 3.0d0, -5.0d0, -7.0d0, 1.0d0, -1.0d0, 1.0d0 ], [3,3]))
    ! ALLOCATE(b, source=reshape([7.0d0, 3.0d0, 5.0d0], [3,1]))

    call solve(A,b,x,info)
    if (info == 0) then
        print *,"x=",x(:,1)
    endif
    stop
end program bem