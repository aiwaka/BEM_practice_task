module subprogram
    use,intrinsic :: iso_fortran_env
    implicit none
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

end module subprogram

program bem
    use,intrinsic :: iso_fortran_env
    use subprogram
    implicit none

    INTEGER(int32) :: i,info
    REAL(real64),ALLOCATABLE :: A(:,:), b(:,:), x(:,:)
    INTEGER(int32),PARAMETER :: DIV_NUM = 100
    REAL(real64),ALLOCATABLE :: points(:)


    ALLOCATE(points(DIV_NUM))

    ALLOCATE(A, source=reshape([ 2.0d0, 4.0d0, 6.0d0, 3.0d0, -5.0d0, -7.0d0, 1.0d0, -1.0d0, 1.0d0 ], [3,3]))
    ALLOCATE(b, source=reshape([7.0d0, 3.0d0, 5.0d0], [3,1]))

    do i = 1, 3
        print *,"[",i,"1 ",i,"2 ",i,"3 ] = ",A(i,:)
    end do
    print *, "b = ",b(:,1)

    call solve(A,b,x,info)
    if (info == 0) then
        print *,"x=",x(:,1)
    endif
    stop
end program bem