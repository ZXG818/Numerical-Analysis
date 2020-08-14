!  此为三对角矩阵求解方程，是书上的例程，不是自己写的。
!
!  FUNCTIONS:
!  Console1 - 书上的例程。
!

!****************************************************************************
!
!  PROGRAM: Console1
!
!  PURPOSE:  三对角矩阵求解
!
!****************************************************************************

    program main
    implicit none
        integer, parameter :: Width = 3;
        integer, parameter :: Row = 5;
        real :: A(Row, Width) = (/ 0, 2, 3, 4, 1, &
                                   1, 3, 4, 5, 2, &
                                   2, 4, 5, 1, 0 /)
        real :: S(Row) = (/ 3, 9, 12, 10, 3 /)
        real :: ans(Row)
        integer :: i
        !equation
        ! a + 2b = 3
        ! 2b + 3c + 4d = 9
        ! 3c + 4d + 5e = 12
        ! 4d + 5e + f = 10
        ! e + 2f = 3
        
        call Gauss_Jordan(A, S, ANS, Row, Width)
        write(*, *) 'Ans:'
        do i=1, Row
            write(*, "(1x, a1, '=', F8.4)") char(96+i), ANS(i)
        end do
    end program main
    
    ! Gauss-Jordan法求解函数
    subroutine Gauss_Jordan(A, S, ANS, Row, Width)
    implicit none
        integer, intent(in) :: Row, Width
        real, intent(inout) :: A(Row, Width)
        real, intent(inout) :: S(Row)
        real, intent(inout) :: ANS(Row)
        real :: B(Row, Width)
        integer :: i
        
        ! 保存原先的矩阵A以及数组S
        B = A
        ANS = S
        
        ! 将B化为对角线矩阵
        call Upper(B, ANS, Row, Width)
        call Lower(B, ANS, Row, Width)
        
        ! 求解
        !forall(i=1:Row) ANS(i) = ANS(i) / B(i, 2)
        do i=1, Row
            ANS(i) = ANS(i) / B(i, 2)
        end do
        return 
    end subroutine 
    
    ! 求上三角矩阵的子程序
    subroutine Upper(M, S, Row, Width)
    implicit none
        integer, intent(in) :: Row, Width
        real, intent(inout) :: M(Row, Width)
        real, intent(inout) :: S(Row)
        integer :: i, j
        
        real :: E
        do i=1, Row-1
            j = i + 1
            E = M(j, 1) / M(i, 2)
            M(j, 1:2) = M(j, 1:2) - M(i, 2:3) * E
            S(j) = S(j) - S(i) * E
        end do
        return
    end subroutine 
    
    ! 求下三角矩阵的子程序
    subroutine Lower(M, S, Row, Width)
    implicit none
        integer, intent(in) :: Row, Width
        real, intent(inout) :: M(Row, Width)
        real, intent(inout) :: S(Row)
        
        integer :: i, j
        real :: E
        
        do i=Row, 2, -1
            j = i-1
            E = M(j,3) / M(i, 2)
            M(j, 3) = M(j, 3) - M(i, 2) * E
            S(j) = S(j) - S(i) * E
        end do
        return 
    end subroutine
