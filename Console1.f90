    MODULE LINEAR_ALGEBRA
    IMPLICIT NONE
        PUBLIC N
        PUBLIC SHOWMAT
        PUBLIC UPPER
        PUBLIC LOWER
        INTEGER, PARAMETER :: N = 3
    CONTAINS
        ! 显示矩阵内容
        SUBROUTINE SHOWMAT(MATRIX)
        IMPLICIT NONE
            REAL, INTENT(IN) :: MATRIX(:, :)
            INTEGER :: I
            DO I=1, SIZE(MATRIX, 2)     ! 遍历矩阵的行
                WRITE(*, "(4F8.2)") MATRIX(:, I)
            END DO
            RETURN 
        END SUBROUTINE
        
        ! 求上三角矩阵的子程序
        SUBROUTINE UPPER(MATRIX)
        IMPLICIT NONE
            REAL, INTENT(INOUT) :: MATRIX(:, :)
            INTEGER :: I, J
            ! 进行矩阵的行变换
            DO I=1, SIZE(MATRIX, 2)-1
                DO J=I+1, SIZE(MATRIX, 2)
                    MATRIX(:, J) = MATRIX(:, J) - MATRIX(I, J) * MATRIX(:, I) / MATRIX(I, I)
                END DO
            END DO
            RETURN
        END SUBROUTINE
        
        ! 求下三角矩阵的子程序
        SUBROUTINE LOWER(MATRIX)
        IMPLICIT NONE
            REAL, INTENT(INOUT) :: MATRIX(:, :)
            INTEGER :: I
            INTEGER :: J
            ! 进行矩阵的行变换
            DO I=SIZE(MATRIX, 2), 2, -1
                DO J=I-1, 1, -1
                    MATRIX(:, J) = MATRIX(:, J) - MATRIX(I, J) &
                                    * (MATRIX(:, I) / MATRIX(I, I))
                END DO
            END DO
        END SUBROUTINE
    END MODULE
    
    PROGRAM MAIN
    USE LINEAR_ALGEBRA
    IMPLICIT NONE
        REAL :: A(N, N) = (/ 1, 2, 1, 3, 2, 3, 2, 3, 4 /)
        REAL :: B(N, N)
        B = A
        CALL SHOWMAT(A)
        CALL UPPER(A)
        WRITE(*, *) "------------------"
        CALL SHOWMAT(A)
        WRITE(*, *) "------------------"
        CALL LOWER(B)
        CALL SHOWMAT(B)
    END PROGRAM
            
