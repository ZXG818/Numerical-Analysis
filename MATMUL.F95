    MODULE MATMUL_TEST
    IMPLICIT NONE
        PRIVATE N
        PUBLIC A
        PUBLIC B
        PUBLIC C
        PUBLIC SHOWMAT
        PUBLIC MATMUL_
        
        ! 此处我定义矩阵，采用第二维代表行，第一维代表列
        ! 与C语言中的正好相反，应当注意！
        INTEGER, PARAMETER :: N = 3  ! 矩阵维数，定义3X3矩阵
        INTEGER :: A(N, N)
        INTEGER :: B(N, N) = (/ 1, 2, 3, 4, 5, 6, 7, 8, 9 /)
        INTEGER :: C(N, N) = (/ 9, 8, 7, 6, 5, 4, 3, 2, 1 /)
        
    CONTAINS
        ! 输出矩阵的子程序
        SUBROUTINE SHOWMAT(MAT)
        IMPLICIT NONE
            INTEGER, INTENT(IN) :: MAT(:, :)
            INTEGER :: I
            DO I=1, SIZE(MAT, 2)   ! 按照行数开始打印
                WRITE(*, "(3I4)") MAT(:, I)
            END DO
            RETURN
        END SUBROUTINE
        
        ! 矩阵乘法子程序
        ! 输出A=B·C
        SUBROUTINE MATMUL_(A, B, C)
        IMPLICIT NONE
            INTEGER, INTENT(OUT) :: A(:, :) ! 结果矩阵
            INTEGER, INTENT(IN) :: B(:, :)
            INTEGER, INTENT(IN) :: C(:, :)
            INTEGER :: I, J, K
            
            ! 检查矩阵维数，矩阵维数不对应的无法进行操作！
            IF(SIZE(B, 1) /= SIZE(C, 2)) THEN
                WRITE(*, *) "DIMENSIONS ARE NOT CORRESPONDING..."
                RETURN
            END IF
            
            ! 进行矩阵相乘
            DO I=1, SIZE(B, 2)      ! 遍历B矩阵的行
                DO J=1, SIZE(C, 1)  ! 遍历C矩阵的列
                    A(J, I) = 0
                    DO K=1, SIZE(B, 1)
                        A(J, I) = A(J, I) + B(K, I) * C(J, K)
                    END DO
                END DO
            END DO
            RETURN
        END SUBROUTINE               
    END MODULE
    
    PROGRAM MAIN
    USE MATMUL_TEST
    IMPLICIT NONE
        CALL SHOWMAT(B)
        WRITE(*, *) "-------------"
        CALL SHOWMAT(C)
        WRITE(*, *) "-------------"
        CALL MATMUL_(A, B, C)  ! 矩阵相乘
        CALL SHOWMAT(A)
    END PROGRAM
