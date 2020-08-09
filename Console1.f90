    ! Gauss-Jordan方法求解联立方程式
    ! 此方法即为线性代数中常用的通过行变换化增广矩阵为行最简型进行求解。
    ! 核心：先对系数矩阵进行上三角化，再进行下三角化，即可得到系数矩阵的行最简型矩阵，
    !      同时，等号右边的列向量也需要进行同步的行变换。
    MODULE NUMERICAL
    IMPLICIT NONE
        PRIVATE UPPER
        PRIVATE LOWER
        PUBLIC SHOWMAT
        PUBLIC SOLVER
    CONTAINS
        ! 打印矩阵
        SUBROUTINE SHOWMAT(MAT)
        IMPLICIT NONE
            REAL, INTENT(IN) :: MAT(:, :)
            INTEGER :: I
            DO I=1, SIZE(MAT, 2)
                WRITE(*, "(3F6.2)") MAT(:, I)
            END DO
            RETURN 
        END SUBROUTINE
        
        ! 上三角化系数矩阵，同时对等号右边列向量做同等行变换
        SUBROUTINE UPPER(MAT, B)
        IMPLICIT NONE
            REAL, INTENT(INOUT) :: MAT(:, :)
            REAL, INTENT(INOUT) :: B(:)
            INTEGER :: I, J
            DO I=1, SIZE(MAT, 2)-1
                DO J=I+1, SIZE(MAT, 2)
                    ! 注意，100语句与200语句的顺序不能调换！
100                 B(J) = B(J) - B(I) * MAT(I, J) / MAT(I, I)                 ! 对等号右边列向量同步行变换
200                 MAT(:, J) = MAT(:, J) - MAT(I, J) * MAT(:, I) / MAT(I, I)  ! 对系数矩阵行变换
                END DO
            END DO
            RETURN
        END SUBROUTINE
        
        ! 下三角化系数矩阵，同时对等号右边列向量做同等行变换
        SUBROUTINE LOWER(MAT, B)
        IMPLICIT NONE
            REAL, INTENT(INOUT) :: MAT(:, :)
            REAL, INTENT(INOUT) :: B(:)
            INTEGER :: I, J
            ! 进行矩阵的行变换
            DO I=SIZE(MAT, 2), 2, -1
                DO J=I-1, 1, -1
                    ! 注意，300语句与400语句的顺序不能调换！
300                 B(J) = B(J) - B(I) * MAT(I, J) / MAT(I, I)
400                 MAT(:, J) = MAT(:, J) - MAT(I, J) &
                                    * MAT(:, I) / MAT(I, I)
                END DO
            END DO
            RETURN
        END SUBROUTINE
        
        ! 求解方程
        SUBROUTINE SOLVER(MAT, X, B)
        IMPLICIT NONE
            REAL, INTENT(INOUT) :: MAT(:, :) ! 系数矩阵
            REAL, INTENT(OUT) :: X(:)        ! 待求解的列向量
            REAL, INTENT(INOUT) :: B(:)      ! 等号右边的列向量
            INTEGER :: I
            ! 进行一次系数矩阵的上三角化，再进行一次下三角化，
            ! 同时等号右边列向量同步进行行变换，
            ! 最终会得到行最简矩阵，非常方便进行求解！
            CALL SHOWMAT(MAT)
            !WRITE(*, *) B
            WRITE(*, *) "------------------"
            CALL UPPER(MAT, B)
            CALL SHOWMAT(MAT)
            !WRITE(*, *) B
            WRITE(*, *) "------------------"
            CALL LOWER(MAT, B)
            CALL SHOWMAT(MAT)
            !WRITE(*, *) B
            ! 求解
            FORALL(I=1:SIZE(MAT, 2)) X(I) = B(I) / MAT(I, I)
            RETURN 
        END SUBROUTINE    
    END MODULE
    
    PROGRAM MAIN
    USE NUMERICAL
    IMPLICIT NONE
        REAL :: MAT(3, 3) = (/ 1, 4, 7, 2, 5, 8, 3, 6, 8 /)
        REAL :: X(3)
        REAL :: B(3) = (/ 12, 15, 17 /)
        CALL SOLVER(MAT, X, B)
        WRITE(*, "(A8, 3F10.4)") "结果为：", X
    END PROGRAM
    
!    测试例子即为方程组：
!    1X + 4Y + 7Z = 12
!    2X + 5Y + 8Z = 15
!    3X + 6Y + 8Z = 17
!    计算结果：X = 1.0000  Y = 1.0000  Z = 1.0000