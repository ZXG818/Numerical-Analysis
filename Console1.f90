    MODULE NUMERICAL
    IMPLICIT NONE
        PRIVATE UPPER
        PRIVATE LOWER
        PUBLIC SHOWMAT
        PUBLIC INVERSE_MAT
    CONTAINS
        ! 打印矩阵
        SUBROUTINE SHOWMAT(MAT)
        IMPLICIT NONE
            REAL, INTENT(IN) :: MAT(:, :)
            INTEGER :: I
            DO I=1, SIZE(MAT, 2)
                WRITE(*, "(3F8.4)") MAT(:, I)
            END DO
            RETURN 
        END SUBROUTINE
        
        ! 上三角化系数矩阵，同时对等号右边列向量做同等行变换
        SUBROUTINE UPPER(MAT, E)
        IMPLICIT NONE
            REAL, INTENT(INOUT) :: MAT(:, :)
            REAL, INTENT(INOUT) :: E(:, :) ! 单位矩阵
            INTEGER :: I, J
            DO I=1, SIZE(MAT, 2)-1
                DO J=I+1, SIZE(MAT, 2)
                    ! 注意，100语句与200语句的顺序不能调换！
100                 E(:, J) = E(:, J) - E(:, I) * MAT(I, J) / MAT(I, I)               
200                 MAT(:, J) = MAT(:, J) - MAT(I, J) * MAT(:, I) / MAT(I, I)
                END DO
            END DO
            RETURN
        END SUBROUTINE
        
        ! 下三角化系数矩阵，同时对等号右边列向量做同等行变换
        SUBROUTINE LOWER(MAT, E)
        IMPLICIT NONE
            REAL, INTENT(INOUT) :: MAT(:, :)
            REAL, INTENT(INOUT) :: E(:, :)
            INTEGER :: I, J
            ! 进行矩阵的行变换
            DO I=SIZE(MAT, 2), 2, -1
                DO J=I-1, 1, -1
                    ! 注意，300语句与400语句的顺序不能调换！
300                 E(:, J) = E(:, J) - E(:, I) * MAT(I, J) / MAT(I, I)
400                 MAT(:, J) = MAT(:, J) - MAT(I, J) * MAT(:, I) / MAT(I, I)
                END DO
            END DO
            RETURN
        END SUBROUTINE
        
        ! 求逆矩阵
        SUBROUTINE INVERSE_MAT(MAT, E)
        IMPLICIT NONE
            REAL, INTENT(INOUT) :: MAT(:, :) ! 系数矩阵
            REAL, INTENT(INOUT) :: E(:, :)
            INTEGER :: I, J
            CALL SHOWMAT(MAT)
            WRITE(*, *) "------------------"
            CALL UPPER(MAT, E)
            CALL SHOWMAT(MAT)
            WRITE(*, *) "------------------"
            CALL LOWER(MAT, E)
            CALL SHOWMAT(MAT)
            WRITE(*, *) "------------------"
            ! 得出MAT的逆矩阵
            FORALL(I=1:SIZE(MAT, 2)) E(:, I) = E(:, I) / MAT(I, I)
            RETURN 
        END SUBROUTINE    
    END MODULE
    
    PROGRAM MAIN
    USE NUMERICAL
    IMPLICIT NONE
        REAL :: MAT(3, 3) = (/ 1, 1, 3, 3, 2, 2, 1, 3, 1 /)
        REAL :: BACKUP(3, 3) = (/ 1, 1, 3, 3, 2, 2, 1, 3, 1 /)
        INTEGER :: I, J
        REAL :: E(3, 3)
        
        ! 生成单位矩阵
        FORALL(I=1:3, J=1:3, I==J) E(I, J) = 1.0
        FORALL(I=1:3, J=1:3, I/=J) E(I, J) = 0.0
        
        CALL INVERSE_MAT(MAT, E)
        CALL SHOWMAT(E)
        WRITE(*, *) MATMUL(BACKUP, E)   ! 使用原矩阵与求得的逆矩阵进行乘积，看看是否是单位矩阵
    END PROGRAM
    