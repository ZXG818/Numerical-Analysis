    ! 求解矩阵的行列式
    ! 将矩阵（方阵）转换为上三角或下三角形式
    ! 主对角线元素的积即为行列式的值
    module LinearAlgebra
    implicit none
        public ShowMat
        public Upper
        public GetDet
    contains
        ! 打印矩阵
        subroutine ShowMat(mat)
        implicit none
            real, intent(in) :: mat(:, :)
            integer :: i
            do i=1, size(mat, 2)
                write(*, "(3F8.2)") mat(:, i)
            end do
            return
        end subroutine
        
        ! 转换为上三角矩阵
        subroutine Upper(mat)
        implicit none
            real, intent(inout) :: mat(:, :)
            integer :: i, j
            do i=1, size(mat, 2)-1
                do j=i+1, size(mat, 2)
                    mat(:, j) = mat(:, j) - mat(i, j) * mat(:, i) / mat(i, i)
                end do
            end do
            return
        end subroutine
        
        ! 对上三角矩阵进行对角线乘积
        subroutine GetDet(mat, det)
        implicit none
            real, intent(inout) :: mat(:, :)
            real, intent(out) :: det
            integer :: i
            
            call Upper(mat)
            det = 1.0
            do i=1, size(mat, 2)
                det = det * mat(i, i)
            end do
            return
        end subroutine
    end module
    
    program main
    use LinearAlgebra
    implicit none
        real :: mat(3, 3) = (/ 1, 2, 1, 3, 2, 3, 2, 3, 4 /)
        real :: det = 0.0
        call ShowMat(mat)
        call Upper(mat)
        write(*, *) "-------------------"
        call ShowMat(mat)
        
        call GetDet(mat, det)
        write(*, *) det
    end program
    