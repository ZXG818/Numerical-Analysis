#! /usr/bin/env python
#! -*- coding:utf-8 -*-
# Author : Xingguang Zhou
# Date   : 2020/9/22
# Program: 矩阵的LU分解
# 更改的内容：2020/9/26在之前的计算方法中，由于上课和书上的例子全部都是对称矩阵（即原矩阵与其转置矩阵相等）
# 所以在计算方法的理解上出现了小的偏差，但是最终结果是正确的，所以在和刘洪权同学的讨论中，
# 他非常细心地给我指正了错误，并帮助我修改了程序，非常感谢他！
# 需要注意的是，LU分解中，L矩阵依旧需要像U矩阵的求解那样，进行二重循环求解。

# Print the matrix
def PrintMatrix(matrix):
    print('===================')
    row = len(matrix)
    for _ in range(row):
        print(matrix[_])
    print('===================')

# LU decomposition matrix
# Input  : An original matrix which need to be decomposed in LU method.
# Return : L matrix and U matrix
def LU(matrix):
    # The row number of the matrix
    row = len(matrix)
    # The coloumn number of the matrix
    col = len(matrix[0])
    
    """ Caution! We assume the row is equal to col, so this matrix is a square matrix. """
    
    # Use the list comprehension
    # Initial the U matrix
    U = [[0 for i in range(row)] for j in range(col)]
    # Initial the L matrix
    L = [[0 for i in range(row)] for j in range(col)]
    
    # The first frame
    for i in range(col):
        U[0][i] = matrix[0][i]
        L[i][0] = matrix[i][0] / matrix[0][0]
    
    # 开始进行第二框及以后的LU分解。
    for i in range(1, row):        
        # 这是我写的方法，运算量比书上的小多了。
        # 开始处理非对角线上的元素
        for j in range(i, col):
            # 开始计算U矩阵的第i行
            summary = 0
            for k in range(i):
                summary = summary + U[k][j] * L[i][k]
            U[i][j] = matrix[i][j] - summary
            # 下面开始计算L矩阵第i列
            summary = 0
            for k in range(j):
                summary = summary + U[k][i] * L[j][k]
            L[j][i] = (matrix[j][i] - summary) / U[i][i]
        
    PrintMatrix(U)
    PrintMatrix(L)                
    


if __name__ == '__main__':
    print("第一次实验：")
    matrix = [[ 4,  -2,   0,  4], 
              [-2,   2,  -3,  1],
              [ 0,  -3,  13, -7],
              [ 4,   1,  -7, 23]]
    LU(matrix)
    
    print("第二次实验：")
    matrix = [[  9,  18,   9, -27],
              [ 18,  45,   0, -45], 
              [  9,   0, 126,   9],
              [-27, -45,   9, 135]]
    LU(matrix)

    print("第三次实验：")
    matrix = [[2, 2, 2],
              [4, 7, 7],
              [6, 18, 22]]
    LU(matrix)

'''
运行结果：
第一次实验：
===================
[4, -2, 0, 4]
[0, 1.0, -3.0, 3.0]
[0, 0, 4.0, 2.0]
[0, 0, 0, 9.0]
===================
===================
[1.0, 0, 0, 0]
[-0.5, 1.0, 0, 0]
[0.0, -3.0, 1.0, 0]
[1.0, 3.0, 0.5, 1.0]
===================
第二次实验：
===================
[9, 18, 9, -27]
[0, 9.0, -18.0, 9.0]
[0, 0, 81.0, 54.0]
[0, 0, 0, 9.0]
===================
===================
[1.0, 0, 0, 0]
[2.0, 1.0, 0, 0]
[1.0, -2.0, 1.0, 0]
[-3.0, 1.0, 0.6666666666666666, 1.0]
===================
第三次实验：
===================
[2, 2, 2]
[0, 3.0, 3.0]
[0, 0, 4.0]
===================
===================
[1.0, 0, 0]
[2.0, 1.0, 0]
[3.0, 4.0, 1.0]
===================

'''
    
