#!/usr/bin/env python
#! -*- coding:UTF-8 -*-
# 日期: 2020/9/18
# 作者: 周星光
# 组织: XJTU NuTHeL
# 备注: 今天是九一八，勿忘国耻，共勉。

# 高斯消元法
# 本程序不使用任何第三方库(如NumPy等)，保持上课学习计算方法的原汁原味，
# 当使用NumPy时，可以直接利用其中的向量操作，所以在高斯消元时可以少写一个内层循环。
# 综上，也就是说当您阅读本程序时，会耗费更多的脑力，如果您今天感觉不适，请立即关闭电脑进行休息。

# 导入数学库，使用其中的fabs绝对值函数
import math

# 定义一个极小值
e = 0.00001

# 打印矩阵
def ShowMatrix(matrix):
    print("=================")
    row = len(matrix)
    # 按行打印
    for i in range(row):
        print(matrix[i][:])
    print("=================")

# 原生的二维矩阵运算函数
def GaussMethod(matrix, b):
    """ 
    matrix为系数矩阵，需要求解的矩阵
    b为方程等号右边的列向量，需要同步求解
    返回值为求解结果x 
    """
    # 矩阵的行数
    row = len(matrix)
    # 矩阵的列数
    col = len(matrix[0])
    # 定义用于存储结果的向量，默认结果全部为0
    x = [0] * col
    # 将matrix系数矩阵转化为上三角矩阵
    for i in range(row-1): # len(matrix)用于求解matrix矩阵的行数
        # 如果程序中主对角线元素近似于0，
        # 则说明该方程组的求解不适合使用高斯消元法，应当考虑列主元高斯消元法。
        if math.fabs(matrix[i][i]) <= e:
            print("Please recheck the triangle element, which maybe near zero...")
            return None
        for j in range(i+1, row):
            # 求解变换系数
            coeff = matrix[j][i] / matrix[i][i]
            # 开始进行行变换
            for k in range(col): # len(matrix[0]) 用于求出矩阵的列数
                matrix[j][k] = matrix[j][k] - coeff * matrix[i][k]
            # 等号右边的列变量同样也需要进行对应的行变换
            b[j] = b[j] - coeff * b[i]
    
    # 打印高斯消元后的临时结果
    ShowMatrix(matrix)
    print(b)
    
    ###############################
    # 接下来进行回代过程
    ###############################
    # 先求列向量x的最后一个元素的值
    x[len(x)-1] = b[col-1] / matrix[col-1][col-1]
    for i in range(row-2, -1, -1):
        for j in range(i, col-1):
            b[i] = b[i] - matrix[i][j+1] * x[j+1]
        x[i] = b[i] / matrix[i][i]
    
    return x

# 由于高斯消元法的应用条件比较苛刻，所以可以考虑更加具有普适性的列主元高斯消元法
def ColMaxGaussMethod(matrix, b):
    """
    列主元高斯消元法
    """
    # 矩阵的行数
    row = len(matrix)
    # 矩阵的列数
    col = len(matrix[0])
    # 存储方程组的解
    x = [0] * col
    # 进行第k次消元
    for k in range(row-1):
        # 开始选取列主元
        max_number = 0
        index = -1
        for i in range(k, row):
            # 注意此处进行比较的时绝对值
            if math.fabs(matrix[i][k]) > math.fabs(max_number):
                max_number = matrix[i][k]
                index = i # 记录是第几行的元素
        # 当这一列的主元为一个接近零的极小值时，退出程序，因为这个方程组无法进行求解。
        if math.fabs(max_number) <= e:
            print("WARNING! Please recheck the formula and identify it can be solved!")
            return None
        # 如果发现列上的最大元素不是当前对角线上的元素，
        # 则进行行交换
        if index != k:
            # 行交换
            for i in range(col):
                matrix[k][i], matrix[index][i] = matrix[index][i], matrix[k][i]
            # 等号右边的列向量同样需要进行行变换操作。
            b[k], b[index] = b[index], b[k]
        # 第k次消元前选列主元，交换列主元行变换完毕。
        # 下面开始进行当前行的高斯消元
        for j in range(k+1, row):
            # 求解变换系数
            coeff = matrix[j][k] / matrix[k][k]
            for i in range(col):
                matrix[j][i] = matrix[j][i] - coeff * matrix[k][i]
            # 等号右边的列变量同样也需要进行对应行变换
            b[j] = b[j] - coeff * b[k]
    
    # 打印高斯消元后的临时结果
    ShowMatrix(matrix)
    print(b)
    # 高斯消元完毕后，进行回代过程
    # 先求列向量x的最后一个元素的值
    x[len(x)-1] = b[col-1] / matrix[col-1][col-1]
    for i in range(row-2, -1, -1):
        for j in range(i, col-1):
            b[i] = b[i] - matrix[i][j+1]*x[j+1]
        x[i] = b[i] / matrix[i][i]
    
    return x



if __name__ == '__main__':
    test_a = [[1, 2, 1], [2, 2, 3], [-1, -3, 0]]
    b = [0, 3, 2]
    result = GaussMethod(test_a, b)
    print('result=', result)
    
    test_b = [[2, 2, 3], [1, 2, 1], [-1, -3, 0]]
    b = [3, 0, 2]
    print('测试列主元高斯消元法')
    result = ColMaxGaussMethod(test_b, b)
    print('result=', result)
    
    test_c = [[1, 2, 1], [2, 2, 3], [-1, -3, 0]]
    b = [0, 3, 2]
    print('测试列主元高斯消元法')
    result = ColMaxGaussMethod(test_c, b)
    print('result=', result)
