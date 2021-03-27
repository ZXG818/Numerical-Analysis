/* 
* 二维矩阵开方计算程序
* 作者：ZXG/watermelon/trace-shadow
* 转载请注明出处：[url]www.0xaa55.com[/url]
* 动手实践，刨根问底，实事求是， 不辞辛苦
* 向弹幕流大佬致敬！
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
// 定义误差
#define myErrno 0.000001
#define in
#define out
 
 
// 求解一元二次方程
// 说明：a为二次项系数，b为一次项系数，c为常数项
void SolveFunction(in double a, in double b, in double c, out double *result)
{
    if (fabs(a - 0) <= myErrno)
    {
        printf("求解一元二次方程，二次项系数为零，只有单根\n");
        result[0] = -c / b;
        result[1] = 0;
        return;
    }
 
    // 二次方程利用求根公式进行求解
    // 先检查b2-4ac
    double delta = b*b - 4 * a*c;
 
    if (delta <= myErrno)
    {
        printf("无解\n");
        result[0] = 0;
        result[1] = 0;
        return;
    }
 
    // 如果有解，则用求根公式进行求解
    result[0] = (-b + sqrt(delta)) / (2 * a);
    result[1] = (-b - sqrt(delta)) / (2 * a);
 
}
 
// 打印二维矩阵各个元素
void PrintMatrix(in int line, in int col, in double matrix[][2])
{
    printf("**************************\n");
    for (int i = 0; i < line; i++)
    {
        for (int j = 0; j < col; j++)
        {
            printf("%f\t", matrix[i][j]);
        }
        printf("\n");
    }
    printf("**************************\n");
}
 
// 两个二维矩阵相乘
// 顺序：matrix1左乘matrix2
void MulMatrix(in double matrix1[][2], in double matrix2[][2], out double result[][2])
{
    result[0][0] = matrix1[0][0] * matrix2[0][0] + matrix1[0][1] * matrix2[1][0];
    result[0][1] = matrix1[0][0] * matrix2[0][1] + matrix1[0][1] * matrix2[1][1];
    result[1][0] = matrix1[1][0] * matrix2[0][0] + matrix1[1][1] * matrix2[1][0];
    result[1][1] = matrix1[1][0] * matrix2[0][1] + matrix1[1][1] * matrix2[1][1];
}
 
 
// 求二维矩阵行列式的值
double GetDet(in double matrix[][2])
{
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
}
 
 
// 由于本程序讨论的为二维方阵，所以在求解特征向量的时候就会方便很多。
// 如果原矩阵可以相似对角化，则在有两个不同特征值的条件下，
// 其特征向量的求解只需要关注矩阵(A-λE)的第一维即可，
// 为什么？因为其(A-λE)的秩需要为1才可以相似对角化。
 
// 获取特征值对应的特征向量
// 这里的特征向量为列向量
void GetEigenvector(in double matrix[][2], in double eigenvalue, out double *result)
{
    result[0] = matrix[0][0] - eigenvalue;
    result[1] = -matrix[0][1];
}
 
// 求解二维矩阵特征值
void GetEigenvalue(in double matrix[][2], out double *eigenvalue)
{
    double a = matrix[0][0];
    double b = matrix[0][1];
    double c = matrix[1][0];
    double d = matrix[1][1];
 
    // 求解一元二次方程
    SolveFunction(1, -a - d, a*d - b*c, eigenvalue);
 
}
 
 
// 二维矩阵求逆
// 二维矩阵求逆有一个非常简便的方法
// 首先求矩阵的行列式的值，取其倒数作为系数放在矩阵前面
// 然后主对角线元素互换位置，副对角线元素取其相反数即可。
void InverseMatrix(in out double matrix[][2])
{
    double det = GetDet(matrix);
 
    double temp[2][2] = { { matrix[0][0], matrix[0][1] },
                          { matrix[1][0], matrix[1][1] } };
 
    matrix[0][0] = (1 / det)*temp[1][1];
    matrix[1][1] = (1 / det)*temp[0][0];
    matrix[1][0] = (1 / det)*(-temp[1][0]);
    matrix[0][1] = (1 / det)*(-temp[0][1]);
}
 
 
 
int main(int argc, char *argv[])
{
    // 举例矩阵
    double matrix[2][2] = { { 4, -2 },
                            { 1,  1 } };
 
    printf("\n初始矩阵为：\n");
    PrintMatrix(2, 2, matrix);
 
 
    // 先求解特征值
    double eigenvalue[2] = { 0 };
    GetEigenvalue(matrix, eigenvalue);
 
    printf("矩阵的特征值：%f\t%f\n", eigenvalue[0], eigenvalue[1]);
 
    if (eigenvalue[0] == eigenvalue[1])
    {
        printf("稍后分析...\n");
        return 0;
    }
 
    // 求特征向量组成的矩阵
    double eigenmatrix[2][2] = { 0 };
    // 先求第一个特征值的特征向量
    double temp[2] = { 0 };
 
    // 特征向量为列向量
    // 这里求特征向量的时候建议在纸上动手算算！！！
    GetEigenvector(matrix, eigenvalue[0], temp);
    eigenmatrix[0][0] = temp[1];    
    eigenmatrix[1][0] = temp[0];
 
    GetEigenvector(matrix, eigenvalue[1], temp);
    eigenmatrix[0][1] = temp[1];
    eigenmatrix[1][1] = temp[0];
 
    printf("\n特征向量组成的矩阵为：\n");
    PrintMatrix(2, 2, eigenmatrix);
 
    // 求特征向量组成的矩阵的逆矩阵
    double inverse_eigenmatrix[2][2] = { { eigenmatrix[0][0], eigenmatrix[0][1] },
                                         { eigenmatrix[1][0], eigenmatrix[1][1] } };
    InverseMatrix(inverse_eigenmatrix);
 
    printf("\n特征向量组成的矩阵的逆矩阵为：\n");
    PrintMatrix(2, 2, inverse_eigenmatrix);
 
 
    // 对已经完成相似对角化后得到的对角矩阵进行开方
    double diag[2][2] = { { sqrt(eigenvalue[0]), 0 },
                          { 0, sqrt(eigenvalue[1]) } };
 
    // 最后求解开方后的矩阵，终于迎来了最终的结果
    double result[2][2] = { 0 };
 
    MulMatrix(eigenmatrix, diag, result);
    double copy_result[2][2] = { 0 };
 
    // 把result的内容复制到temp中，方便接下来的矩阵相乘
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            copy_result[i][j] = result[i][j];
        }
    }
 
    MulMatrix(copy_result, inverse_eigenmatrix, result);
 
    printf("\n最终的结果矩阵为：\n");
    PrintMatrix(2, 2, result);
 
 
    // 验证结果矩阵
    double test[2][2] = { 0 };
    MulMatrix(result, result, test);
    printf("\n验证结果矩阵，两个结果矩阵相乘得到的矩阵为：\n");
    PrintMatrix(2, 2, test);
    return 0;
}