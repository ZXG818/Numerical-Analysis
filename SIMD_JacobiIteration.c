/* 
* Use the SIMD to accelerate the calculation speed of the 4X4 matrix in Jacobi Iteration
* Author       : Zhou Xingguang
* Organization : Xi'an Jiaotong University - NuTHeL
* Date         : 2020/11/13
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <xmmintrin.h>
#include <smmintrin.h>
#include <Windows.h>

#define STOP_STEP 50
#define EPSILON 1.0E-3


// 4x4 matrix
void JacobiIter(float (*matrix)[4], float *b, float *x)
{
	int i;
	float *x_front = NULL;
	__m128 _matrix_0 = _mm_loadu_ps(matrix[0]);
	__m128 _matrix_1 = _mm_loadu_ps(matrix[1]);
	__m128 _matrix_2 = _mm_loadu_ps(matrix[2]);
	__m128 _matrix_3 = _mm_loadu_ps(matrix[3]);
	__m128 _x = _mm_set1_ps(0);
	__m128 _b = _mm_loadu_ps(b);
	__m128 _x_front = _mm_set1_ps(0);
	__m128 _diagonal = _mm_set_ps(matrix[3][3],
		matrix[2][2],
		matrix[1][1],
		matrix[0][0]);

	// DEBUG USE
	__m128 temp;
	__m128 _sum0, _sum1, _sum2, _sum3;
	__m128 _sum_final;

	x_front = (float *)malloc(sizeof(float)*4);
	memset(x_front, 0, sizeof(float)* 4);

	_x = _mm_loadu_ps(x);
	_x_front = _mm_loadu_ps(x_front);

	for (i = 0; i < STOP_STEP; i++)
	{
		// only the SSE4.1 can use the dot multiply between two _m128 vectors
		temp = _mm_mul_ps(_matrix_0, _x_front);
		_sum0 = _mm_dp_ps(temp, _mm_set1_ps(1.0), 0xFF);

		temp = _mm_mul_ps(_matrix_1, _x_front);
		_sum1 = _mm_dp_ps(temp, _mm_set1_ps(1.0), 0xFF);

		temp = _mm_mul_ps(_matrix_2, _x_front);
		_sum2 = _mm_dp_ps(temp, _mm_set1_ps(1.0), 0xFF);

		temp = _mm_mul_ps(_matrix_3, _x_front);
		_sum3 = _mm_dp_ps(temp, _mm_set1_ps(1.0), 0xFF);

		// make the multiple result.
		__m128 m000 = _mm_move_ss(_sum1, _sum0);
		__m128 m111 = _mm_movelh_ps(m000, _sum2);
		__m128 m222 = _mm_shuffle_ps(m111, m111, _MM_SHUFFLE(0, 1, 2, 3));
		_sum_final = _mm_move_ss(m222, _sum3);
		_sum_final = _mm_shuffle_ps(_sum_final, _sum_final, _MM_SHUFFLE(0, 1, 2, 3));

		temp = _mm_add_ps(_mm_sub_ps(_b, _sum_final), _mm_mul_ps(_diagonal, _x_front));
		_x = _mm_div_ps(temp, _diagonal);

		// begin to judge the 2-Norm of the subtract vector between the solution vector and its front vector
		__m128 square = _mm_mul_ps(_mm_sub_ps(_x, _x_front), _mm_sub_ps(_x, _x_front));
		temp = _mm_dp_ps(square, _mm_set1_ps(1.0), 0xFF);
		__m128 sqrt = _mm_sqrt_ps(temp);
		// if the 2-Norm is less than 1E-3, then we get the solution, return the solution.
		if (sqrt.m128_f32[0] < EPSILON)
		{
			//printf("[*] 2-Norm is less than 1E-3...\n");
			goto end;
		}

		// record the front value
		_x_front = _x;

		// print the result of each iteration
		/*printf("Iteration : %d\n", i + 1);
		printf("x = %f, %f, %f, %f\n", _x.m128_f32[0],
			_x.m128_f32[1],
			_x.m128_f32[2],
			_x.m128_f32[3]);*/
	}
end:
	_mm_storeu_ps(x, _x);
	_mm_sfence();
	_ReadWriteBarrier();
	free(x_front);
	x_front = NULL;
	return;
}


int main(void)
{
	int i, j;
	float matrix[4][4] = { { 2.52, 0.95, 1.25, -0.85 },
	{ 0.39, 1.69, -0.45, 0.49 },
	{ 0.55, -1.25, 1.96, -0.98 },
	{ 0.23, -1.15, -0.45, 2.31 } };

	float b[4] = { 1.38, -0.34, 0.67, 1.52 };
	float x[4] = { 0 };   // make the initial value to zero

	// time consumption test
	SetThreadAffinityMask(GetCurrentThread(), 1);
	LARGE_INTEGER begin_JacobiIter;
	LARGE_INTEGER end_JacobiIter;
	LARGE_INTEGER freq;

	QueryPerformanceFrequency(&freq);
	QueryPerformanceCounter(&begin_JacobiIter);
	
	for (i = 0; i < 100; i++)
	{
		for (j = 0; j < 1048576; j++)
		{
			JacobiIter(matrix, b, x);
		}
	}
	
	QueryPerformanceCounter(&end_JacobiIter);
	
	// print the result
	printf("The final result: x = %f, %f, %f, %f\n", x[0], x[1], x[2], x[3]);

	printf("Time consumption:%f\n",
			(float)(end_JacobiIter.QuadPart - begin_JacobiIter.QuadPart) / (float)freq.QuadPart);

	return 0;
}