# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 14:54:24 2020
@copyright: Xi'an Jiaotong University - NuTHeL
@author: Zhou Xingguang
"""

import numpy as np
import copy

class ConjugateGradient(object):
    def __init__(self, matrix, b, x):
        # define the little number
        self.__ERR = 1E-6  # private parameter
        # define the MAX_STEP
        self.__MAX_STEP = 1000 # private parameter
        self.matrix = np.array(copy.deepcopy(matrix))
        self.b = np.array(copy.deepcopy(b)).reshape(np.size(b), 1)
        self.x = np.array(copy.deepcopy(x)).reshape(np.size(b), 1)
        self.last_x = self.x[:, :]
        
        self.r = np.array(np.zeros((np.size(b), 1), dtype=float))
        self.last_r = self.r[:, :]
        self.d = np.array(np.zeros((np.size(b), 1), dtype=float))
        self.last_d = self.d[:, :]
        
        self.alpha = 0
        self.beta = 0
        self.time = 0
        
    
    def CG(self):
        self.r = self.b - np.matmul(self.matrix, self.x)
        self.last_r = self.r[:, :]
        self.d = self.r[:, :]
        for time in range(self.__MAX_STEP):
            self.alpha = np.linalg.norm(self.r, ord=2)**2 / np.matmul(np.matmul(np.transpose(self.d), self.matrix), self.d)
            self.last_x = self.x[:, :]
            self.x = self.x + self.alpha * self.d
            self.last_r = self.r[:, :]
            self.r = self.b - np.matmul(self.matrix, self.x)
            self.beta = (np.linalg.norm(self.r, ord=2)**2) / (np.linalg.norm(self.last_r, ord=2)**2)
            self.d = self.r + self.beta*self.d
            self.alpha = (np.linalg.norm(self.r)**2) / np.matmul(np.matmul(np.transpose(self.d), self.matrix), self.d)
            if np.linalg.norm(self.r) <= self.__ERR:
                print("iteration has been converged...\n")
                self.time = time + 1
                return self.x
        print("The result is not converged, please recheck the method")
        return None

def GenerateExercise(matrix, b):
    cnt = np.size(b)
    b[0, 0] = -1
    b[cnt-1, 0] = -1
    
    for i in range(cnt):
        if i == 0:
            matrix[i, i] = -2
            matrix[i, i+1] = 1
        elif i == cnt-1:
            matrix[i, i] = -2
            matrix[i, i-1] = 1
        else:
            matrix[i, i] = -2
            matrix[i, i-1] = 1
            matrix[i, i+1] = 1
    return matrix, b


if __name__ == '__main__':
    cnt = 400
    matrix = np.zeros((cnt, cnt), dtype=float)
    x = np.zeros((cnt, 1), dtype=float)
    b = np.zeros((cnt, 1), dtype=float)
    
    matrix, b = GenerateExercise(matrix, b)
    
    CG = ConjugateGradient(matrix, b, x)
    result = CG.CG()
    print(result)
    print(CG.time)
    
    