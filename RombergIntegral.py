# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 19:42:39 2020
Copyright@ Xi'an Jiaotong University - NuTHeL
@author: Zhou Xingguang
"""

import numpy as np

class Romberg(object):
    def __init__(self, begin, end, N, func):
        # critera condition
        self.__MAX_STEP = 100
        self.__ERR = 1E-12
        
        self.begin = begin
        self.end = end
        self.N = N
        self.func = func
    
    # private method
    def __GetPoint(self):
        self.x = np.linspace(self.begin, self.end, self.N)
        self.y = self.func(self.x)
    
    def TrapeIntegral(self):
        self.__GetPoint()
        result = 0.0
        for i in range(self.N-1):
            result += (self.x[i+1]-self.x[i])*(self.y[i+1]+self.y[i]) / 2.0
        return result
    
    def RombergIntegral(self):
        R = 0
        T = []
        ''' Initial Romberg ''' 
        T.append(self.TrapeIntegral())
        self.N = 2*self.N-1
        T.append(self.TrapeIntegral())
        self.N = 2*self.N-1
        T.append(self.TrapeIntegral())
        self.N = 2*self.N-1
        T.append(self.TrapeIntegral())
        
        S = []
        S.append(T[1] + (T[1]-T[0])/3)
        S.append(T[2] + (T[2]-T[1])/3)
        S.append(T[3] + (T[3]-T[2])/3)
        
        C = []
        C.append(S[1] + (S[1]-S[0])/15)
        C.append(S[2] + (S[2]-S[1])/15)
        
        R = []
        print(T)
        print(self.N)
        R.append(C[1]+(C[1]-C[0])/63)
        # reserve the memory
        for i in range(self.__MAX_STEP):
            self.N = 2*self.N-1
            T[0:3] = T[1:4]
            T[3] = self.TrapeIntegral()
            S[0:2] = S[1:3]
            S[2] = T[3] + (T[3]-T[2])/3
            C[0] = C[1]
            C[1] = S[2] + (S[2]-S[1])/15
            R.append(C[1] + (C[1]-C[0])/63)
            if np.fabs(R[i+1] - R[i]) <= self.__ERR:
                print("Iteration has converged...\n")
                break
            print(T[3])
            print(self.N)
        return R
    

def func(x):
    return np.sin(x)/x


if __name__ == '__main__':
    r = Romberg(1e-16, 1, 2, func)
    result = r.RombergIntegral()
    print(result)
