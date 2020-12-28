# -*- coding: utf-8 -*-
"""
Created on Sun Dec 27 19:30:57 2020
Copyright@Xi'an Jiaotong University - NuTHeL
@author: Zhou Xingguang
"""

import numpy as np
import matplotlib.pyplot as plt

class CubicSpline(object):
    def __init__(self, x, y):
        self.x = x[:]
        self.y = y[:]
    
    @staticmethod
    def SecondOrderDiff(x, y):
        return np.array([((y[i+2]-y[i+1])/(x[i+2]-x[i+1])-(y[i+1]-y[i])/(x[i+1]-x[i]))/(x[i+2]-x[i]) for i in range(len(x)-2)])
    
    def Interpolation(self):
        h = np.array([self.x[i+1]-self.x[i] for i in range(len(self.x)-1)])
        u = np.array([h[i] / (h[i]+h[i+1]) for i in range(len(h)-1)])
        my_lambda = 1 - u
        d = 6 * CubicSpline.SecondOrderDiff(self.x, self.y)
        # generate the M matrix
        Matrix = 2*np.identity(np.size(d), dtype=float)
        Matrix[0, 1] = my_lambda[0]
        Matrix[np.size(d)-1, np.size(d)-2] = u[np.size(d)-1]
        for i in range(1, np.size(d)-1):
            Matrix[i, i-1] = u[i]
            Matrix[i, i+1] = my_lambda[i]
        # solve the function group
        M = np.linspace(0, 1, len(self.x))
        M[0] = 0
        M[1:len(self.x)-1] = np.linalg.solve(Matrix, d)
        M[len(self.x)-1] = 0
        return M, h
    
    def GenerateFunc(self, x):
        M, h = self.Interpolation()
        # generate the coefficient matrix
        coeff = np.zeros((np.size(self.x)-1, 4), dtype=float)
        for i in range(0, np.size(self.x)-1):
            coeff[i][0] = M[i] / (6*h[i])
            coeff[i][1] = M[i+1] / (6*h[i])
            coeff[i][2] = (self.y[i] - h[i]**2*M[i]/6)/h[i]
            coeff[i][3] = (self.y[i+1] - h[i]**2*M[i+1]/6)/h[i]
        
        # get the x in which section
        location = 0
        for i in range(len(self.x)-1):
            if self.x[i] <= x <= self.x[i+1]:
                break
            location += 1
        
        result = (self.x[location+1]-x)**3*coeff[location][0]
        result += (x-self.x[location])**3*coeff[location][1]
        result += (self.x[location+1]-x)*coeff[location][2]
        result += (x-self.x[location]) * coeff[location][3]
        return result


if __name__ == '__main__':
    x = np.linspace(-1, 1, 11)
    y = 1/(1+25*x**2)
    C = CubicSpline(x, y)
    for i in range(10):
        test_x = np.linspace(x[i], x[i+1], 5)
        test_y = []
        for j in range(0, 5):
            test_y.append(C.GenerateFunc(test_x[j]))
        plt.scatter(test_x[0], test_y[0])
        plt.scatter(test_x[1:5], test_y[1:5], marker='x')
    
    # draw the analytic solution
    plt.plot(test_x, test_y)
    
    Graphic_x = np.linspace(-1, 1, 100)
    Graphic_y = 1/(1+25*Graphic_x**2)
    plt.plot(Graphic_x, Graphic_y)
    plt.show()