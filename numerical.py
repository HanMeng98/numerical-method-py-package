# This is a package for numerical methods I have learnt in ME581 in Purdue
'''''''''
@author: Han Meng
@contact: HanMeng98@outlook.com
@file: numerical.py 
@time: 11/07/2019
@update: 11/07/2019
'''''''''
import numpy as np
import sympy

'''''''''
Newton's iteration method is used to find the root of f(x)=0
coeffs: 
    - start_point - start value
    - var         - the symbol variable, like x in f(x), t in f(t) etc
    - f           - sympy funciton f(x)
    - epoch       - numbers of iteration
    - TOL         - tolerance for the approximate root (>0)
returns:
    - p           - list of iteration output. p[-1] is the final approximate value
'''''''''
def newtonMethod(start_point, var, f, epoch = 0, TOL = -1):
    p = []
    p.append(start_point)
    df = sympy.diff(f, var)
    if epoch != 0 and TOL != -1:
        case = 0
    elif epoch != 0 and TOL == -1:
        case = 1
    elif epoch == 0 and TOL != -1:
        case = 2
    else:
        print("You need to set an ending condition!")

    # Case 0, choose epoch and tolerance to be the ending condition
    if case == 0:
        for i in range(epoch):
            assert df.subs(var,p[i]) != 0, "DENOMINATOR CANNOT BE 0!"
            p_tmp = p[i] - f.subs(var, p[i]) / df.subs(var, p[i])
            p.append(p_tmp)
            if abs(p[-1]-p[-2]) < TOL:
                print("Iteration finished ahead of epoches.")
                break
            
    # Case 1, choose epoch to be the ending condition
    if case == 1:
        for i in range(epoch):
            assert df.subs(var,p[i]) != 0, "DENOMINATOR CANNOT BE 0!"
            p_tmp = p[i] - f.subs(var, p[i]) / df.subs(var, p[i])
            p.append(p_tmp)
    
    # Case 2, choose tolerance to be the ending condition
    if case == 2:
        assert df.subs(var, p[-1]) != 0, "DENOMINATOR CANNOT BE 0!"
        p.append(p[-1] - f.subs(var, p[-1]) / df.subs(var, p[-1]))
        while(abs(p[-1] - p[-2]) > TOL):
            assert df.subs(var, p[-1]) != 0, "DENOMINATOR CANNOT BE 0!"
            p_tmp = p[-1] - f.subs(var, p[-1]) / df.subs(var, p[-1])
            p.append(p_tmp)
    return p

'''''''''
Secant Method is used to find the root of f(x)=0
coeffs: 
    - start_point1 - first iteration value
    - start_point2 - second iteration value
    - var          - the symbol variable, like x in f(x), t in f(t) etc
    - f            - sympy funciton f(x)
    - epoch        - numbers of iteration
    - TOL          - tolerance for the approximate root (>0)
returns:
    - x            - list of iteration output. x[-1] is the final approximate value
'''''''''
def secantMethod(start_point1, start_point2, var, f, epoch = 0, TOL = -1):
    x = []
    x.append(start_point1)
    x.append(start_point2)
    if epoch != 0 and TOL != -1:
        case = 0
    elif epoch != 0 and TOL == -1:
        case = 1
    elif epoch == 0 and TOL != -1:
        case = 2
    else:
        print("You need to set an ending condition!")

    # Case 0, choose epoch and tolerance to be the ending condition
    if case == 0:
        for i in range(epoch):
            x_tmp = x[-1] - f.subs(var, x[-1])*(x[-1]-x[-2])/(f.subs(var, x[-1])-f.subs(var, x[-2]))
            x.append(x_tmp)
            if abs(x[-1]-x[-2]) < TOL:
                print("Iteration finished ahead of epoches.")
                break
    # Case 1, choose epoch to be the ending condition
    if case == 1:
        for i in range(epoch):
            x_tmp = x[-1] - f.subs(var, x[-1])*(x[-1]-x[-2])/(f.subs(var, x[-1])-f.subs(var, x[-2]))
            x.append(x_tmp)
    # Case 2, choose tolerance to be the ending condition
    if case == 2:
        count = 0
        while(abs(x[-1]-x[-2]) > TOL):
            x_tmp = x[-1] - f.subs(var, x[-1])*(x[-1]-x[-2])/(f.subs(var, x[-1])-f.subs(var, x[-2]))
            x.append(x_tmp)
            count += 1
            if count > 1000:
                print("It have iterated over 1000 epoches.\n\r You may need to check the params.\n")
    return x