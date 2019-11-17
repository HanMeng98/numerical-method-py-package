# This is a package for numerical methods I have learnt in ME581 in Purdue
'''''''''
@author: Han Meng
@contact: HanMeng98@outlook.com
@file: numerical.py 
@time: 11/07/2019
@update: 11/16/2019
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
        count = 0
        for i in range(epoch):
            count += 1
            assert float(df.subs(var,p[i])) != 0, "DENOMINATOR CANNOT BE 0!"
            p_tmp = p[i] - float(f.subs(var, p[i])) / float(df.subs(var, p[i]))
            p.append(p_tmp)
            if abs(p[-1]-p[-2]) < TOL:
                print("Iteration finished ahead of epoches.")
                break
            
    # Case 1, choose epoch to be the ending condition
    if case == 1:
        count = 0
        for i in range(epoch):
            assert float(df.subs(var,p[i])) != 0, "DENOMINATOR CANNOT BE 0!"
            count += 1
            p_tmp = p[i] - float(f.subs(var, p[i]) / df.subs(var, p[i]))
            p.append(p_tmp)
    
    # Case 2, choose tolerance to be the ending condition
    if case == 2:
        count = 0
        assert float(df.subs(var, p[-1])) != 0, "DENOMINATOR CANNOT BE 0!"
        p.append(p[-1] - float(f.subs(var, p[-1]) / df.subs(var, p[-1])))
        while(abs(p[-1] - p[-2]) > TOL):
            assert float(df.subs(var, p[-1])) != 0, "DENOMINATOR CANNOT BE 0!"
            count += 1
            p_tmp = p[-1] - float(f.subs(var, p[-1]) / df.subs(var, p[-1]))
            p.append(p_tmp)
            if count>100: break
    
    if count <= 100 or count == epoch:
        return p[-1]
    else:
        return None
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
    
'''''''''
Euler Backward Method is used to solve the Inital value problem
It is an implicit method using Newton Method. It is more accurate but cost more time.
coeffs: 
    - func         - dy/dx = func
    - step         - iteration step
    - start_point  - x interval left side
    - final_point  - x interval right side
    - y0           - initial value y(0)
returns:
    - y_list       - list of approximate values evaluated at each step
'''''''''

def eulerBackwardIVP(func, step, start_point, final_point, y0):
    n = round((final_point - start_point)/step)
    x_list = list(start_point+step*i for i in range(n+1))
    y_list = [y0]
    y = sympy.symbols("y")
    for i in range(n):
        f = y - (y_list[i]+step*func(x_list[i+1],y))
        y_next = newtonMethod(0,y,f,TOL=1e-2)
        if y_next == None:
            break
        else:
            y_list.append(y_next)
    return y_list
'''''''''
Euler Method is used to solve the Inital value problem
It is an explicit method with lower accuracy.
coeffs: 
    - func         - dy/dx = func
    - step         - iteration step
    - start_point  - x interval left side
    - final_point  - x interval right side
    - y0           - initial value y(0)
returns:
    - y_list       - list of approximate values evaluated at each step
'''''''''
def eulerMethodIVP(func, step, start_point, final_point, y0):
    n = round((final_point - start_point)/step)
    x_list = list(start_point+step*i for i in range(n+1))
    y_list = [y0]
    for i in range(n):
        y_next = y_list[i] + step*func(x_list[i], y_list[i])
        y_list.append(y_next)
    return y_list