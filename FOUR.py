import matplotlib.pyplot as plt
import scipy as sp
import math
import random
from MAIN import p,Lam,NLL,para_theta,theta_error,m_error
data = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_ONE.txt",encoding='utf-8-sig') #data one is the real data
data_2 = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_TWO.txt",encoding='utf-8-sig') #data two is without the oscillations

'''
parabolic minimiser for m
'''

def para_m(y,theta,function):
    x = sp.copy(y)
    x_3_old = x[1]

    y = function(theta,x,data)
    y = list(y)

    def x_3(x,y):
        numerator = (x[2]**2-x[1]**2)*y[0] + (x[0]**2 - x[2]**2)*y[1] + (x[1]**2 - x[0]**2)*y[2]
        denom = (x[2] - x[1])*y[0] + (x[0] - x[2])*y[1] + (x[1] - x[0])*y[2]
        return (0.5*numerator)/denom   
    
    def x_replace(x_3_new,x,y):
        x[y.index(max(y))] = x_3_new
        return x
    
    x_3_new = x_3(x,y)
 
    x = x_replace(x_3(x,y),x,y)

    
    while abs(x_3_new - x_3_old) > 1e-5:
        x_3_old = x_3_new

        y = list(function(theta,x,data))
        x_3_new = x_3(x,y)
        x = x_replace(x_3_new,x,y)
    x.sort()
    return x

'''
univariate function
'''

def univariate(theta_guess,m_guess,data,function):
    print("start univariate")
    looper = 0
    
    theta_plot = []
    m_plot = []
    
    theta_plot.append(theta_guess[1])
    m_plot.append(m_guess[1])
    
    theta_old = sp.copy(theta_guess)
    m_old = sp.copy(m_guess)
   
    theta_new = para_theta(theta_old,m_old[1],function) #theta is minimised first
    m_new = para_m(m_old, theta_new[1], function)
 
    theta_plot.append(theta_new[1])
    m_plot.append(m_new[1])
  
    def difference(function,data, theta_old, theta_new, m_old, m_new): #difference function for difference in NLL between loops
        return abs(function(theta_old[1],m_old[1],data) - function(theta_new[1],m_new[1],data))
    
    
    while difference(function,data, theta_old, theta_new, m_old, m_new) > 0.00000000001: #condition to stop
        looper += 1
        #print("this is loop number", looper)
        
        theta_old = sp.copy(theta_new)
        m_old = sp.copy(m_new)
        
        theta_new = para_theta(theta_old,m_old[1],function)
        m_new = para_m(m_old, theta_new[1], function)
        
        theta_plot.append(theta_new[1])
        m_plot.append(m_new[1])
        
    
    print("number of times looped peforomed is ",looper)
    return theta_new,m_new,theta_plot,m_plot

'''
Annleaing Function
'''
def Annealing(theta_guess,m_guess,data,function, theta_range,m_range):
    print("Start Annealing")
    def bolt(E,T):
        K = 1.38064852e-2
        return sp.exp(-E/(K*T))

    def pacc(E,T):
        if E <= 0:
            return 1
        if E > 0:
            return bolt(E,T)
        
    def change(function,theta_old,theta_new,m_old,m_mew):
        return function(theta_new,m_new,data) - function(theta_old,m_old,data)
    
    theta_plot = []
    m_plot = []
       
    theta_old = theta_guess
    m_old = m_guess
    
    theta_new = random.random()*theta_range
    m_new = random.random()*m_range
    
    looper = 0
    
    for T in range(10000,0,-1):
        if pacc(change(function,theta_old,theta_new,m_old,m_new),T/100) > 0.5:
            theta_old = theta_new
            m_old = m_new
            #print("new value taken at loop number ",looper)
            theta_plot.append(theta_old)
            m_plot.append(m_old)
            
        
        theta_new = random.random()*0.1 + theta_old 
        m_new = random.random()*0.1 + theta_new
        looper+=1


    return theta_old,m_old,theta_plot,m_plot
        


'''
Plot of the univariate method
'''    
theta_guess = sp.array([0,0.5,0.7])
m_guess = sp.array([0.001,0.002,0.004])

x = univariate(theta_guess,m_guess,data,NLL)

theta = x[0][1]
m = x[1][1]
theta_plot = x[2]
m_plot = x[3]
plt.plot(theta_plot,m_plot)
plt.xlabel("$\Theta_{23}$")
plt.ylabel("$\Delta m_{23}^2$")
plt.title("Univariate Method for Minimising $\Delta m_{23}^2$ and $\Theta_{23}$ ")

print("Univariate produces a value of", theta,"for theta")
print("Univariate produces a value of'", m," for m")


'''
Plot of the annleaing method

'''
theta_guess = 0.7
m_guess = 0.2
x = Annealing(theta_guess,m_guess,data,NLL,sp.pi/4,0.5)

plt.figure()
theta = x[0]
m = x[1]
theta_plot = x[2]
m_plot = x[3]
plt.plot(theta_plot,m_plot)
plt.xlabel("$\Theta_{23}$")
plt.ylabel("$\Delta m_{23}^2$")
plt.title("Annealing Method for Minimising $\Delta m_{23}^2$ and $\Theta_{23}$ ")

theta_range = sp.linspace(0,sp.pi/4,500)
m_range = sp.linspace(0,0.004,500)


error  = theta_error(theta,list(theta_range),m,NLL,data), m_error(m,list(m_range),theta,NLL,data)

print("Annealing produces a value of", theta,"for theta with an error of ", error[0])
print("Annealing produces a value of", m," for m with an an error of ", error[1])






'''
Validation of the anneleaing
'''
def f(x,y,data):
    return x*sp.sin(x-3*sp.pi/2)*sp.sin(y-3*sp.pi/2) #function which to find minumum

'''
Heatmap of function
Black region shows region of minimum
'''
heat = sp.zeros([100,100])
x_plot = sp.linspace(0,2*sp.pi,100)
y_plot = sp.linspace(0,2*sp.pi,100)

for i in range(100):
    for j in range(100):
        if f(x_plot[i],y_plot[j],data) < 0:
            heat[i][j] = -abs(f(x_plot[i],y_plot[j],data))
        else:
            heat[i][j] = 0
plt.figure()
plt.yticks([0,50,99], ["0", "$\pi$", "2$\pi$"])
plt.xticks([0,50,99], ["0", "$\pi$", "2$\pi$"])
a = heat.transpose()
plt.imshow(a, cmap='hot', interpolation='nearest')
plt.title("Heat Map of function f(x,y)")
plt.ylabel("y")
plt.xlabel('X')
plt.show()

x =  Annealing(sp.pi/2,sp.pi/2,data,f,2*sp.pi,2*sp.pi)

plt.figure()
x_end = x[0]
y_end = x[1]
x_anneal = x[2]
y_anneal = x[3]

plt.plot(x_anneal,y_anneal)
plt.xlabel("x")
plt.ylabel("y")
plt.xlim(0,2*sp.pi)
plt.ylim(0,2*sp.pi)
plt.title("Annealing Method for Minimising f(x,y)")

print("x is", x_end)
print("y is", y_end)
print("Z is", f(x_end,y_end,data))





        

