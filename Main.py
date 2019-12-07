import matplotlib.pyplot as plt
import scipy as sp
import math
import random
data = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_ONE.txt",encoding='utf-8-sig')
data_2 = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_TWO.txt",encoding='utf-8-sig')

#%%
R = sum(data)/sum(data_2) #ratio matching

#function oscilation
def p(E,theta, m, L):
    return 1 - (sp.sin(2*theta)**2)*sp.sin((1.267*(m**2)*L)/E)**2

#plots
E = sp.linspace(0.05,10.0,len(data))
plt.figure()
plt.xlabel("energy (GeV)")
plt.ylabel("count")
plt.title("Unoscillated")
plt.bar(E,data,width = 0.05)


plt.figure()
plt.title("Unoscillated flux")
plt.xlabel("energy (GeV)")
plt.ylabel("count")
plt.bar(sp.linspace(0,10,len(data_2)),data_2,width=0.05)

plt.figure()
plt.plot(E,p(E,sp.pi/4,sp.sqrt(2.4e-3),295))
plt.xlabel("Energy (GeV)")
plt.ylabel("probablity")
plt.title("model")

plt.figure()

def Lam(Energy,theta,m):
    index = math.floor(Energy/0.05) - 1
    return data_2[index]*p(Energy,theta,sp.sqrt(m),295)
    
prediction = []
for i in range(len(data_2)):
    prediction.append(Lam(E[i],sp.pi/4,2.4e-3))

plt.title("simulation")
plt.bar(E,prediction,width = 0.05)

#%%
def NLL(theta,m,data):
    x = []
    for i in range(len(data)):
        if data[i] != 0:
            x.append(Lam(E[i],theta,m) - data[i] + data[i]*sp.log(data[i]/Lam(E[i],theta,m)))
    return sum(x)

theta = sp.linspace(0,sp.pi/4,500)
m = sp.linspace(0.001,0.004,500)


plt.plot(theta, NLL(theta,2.4e-3,data))
plt.xlabel("theta")
plt.ylabel("NLL")


plt.figure()
plt.plot(m, NLL(0.7,m,data))
plt.xlabel("m")
plt.ylabel("NLL")

#%% 
guess = sp.array([0,0.5,0.7])
guess_m = sp.array([0.001,0.002,0.004])
def para_theta(y,m,function):
    x = sp.copy(y)
    #print("x is", x, "m is", m)
    x_3_old = x[1]
    y = function(x,2.4e-3,data)

    y = list(y)
    #print("y is" , y)

    def x_3(x,y):
        numerator = (x[2]**2-x[1]**2)*y[0] + (x[0]**2 - x[2]**2)*y[1] + (x[1]**2 - x[0]**2)*y[2]
        denom = (x[2] - x[1])*y[0] + (x[0] - x[2])*y[1] + (x[1] - x[0])*y[2]
        return (0.5*numerator)/denom
    
    
    def x_replace(x_3_new,x,y):
    
        x[y.index(max(y))] = x_3_new
        return x
    
    #print("marker 7")
    x_3_new = x_3(x,y)
 
    x = x_replace(x_3(x,y),x,y)

   
    while abs(x_3_new - x_3_old) > 1e-5:
        x_3_old = x_3_new

        y = list(function(x,m,data))

        x_3_new = x_3(x,y)
        x = x_replace(x_3_new,x,y)
    x.sort()

    return x



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

para_theta(guess,2.4e-3,NLL)


para_m(guess_m, 0.7,NLL)


#%%
def univariate(theta_guess,m_guess,data,function):
    looper = 0
    
    theta_plot = []
    m_plot = []
    
    theta_plot.append(theta_guess[1])
    m_plot.append(m_guess[1])
    

    theta_old = sp.copy(theta_guess)

    m_old = sp.copy(m_guess)
   
    theta_new = para_theta(theta_old,m_old[1],function)

    m_new = para_m(m_old, theta_new[1], function)
 
    theta_plot.append(theta_new[1])
    m_plot.append(m_new[1])
  
    def difference(function,data, theta_old, theta_new, m_old, m_new):
        return abs(function(theta_old[1],m_old[1],data) - function(theta_new[1],m_new[1],data))
    
    
    
    while difference(function,data, theta_old, theta_new, m_old, m_new) > 0.00000000001:
        looper += 1
        print("this is loop number", looper)
        
        theta_old = sp.copy(theta_new)
        m_old = sp.copy(m_new)
        
        theta_new = para_theta(theta_old,m_old[1],function)
        m_new = para_m(m_old, theta_new[1], function)
        
        theta_plot.append(theta_new[1])
        m_plot.append(m_new[1])
        
    plt.figure()
    plt.plot(theta_plot,m_plot,'x')
    plt.xlabel("theta")
    plt.ylabel("m")
    
    plt.figure()
    plt.plot(sp.linspace(0,looper+2,looper+2), theta_plot,'x')
    #print(m_plot,theta_plot)
    print("number of times looped peforomed",looper)
    return theta_new,m_new
    
theta_guess = sp.array([0,0.5,0.7])
m_guess = sp.array([0.001,0.002,0.004])
univariate(theta_guess,m_guess,data,NLL)
        
#%% monte carlo

def Annealing(theta_guess,m_guess,data,function, theta_range,m_range):
    
    def bolt(E,T):
        K = 1.38064852e-23
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
            print("new value taken",looper)
            theta_plot.append(theta_old)
            m_plot.append(m_old)
            
        theta_new = random.random()*theta_range
        m_new = random.random()*m_range
        looper+=1
    plt.plot(theta_plot,m_plot)
    
    return theta_old,m_old
        
Annealing(0.7,0.2,data,NLL,sp.pi/4,0.5)

