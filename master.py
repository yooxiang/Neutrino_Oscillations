import matplotlib.pyplot as plt
import scipy as sp
import math
data = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_ONE.txt",encoding='utf-8-sig')
data_2 = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_TWO.txt",encoding='utf-8-sig')
#%%
R = sum(data)/sum(data_2)

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

theta = sp.linspace(0,sp.pi/2,500)
plt.plot(theta, NLL(theta,2.4e-3,data))
plt.xlabel("theta")
plt.ylabel("NLL")
#parabolic miinimeser

#%% 
guess = sp.array([0,0.5,0.7])
def para_theta(x,m,function):
    x_3_old = x[1]
    y = function(x,2.4e-3,data)

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

        y = list(function(x,2.4e-3,data))
        x_3_new = x_3(x,y)
        x = x_replace(x_3_new,x,y)
    x.sort()
    return x

para_theta(guess,NLL)


def para_m(x,theta,function):
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
#%% doing the error
def error(x,function, variable_min):
    Min = para_theta(x,function)[0]
    Min = min(variable_min, key=lambda y:abs(y - Min)) # find in list of theta
    
    plus = variable_min[variable_min.index(Min) +1] #


    while abs(NLL(Min,data) - NLL(plus,data)) < 0.5 :
        print("index is", variable_min.index(plus), "plus is", plus, "difference", NLL(Min,data) - NLL(plus,data))
        plus = variable_min[variable_min.index(plus)+1] ##
  

    minus =variable_min[variable_min.index(Min) - 1] 
    
    while abs(NLL(Min,data) - NLL(minus,data)) < 0.5 :
        print("index is", theta.index(minus), "minus is", minus, "difference", NLL(Min,data) - NLL(minus,data))
        minus = variable_min[variable_min.index(minus) -1 ] ##
        
    return (abs(Min - minus),abs(Min-plus))

theta = list(sp.linspace(0,sp.pi/2,500))
error(guess,NLL,theta)
#%% The univariate method
def univariate(theta_guess,m_guess,data,function):
    x = theta_guess
    y = function(theta,x,data)
    y = list(y)
    
    def x_3(x,y):
        numerator = (x[2]**2-x[1]**2)*y[0] + (x[0]**2 - x[2]**2)*y[1] + (x[1]**2 - x[0]**2)*y[2]
        denom = (x[2] - x[1])*y[0] + (x[0] - x[2])*y[1] + (x[1] - x[0])*y[2]
        return (0.5*numerator)/denom   
    
    def x_replace(x_3_new,x,y):
        x[y.index(max(y))] = x_3_new
        return x
 
    x = x_replace(x_3(x,y),x,y)
    
    theta_new = x
    theta_old = theta_guess
    m_new,m_old = m_guess,m_guess
    
    while abs(function(theta_new[0], m_new[0], data) - function(theta_old[0],m_old[0],data)) < 0.0001:
        print("new value is", function(theta_new[0], m_new[0], data))
        theta_old,m_old = theta_new,m_new
        
        theta_new = para_theta(theta_new,m_new[0],function)
        m_new = para_m(m_new,theta[0],function)

univariate()
        
    
        


    
