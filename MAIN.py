import scipy as sp
import math
data = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_ONE.txt",encoding='utf-8-sig') #data one is the real data
data_2 = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_TWO.txt",encoding='utf-8-sig') #data two is without the oscillation
E = sp.linspace(0.05,10.0,len(data))

def p(E,theta, m, L):
    return 1 - (sp.sin(2*theta)**2)*sp.sin((1.267*(m**2)*L)/E)**2

def Lam(Energy,theta,m): #lamda function
    index = math.floor(Energy/0.05) - 1
    return data_2[index]*p(Energy,theta,sp.sqrt(m),295)

def NLL(theta,m,data): #NLL function
    x = []
    for i in range(len(data)):
        if data[i] != 0: #only appends non zero values
            x.append(Lam(E[i],theta,m) - data[i] + data[i]*sp.log(data[i]/Lam(E[i],theta,m)))
    return sum(x)

def para_theta(y,m,function): #parabolic minimiser for theta
    
    x = sp.copy(y)
    x_3_old = x[1] #intialise
    y = function(x,2.4e-3,data)
    y = list(y)

    def x_3(x,y): # function to return new third point
        numerator = (x[2]**2-x[1]**2)*y[0] + (x[0]**2 - x[2]**2)*y[1] + (x[1]**2 - x[0]**2)*y[2]
        denom = (x[2] - x[1])*y[0] + (x[0] - x[2])*y[1] + (x[1] - x[0])*y[2]
        return (0.5*numerator)/denom
    
    def x_replace(x_3_new,x,y): # function that replaces the previous biggest point with the new point
    
        x[y.index(max(y))] = x_3_new
        return x
    
    x_3_new = x_3(x,y)
    x = x_replace(x_3(x,y),x,y)

    while abs(x_3_new - x_3_old) > 1e-5: #condition to stop for difference between sucessive points
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

def theta_error(minimum,Range,m,function,data):
    Min = min(Range, key=lambda y:abs(y - minimum)) # find in list of theta
    
    plus = Range[Range.index(Min)+1] #initialise

    while abs(function(Min,m,data) - function(plus,m,data)) < 0.5 : #condition to find where points differ by 0.5
        plus = Range[Range.index(plus)+1] 
  
    minus =Range[Range.index(Min) - 1] 
    
    while abs(function(Min,m,data) - function(minus,m,data)) < 0.5 :
        #print("index is", theta.index(minus), "minus is", minus, "difference", NLL(Min,data) - NLL(minus,data))
        minus = Range[Range.index(minus) -1 ] ##
    
    if abs(minimum - minus) > abs(minimum - plus): #choose the larger error
        return abs(minimum - minus)
    else:
        return abs(minimum - plus)
    
def m_error(minimum,Range,theta,function,data):
    Min = min(Range, key=lambda y:abs(y - minimum)) # find in list of theta
    
    plus = Range[Range.index(Min)+1] #initialise

    while abs(function(theta,Min,data) - function(theta,plus,data)) < 0.5 : #condition to find where points differ by 0.5
        plus = Range[Range.index(plus)+1] 
  
    minus =Range[Range.index(Min) - 1] 
    
    while abs(function(theta,Min,data) - function(theta,minus,data)) < 0.5 :
        #print("index is", theta.index(minus), "minus is", minus, "difference", NLL(Min,data) - NLL(minus,data))
        minus = Range[Range.index(minus) -1 ] ##
    
    if abs(minimum - minus) > abs(minimum - plus): #choose the larger error
        return abs(minimum - minus)
    else:
        return abs(minimum - plus)
    
    