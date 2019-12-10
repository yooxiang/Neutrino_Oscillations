import matplotlib.pyplot as plt
import scipy as sp
import math
data = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_ONE.txt",encoding='utf-8-sig') #data one is the real data
data_2 = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_TWO.txt",encoding='utf-8-sig') #data two is without the oscillation

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


'''
Plots showing the data read in from the data given
'''
R = sum(data)/sum(data_2) #ratio matching because the data is the not the same length
E = sp.linspace(0.05,10.0,len(data)) #starting from 0.05, not 0
plt.figure()
plt.xlabel("energy (GeV)")
plt.ylabel("count")
plt.title("Date to fit")
plt.bar(E,data,width = 0.05)

plt.figure()
plt.title("Unoscillated flux")
plt.xlabel("Energy (GeV)")
plt.ylabel("count")
plt.bar(sp.linspace(0,10,len(data_2)),data_2,width=0.05) #data without oscillation

plt.figure()
plt.plot(E,p(E,sp.pi/4,sp.sqrt(2.4e-3),295))
plt.xlabel("Energy (GeV)")
plt.ylabel("probablity")
plt.title("model") #oscillation probabiity from Eqn. 1

'''
The convolution of Eqn. 1 and the unoscillated data
'''
prediction = []
for i in range(len(data_2)):
    prediction.append(Lam(E[i],sp.pi/4,2.4e-3))
plt.figure()
plt.title("Oscillated Event Rate Prediction")
plt.xlabel("Energy (GeV)")
plt.ylabel("count")
plt.bar(E,prediction,width = 0.05)

'''
Plot of the negative log likelihood 
'''
plt.figure()
theta = sp.linspace(0,sp.pi/4,500)
plt.plot(theta, NLL(theta,2.4e-3,data))
plt.title("Negative Log Likelihood for range of Theta")
plt.xlabel("theta")
plt.ylabel("NLL")

'''
The parabolic minimiser
'''
theta_guess = sp.array([0,0.5,0.7])
minimum = para_theta(theta_guess,2.4e-3,NLL)[1]
print("Parabolic Minimiser produces a value of ",minimum)

'''
Finding accuracy of results
'''
print("error in theta is", theta_error(minimum,list(theta),2.4e-3,NLL,data))
