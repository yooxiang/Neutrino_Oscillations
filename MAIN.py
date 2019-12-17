import scipy as sp
import math
import random
data = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_ONE.txt",encoding='utf-8-sig') #data one is the real data
data_2 = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_TWO.txt",encoding='utf-8-sig') #data two is without the oscillation
E = sp.linspace(0.05,10.0,len(data))

def p(E,theta, m, L): #probablity function
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

    
    y = [] 
    for i in range(len(x)):
        y.append(function(x[i],2.4e-3,data))
    
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
        y = [] 
        for i in range(len(x)):
            y.append(function(x[i],m,data))
        y = list(y)
       
        
        x_3_new = x_3(x,y)
        x = x_replace(x_3_new,x,y)
    x.sort()
    return x

def theta_error(minimum,Range,m,function,data): #find error in theta
    Min = min(Range, key=lambda y:abs(y - minimum)) # find in list of theta
    
    plus = Range[Range.index(Min)+1] #initialise

    while abs(function(Min,m,data) - function(plus,m,data)) < 0.5 : #condition to find where points differ by 0.5
        plus = Range[Range.index(plus)+1] 
  
    minus =Range[Range.index(Min) - 1] 
    
    while abs(function(Min,m,data) - function(minus,m,data)) < 0.5 :
        
        minus = Range[Range.index(minus) -1 ] ##
    
    if abs(minimum - minus) > abs(minimum - plus): #choose the larger error
        return abs(minimum - minus)
    else:
        return abs(minimum - plus)


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
        theta_plot.append(theta_new[1])
        m_plot.append(m_new[1])
        
        m_new = para_m(m_old, theta_new[1], function)
        
        theta_plot.append(theta_new[1])
        m_plot.append(m_new[1])
        

    print("number of times looped peformed is ",looper)
    return theta_new,m_new,theta_plot,m_plot


def Annealing(theta_guess,m_guess,data,function, theta_range,m_range,step_theta,step_m):

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
   

    for T in range(3000,0,-1): # temperature starts at 3000
     
        if pacc(change(function,theta_old,theta_new,m_old,m_new),T/100) > 0.5: #probablity condition
            
            if theta_new < theta_range and theta_new > 0 and m_new < m_range and m_new > 0: #ensures in range
                theta_old = theta_new
                m_old = m_new
                #print("new value taken at loop number ",looper)
                theta_plot.append(theta_old)
                m_plot.append(m_old)
            

        theta_new = (-1 + random.random()*2)*step_theta + theta_old #new point 

        m_new = (-1 + random.random()*2)*step_m  + m_old # new point

        looper+=1


    return theta_old,m_old,theta_plot,m_plot
        
def cross_section(theta_guess,m_guess,data,cross_guess):
    #previous functions with cross parameter added
    def p_cross(E,theta, m, L, cross):
        return (cross*E)*(1 - (sp.sin(2*theta)**2)*sp.sin((1.267*(m**2)*L)/E)**2)
    

    def Lam_cross(Energy,theta,m,cross):
        index = math.floor(Energy/0.05) - 1
        #print("m is", m)
        if sp.isreal((p_cross(Energy,theta,sp.sqrt(m),295,cross))) == False:
            print("lam is fake")
        #print(p_cross(Energy,theta,sp.sqrt(m),295,cross))
        #print("new")
        return data_2[index]*p_cross(Energy,theta,sp.sqrt(m),295,cross)

    def NLL_cross(theta,m,data,cross):
        x = []
        for i in range(len(data)):
            if data[i] != 0:
                x.append(Lam_cross(E[i],theta,m,cross) - data[i] + data[i]*sp.log(data[i]/Lam_cross(E[i],theta,m,cross)))
                if sp.isreal(sp.log(data[i]/Lam_cross(E[i],theta,m,cross))) == False:
                    print("the log is fake",data[i],Lam_cross(E[i],theta,m,cross),data[i]/Lam_cross(E[i],theta,m,cross))
        if sp.isreal(sum(x)) == False:
            print("FAKE")
            print("values are ",theta,m,cross)
        return sum(x)
    
    function = NLL_cross

    def para_theta_cross(z,m,function,cross):
        x = sp.copy(z)
        #print("x is", x, "m is", m)
        x_3_old = x[1]
        y =[]
        for i in range(len(x)):
            y.append(function(x[i],m,data,cross))
        #print("x is",x,cross)
        #print("y is", y)
        y = list(y) ##
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
            y =[]
            for i in range(len(x)):
                y.append(function(x[i],m,data,cross))
        
            y = list(y)
    
            x_3_new = x_3(x,y)
            x = x_replace(x_3_new,x,y)
        x.sort()
    
        return x



    def para_m_cross(y,theta,function,cross):
        x = sp.copy(y)
        x_3_old = x[1]
    
        
        y =[]
        for i in range(len(x)):
            y.append(function(theta,x[i],data,cross))
        y = list(y) ##
    
        
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
    
            y =[]
            for i in range(len(x)):
               y.append(function(theta,x[i],data,cross))
            y = list(y) ##
            
            x_3_new = x_3(x,y)
            x = x_replace(x_3_new,x,y)
        x.sort()

        return x
    
    def para_cross_cross(theta,m,function,cross):
       
        x = sp.copy(cross)
        #print("x is", x, "m is", m)
        x_3_old = x[1]
        
        y =[]
        for i in range(len(x)):
            y.append(function(theta,m,data,x[i]))
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
    
            y =[]
            for i in range(len(x)):
                y.append(function(theta,m,data,x[i]))
            y = list(y) ##
    
            x_3_new = x_3(x,y)
            x = x_replace(x_3_new,x,y)
        x.sort()

        return x
    

    def univariate_cross(theta_guess,m_guess,data,function,cross_guess):
    
        looper = 0
        
        theta_plot = []
        m_plot = []
        cross_plot = []
        
        theta_plot.append(theta_guess[1])
        m_plot.append(m_guess[1])
        cross_plot.append(cross_guess[1])
        

        theta_old = sp.copy(theta_guess)

        m_old = sp.copy(m_guess)
        
        cross_old = sp.copy(cross_guess)
       
        theta_new = para_theta_cross(theta_old,m_old[1],function,cross_guess[1])

        m_new = para_m_cross(m_old, theta_new[1], function,cross_guess[1])

        cross_new = para_cross_cross(theta_new[1],m_new[1],function,cross_old)
        

     
        theta_plot.append(theta_new[1])
        m_plot.append(m_new[1])
        cross_plot.append(cross_new[1])
        

        
 
        def difference(function,data, theta_old, theta_new, m_old, m_new, cross_old,cross_new):
            return abs(function(theta_old[1],m_old[1],data,cross_old[1]) - function(theta_new[1],m_new[1],data,cross_new[1]))
        
        
        while difference(function,data, theta_old, theta_new, m_old, m_new,cross_old,cross_new) > 0.0000001:
 
            looper += 1
         
            
            theta_old = sp.copy(theta_new)
            m_old = sp.copy(m_new)
            cross_old = sp.copy(cross_new)
            
            theta_new = para_theta_cross(theta_old,m_old[1],function,cross_old[1])
            m_new = para_m_cross(m_old, theta_new[1], function,cross_old[1])
            cross_new = para_cross_cross(theta_new[1],m_new[1],function,cross_old)
            
            cross_plot.append(cross_new[1])
            theta_plot.append(theta_new[1])
            m_plot.append(m_new[1])
      
            
       
        print("number of times looped peforomed",looper)
        
        
        prediction  = []
        for i in range(len(data_2)):
            prediction.append(Lam_cross(E[i],theta_new[1],m_new[1],cross_new[1]))
        plt.title("Oscillated Event Rate Prediction with Cross Section")
        plt.xlabel("Energy (GeV)")
        plt.ylabel("Count")
        plt.bar(E,prediction,width = 0.05)
      
     
        return theta_new,m_new,cross_new
    
    
   
    return univariate_cross(theta_guess,m_guess,data,function,cross_guess)
    