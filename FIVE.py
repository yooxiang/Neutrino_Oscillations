import matplotlib.pyplot as plt
import scipy as sp
import math
import random
from MAIN import p,Lam,NLL,para_theta,theta_error,m_error
data = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_ONE.txt",encoding='utf-8-sig') #data one is the real data
data_2 = sp.loadtxt("/Users/yuxiang/Documents/COMPUTING/PROJECT/DATA_TWO.txt",encoding='utf-8-sig') #data two is without the oscillations
E = sp.linspace(0.05,10.0,len(data)) #starting from 0.05, not 0
theta = sp.linspace(0,sp.pi/4,500)

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

x = cross_section(sp.array([0,0.5,0.7]),sp.array([0.002,0.003,0.004]),data,sp.array([0.1,0.5,1]))
theta_range = sp.linspace(0,sp.pi/4,500)
m_range = sp.linspace(0,0.004,500)
error_theta= theta_error(x[0][1],list(theta_range),x[1][1],NLL,data)
error_m =  m_error(x[1][1],list(m_range),x[0][1],NLL,data)


print("Univariate produces a value of", x[0][1],"for theta with an error of ", error_theta)
print("Univariate produces a value of", x[1][1]," for m with an an error of ", error_m)
print("univariate prodcues a value of", x[2][1], "for cross section")
