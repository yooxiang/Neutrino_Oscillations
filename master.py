import matplotlib.pyplot as plt
import scipy as sp
data = sp.loadtxt("DATA_ONE.txt",encoding='utf-8-sig')
data_2 = sp.loadtxt("DATA_TWO.txt",encoding='utf-8-sig')

R = sum(data)/sum(data_2)

#function
def p(E,theta, m, L):
    return 1 - (sp.sin(2*theta)**2)*sp.sin((1.267*(m**2)*L)/E)**2

#plots
E = sp.linspace(0,10,len(data))
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
def lamda(E,theta): #returns the prediction as a function of E
    return R*(E*p(0.05*i,theta,sp.sqrt(2.4e-3),295))
    
prediction = []
for i in range(len(data_2)):
    prediction.append(lamda(data_2[i],sp.pi/4))
plt.bar(E,prediction,width= 0.05)

def NLL(u,data):
    x = 0
    for i in range(len(data)):
        x += lamda(i*0.05,u) - data[i] - data[i]*sp.log(data[i]/lamda(i*0.05,u))
    
plt.figure()    
theta = sp.linspace(0,sp.pi,100)
plt.bar(theta,NLL(theta,data))
        
