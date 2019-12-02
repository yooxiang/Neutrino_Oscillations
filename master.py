import matplotlib.pyplot as plt
import scipy as sp
import math
data = sp.loadtxt("DATA_ONE.txt",encoding='utf-8-sig')
data_2 = sp.loadtxt("DATA_TWO.txt",encoding='utf-8-sig')

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

def Lam(Energy,u):
    index = math.floor(Energy/0.05) - 1
    return data_2[index]*p(Energy,u,sp.sqrt(2.4e-3),295)
    
prediction = []
for i in range(len(data_2)):
    prediction.append(Lam(E[i],sp.pi/4))

plt.title("simulation")
plt.bar(E,prediction,width = 0.05)

#%%
def NLL(u,data):
    x = 0
    for i in range(len(data)):
        x += Lam((0.05)+i*0.05,u) - data[i] - data[i]*sp.log(data[i]/Lam((0.05)+i*0.05,u))
    return x
    
plt.figure()    
theta = sp.linspace(0,sp.pi/4,500)

plt.plot(theta,NLL(theta,data))
        
