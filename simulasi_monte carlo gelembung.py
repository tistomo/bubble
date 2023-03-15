# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 10:23:31 2023

@author: arfan
"""

import numpy as np
import matplotlib.pyplot as plt
import random

#simulasi untuk 0.1 lpm
n=10000 #jumlah sampel
a=[0 for x in range(n)]    #parameter 1
b=[0 for x in range(n)]    #parameter 2
d=[0 for x in range(n)]    #parameter 3 /diameter 
for i in range(n):
    a[i]=random.uniform(1,2)
    b[i]=random.uniform(0,1)
    d[i]=-np.log(a[i]*b[i])/4+0.249
    
plt.style.use('ggplot')
plt.hist(d)
plt.xlabel('diameter, mm')
plt.ylabel('frequency')
plt.show()  
mean=np.mean(d)
plt.close

#mencari panjang heat exchanger
#t_amb=[x for x in range(n)]   # suhu ruangan, celcius
ti=[x for x in range(n)]        # suhu bak air
tsat=[x for x in range(n)]      #suhu saturator / target
l=[x for x in range(n)]      #panjang tubing
r_out=9.5/2 /1000             #jari jari tubing luar, m
r_in=r_out-0.001                #jari jari tubing dalam, m
delta_tin=[x for x in range(n)] # fungsi logaritmik
alpha_in=10 #convective heat transfer coefficient of air, W/m^2/K
alpha_out=2000 #convective heat transfer coefficient of water, W/m^2/K                      
flow_rate= 5 #flow rate udara, lpm
m_flowrate=1.204*flow_rate/60/1000   #mass flow rate, kg/s
cp_a=1.005*1000  #specific heat of gas, J/kg K
l_amda=15 #thermal conductivity of stainless steel, W/K/m
ho=1/(alpha_out*2*np.pi*r_out) #hambatan panas dari OD
hi=1/(alpha_in*2*np.pi*r_in) #hambatan panas dari ID
k=np.log(r_out/r_in)/(l_amda*2*np.pi) #hambatan panas dari termal konductivity stainless steel
A=ho+hi+k
for i in range(n):
  t_amb[i]=random.uniform(23,24)      #fluktuasi 3  deg C dengan rectangular distribusi
  ti[i]=random.uniform(-0.02,0.02)+ random.uniform(-0.05,0.05)+random.normalvariate(24.89,0.05)
  tsat[i]=ti[i]-0.001
  delta_tin[i]=(t_amb[i]-tsat[i])/np.log((t_amb[i]-ti[i])/(tsat[i]-ti[i]))
  l[i]=m_flowrate*cp_a*(t_amb[i]-tsat[i])/delta_tin[i]*A# panjang pipa stainless steel yang diperlukan
#plt.hist(l)
#plt.xlabel('length, m')
#plt.ylabel('frequency')
#plt.show

#mencari kecepatan gelembung
r=[x for x in range(n)]    #parameter 4 /jari jari
Fa=[x for x in range(n)]   #gaya archimedes
V=[x for x in range(n)]     #volume gelembung
Psat=[x for x in range(n)]  #tekanan saturator
PL=[x for x in range(n)]  #tekanan air pada kedalaman tertentu
Pb=[x for x in range(n)]  #tekanan gelembung
rho_g=[x for x in range(n)]  #densitas gelembung
W=[x for x in range(n)]   #gaya berat
F=[x for x in range(n)]   #gaya stokes
mb=[x for x in range(n)]   #massa benda
C=[x for x in range(n)]   #parameter pengganti gaya arcimedes dan gaya berat
D=[x for x in range(n)]   #parameter pengganti gaya stoke
velocity=[x for x in range(n)]   #kecepatan gelembung
#depth=[x for x in range(n)] #kedalaman air
rho_w=997.13                # densitas air pada suhu 25 deg C, kg/m3
g=9.8                       # gravitasi bumi, m/s2
Rg=8.314463                  #ketetapan gas, JK-1mol-1
s=71.99/1000                 #surface tension water-air, N/m
h=7/100                     # kedalaman bubble aerator, m
psat=100090.8089            #tekanan saturator pada 0.1 lpm
sd_psat=50                  #standar deviasi tekanan saturator
eta= 0.891/1000          #viscosity of water kg/m s
for i in range(n):
    r[i]=d[i]/2/1000        #jari jari gelembung, m
    Psat[i]=random.normalvariate(psat,sd_psat)  #tekanan atmosfer di saturator
    PL[i]=Psat[i]+rho_w*g*h 
    Pb[i]=PL[i]+2*s/r[i]
    rho_g[i]=Pb[i]*28/(Rg*(ti[i]+273.15)*1000)   
    V[i]=4/3*np.pi*r[i]**3 # volume gelembung, m3
    mb[i]=rho_g[i]*V[i]     #massa benda
    Fa[i]=rho_w*g*V[i]      #gaya archimedes ke atas
    W[i]=rho_g[i]*g*V[i]    #gaya berat
    F[i]=6*np.pi*eta*r[i]   #parameter gaya stokes
    C[i]=(Fa[i]-W[i])/mb[i] #parameter pengganti Fa dan W
    D[i]=F[i]/mb[i]         #parameter pengganti gaya stokes
    t = 0.5              #waktu, s
    velocity[i]=-C[i]/D[i]*np.exp(-D[i]*t)+C[i]/D[i]
    
print("densitas gelembung", np.mean(rho_g))    
print("diameter",np.max(d))   
print("gaya archimedes",np.max(Fa))
print("gaya berat",np.max(W))             
print("kecepatan",np.mean(velocity)) 
print("kecepatan minimal",np.min(velocity))
print("kecepatan maksimal",np.max(velocity))
plt.hist(velocity)
plt.xlabel('velocity, m/s')
plt.ylabel('frequency')
plt.show()
#print(velocity)    

#mencari koefisien convection oleh air
cp_w=4.18*1000 #specific heat of water, J/kg K
k_w=0.598 #thermal conductivity pada suhu 25 C
Pr=cp_w*eta/k_w #Prandtl number
Re=[x for x in range(n)] 
Nu=[x for x in range(n)] 
hconv=[x for x in range(n)]     #parameter konveksi 
As=[x for x in range(n)]     #luas area
j=[x for x in range(n)]     #parameter exponensial
suhu_akhir=[x for x in range(n)] #suhu gelembung lepas ke permukaan
for i in range(n):
    Re[i]=rho_w*velocity[i]*d[i]/eta   #Reynolds number
    Nu[i]=2+(0.4*Re[i]**0.5+0.06**(2/3))*Pr**0.4 #Nusselt number
    hconv[i]=Nu[i]*k_w/d[i]
    t=0.1
    As[i]=4*np.pi*r[i]**2
    j[i]=hconv[i]*As[i]/(rho_g[i]*V[i]*cp_a)
    suhu_akhir[i]=np.exp(-j[i]*t)*(ti[i]-tsat[i])+ti[i]
y =np.random.normal(np.mean(suhu_akhir),np.std(suhu_akhir),n)
plt.style.use('seaborn-deep')
plt.hist(suhu_akhir,label="simulasi")
plt.hist(y,label="gaussian")
plt.xlabel('suhu, degC')
plt.ylabel('frequency')
plt.title("suhu dew point")
plt.legend(loc="upper right")
plt.show()
print("suhu akhir",np.mean(suhu_akhir))   
print("simpangan suhu akhir",np.std(suhu_akhir))   

