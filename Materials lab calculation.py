#imports
from scipy import stats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import trapz
import os

os.getcwd()
 
###***meaning of variable names:
    #all variable starting/contain with 'a' meaning its relevent to aluminium
    #all variable starting/contain with 'l' meaning its relevent to low steel carbon
    #all variable starting/contain with 'h' meaning its relevent to high steel carbon
    
##Tensile testing
# 1) Calculate and plot a stress-strain graph //
# 2) Calculate and print the Young’s modulus (E) (linear regression) 
# 3) Calculate and print the 0.2% proof stress, ultimate tensile strength and failure strain
################################################################################################
 
#read data from csv file
Hdata=pd.read_csv("HC_tensile.csv")
Adata=pd.read_csv("Aluminium_tensile.csv")
Ldata=pd.read_csv("LC_tensile.csv")

#extract the force values of each material from the csv file
hf=Hdata['Force']
lf=Ldata['Force']
af=Adata['Force']

#use the  tensile strain of each material from the csv file
H= Hdata['Tensile strain (Strain 1)']
L= Ldata['Tensile strain (Strain 1)']
A= Adata['Tensile strain (Strain 1)']

###calculations of stresses

H_area=(0.00492**2)*3.14 # area of high carbon steel specimen
Hstress=(hf/H_area)/10**3# Stress of high carbon steel specimen(Pa)

L_area=(0.00508**2)*3.14# area of low carbon steel specimen
Lstress=(lf/L_area)/10**3# Stress of low carbon steel specimen(Pa)

A_area=(0.00475**2)*3.14# area of aluminium specimen
Astress=(af/A_area)/10**3# Stress of aluminium specimen(Pa)

# stress strain of three materials on same graph(overall view )

fig,ax = plt.subplots()# plot
hplt = ax.plot(H,Hstress,linewidth=1,color='orange',linestyle='-',label='High carbon steel')# line graph for high carbon steel 
lplt = ax.plot(L,Lstress,linewidth=1,color='pink',linestyle='-',label='Low carbon steel')# line graph for low carbon steel 
aplt = ax.plot(A,Astress,linewidth=1,color='blue',linestyle='-',label='Aluminium')# line graph for aluminium
ax.set_title('stress strain graph (overall view)',fontsize=14)# graph title 
ax.set_xlabel('strain(%)',fontsize=14)# x axis label 
ax.set_ylabel('Engineering stress(MPa)',fontsize=14)# y axis label 
ax.legend(fontsize = 7)#legend of the graph
plt.grid()#grid

#aluminium graph 
##########################
fig,ax3 = plt.subplots()
plt3 = ax3.plot(A,Astress,linewidth=1,color='blue',linestyle='-',label='Aluminium')
ax3.set_title('Aluminium stress strain graph',fontsize=14)
ax3.set_xlabel('strain(%)',fontsize=14)
ax3.set_ylabel('Engineering stress(MPa)',fontsize=14)


# #young's mudulus of aluminium
alrange=A.head(90)# only use data before plastic deformation
Alreg=stats.linregress(alrange,Astress[0:90]) # line regression info
ya=Alreg.slope
print("the young's modulus of aluminium is :" ,ya/10,"Gpa")#unit conversion from Mpa to Gpa divide by 10 as strain is in percentage   

#0.2% proof stress of aluminium
ax=np.linspace(0,1)
ay=ya*ax
aproof=plt.plot((ax+0.2)[0:23],ay[0:23],linestyle="--",color="green")

#equation of aluminium parallel equation
az=ya*A-0.2*ya

#difference between stress and y values of 0.2% stress parallel line
adiff=np.abs(Astress-az)

#the minimum value form the differnece 
aindex=np.argmin(adiff)
aPstress= Astress[aindex]
aintersect=ax3.plot(A[aindex],aPstress,marker="o",label="intersection")
print("the 0.2% stress of aluminium is ",aPstress,"Gpa")
#ultimate tensile strength
aultimate=np.max(Astress)
apointmax=ax3.plot(A[Astress.argmax()],aultimate,marker="o",label="max stress point")#[Astress.argmax()the element that Astress is maxmimum
print("The ultimate tensile strength of aluminium is",aultimate,"Mpa")

#Faliure strain
astrainmax=np.max(A)
apointfail=ax3.plot(astrainmax,Astress[267],marker="o",label="fail point")
print("The failure strain is of aluminium is",astrainmax,"%")
print(" ")
ax3.legend(fontsize = 8)#legend
plt.grid()#grid
plt.show()#show graph and points 


#high carbon steel graph 
#################################
fig,ax1 = plt.subplots()
plt1 = ax1.plot(H,Hstress,linewidth=1,color='orange',linestyle='-',label='High carbon steel')
ax1.set_title('high carbon steel stress strain graph',fontsize=14)
ax1.set_xlabel('strain(%)',fontsize=14)
ax1.set_ylabel('Engineering stress(MPa)',fontsize=14)

# #young's mudulus of high carbon steel
hrange=H.head(200)# only use data before plastic deformation
Hreg=stats.linregress(hrange,Hstress[0:200]) # line regression info
yh=Hreg.slope
print("the young's modulus of High carbon steel is :" ,yh/10,"Gpa")

#0.2% proof stress of high carbon steel
hx=np.linspace(0,0.4)
hy=yh*hx
hproof=plt.plot((hx+0.2)[0:40],ay[0:40],linestyle="--",color="green")

#equation of high carbon parallel equation
hz=yh*H-0.2*yh

#difference between stress and y values of 0.2% stress parallel line
hdiff=np.abs(Hstress-hz)

#the minimum value form the differnece 
hindex=np.argmin(hdiff)
hPstress= Hstress[hindex]
hintersect=ax1.plot(H[hindex],hPstress,marker="o",label="intersection")
print("the 0.2% stress of high carbon steel is ",hPstress,"Mpa")
#ultimate tensile strength
hultimate=np.max(Hstress)
hpointmax=ax1.plot(H[Hstress.argmax()], hultimate ,marker="o",label="max stress point")#[Hstress.argmax()the element that Astress is maxmimum
print("The ultimate tensile strength of high carbon steel is",hultimate,"Mpa")

#Faliure strain
hstrainmax=np.max(H)
hpointfail=ax1.plot(hstrainmax,Hstress[220],marker="o",label="fail point")
print("The failure strain is of high carbon steel is",hstrainmax,"%")
print(" ")
ax1.legend(fontsize = 8)#legend
plt.grid()#grid
plt.show()#show graph and points 


# low carbon steel graph
####################################
fig,ax2 = plt.subplots()
plt2 = ax2.plot(L,Lstress,linewidth=1,color='pink',linestyle='-',label='Low carbon steel')
ax2.set_title('low carbon steel stress strain graph',fontsize=14)
ax2.set_xlabel('strain(%)',fontsize=14)
ax2.set_ylabel('Engineering stress(MPa)',fontsize=14)


# #young's mudulus of low carbon steel
Lrange=L.head(100)# only use data before plastic deformation
Lreg=stats.linregress(Lrange,Lstress[0:100]) # line regression info
yl=Lreg.slope
print("the young's modulus of Low carbon steel is :" ,yl/10,"Gpa")

#0.2% proof stress of low carbon steel
lx=np.linspace(0,0.18)
ly=yl*lx
lproof=plt.plot((lx+0.2)[0:40],ly[0:40],linestyle="--",color="green")

#equation of low carbon steel parallel equation
lz=yl*L-0.2*yl

#difference between stress and y values of 0.2% stress parallel line
ldiff=np.abs(Lstress-lz)

#the minimum value form the differnece 
lindex=np.argmin(ldiff)
lPstress= Lstress[lindex]
lintersect=ax2.plot(L[lindex],lPstress,marker="o",label="intersection")
print("the 0.2% stress of low carbon steel is ",lPstress,"Mpa")
#ultimate tensile strength
lultimate=np.max(Lstress)
lpointmax=ax2.plot(L[Lstress.argmax()],lultimate,marker="o",label="max stress point")#[Hstress.argmax()the element that Astress is maxmimum
print("The ultimate tensile strength of low carbon steel is",lultimate,"Mpa")

#Faliure strain
lstrainmax=np.max(L)
lpointfail=ax2.plot(lstrainmax,Lstress[173],marker="o",label="fail point")
print("The failure strain is of aluminium is",lstrainmax,"%")
ax2.legend(fontsize = 8)#legend
plt.grid()#grid
plt.show()#show graph and points 


## torsion testing
# 1) Calculate and plot a stress-strain graph
# 2) Calculate and print the Shear modulus (G)
# 3) Toughness (trapezium rule of each material)

#####################################################################################################
# extract torsion result from csv 

hdata=pd.read_csv("Materials torsion H.csv")
ldata=pd.read_csv("Materials torsion L.csv")
adata=pd.read_csv("Materials torsion al.csv")
#Angle of rotation for each material(degree)
a=adata['AL Rotaion']
h=hdata['HC Rotaion']
l=ldata['LC Rotaion']

# convert degrees into radiance
arot=(a*3.1415926)/360
hrot=(h*3.1415926)/360
lrot=(l*3.1415926)/360
#Torque of each material(Nm)
at= adata['AL Torque']#.to_numpy()
ht= hdata['HC Torque']#.to_numpy()
lt= ldata['LC Torque']#.to_numpy()

##Shear stress and strain of each material
#initial radii(m)
iar=(5.88/2)/1000
ilr=(5.94/2)/1000
ihr=(5.84/2)/1000
#gauge length(m)
agl=76.6/1000
hgl=77.5/1000
lgl=77.5/1000

#polar moment of inertia (J)
aj=(3.14*iar**4)/2
hj=(3.14*ilr**4)/2
lj=(3.14*ihr**4)/2

#shear stress
asstress=(at*iar)/aj
hsstress=(ht*ihr)/hj
lsstress=(lt*ilr)/lj

# Shear strain 
asstrain=(arot*iar)/agl
hsstrain=(hrot*ihr)/hgl
lsstrain=(lrot*ilr)/lgl

#graph of shear stress/strain

fig,ax1 = plt.subplots()
plt1 = ax1.plot(hsstrain,hsstress,linewidth=1,color='orange',linestyle='-',label='High carbon steel')
plt2 = ax1.plot(lsstrain,lsstress,linewidth=1,color='pink',linestyle='-',label='Low carbon steel')
plt3 = ax1.plot(asstrain,asstress,linewidth=1,color='blue',linestyle='-',label='Aluminium')
ax1.set_title('shear stress strain graph',fontsize=14)
ax1.set_xlabel('Shear strain',fontsize=14)
ax1.set_ylabel('Shear stress(Pa)',fontsize=14)
ax1.legend(fontsize = 8)
plt.grid()#grid

# ## Shear modulus(G) gradient of graph with value up to where the line is straight

Alreg=stats.linregress(asstrain[0:4],asstress[0:4]) # range select 
Hlreg=stats.linregress(hsstrain[0:6],hsstress[0:6])
Ilreg=stats.linregress(lsstrain[0:6],lsstress[0:6])
#gradient 
aG=Alreg.slope
hG=Hlreg.slope
lG=Ilreg.slope

#shear modulus in Mpa
print(' ')
print("the shear modulus of aluminium is ",aG/10**6,"(MPa)\nthe shear modulus of low carbon steel is ",lG/10**6,"(MPa)\nthe shear modulus of high carbon steel is ",hG/10**6,'(MPa)')
print(' ')
# #Toughness[Pa] (trapezium rule of each material)
atoughness=trapz(asstress,dx=1)
htoughness=trapz(hsstress,dx=1)
ltoughness=trapz(lsstress,dx=1)

print("the toughness of aluminium is ",atoughness,"(Pa)\nthe toughness of low carbon steel is ",ltoughness,"(Pa)\nthe toughness of high carbon steel is ",htoughness,'(Pa)')
print(' ')

## impact testing
###############################################################################

#1) Calculate and print the energy absorbed for each material at high and low temperature
#Energy lost to fracture(J)
temp=np.array([20,-48,20,-48,20,-47])
energyf=np.array([0.25,0.65,0.20,0.61,0.21,0.20])#aluminium,High carbon steel,Low carbon steel
im= energyf-0.04 # energy lost through friction and resistance=0.04(J)
print("for aluminium:")
print("Energy absorbed when the temperature is " ,temp[0], " centigrade is ",im[0], " J")
print("Energy absorbed when the temperature is ",temp[1]," centigrade is ", energyf[1], " J")
print(" ")

print("for high carbon:")
print("Energy absorbed when the temperature is " ,temp[2], " centigrade is ",im[2], " J")
print("Energy absorbed when the temperature is ",temp[3]," centigrade is ", im[3], " J")
print(" ")

print("for low carbon:")
print("Energy absorbed when the temperature is " ,temp[4], " centigrade is ",im[4], " J")
print("Energy absorbed when the temperature is ",temp[5]," centigrade is ", im[5], " J")
