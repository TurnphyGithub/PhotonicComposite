# -*- coding: utf-8 -*-
"""
@author: Turnphy
"""

import numpy as np
from matplotlib import pyplot as plt
from tkinter import *
from TransferMatrixEngine import (coh_tmm, position_resolved, find_in_structure_with_inf)
from numpy import pi, linspace, inf, array
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

try:
    import colorpy.illuminants
    import colorpy.colormodels
    from . import color
    colors_were_imported = True
except ImportError:
    # without colorpy, you can't run sample5(), but everything else is fine.
    colors_were_imported = False


# "5 * degree" is 5 degrees expressed in radians
# "1.2 / degree" is 1.2 radians expressed in degrees
degree = pi/180





root=Tk()
root.title("Calculate reflectance")


# Define the buttons
def saveStruc():
    global index1,thickness1,index2,thickness2,repeats
    index1,thickness1,repeats=float(n1.get()),float(d1.get()), float(period1.get())
    index2,thickness2=float(n2.get()),float(d2.get())
    print("n2 = " + str(index2))
    return

def AngularRef(wavelength0, thetal1, thetal2):
    """
    This code use TMM
    """
    global index1,thickness1,index2,thickness2,repeats
    d_list = [inf, 300, thickness1, thickness2, inf] #in nm
    n_list = [2.2, 1, index1+0.01j, index2+0.001j, 3.46]    
    
    
    # list of theta to plot in deg
    thetas=linspace(thetal1, thetal2, num=1000)
    # list of wavenumbers to plot in nm^-1
    # ks = linspace(0.0001, .01, num=400)
    # initialize lists of y-values to plot
    Rangle = []
    for theta in thetas:
		# For normal incidence, s and p polarizations are identical.
		# I arbitrarily decided to use 's'.        
        Rangle.append(coh_tmm('s', n_list, d_list, theta*pi/180, wavelength0)['R'])
        # R45.append(unpolarized_RT(n_list, d_list, 45*degree, 1/k)['R'])
    # kcm = ks * 1e7 #ks in cm^-1 rather than nm^-1
    plt.figure()
    plt.plot(thetas, Rangle, 'blue')
    plt.xlabel('Incident Angle (degree)')
    plt.ylabel('Reflectance')
    # plt.title('Reflection of unpolarized light at 0$^\circ$ incidence (blue), '
    #           '45$^\circ$ (purple)')
    plt.show()
    
def Spectral(angle0, wavelength1, wavelength2):
    """
    This code use TMM
    """
    global index1,thickness1,index2,thickness2,repeats
    d_list = [inf, 300, thickness1, thickness2, inf] #in nm
    n_list = [2.2, 1, index1+0.01j, index2+0.001j, 3.46]    
    
    
    # list of theta to plot in deg
    wavelengths=linspace(wavelength1, wavelength2, num=1000)
    # thetas=linspace(thetal1, thetal2, num=1000)
    # list of wavenumbers to plot in nm^-1
    # ks = linspace(0.0001, .01, num=400)
    # initialize lists of y-values to plot
    Rw = []
    for wavelength in wavelengths:
		# For normal incidence, s and p polarizations are identical.
		# I arbitrarily decided to use 's'.
        Rw.append(coh_tmm('s', n_list, d_list, angle0, wavelength)['R'])
        # R45.append(unpolarized_RT(n_list, d_list, 45*degree, 1/k)['R'])
    # kcm = ks * 1e7 #ks in cm^-1 rather than nm^-1
    plt.figure()
    plt.plot(wavelengths, Rw, 'red')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Reflectance')
    # plt.title('Reflection of unpolarized light at 0$^\circ$ incidence (blue), '
    #           '45$^\circ$ (purple)')
    plt.show()
    
def CalRef():
    global clicked
    if clicked.get()=='Angle (deg)':
        print(1)
        thetal1, thetal2=float(lim1.get()),float(lim2.get())
        wavelength0=float(lim0.get())
        print('CalRef.theta1')
        print(thetal1)
        print('CalRef.theta1')
        AngularRef(wavelength0, thetal1, thetal2)
    else:
        wavelength1, wavelength2=float(lim1.get()),float(lim2.get())
        angle0=float(lim0.get())
        Spectral(angle0, wavelength1, wavelength2)
        print('Wait...')
    return

# Calculate the electric field
def E_y(angle0, wavelength0, depth=600):
    """
    Here is an example where we plot absorption and Poynting vector
    as a function of depth.
    """

    global index1,thickness1,index2,thickness2,repeats
    d_list = [inf, 300, thickness1, thickness2, inf] #in nm
    n_list = [2.2, 1, index1+0.01j, index2+0.001j, 3.46] 
    # d_list = [inf, 300, 310,1330, inf] #in nm
    # n_list = [2.2, 1, 2.08+0.01j, 1.41+0.001j, 3.46]
    depth=sum(d_list[1:-2])+depth
    th_0 = pi*angle0/180
    lam_vac =wavelength0
    pol = 's'
    coh_tmm_data = coh_tmm(pol, n_list, d_list, th_0, lam_vac)

    ds = linspace(-50, depth, num=1000) #position in structure
    poyn = []
    Ey = []
    for d in ds:
        layer, d_in_layer = find_in_structure_with_inf(d_list, d)
        data = position_resolved(layer, d_in_layer, coh_tmm_data)
        poyn.append(data['poyn'])
        Ey.append(data['Ey'])
    # convert data to numpy arrays for easy scaling in the plot
    poyn = array(poyn)
    Ey = array(Ey)
    plt.figure(figsize=(4, 2.5), dpi=200, facecolor='w', edgecolor='k')
    plt.plot(ds, abs(Ey)**2, 'purple')
    zintface=d_list[1]
    plt.plot([zintface,zintface],[0, max(abs(Ey)**2)], 'k--')
    plt.xlabel('depth (nm)',fontsize=14)
    plt.ylabel('E-squared (a.u.)',fontsize=14)
    plt.xticks(fontsize= 12) 
    plt.yticks(fontsize= 12)
    plt.tight_layout()
    # plt.title('E-squared')
    plt.show()
    
def ShowField():
    global theta2f,lambda2f, z2f
    theta2f,lambda2f, z2f=float(angle2f.get()),float(wavelength2f.get()), float(depth2f.get())
    E_y(theta2f,lambda2f, z2f)
    
    return
def CalChoice(modes):
    global lim0, lim1, lim2
    global l0,l1,l2,l3
    if modes=='Angle (deg)':
        l0.destroy()
        l1.destroy()
        l2.destroy()
        l3.destroy()        
        
        l0=Label(frame01, text="Setting the working wavelength (nm) and angle limits (deg) ")
        l0.grid(row=1, column=0, columnspan=2) 
        l1=Label(frame01, text="wavelength ")
        l1.grid(row=2, column=0)
        lim0=Entry(frame01, width = 10, borderwidth=5)        
        lim0.grid(row=2, column=1, columnspan=1, padx=5, pady=5)            
        l2=Label(frame01, text="theta1 ")
        l2.grid(row=3, column=0)
        lim1=Entry(frame01, width = 10, borderwidth=5)
        lim1.grid(row=3, column=1, columnspan=1, padx=5, pady=5)
        l3=Label(frame01, text="theta2 ")
        l3.grid(row=4, column=0)
        lim2=Entry(frame01, width = 10, borderwidth=5)
        lim2.grid(row=4, column=1, columnspan=1, padx=5, pady=5)
        
        lim0.insert(0, 1550)
        lim1.insert(0, 30)
        lim2.insert(0, 70)
        
    if modes=="Wavelength (nm)":
        l0=Label(frame01, text="Setting the working angle (deg) and wavelength limits (nm) ")
        l0.grid(row=1, column=0, columnspan=2)
        l0.destroy()
        l1.destroy()
        l2.destroy()
        l3.destroy()
        l1=Label(frame01, text="angle ")
        l1.grid(row=2, column=0)
        lim0=Entry(frame01, width = 10, borderwidth=5)
        lim0.grid(row=2, column=1, columnspan=1, padx=5, pady=5)   
        l2=Label(frame01, text="wavelength1 ")
        l2.grid(row=3, column=0)
        lim1=Entry(frame01, width = 10, borderwidth=5)
        lim1.grid(row=3, column=1, columnspan=1, padx=5, pady=5)
        l3=Label(frame01, text="wavelength2 ")
        l3.grid(row=4, column=0)
        lim2=Entry(frame01, width = 10, borderwidth=5)
        lim2.grid(row=4, column=1, columnspan=1, padx=5, pady=5)
        
        lim0.insert(0, 0)        
        lim1.insert(0, 1250)
        lim2.insert(0, 1750)        
        
        

def CleTable():
    n1.delete(0,END)
    n2.delete(0,END)
    d1.delete(0,END)
    d2.delete(0,END)
    period1.delete(0,END)
    return

DefaultValues=[2.08,310,1.41,1330,1]
def SetDefault():
    n1.delete(0,END)
    n2.delete(0,END)
    d1.delete(0,END)
    d2.delete(0,END)
    period1.delete(0,END)
    global index1,thickness1,index2,thickness2,repeats
    index1,thickness1,index2,thickness2,repeats=DefaultValues
    n1.insert(0,index1)
    d1.insert(0,thickness1)
    n2.insert(0,index2)
    d2.insert(0,thickness2)
    period1.insert(0,repeats)
    
    return


# Layout
frame00=LabelFrame(root,text="Define the structure here", padx=50, pady=50)
frame00.grid(row=0, column=0, columnspan=2, padx=100,pady=10)
frame01=LabelFrame(root,text="Calculate the reflectance", padx=10, pady=20)
frame01.grid(row=1, column=0, padx=10,pady=10)
frame02=LabelFrame(root,text="Show field profile", padx=10, pady=20)
frame02.grid(row=1, column=1, padx=10,pady=10)
# Frame for the structures Frame00 blanks
Position=[0,0]
i,j= Position   
DefineStructureLabel=Label(frame00, text="From Top to Bottom, Units(RIU, nm, #)", anchor=CENTER)
DefineStructureLabel.grid(row=0, column=0, columnspan=10, padx=10, pady=5)

Label(frame00, text="n1 ").grid(row=1, column=0)
n1=Entry(frame00, width = 10, borderwidth=5)
n1.grid(row=1, column=1, columnspan=1, padx=5, pady=5)
Label(frame00, text="d1 ").grid(row=1, column=2)
d1=Entry(frame00, width = 10, borderwidth=5)
d1.grid(row=1, column=3, columnspan=1, padx=5, pady=5)
Label(frame00, text="n2 ").grid(row=1, column=4)
n2=Entry(frame00, width = 10, borderwidth=5)
n2.grid(row=1, column=5, columnspan=1, padx=5, pady=5)
Label(frame00, text="d2 ").grid(row=1, column=6)
d2=Entry(frame00, width = 10, borderwidth=5)
d2.grid(row=1, column=7, columnspan=1, padx=5, pady=5)
Label(frame00, text="repeats ").grid(row=1, column=8)
period1=Entry(frame00, width = 10, borderwidth=5)
period1.grid(row=1, column=9, columnspan=1, padx=5, pady=5)

iloc=1
Label(frame00, text="n1 ").grid(row=1+iloc, column=0)
n11=Entry(frame00, width = 10, borderwidth=5)
n11.grid(row=1+iloc, column=1, columnspan=1, padx=5, pady=5)
Label(frame00, text="d1 ").grid(row=1+iloc, column=2)
d11=Entry(frame00, width = 10, borderwidth=5)
d11.grid(row=1+iloc, column=3, columnspan=1, padx=5, pady=5)
Label(frame00, text="n2 ").grid(row=1+iloc, column=4)
n21=Entry(frame00, width = 10, borderwidth=5)
n21.grid(row=1+iloc, column=5, columnspan=1, padx=5, pady=5)
Label(frame00, text="d2 ").grid(row=1+iloc, column=6)
d21=Entry(frame00, width = 10, borderwidth=5)
d21.grid(row=1+iloc, column=7, columnspan=1, padx=5, pady=5)
Label(frame00, text="repeats ").grid(row=1+iloc, column=8)
period11=Entry(frame00, width = 10, borderwidth=5)
period11.grid(row=1+iloc, column=9, columnspan=1, padx=5, pady=5)

iloc=2
Label(frame00, text="n1 ").grid(row=1+iloc, column=0)
n12=Entry(frame00, width = 10, borderwidth=5)
n12.grid(row=1+iloc, column=1, columnspan=1, padx=5, pady=5)
Label(frame00, text="d1 ").grid(row=1+iloc, column=2)
d12=Entry(frame00, width = 10, borderwidth=5)
d12.grid(row=1+iloc, column=3, columnspan=1, padx=5, pady=5)
Label(frame00, text="n2 ").grid(row=1+iloc, column=4)
n22=Entry(frame00, width = 10, borderwidth=5)
n22.grid(row=1+iloc, column=5, columnspan=1, padx=5, pady=5)
Label(frame00, text="d2 ").grid(row=1+iloc, column=6)
d22=Entry(frame00, width = 10, borderwidth=5)
d22.grid(row=1+iloc, column=7, columnspan=1, padx=5, pady=5)
Label(frame00, text="repeats ").grid(row=1+iloc, column=8)
period12=Entry(frame00, width = 10, borderwidth=5)
period12.grid(row=1+iloc, column=9, columnspan=1, padx=5, pady=5)

# Default photonic structure
n1.insert(0,2.08)
d1.insert(0,310)
n2.insert(0,1.41)
d2.insert(0,1330)
period1.insert(0,1)
index1,thickness1,index2,thickness2,repeats=2.08,310,1.41,1330,1
# Set the buttons in Frame00
SetDefaultButton =Button(frame00, text= "Default", command=SetDefault)
SetDefaultButton.grid(row=10, column=0, columnspan=5)


CleTableButton =Button(frame00, text= "Clear", command=CleTable)
CleTableButton.grid(row=10, column=3, columnspan=5)


OKButton =Button(frame00, text= "OK", command=saveStruc)
OKButton.grid(row=10, column=7, columnspan=5)

# Options in Frame01
Label(frame01, text="Select mode").grid(row=0, column=0)
options=['Angle (deg)', 'Wavelength (nm)']
clicked=StringVar()
clicked.set(options[1]) # default option

## defaults

l0=Label(frame01, text="Setting the working angle (deg) and wavelength limits (nm) ")
l0.grid(row=1, column=0, columnspan=2)
l1=Label(frame01, text="angle ")
l1.grid(row=2, column=0)
lim0=Entry(frame01, width = 10, borderwidth=5)
lim0.grid(row=2, column=1, columnspan=1, padx=5, pady=5)   
l2=Label(frame01, text="wavelength1 ")
l2.grid(row=3, column=0)
lim1=Entry(frame01, width = 10, borderwidth=5)
lim1.grid(row=3, column=1, columnspan=1, padx=5, pady=5)
l3=Label(frame01, text="wavelength2 ")
l3.grid(row=4, column=0)
lim2=Entry(frame01, width = 10, borderwidth=5)
lim2.grid(row=4, column=1, columnspan=1, padx=5, pady=5)
lim0.insert(0, 0)        
lim1.insert(0, 1250)
lim2.insert(0, 1750)  

# Dropdown 
drop = OptionMenu(frame01, clicked,  *options, command=CalChoice)
drop.grid(row=0, column=1)

# Set the buttons in Frame01
CalRefButton =Button(frame01, text= "Calculate", command=CalRef)
CalRefButton.grid(row=10, column=0,columnspan=2)


# Frame for the structures Frame02 blanks
Label(frame02, text="Setting the working conditions:\
      ").grid(row=1, column=0, columnspan=2)
Label(frame02, text="angle ").grid(row=2, column=0)
angle2f=Entry(frame02, width = 10, borderwidth=5)
angle2f.grid(row=2, column=1, columnspan=1, padx=5, pady=5)   
Label(frame02, text="wavelength ").grid(row=3, column=0)
wavelength2f=Entry(frame02, width = 10, borderwidth=5)
wavelength2f.grid(row=3, column=1, columnspan=1, padx=5, pady=5)
Label(frame02, text="+ Depth(nm) ").grid(row=4, column=0)
depth2f=Entry(frame02, width = 10, borderwidth=5)
depth2f.grid(row=4, column=1, columnspan=1, padx=5, pady=5)
angle2f.insert(0, 50.42)        
wavelength2f.insert(0, 1550)
depth2f.insert(0, 1000)  
theta2f,lambda2f, z2f=50.42,1550,1000


# Set the buttons in Frame02

ShowFieldButton=Button(frame02, text= "Show Field Profile",command=ShowField)
ShowFieldButton.grid(row=10, column=0,columnspan=2)

root.mainloop()