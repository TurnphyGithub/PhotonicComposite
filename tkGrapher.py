# -*- coding: utf-8 -*-
"""
@author: Turnphy
"""
import numpy as np
from matplotlib import pyplot as plt
from tkinter import *
from tkinter import filedialog
from PIL import ImageTk
import PIL
from scipy.fftpack import fft
from scipy import interpolate
from detect_peaks_py import detect_peaks
# from tkinter import RadioButton

def QuickLook():
    filenameQL=filedialog.askopenfilename(title='select the data',\
filetypes=(("txt files", '*.txt'),("csv files", '*.csv')))
    dat=np.loadtxt(filenameQL)
    xx,yy=dat.T, dat.T
    x=xx[0]
    y=yy[1]
    ylimt=5*np.median(y)
    plt.figure()
    plt.plot(x,y)
    plt.ylim((0,ylimt))
    plt.xlabel('xData')
    plt.ylabel('yData')    

def LoadDataFile():
    global DP_
    DP_.filename=filedialog.askopenfilename(title='select the data',\
filetypes=(("txt files", '*.txt'),("csv files", '*.csv')))
    print(DP_.filename)
    return

def nextpow2(n):
    """ Function to calculate next power of 2 Value,
    it returns 2^n"""
    x=1
    while x<n:
        x=x*2
    return x

def PlotGraph():
    try:
        da=np.loadtxt(DP_.filename)
        if DP_modes.get()==1:
            xx,yy=da.T, da.T
            x=xx[0][985:2800]
            y=yy[1][985:2800]
            plt.figure()
            plt.plot(x,y)
        if DP_modes.get()==2:
            wavelength=(da.T)[0][985:3000]
            reflection=(da.T)[1][985:3000]
            plt.figure()
            plt.subplot(311)
            plt.plot(wavelength, reflection)
            k=np.flipud(1000/wavelength)
            reflection=np.flipud(reflection)
        
        
            plt.subplot(312)
            plt.plot(k, reflection)     
            
            # sample spacing
            T = 0.0002
            tck=interpolate.splrep(k, reflection, s=0)
            kk= np.arange(1.3, 1.8, T)
            RK=interpolate.splev(kk, tck, der=0)
            
            RK = RK - np.mean(RK)
            L = len(RK)         
            RK=RK*np.hanning(L)
            NFFT=nextpow2(L)*128
            yf = fft(RK,NFFT)/L
            xf = np.linspace(0.0, 1000.0/(2.0*T), NFFT//2)
            fR =2.0 * np.abs(yf[0:NFFT//2])
            ffb=xf[1:4000]
            fFRb=fR[1:4000]
            tcf=interpolate.splrep(ffb, fFRb, s=0)
            ff=np.linspace(ffb[0],ffb[-1], num=400000)
            fFR=interpolate.splev(ff, tcf, der=0)
            peaks= detect_peaks(fFR)
            maxind=np.argmax(fFR[peaks], axis = None, out = None) 
            print(ff[peaks[maxind]])
            plt.subplot(313)
            plt.plot(ff, fFR)
            plt.plot(ff[peaks], fFR[peaks], "x")
            for i in range(3):
                print(str(ff[peaks[i]])+','+str(fFR[peaks[i]]))
                print('\n')
            #plt.plot(k,reflection)
            plt.grid()
            plt.show()
            
            
    except:
        messagebox.showinfo('Operation error','Please load the file first')
        # Label(DP_, text='Please load the file first').grid(row=10, column=0)
    return


    

def AddRadioButton():
    """This is a helper function of DP_create;
    It's to create the selections
    """
    global DP_modes, DP_
    DP_modes= IntVar()
    DP_modes.set(1)
    Radiobutton(DP_,text='Reflectance Spectrum', variable=DP_modes, value=1).grid(row=1, column=0)
    Radiobutton(DP_,text='FFT Spectrum', variable=DP_modes, value=2).grid(row=2, column=0)
    return


    
def DP_create():
    global DP_    
    DP_= Toplevel()
    DP_.title('Data Processor')
    LoadData=Button(DP_, text='Load the data', padx=50, pady=10,command= LoadDataFile)
    LoadData.grid(row=0, column=0, padx=50, pady=50)
    AddRadioButton()
    Button(DP_,text='Plot the data', padx=50, pady=10,command= PlotGraph).grid(row=3, column=0)
    return

root=Tk()
root.title("Optical analyzer v.0.1")
root.iconbitmap("light.ico")
def donothing():
    return
menu=Menu(root)
menubar=Menu(root)
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label="Quick view", command=QuickLook)
filemenu.add_command(label="Data Processing ", command=DP_create)
filemenu.add_command(label="Multilayer Simulation", command=DP_create)
filemenu.add_separator()
filemenu.add_command(label="Exit", command=root.quit)

menubar.add_cascade(label="File", menu=filemenu)

helpmenu = Menu(menubar, tearoff=0)
helpmenu.add_command(label="Help Index", command=donothing)
helpmenu.add_command(label="About...", command=donothing)
menubar.add_cascade(label="Help", menu=helpmenu)

root.config(menu=menubar)


MainIcon_img= ImageTk.PhotoImage(PIL.Image.open("Images/SpetralAnalysis0.png"))
MainIcon=Label(image=MainIcon_img)
MainIcon.grid(row=0,column=0, columnspan=2)

QLIcon_img= ImageTk.PhotoImage(PIL.Image.open("Images/magnifyingGlass.png"))
QL=Button(root,text="        Quick View       ", padx=40, pady=10, image=QLIcon_img,\
compound = LEFT, command=QuickLook, font='Arial 10 bold')
QL.grid(row=1, column=0, columnspan=2, padx=50,pady=10)

DPIcon_img =  ImageTk.PhotoImage(PIL.Image.open("Images/DataProcessor.png"))
DP=Button(root,text="    Data Processing   ", padx=40, pady=10, image=DPIcon_img,\
compound = LEFT, command=DP_create,font='Arial 10 bold')
DP.grid(row=2, column=0, columnspan=2, padx=50,pady=10)

MSIcon_img =  ImageTk.PhotoImage(PIL.Image.open("Images/Simulation.png"))
MS=Button(root,text="Multilayer Simulation", padx=40, pady=10, image=MSIcon_img,\
compound = LEFT, command=DP_create, font='Arial 10 bold')
MS.grid(row=3, column=0, columnspan=2, padx=50,pady=10)

root.mainloop()


