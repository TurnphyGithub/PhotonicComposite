# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 17:19:29 2019

@author: Turnphy
"""

from tkinter import filedialog
from PIL import ImageTk
import numpy as np
from matplotlib import pyplot as plt
from tkinter import *
from TransferMatrixEngine import (coh_tmm, position_resolved, find_in_structure_with_inf)
from numpy import pi, linspace, inf, array
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def helpManual():
    global HP_
    HP_= Toplevel()
    text1='This software is used for optical signal analysis:'
    Label(HP_, text=text1, font='Arial 10 bold').grid(row=0, column=0)
    text2='Quick view can give you a quick snapshot of the data.'
    Label(HP_, text=text2, font='Arial 10 normal').grid(row=1, column=0)
    text3='Data processing allows direct plot and FFT analysis.'
    Label(HP_, text=text3, font='Arial 10 normal').grid(row=2, column=0)
    text4='Multilayer simulation uses Transfer matrix, both \n\
reflectance and field profile can be shown'
    Label(HP_, text=text4, font='Arial 10 normal').grid(row=3, column=0)    
    return

def aboutMe():
    global AM_
    AM_= Toplevel()
    text1='Tengfei Cao, graduate student at Weiss Lab, Vanderbilt University.'
    Label(AM_, text=text1, font='Arial 10 bold').grid(row=0, column=0)   
    text2="Note: The transfer matrix internal engine was adapted from\n\
Steven Byrnes' transfer matrix code, for more please visit\n\
https://sjbyrnes.com/ , https://github.com/sbyrnes321/tmm/."
    Label(AM_, text=text2, font='Arial 10 italic').grid(row=1, column=0)  
    
    return