# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 16:37:35 2021
updated 18.04.2024

integration of the spectral line absorbance and determination of half width (FWHM)
-------------------------------------------------------------------------------------------------------------
before start:
-------------
set variable noise, which describes signal level at saturation and should be excluded from analysis of the line
default:
noise = 0.046 - determine the level where line is saturated

Interaction:
-------------
Select spectral line with SPACE and mouse cursor and read the data at the screen output
 
@author: Pipa
"""
import os
import sys
from numpy import zeros # used to load the asci file
from matplotlib import pyplot as plt
from matplotlib.widgets import Cursor
import numpy as np
import math
from math import sqrt
from PyQt6 import QtWidgets
from scipy.optimize import curve_fit

noise = 0.07 # for determination of the saturation level, 0.01 = 1%
sel_x = []
data_sig = []
data_sig_mm = []

def get_data_folder():
    path = os.path.dirname(os.path.abspath(__file__))
    # print(path)
    os.chdir(path)
    os.chdir("..\\")
    # print(os.getcwd())
    with open("set_dir.txt", "rt") as f:
        rpath = f.read()
        print(f"{rpath=}")
        os.chdir(rpath)


def int_Folder():
    get_data_folder()
    global noise
    #  data = [finfo, data_x, data_y, data_y1, data_y2]
    #  finfo = [filename, Header]
    #  file_matrix = [data[i]] 
    #### finfo = file_matrix[i][0]
    #### data_x = file_matrix[i][1]
    tmp = open_folder()
    file_matrix = tmp[0]
    
    fig_abs = plt.figure(num='abs Figure', figsize=(15, 7), 
                              facecolor='#E0EDCA' ,clear=True,)
    ax_abs = fig_abs.add_subplot(111, facecolor='#DAE1CB')
    
    L_names = []
    for i in range(len(file_matrix)):
        for ii in range(len(file_matrix[i][2])):
            # mark seturated points in signal data (data_y)
            if file_matrix[i][2][ii]<=noise :
                file_matrix[i][2][ii]= noise
            elif file_matrix[i][2][ii]> noise :
                file_matrix[i][2][ii]= file_matrix[i][2][ii]
            else : print('else condition int_Folder file_matrix[i][2] ')
        file_matrix[i][2]= [math.log(1/num) for num in file_matrix[i][2]]
        # mark seturated points in data_y1
        for ii in range(len(file_matrix[i][3])):
            if file_matrix[i][3][ii]<=noise :
                file_matrix[i][3][ii]= noise
            elif file_matrix[i][3][ii]> noise :
                file_matrix[i][3][ii]= file_matrix[i][3][ii]
            else : print('else condition int_Folder file_matrix[i][3] ')
        file_matrix[i][3]= [math.log(1/num) for num in file_matrix[i][3]]
#        print('type(file_matrix[i][3]) =', type(file_matrix[i][3]))
        # mark seturated points in data_y2
        for ii in range(len(file_matrix[i][4])):
            if file_matrix[i][4][ii]<=noise :
                file_matrix[i][4][ii]= noise
            elif file_matrix[i][4][ii]>noise :
                file_matrix[i][4][ii]= file_matrix[i][4][ii]
            else : print('else condition int_Folder file_matrix[i][4] ')
            
        file_matrix[i][4]= [math.log(1/num) for num in file_matrix[i][4]]
        # plot min and max 
        ax_abs.plot(file_matrix[i][1], file_matrix[i][3], color='xkcd:grey', marker='.', linewidth=0)
        ax_abs.plot(file_matrix[i][1], file_matrix[i][4],  color='xkcd:grey', marker='.', linewidth=0)
        # plot main signal 
        ax_abs.plot(file_matrix[i][1], file_matrix[i][2], '.')
#        print('type(l) =', type(l))

        line_name = 'L'+ str(i)
        nul = zeros(len(file_matrix[i][1]))
        line_name, = ax_abs.plot(file_matrix[i][1], nul, '-')
        L_names.append(line_name)

        
    cursor = Cursor(ax_abs, useblit=True, color='red', linewidth=1)
    cid_key = fig_abs.canvas.mpl_connect('key_press_event', lambda event: 
                                         on_key_2(event, fig_abs, ax_abs, 
                                         file_matrix, L_names)
                                        )
    plt.show()
    return (cid_key, cursor)

def on_key_2(event, fig_abs, ax_abs, file_matrix, L_names):
    print('you pressed', event.key, event.xdata, event.ydata)
    global noise
    
    if event.key == ' ' :
        sel_x.append(event.xdata)
        if len(sel_x)== 1:
            print('space is pressed first time')
            
        elif len(sel_x)== 2:
            print('space is pressed second time')
            # search for index of sel_x in x = file_matrix[j][1]
            for j in range(len(file_matrix)):
                i = 0
                while file_matrix[j][1][i] < min(sel_x[0], sel_x[1]):
                    i=i+1
                start_index = i
                i = 0
                while file_matrix[j][1][i] < max(sel_x[0], sel_x[1]):
                    i = i+1
                stop_index = i
#            print('start/stop_index =', start_index, '/', stop_index)
#                sel_index.append([start_index, stop_index])

                ##############################################################
                ### start with SIG data ### ### start with SIG data ### 
                x = file_matrix[j][1][start_index:stop_index]
                y = file_matrix[j][2][start_index:stop_index]
                i_peak = y.index(max(y))
#                print('x_0 [', j,'] = ', x[i_peak])
#                print('y_0 = ', y[i_peak])
                
                fwhm_sig = get_fwhm(x, y)
                sig_x1 = fwhm_sig[0][0]
                sig_y1 = fwhm_sig[0][1]
                sig_x2 = fwhm_sig[1][0]
                sig_y2 = fwhm_sig[1][1]
             
                ax_abs.plot(x[i_peak], y[i_peak], 'ro' )
                ax_abs.plot(sig_x1, sig_y1, 'ro' )  
                ax_abs.plot(sig_x2, sig_y2, 'ro' )
                
#                plt.draw()
                
                sig_A0 = y[i_peak]
                sig_fwhm0 = sig_x2 - sig_x1
                sig_w0 = sig_fwhm0/2.35482 # transfer to the factor in Gaussian curve
                sig_area0 = integral_trapez(x, y)
#                print('A0 = ', sig_A0)
#                print('fwhm0 = ', sig_fwhm0)
#                print('area0 = ', sig_area0)
                
                ######################################### 
                ############# fit Gaus ##################
                y_tmp = []
                x_tmp = []
                for i in range(len(y)):
                    if y[i] >= math.log(1/noise): # remoove seturated points
                        pass
                    elif y[i] < math.log(1/noise):
                        y_tmp.append(y[i])
                        x_tmp.append(x[i])
                
                popt, pcov = curve_fit(my_Gauss, x_tmp, y_tmp, p0 = [sig_A0, x[i_peak], sig_w0])
#               print('popt = ', *popt)                  #   
                gauss_y = my_Gauss(x, *popt)             #
                fwhm_G = popt[2]*2.35482                 #
                area_G = sqrt(2*math.pi)*popt[0]*popt[2] #
                peak_G = popt[1]
#               print('fwhm_G = ', fwhm_G)               #
#               print('area_G = ', area_G)               #
                #plot selected eta data                  #
                
                L_names[j].set_xdata(x)                  #
                L_names[j].set_ydata(gauss_y)            #
                
#                plt.draw()                               # 
                #####################################    #
                ### collect signal data #####
                area = sig_area0
                area_err = sqrt( (sig_area0 - area_G)**2 )
                fwhm = sig_fwhm0
                fwhm_err = sqrt( (sig_fwhm0 - fwhm_G)**2 )
                data_sig.append([[area, area_err, area_G, peak_G], [fwhm, fwhm_err, fwhm_G]]) 
                
                #############################################
                ### MIN data ### ### MIN data ###############
                y_tmp = []
                x_tmp = [] 
                y_min = file_matrix[j][3][start_index:stop_index]
                for i in range(len(y_min)):
                    if y_min[i] >= math.log(1/noise): 
                       pass
                    elif y_min[i] < math.log(1/noise):
                        y_tmp.append(y_min[i])
                        x_tmp.append(x[i])
                    else: print('else condition for selection of min data')
                
                popt, pcov = curve_fit(my_Gauss, x_tmp, y_tmp, p0 = [sig_A0, x[i_peak], sig_w0])
#               print('popt = ', *popt)                     #   
                gauss_y = my_Gauss(x, *popt)                #
                fwhm_Gmin = popt[2]*2.35482                 #
                area_Gmin = sqrt(2*math.pi)*popt[0]*popt[2] #
#               print('fwhm_Gmin = ', fwhm_Gmin)            #
#               print('area_Gmin = ', area_Gmin)            #
                #plot selected eta data                     #
                plt.plot(x, gauss_y, 'k')
                
#                plt.draw()                                  #
                #############################################
                ### MAX data ### ### MAX data ###############
                y_tmp = []
                x_tmp = []    
                y_max = file_matrix[j][4][start_index:stop_index]
                for i in range(len(y_max)):
                    if y_max[i] >= math.log(1/noise): 
                       pass
                    elif y_max[i] < math.log(1/noise):
                        y_tmp.append(y_max[i])
                        x_tmp.append(x[i])
                    else: print('else condition for selection of max data')
                
                popt, pcov = curve_fit(my_Gauss, x_tmp, y_tmp, p0 = [sig_A0, x[i_peak], sig_w0])
#               print('popt = ', *popt)                     #   
                gauss_y = my_Gauss(x, *popt)                #
                fwhm_Gmax = popt[2]*2.35482                 #
                area_Gmax = sqrt(2*math.pi)*popt[0]*popt[2] #
#               print('fwhm_Gmax = ', fwhm_Gmax)            #
#               print('area_Gmax = ', area_Gmax)            #
                #plot selected eta data                     #
                plt.plot(x, gauss_y, 'k')                   #
                plt.autoscale(False)
                
                plt.draw()                                  # 
                #####################################       #   
                ### collect min/max data #####
                area = (area_Gmin+area_Gmax)/2
                area_err = sqrt( ((area_Gmin-area_Gmax)/2)**2 )
                fwhm = (fwhm_Gmin + fwhm_Gmax)/2
                #fwhm_err = sqrt( ((fwhm_Gmin - fwhm_Gmax)/2)**2 )
                fwhm_err = sqrt( ((fwhm_Gmin - fwhm_Gmax))**2 ) # larger error bar
                data_sig_mm.append([[area, area_err], [fwhm, fwhm_err]])
                
                

            ##################################################################   
            ### calculation of the error bar and averaging ###################           
            area = []
            area_G = []
            area_err_fitR = 0
            area_err_offsetR = 0
            area_errR = 0
            for i in range(len(data_sig)):
 #               print('sig_area    = ', data_sig[i][0][0], '+/-', data_sig[i][0][1], 'area_G = ', data_sig[i][0][2] )
                print('peak_G', data_sig[i][0][3])
                area.append(data_sig[i][0][0])
                area_G.append(data_sig[i][0][2])
                area_err_fitR += (data_sig[i][0][1]/data_sig[i][0][0])/len(data_sig)
 #               print('sig_area_mm = ', data_sig_mm[i][0][0], '+/-', data_sig_mm[i][0][1] )
 #               print(data_sig_mm[i][0][1]/data_sig_mm[i][0][0])
                area_err_offsetR += (data_sig_mm[i][0][1]/data_sig_mm[i][0][0])/len(data_sig_mm)
            area_av = average(area)
            area_G_av = average(area_G)
            
            print( ' ' )
            area_err_sctrR = area_av[2]/area_av[0]
  #          print('area_err_sctrR = ', area_err_sctrR*100 )
  #          print('area_err_fitR = ', area_err_fitR*100)
  #          print('area_err_offsetR = ', area_err_offsetR*100)
            area_errR = sqrt(area_err_sctrR**2 + area_err_fitR**2 + area_err_offsetR**2)
  #          print('area_errR', area_errR*100)
            print('aera relative errors')
            print('total = scattering + fit + 0_offset ')
            print("%3.1f" %(area_errR*100),'%  =',
                  "%3.1f" %(area_err_sctrR*100),'% + ', 
                  "%3.1f" %(area_err_fitR*100) ,'% + ',
                  "%3.1f" %(area_err_offsetR*100) ,'%')                                    
            
            fwhm = []
            fwhm_G = []
            fwhm_err_fitR = 0
            fwhm_err_offsetR = 0
            fwhm_errR = 0
            for i in range(len(data_sig)):
#                print('sig_fwhm    = ', data_sig[i][1][0], '+/-', data_sig[i][1][1] )
#                print('sig_fwhm_mm = ', data_sig_mm[i][1][0], '+/-', data_sig_mm[i][1][1] )
                fwhm.append(data_sig[i][1][0])
                fwhm_G.append(data_sig[i][1][2])
                fwhm_err_fitR += (data_sig[i][1][1]/data_sig[i][1][0])/len(data_sig)
                fwhm_err_offsetR += (data_sig_mm[i][1][1]/data_sig_mm[i][1][0])/len(data_sig_mm)
            fwhm_av = average(fwhm)
            fwhm_G_av = average(fwhm_G)
            fwhm_err_sctrR = fwhm_av[2]/fwhm_av[0]
            fwhm_errR = sqrt(fwhm_err_sctrR**2 + fwhm_err_fitR**2 + fwhm_err_offsetR**2)
            print('FWHM relative errors')
            print('total = scattering + fit + 0_offset ')
            print("%3.1f" %(fwhm_errR*100),'%  =',
                  "%3.1f" %(fwhm_err_sctrR*100),'% + ', 
                  "%3.1f" %(fwhm_err_fitR*100) ,'% + ',
                  "%3.1f" %(fwhm_err_offsetR*100) ,'%')
            
            print('         Area_G           /           FWHM_G')
            print('         Area             /           FWHM')
            
#            print(sig_area0, '+/-', area_err, '  ///   ', sig_fwhm0, '+/-', fwhmErr)
            print(
                  "{:.3e}".format(area_G_av[0]), '+/-', "{:.1e}".format(area_av[0]*area_errR), 
                  '    ', 
                  "{:.3e}".format(fwhm_G_av[0]), '+/-', "{:.1e}".format(fwhm_av[0]*fwhm_errR),
                  '  //'
                  )
            
            
            
#            print(sig_area0, '+/-', area_err, '  ///   ', sig_fwhm0, '+/-', fwhmErr)
            print(
                  "{:.3e}".format(area_av[0]), '+/-', "{:.1e}".format(area_av[0]*area_errR), 
                  '    ', 
                  "{:.3e}".format(fwhm_av[0]), '+/-', "{:.1e}".format(fwhm_av[0]*fwhm_errR),
                  '  //'
                  )
            
            sel_x.clear()
            data_sig.clear()
            data_sig_mm.clear()
            
    if event.key == 'escape' :
        print('##############################################')
        print('clear variables on escape key ')
        sel_x.clear()
        data_sig.clear()
        data_sig_mm.clear()
        
        
        for j in range(len(file_matrix)):
            nul = zeros(len(file_matrix[j][1]))
            L_names[j].set_xdata(file_matrix[j][1])                  #
            L_names[j].set_ydata(nul)
        plt.draw()

def average(A):      # averaging elements of the arrey 
    
    sumA = 0         # sum of the all elements
    avA = 0          # return value B[0], average of the all elements
    deltaAsum = 0    # sum of the squared deviations
    deltaA = 0       # return value B[2], deviation from avarage
    n = len(A)
    for i in range (0, n):
        sumA += A[i]
        
    avA = sumA/n
    
    for i in range (0, n):
        deltaAsum += ( A[i] - avA )**2
        
    deltaA = (deltaAsum/(n*(n-1)))**0.5
#    print(avA, "+/-", deltaA) # for debugging 
    B = [avA, "+/-", deltaA]
    return B
                
    
            
def get_fwhm(x, y):
    
    A0 = max(y)
    i_cross = [] # index of the point befor y = max(y)/2

    for i in range(len(x)-1) :
#       print('y[i]-max(y)/2 = ', y[i]-max(y)/2)
        if (y[i]-A0/2)*(y[i+1]-A0/2) < 0 :
            i_cross.append(i)
#                    print('i_width [i] = ', i_width)
        else : pass
    
    i_width = []
    i_width.append(i_cross[0])
    i_width.append(i_cross[-1])
    
#    print('i_width = ', i_width)
    # determination of the points at y = max(y)/2 in linear interpolation
    a1 = (y[i_width[0]+1]-y[i_width[0]])/(x[i_width[0]+1]-x[i_width[0]])
    b1 = (-1)*a1*x[i_width[0]]+y[i_width[0]]
    x1 = 1/a1 * (A0/2-b1)
    y1 = a1*x1 + b1
    
    a2 = (y[i_width[1]+1]-y[i_width[1]])/(x[i_width[1]+1]-x[i_width[1]])
    b2 = (-1)*a2*x[i_width[1]]+y[i_width[1]]
    x2 = 1/a2 * (A0/2-b2)
    y2 = a1*x1 + b1
    
    point1 = [x1, y1]
    point2 = [x2, y2]
    return (point1, point2)
    
            
def integral_trapez(x, y):
    I = 0
    for i in range(0, len(x)-1):
        I += (x[i+1]-x[i])*(y[i+1]+y[i])/2
#    print('final I = ', I)
    return I   
    
def my_Gauss(x, A, x0, w):
    y = A * np.exp(-(x-x0)**2/(2*w**2))    
    return y

def open_folder():
    app = QtWidgets.QApplication(sys.argv)
    dname = QtWidgets.QFileDialog.getExistingDirectory()
    path = dname + '/'
#    print(path)
#    print('dname =', dname)
#    listOffiles = [f for f in os.listdir(dname) if f.endswith("n_shft.txt")]
    listOffiles = [f for f in os.listdir(dname)]
#    print('riding file names in the folder')
#    print('listOffiles =', listOffiles)
    i=0
    file_matrix = []
    for i in range(len(listOffiles)):
#        print('listOffiles[f] = ', listOffiles[i])
        filename = listOffiles[i]
        data = read_shft_file(path, filename)
#        plt.plot(data[1], data[2])
        file_matrix.append(data)
    return (file_matrix, path)
    
    
def read_shft_file(path, filename):
    fullname = path + filename
    h = 4 # number of header elements
    if os.path.isfile(fullname):
        f = open(fullname,'r')       # Open the file for read
        line = f.readlines()   # read all lines of the file
        f.close()
        l = int(len(line))  # number of lines elements
#        print('number of lines =', l)       # ptint number for debaging

        Header = list()    
        for i in range(h):
            Htmp = line[i]
            Header.append(Htmp)
#            print( 'Header line[i] = ', line[i])
#        print ('read header = ', Header)
        
        n = l - h
#         print('number of data lines =', n) # for debaging
        data_x = zeros(n)   # relative wavenumbers 
        data_y = zeros(n) # main data
        data_y1 = zeros(n) # estimation of the Io jitter max
        data_y2 = zeros(n) # estimation of the Io jitter min
        for i in range(h, l):
            # change comma to point two times:
#            line[i] =  line[i].replace(',','.',2)
            data_tmp = line[i].split('\t')
            data_x[i-h] = data_tmp[0]
            data_y[i-h] = data_tmp[1]
            data_y1[i-h] = data_tmp[2]
            data_y2[i-h] = data_tmp[3]

    else: 
        print('This file does not exist!')
        data_x = []
        data_y = []
        data_y1 = []
        data_y2 = []
    finfo = [filename, Header]
    data = [finfo, data_x, data_y, data_y1, data_y2]
#    print (data)
    return data


if __name__ == "__main__":
    int_Folder()