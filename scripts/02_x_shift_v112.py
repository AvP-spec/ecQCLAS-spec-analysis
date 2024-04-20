# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 10:44:23 2021
updated 18.04.2024

alignment of relative spectral positions for several spectral records
-------------------------------------------------------------------------------------------------------------
Input:
-------------
on start select folder with Qt open file dialog which was processed in previous step (sig_n or ref_n)

Interaction:
-------------
select spectral line to which wavenamber at maximum absorbance all records will be aligned (SPACE + maus cursor)
wavenamber at maximum = 0

Output:
-------------
aligned spectra will be saved in a folder with "selected_folder" + "_shft", alike "sig_n_shft"

@author: Pipa
"""
import os
import sys
from numpy import zeros # used to load the asci file
from matplotlib import pyplot as plt
from matplotlib.widgets import Cursor
import numpy as np
import math
from PyQt6 import QtWidgets
from scipy.optimize import curve_fit

sel_x = [] # to transfer selected x to index 
index = []

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

def x_shift():
    get_data_folder()
    #  data = [finfo, data_x, data_y, data_y1, data_y2]
    #  finfo = [filename, Header]
    #  file_matrix = [data[i]]
    tmp = open_folder()
    file_matrix = tmp[0]
    path = tmp[1]
#    print('len(file_matrix) = ', len(file_matrix))
    fig_tran = plt.figure(num='transmittance Figure', figsize=(15, 7), 
                              facecolor='#E0EDCA' ,clear=True,)
    ax_tran = fig_tran.add_subplot(111, facecolor='#DAE1CB')
    for i in range(len(file_matrix)):
        ax_tran.plot(file_matrix[i][1], file_matrix[i][2], label=i)
#        plt.plot(file_matrix[4][1], file_matrix[4][3])
#        plt.plot(file_matrix[4][1], file_matrix[4][4])
    plt.legend(loc="upper left")
#    plt.ylim(-1.5, 2.0)

    cursor = Cursor(ax_tran, useblit=True, color='red', linewidth=1)
    cid_key = fig_tran.canvas.mpl_connect('key_press_event', lambda event: 
                                         on_key_1(event, fig_tran, ax_tran, 
                                         file_matrix, path)
                                        )
    plt.show()
    return (cid_key, cursor)
              
    
def on_key_1(event, fig_tran, ax_tran, file_matrix, path):
 #   print('you pressed', event.key, event.xdata, event.ydata)
    print('len(file_matrix) = ', len(file_matrix))
    #data = [finfo, data_x, data_y, data_y1, data_y2]
    #file_matrix = [data[i]]
    start_index = 0
    stop_index = 0
    if event.key == ' ' :
        sel_x.append(event.xdata)    
        if len(sel_x)== 1:
            print('space is pressed first time')
        elif len(sel_x)== 2:
            print('space is pressed second time')
            i=0
            for i in range(len(file_matrix)):
                ii=0
                while file_matrix[i][1][ii] < min(sel_x[0], sel_x[1]):
                    ii+=1
                start_index = ii
                ii = 0
                while file_matrix[i][1][ii] < max(sel_x[0], sel_x[1]):
                    ii+=1
                stop_index = ii
                sel_index = [start_index, stop_index]
                index.append(sel_index)
                
 #           print('index = ', index)
 
            print('path = ', path )
            old_path_splt = path.split('/')
  #            print('path_splt = ',  old_path_splt)
            old_path_splt[-2] = old_path_splt[-2] + '_shft'
            for i in range (0, len(old_path_splt)-1):
                new_path = '/'.join(old_path_splt) 
            
            if os.path.exists(new_path):
                print("folder alrady exist") 
                sys.exit(0)
            else:
                os.mkdir(new_path) 
 #          pass
                print('new_path = ', new_path )

            i=0
            for i in range(len(file_matrix)):
                x = list( file_matrix[i][1][index[i][0]:index[i][1]] )
                y = list( file_matrix[i][2][index[i][0]:index[i][1]] )
#                plt.plot(x, y, 'r')
#                peak_y.append(min(y[i]))
                ab_y = [math.log(1/num) for num in y]
                fwhm_peak = get_fwhm(x, ab_y)
                fwhm_x1 = fwhm_peak[0][0]
 #               fwhm_y1 = fwhm_peak[0][1]
                fwhm_x2 = fwhm_peak[1][0]
 #               fwhm_y2 = fwhm_peak[1][1]
                
                A0 = max(y)
                w0 = (fwhm_x2 - fwhm_x1)/2.35482
                peak_i = y.index(min(y))
                x0 = x[peak_i]
  #              print('x0 = ', x0)
                popt, pcov = curve_fit(my_Gauss, x, ab_y, p0 = [A0, x0, w0])
                # select x at Gaus maximum 
                x0 = popt[1]            
                x_new = [ num - x0 for num in file_matrix[i][1] ]
                plt.close('transmittance Figure')
                plt.plot(x_new, file_matrix[i][2], label=i)
                finfo = file_matrix[i][0]
                fname_old = finfo[0] 
                fname_tmp = fname_old.split('.')
                fname = fname_tmp[0]+'_shft.txt'
                Header = finfo[1]
                columns = [x_new, file_matrix[i][2], file_matrix[i][3], file_matrix[i][4]]
                WriteFile(new_path, fname, Header, columns)
#                print('path for writing =', path)
                
            sel_x.clear()
            index.clear() 
            plt.legend(loc="upper left")
#            plt.show()
#            plt.draw()
#            fig_tran.canvas.draw()
#            plt.close('transmittance Figure')
            ##############################################
            # EXIT # EXIT # EXIT # EXIT # EXIT #
            ##############################################
            
        else: print('else condition for space counter ')
        
    elif event.key == 'escape' :
        sel_x.clear()
        index.clear()
        
def WriteFile(path, fname, Header, columns):
    
#    namenew = name + ".txt"
    writename = os.path.join(path, fname)
    
    if os.path.isfile(writename):
         sys.exit('WriteFile = File already exist!')
    
    print('recording in file  ' , writename )
    filestream = open(writename,'w')
    
    for i in range(0, len(Header)):
#        filestream.writelines(Header[i]+'\n') # '\n' is included during readout
        filestream.writelines(Header[i])
        
    line = ""  
    for i in range(0, len(columns[0])):        # for every line i  
        line = ""                              # start with empty line
        for j in range(0, len(columns)):       # write number of collumns 
            line += str(columns[j][i]) +'\t'
        filestream.writelines(line+'\n')       # write to the file 
#        print(line)                           # for debugging 
        
    filestream.close()    
        
def get_fwhm(x, y):
    
    A0 = max(y)
    i_width = [] # index of the point befor y = max(y)/2

    for i in range(len(x)-1) :
#       print('y[i]-max(y)/2 = ', y[i]-max(y)/2)
        if (y[i]-A0/2)*(y[i+1]-A0/2) < 0 :
            i_width.append(i)
#                    print('i_width [i] = ', i_width)
        else : pass
#    print('i_width get_fwhm = ', i_width)
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
        
def my_Gauss(x, A, x0, w):
    y = A * np.exp(-(x-x0)**2/(2*w**2))    
    return y
    
    
def open_folder():
    app = QtWidgets.QApplication(sys.argv)
    dname = QtWidgets.QFileDialog.getExistingDirectory()
    path = dname + '/'
#    print(path)
#    print('dname =', dname)
    listOffiles = [f for f in os.listdir(dname) if f.endswith(".txt")]
#    print('riding file names in the folder')
#    print('listOffiles =', listOffiles)
    i=0
    file_matrix = []
    for i in range(len(listOffiles)):
#        print('listOffiles[f] = ', listOffiles[i])
        filename = listOffiles[i]
        data = read_nfile(path, filename)
#        plt.plot(data[1], data[2])
        file_matrix.append(data)
        
    
    return (file_matrix, path)
    
    
def read_nfile(path, filename):
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
    x_shift()