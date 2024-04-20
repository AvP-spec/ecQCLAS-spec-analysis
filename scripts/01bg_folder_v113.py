# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 19:32:49 2021
updated 18.04.2024
background (bg) normalization
-------------------------------------------------------------------------------------------------------------
before start:
-------------
set offset_err = 0.05  # default: 0.05 = 5% from max signal
The uncertainty of the zero line (absence of the signal on the detector) determined from figures of the previous step
Also specific measurements should be made for this estimation

Input:
-------------
on start select folder with Qt open file dialog which was processed in previous step (sig or ref)

Interaction:
-------------
1) select bg region with mouse cursor and SPACE key (laser intensity which is not disturbed by absorption )
2) press UP and DOWN keys to change polynomial order of the fit
2a) ESCAPE - start selection again
3) ENTER to save results
4) script can be repeated on the created file

Output:
-------------
data normalized on background saved in folder with name  "selected_folder" + "_n"
 
@author: Pipa
"""

import os
import sys
from numpy import zeros # used to load the asci file
from matplotlib import pyplot as plt
from matplotlib.widgets import Cursor
import numpy as np
from PyQt6 import QtWidgets
#from scipy.optimize import curve_fit


ofset_err = 0.05  # 0.05 = 5% uncertanty of the of set level 
polynomial_order = 3 
#sin_counter = 1
sel_x = [] 
sel_y = []
space_counter = [] 
mtrx_poly = []

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


def bg_folder():
    get_data_folder()
    data = open_folder() # data = file_matrix, path
    file_matrix = data[0]
#    print('file_matrix =', file_matrix)
    fig_trm = plt.figure(num='transmittion Figure', figsize=(15, 7), 
                              facecolor='#E0EDCA' ,clear=True,)
    ax_trm = fig_trm.add_subplot(111, facecolor='#DAE1CB')
    
    L_names = []
    P_names = []
    for i in range(len(file_matrix)):
 #       print('len(file_matrix[i][1])= ', len(file_matrix[i][1]) )
        # plot data 
        ax_trm.plot(file_matrix[i][1], file_matrix[i][2])
        # plot for points selection
        points_name = 'P'+ str(i)
        points_name, = ax_trm.plot(file_matrix[i][1], file_matrix[i][2], '.')
        P_names.append(points_name)
        # plot lines for polynomial fit
        line_name = 'L'+ str(i)
        line_name, = ax_trm.plot(file_matrix[i][1], file_matrix[i][2], 'w')
        L_names.append(line_name)
        # format the list of the selected points for bg polynomial fit
        x = []
        y = []
        sel_x.append(x)
        sel_y.append(y)
        # format the list for bg polynomial fit
        poly = []
        mtrx_poly.append(poly)
        
    cursor = Cursor(ax_trm, useblit=True, color='red', linewidth=2)
    cid_key = fig_trm.canvas.mpl_connect('key_press_event', lambda event: 
                                         on_key_bg(event, fig_trm, ax_trm, data,
                                                    L_names, P_names)
                                         )
    plt.show()
    return (cursor, cid_key)


def on_key_bg(event, fig_trm, ax_trm, data, L_names, P_names):    
#    print('you pressed', event.key, event.xdata, event.ydata)
    # data = file_matrix, path
    # file_matrix = file[i]
    #file = [finfo, data_x, data_y]
    
    file_matrix = data[0]
    
    global polynomial_order
    poly = []
    
    if event.key == ' ' :
#        print('event.xdata =',  event.xdata)
        space_counter.append(event.xdata)
        if len(space_counter)== 1:
            print('space is pressed first time')
                
        elif len(space_counter)== 2:
            print('space is pressed second time')
            min_nu = min(space_counter[0], space_counter[1])
            max_nu = max(space_counter[0], space_counter[1])
            
            iii = 0
            for iii in range( len(file_matrix) ):
                   
                    i=0
                    for i in range( len(file_matrix[iii][1]) ):
                        if file_matrix[iii][1][i]< min_nu :
                            pass
                        elif file_matrix[iii][1][i] >= min_nu :
                            if file_matrix[iii][1][i] <= max_nu:
                                sel_x[iii].append(file_matrix[iii][1][i])
                                sel_y[iii].append(file_matrix[iii][2][i])
                        elif file_matrix[iii][1][i] > max_nu: 
                            pass
                        else: print('else condition in selection of data points')
            
#            print('len(sel_x) = ', len(sel_x))
            iii = 0
            for iii in range( len(file_matrix) ):
 #               ax_trm.plot(sel_x[iii],  sel_y[iii], 'b.')
                P_names[iii].set_xdata(sel_x[iii])
                P_names[iii].set_ydata(sel_y[iii])
                
            iii = 0
            for iii in range( len(file_matrix) ):
                z = np.polyfit(sel_x[iii], sel_y[iii], polynomial_order)
                poly = np.polyval(z, file_matrix[iii][1])
                L_names[iii].set_xdata(file_matrix[iii][1])
                L_names[iii].set_ydata(poly)
                mtrx_poly[iii] = list(poly.copy())
                
            space_counter.clear()
            fig_trm.canvas.draw()    

        else : 
            print('else condition awoked in bg selection on key 1')
            space_counter.clear()
            sel_x.clear()
            sel_y.clear()
            poly.clear()
            mtrx_poly.clear()
            data.clear()

    if event.key == 'up' :
#       print('if condition for key up')
        polynomial_order += 1
        print('polynomial_order =', polynomial_order)
        for iii in range( len(file_matrix) ):
                z = np.polyfit(sel_x[iii], sel_y[iii], polynomial_order)
                poly = np.polyval(z, file_matrix[iii][1])
                L_names[iii].set_xdata(file_matrix[iii][1])
                L_names[iii].set_ydata(poly)
                mtrx_poly[iii] = list(poly.copy())

        plt.draw()         
 #       print('key up len(mtrx_poly) = ', len(mtrx_poly))
 #       print('key up len(mtrx_poly[0]) = ', len(mtrx_poly[0]))

    elif event.key == 'down' :
#       print('if condition for key down')
        polynomial_order -= 1
        print('polynomial_order =', polynomial_order)
        for iii in range( len(file_matrix) ):
                z = np.polyfit(sel_x[iii], sel_y[iii], polynomial_order)
                poly = np.polyval(z, file_matrix[iii][1])
                L_names[iii].set_xdata(file_matrix[iii][1])
                L_names[iii].set_ydata(poly) 
                mtrx_poly[iii] = list(poly.copy())

        plt.draw()
 #       print('key up len(mtrx_poly) = ', len(mtrx_poly))
        
    elif event.key == 'escape' :
        polynomial_order = 3
        space_counter.clear()
        poly.clear()        
        sel_x.clear()
        sel_y.clear()
        mtrx_poly.clear()
        iii = 0
        for iii in range( len(file_matrix) ):
            x = []
            y = []
            poly = []
            sel_x.append(x)
            sel_y.append(y)
            mtrx_poly.append(poly)
            P_names[iii].set_xdata(file_matrix[iii][1])
            P_names[iii].set_ydata(file_matrix[iii][2])
            L_names[iii].set_xdata(file_matrix[iii][1])
            L_names[iii].set_ydata(file_matrix[iii][2])

        fig_trm.canvas.draw()
                
        print('escape polynomial_order = ', polynomial_order)
    

    
    # ENTER # ENTER # ENTER # ENTER # ENTER # ENTER # ENTER #
    elif event.key == 'enter' :
        plt.close('transmittion Figure')
        global ofset_err

        mtrx_sig_n = []
        mtrx_sig_n_min = []
        mtrx_sig_n_max = []
        
        iii = 0
        for iii in range( len(file_matrix) ):
            sig = file_matrix[iii][2]            
            sig_n = sig / mtrx_poly[iii]
            mtrx_sig_n.append(sig_n)
            
            sig_min = np.asarray([number - ofset_err*max(sig) for number in sig ])
            y_fit_min = np.asarray([number - ofset_err*max(sig) for number in mtrx_poly[iii]])
            sig_n_min = sig_min / y_fit_min 
            mtrx_sig_n_min.append(sig_n_min)
            
            sig_max = np.asarray([number + ofset_err*max(sig) for number in sig ])
            y_fit_max = np.asarray([number + ofset_err*max(sig) for number in mtrx_poly[iii]])
            sig_n_max = sig_max / y_fit_max         
            mtrx_sig_n_max.append(sig_n_max)
            
#            plt.plot(file_matrix[iii][1], sig_n_min, color='xkcd:silver' )
#            plt.plot(file_matrix[iii][1], sig_n_max, color='xkcd:silver' )
 #           plt.plot(file_matrix[iii][1], sig_n )
        
  #      plt.draw()
        
        ##########################
        # write to file
        ######################
        
        # data = file_matrix, path
        # file_matrix = file[i]
        #file = [finfo, data_x, data_y]
        # finfo = [filename, Header]
        
  #      print('old_path = ', data[1] )
        old_path_splt = data[1].split('/')
  #      print('path_splt = ',  old_path_splt)
        old_path_splt[-2] = old_path_splt[-2] + '_n'
        
        for i in range (0, len(old_path_splt)-1):
            path = '/'.join(old_path_splt) 
            
        if os.path.exists(path):
           print("folder alrady exist") 
           sys.exit(0)
        else:
           os.mkdir(path) 
 #          pass
#        print('path = ', path )
        
        iii = 0        
        for iii in range(0, len(file_matrix) ):
#            print('old_fname = ', data[0][iii][0][0])
            old_fname_spit = data[0][iii][0][0].split('.')
            fname = old_fname_spit[0] + '_n' + '.txt'
            print('fname =', fname)           
#            print('old_header = ', data[0][iii][0][1])
            etalon_cnst = data[0][iii][0][1][2]
            print(etalon_cnst)
            col_name = 'x_rel_nu          sig_n       sig_min         sig_max '
            Header = [path +'\n', fname+'\n', etalon_cnst, col_name+'\n']
            rel_nu = file_matrix[iii][1]
            sig_n = mtrx_sig_n[iii]
            
            sig_n_min = mtrx_sig_n_min[iii]
            plt.plot(rel_nu, sig_n_min, color='xkcd:silver' )
            sig_n_max = mtrx_sig_n_max[iii]
            plt.plot(rel_nu, sig_n_max, color='xkcd:silver' )
            columns = [rel_nu, sig_n, sig_n_min, sig_n_max]
            plt.plot(rel_nu, sig_n )
            
  #          plt.draw()
            WriteFile(path, fname, Header, columns)
                 
        plt.draw()
        space_counter.clear()
        sel_x.clear()
        sel_y.clear()
        mtrx_poly.clear()
        poly.clear()
        mtrx_sig_n.clear()
        mtrx_sig_n_min.clear()
        mtrx_sig_n_max.clear()
        
        

def WriteFile(path, fname, Header, columns):
    
#    namenew = name + ".txt"
    writename = os.path.join(path, fname)
    
    if os.path.isfile(writename):
         sys.exit('WriteFile = File already exist!')
    
    print('recording in file  ' , writename )
    filestream = open(writename,'w')
    
    for i in range(0, len(Header)):
 #       filestream.writelines(Header[i]+'\n') # '\n' is included during readout
        filestream.writelines(Header[i])
        
    line = ""  
    for i in range(0, len(columns[0])):        # for every line i  
        line = ""                              # start with empty line
        for j in range(0, len(columns)):       # write number of collumns 
            line += str(columns[j][i]) +'\t'
        filestream.writelines(line+'\n')       # write to the file 
#        print(line)                           # for debugging 
        
    filestream.close()    



def open_folder():
    app = QtWidgets.QApplication(sys.argv)
    dname = QtWidgets.QFileDialog.getExistingDirectory()
    path = dname + '/'
#    print(path)
#    print('dname =', dname)
    listOffiles = [f for f in os.listdir(dname)  if f.endswith(".txt")]
#    print('riding file names in the folder')
    print('listOffiles =', listOffiles)
    i=0
    file_matrix = []
    hlines = 4 # number of header lines - 1 
    for i in range(len(listOffiles)):
#        print('listOffiles[f] = ', listOffiles[i])
        filename = listOffiles[i]
        data = read_file(path, filename, hlines)
#        plt.plot(data[1], data[2])
#        print("header =", data[0])
        file_matrix.append(data)
    
    return (file_matrix, path)
    
    
def read_file(path, filename, hlines):
    fullname = path + filename
    h = hlines # number of the header lines
    Header = list() 
    if os.path.isfile(fullname):
        f = open(fullname,'r')       # Open the file for read
        line = f.readlines()   # read all lines of the file
        f.close()
        l = int(len(line))  # number of lines elements
#        print('number of lines =', l)       # ptint number for debaging

#        Header = list()    
        for i in range(h):
            Htmp = line[i]
            Header.append(Htmp)
#            print( 'Header line[i] = ', line[i])
#        print ('read header = ', Header)
        
        n = l - h
#         print('number of data lines =', n) # for debaging
        data_x = zeros(n)   # relative wavenumbers 
        data_y = zeros(n) # main data

        for i in range(h, l):
            # change comma to point two times:
#            line[i] =  line[i].replace(',','.',2)
            data_tmp = line[i].split('\t')
            data_x[i-h] = data_tmp[0]
            data_y[i-h] = data_tmp[1]

    else: 
        print('This file does not exist!')
        data_x = []
        data_y = []

    finfo = [filename, Header]
    data = [finfo, data_x, data_y]
#    print (data)
    return data

if __name__ == "__main__":
    bg_folder()