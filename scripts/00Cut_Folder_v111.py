# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 12:38:54 2021
updated 18.04.2024
Select region of interest and make x calibration
-------------------------------------------------------------------------------------------------------------
before start:
-------------
set etalon_constant
etalon_constant = 1 if it is unknown
examples used for 3 inch Germanium etalon:
etalon_constant = 0.0163 (window 1a, 1b: 1300-1400 cm^-1)
etalon_constant = 0.01617 (window 2a: 2100 cm^-1)

Input:
-------------
on start select folder with LAS data with Qt open file dialog
the folder should not contain folders "sig" and "ref"

Interaction:
-------------
1) the script reads three channels of LAS data which ends with '_etalon.asc', '_reference.asc' '_signal.asc'
2) plot the channels and suggest to select working range of the data (press SPACE at the cursor position)
3) select polynomial order for fitting middle line of etalon channel with up and down arrows and press enter
4) an etalon calibration graph will be shown. select polynomial order >= 9 with UP and DOWN arrows and press ENTER
5) control that min and max on etalon data was determined correctly and close all figures.

Output:
-------------
data calibrated on relative wavenumbers are saved in folders ref and sig

@author: Pipa
"""

import os
import sys
from numpy import zeros # used to load the asci file
from matplotlib import pyplot as plt
from matplotlib.widgets import MultiCursor
import numpy as np
# from PyQt6.QtWidgets import QApplication
# from PyQt6.QtWidgets import QFileDialog
from PyQt6 import QtWidgets

etalon_constant = 0.0508   # for measurements of the etalon constant
#etalon_constant = 0.016317 # (window 1a, 1b: 1300-1400 cm^-1)
#etalon_constant = 0.01617  # window 2a: 2100 cm^-1
polynomial_order = 6
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

def set_data_folder(dname):
    path = os.path.dirname(os.path.abspath(__file__))
    # print(path)
    os.chdir(path)
    os.chdir("..\\")
    # print(os.getcwd())
    with open("set_dir.txt", "wt") as f:
        f.write(dname)


def cut_folder(step_data):
    get_data_folder()
    if step_data == 0:
        app = QtWidgets.QApplication(sys.argv)
        dname = QtWidgets.QFileDialog.getExistingDirectory()
        set_data_folder(dname)

        # fname = QFileDialog.getOpenFileName()
        # print(fname)
        # dname = fname[0].split("/")
        # dname.pop(-1)
        # dname = "/".join(dname)
        # # app.exec() # not required but on the logic of PyQt 6

        print(dname)


        if os.path.exists(dname + '/ref'):
            print("ref folder alrady exist, 00Cut_Folder script")
            sys.exit(0)
        else:
            os.mkdir(dname + '/ref') #for saving new data
        if os.path.exists(dname + '/sig'):
            print("sig folder alrady exist, 00Cut_Folder script")
            sys.exit(0)
        else:
            os.mkdir(dname + '/sig') #for saving new data

        data_matrix = open_LAS_data(dname)  #data = [f_info, etaData[1], etaData[2], sigData[2] ,refData[2]]
        print('len(data_matrix) = ', len(data_matrix))
        file_matrix = []
        # substruct zero level
        for data in data_matrix :
            f_info = data[0]       #  path = f_info[0], meas_id = f_info[1]
            print('meas_id=', f_info[1])
            data_x = data[1]
            eta = data[2]*(-1) - av_min(list(data[2]*(-1)))
            sig = data[3] - av_min(list(data[3]))
            ref = data[4]*(-1) - av_min(list(data[4]*(-1)))
            data = [f_info, data_x, eta, sig, ref]
            file_matrix.append(data)

#           plt.plot(data_x, sig)
#       plt.plot(file_matrix[0][1], file_matrix[0][3] )
        fig_cut = plt.figure(num='cut Figure', figsize=(15, 9), facecolor='#E0EDCA' ,clear=True,)
        fig_cut.subplots_adjust(hspace=0)
        ax_eta = fig_cut.add_subplot(311,  facecolor='#F3E9BE')
        ax_sig = fig_cut.add_subplot(312, sharex=ax_eta, facecolor='#DAE1CB')
        ax_ref = fig_cut.add_subplot(313, sharex=ax_eta, facecolor='#DAE1CB')
        nul = zeros(len(data_x))
        for data in file_matrix :
            data_x = data[1]
            eta = data[2]
            sig = data[3]
            ref = data[4]
            ax_eta.plot(data_x, eta)
            ax_sig.plot(data_x, sig)
            ax_ref.plot(data_x, ref)
        ax_eta.plot(data_x, nul, 'w')
        ax_sig.plot(data_x, nul, 'w')
        ax_ref.plot(data_x, nul, 'w')

        eta_sel_line, = ax_eta.plot(data_x, nul, 'w')
        Poly, = ax_eta.plot(data_x, nul, 'r')

        multi = MultiCursor(fig_cut.canvas, (ax_eta, ax_sig, ax_ref), color='r',
                                            lw=1,  horizOn=False, vertOn=True)

        cid_key = fig_cut.canvas.mpl_connect('key_press_event', lambda event:
                                               on_key_cut(event, fig_cut, ax_eta,
                                                      file_matrix, eta_sel_line, Poly)
                                             )
        plt.show()
        return(multi, cid_key)

    elif step_data[0] == 1 :
        print ('second step started')
        #c_data = [f_info, c_x, c_eta, c_sig, c_ref]
        # c_matrix.append(c_data)            
        # global polynomial_order
        step = step_data[0]
        print('step = ', step)
        c_matrix = step_data[1]
        
        cn_matrix = []
        ######################################################################
        #normalization of the etalon measuremtns, preparatin for peak-finder
        for data in c_matrix :
            data_x = data[1]
            eta_meas = data[2]
            y_poply = Poly_calc(data_x, eta_meas, polynomial_order)
            eta = eta_meas / y_poply
            
            fig_eta = plt.figure(figsize=(15, 7), facecolor='#E0EDCA' ,clear=True,)
            ax = fig_eta.add_subplot(111, facecolor='#DAE1CB')
            ax.plot(data_x, eta, 'g')
            selected_eta, = ax.plot(data_x, eta, 'g.')
            
            # searching for data before crossing y=1 to define 'peaks' of the etalon
            cross_data_eta = []
            cross_data_x = []
            sel_index = []
            
            for i in range(0, len(data_x)-2):
                if (eta[i]-1)*(eta[i+1]-1) < 0:
                    # check that next point above 1 (noise filtration)
                    if (eta[i+1]-1)*(eta[i+2]-1) > 0:
                        sel_index.append(i)
                        cross_data_x.append(data_x[i])
                        cross_data_eta.append(eta[i])
                    elif (eta[i+1]-1)*(eta[i+2]-1) < 0:
                        print('noisy etalon signal, one point was skipt')
                    else: print('else condition in searching for data before crossing y=1')
                    
                else: pass
            
            ax.plot(cross_data_x, cross_data_eta, 'ro')
            plt.show()
            
            ######################################################################
            ######################################################################
            # fit 'peaks' of the etalon by 4-order polinomial end search for extremums
            peak_data_x = []
            peak_eta = []
            number_of_peak_points = []
            extr_eta = []
            extr_x = []
            poly_x = []
            y_poly = []
            
            for i in range(0, len(cross_data_x)-1):
#                print('index j = ', sel_index[i], 'index j+1 = ', sel_index[i+1])
                number_of_peak_points.append(sel_index[i+1] - sel_index[i])
                
                for j in range(sel_index[i], sel_index[i+1]):            
                    peak_data_x.append(data_x[j])
                    peak_eta.append(eta[j])
                    delta = (data_x[j+1]-data_x[j])/4
                    poly_x.append(data_x[j])
                    poly_x.append(data_x[j]+delta)
                    poly_x.append(data_x[j]+2*delta)
                    poly_x.append(data_x[j]+3*delta)
                
                z = np.polyfit(peak_data_x, peak_eta, 4)
                
                for ii in range (len(poly_x)):
                    p = np.polyval(z, poly_x[ii])
                    y_poly.append(p)

                ax.plot(poly_x, y_poly, 'r')
#                print('eta peak plotted')
                if y_poly[len(poly_x)//2]-y_poly[0] > 0:
                    index = y_poly.index(max(y_poly))
#                    print('max index = ', index)
                    extr_x.append(poly_x[index])
                    extr_eta.append(y_poly[index])
                    
                elif y_poly[len(poly_x)//2]-y_poly[0] < 0:
                    index = y_poly.index(min(y_poly))
#                    print('min index = ', index)
                    extr_x.append(poly_x[index])
                    extr_eta.append(y_poly[index]) 
                    
                else : print('else condition 2: search eta min and max ')  
                
                poly_x.clear()
                y_poly.clear()
                peak_data_x.clear()
                peak_eta.clear()
                
            ax.plot(extr_x, extr_eta, 'x', c='black',)     
            plt.draw()

            
            # calculate number of min and max for relative wave numbers (rel_nu)
            rel_nu = []
            global etalon_constant
            for i in range(0, len(extr_x)):
                rel_nu.append(i*etalon_constant/2)

            # park all the data 
            n_data = [data[0], data[1], eta, data[3], data[4], 
                        extr_x, rel_nu]
            cn_matrix.append(n_data)
                
            
        ######################################################################
        # fit callibration curve
        fig_callib = plt.figure(num='callibration', figsize=(15, 7), facecolor='#E0EDCA' ,clear=True,)
        ax_callib = fig_callib.add_subplot(111, facecolor='#000000')
        
#        global polynomial_order
        print('polynomial_order =', polynomial_order)
        
        # plot all callibration data 
        for data in cn_matrix :
            extr_x = data[5]
            rel_nu = data[6]
#            ax_callib.plot(extr_x, rel_nu, '-r.') 
            ax_callib.plot(extr_x, rel_nu, '-.')                          
 #           z = np.polyfit(extr_x, rel_nu, polynomial_order)
 #           p = np.polyval(z, extr_x)
 #           Poly, = ax_callib.plot(extr_x, p, '-y')
 #           ax_callib.plot(extr_x, p, '-y')
        # repeat first data with yellow color
        extr_x0 = cn_matrix[0][5]
        rel_nu0 = cn_matrix[0][6]
        ax_callib.plot(extr_x0, rel_nu0, '-w.')
        z = np.polyfit(extr_x0, rel_nu0, polynomial_order)
        p = np.polyval(z, extr_x0)
        Poly, = ax_callib.plot(extr_x0, p, '-y')

        #c_data = [f_info, c_x, c_eta, c_sig, c_ref]  
        # transfer data
        tr_data = cn_matrix
        cid_key = fig_callib.canvas.mpl_connect('key_press_event', lambda event: 
                                            on_key_cal(event, fig_callib, ax_callib,
                                                    extr_x0, rel_nu0, Poly, tr_data)
                                                )
        plt.draw()
        plt.show()
        return (cid_key)    
            

    else: print('else condition in cut_folder function, step_data non-identified')
    
# on_key_cal - function on_key for etalon callibration
def on_key_cal(event, fig_callib, ax_callib, extr_x, rel_nu, Poly, cn_matrix):
    print('you pressed', event.key, event.xdata, event.ydata)    
    global polynomial_order

    if event.key == 'up' :
#       print('if condition for key right')
        polynomial_order += 1
        print('polynomial_order =', polynomial_order)
        z = np.polyfit(extr_x, rel_nu, polynomial_order)
        p = np.polyval(z, extr_x)
        Poly.set_ydata(p)
        plt.draw()        

    elif event.key == 'down' :
#       print('if condition for key right')
        polynomial_order -= 1
        print('polynomial_order =', polynomial_order)
        z = np.polyfit(extr_x, rel_nu, polynomial_order)
        p = np.polyval(z, extr_x)
        Poly.set_ydata(p)
        plt.draw()
    # ENTER # ENTER # ENTER # ENTER # ENTER # ENTER # ENTER #
    elif event.key == 'enter' :
        print('Enter is pressed')
        plt.close('callibration')
        print('polynomial_order =', polynomial_order)
        #  n_data = [f_info[0], data_x[1], eta[2], sig[3], ref[4], extr_x[5], rel_nu[6] ]
        #  cn_matrix.append(n_data)
        for data in cn_matrix :
            f_info = data[0]
            data_x = data[1]            
            sig = data[3]
            ref = data[4]
            extr_x = data[5]
            rel_nu = data[6]
            
            z = np.polyfit(extr_x, rel_nu, polynomial_order)
            x_rel_nu = np.polyval(z, data_x)
            
            # Write to file

            meas_id = f_info[1]
            #first header line
            h0 = 'etalon_constant = ' + str(etalon_constant)
            
            # one file with all data 
#            path = f_info[0]
#            eta = data[2]
#            columns = [x_rel_nu, data_x, eta, sig, ref]                        
#            fname = meas_id + '_cut_and_cal.txt'            
#            h1 = '      x_rel_nu           data_x          eta           sig      ref'
#            Header = [path, fname, h0, h1]
#            WriteFile(path, fname, Header, columns)
            
            columns = [x_rel_nu, sig]
            path = f_info[0]+ 'sig/'
            fname = meas_id + '_sig.txt'
            h1 = '        x_rel_nu        sig'
            Header = [path, fname, h0, h1]
            WriteFile(path, fname, Header, columns)
        
            columns = [x_rel_nu, ref]
            path = f_info[0]+ 'ref/'
            fname = meas_id + '_ref.txt'
            h1 = '        x_rel_nu        ref'
            Header = [path, fname, h0, h1]
            WriteFile(path, fname, Header, columns)

def WriteFile(path, fname, Header, columns):
    
#    namenew = name + ".txt"
    writename = os.path.join(path, fname)
    
    if os.path.isfile(writename):
        sys.exit('WriteFile = File already exist!')
    
    print('recording in file  ' , writename )
    filestream = open(writename,'w')
    
    for i in range(0, len(Header)):
        filestream.writelines(Header[i]+'\n')
        
    line = ""  
    for i in range(0, len(columns[0])):        # for every line i  
        line = ""                              # start with empty line
        for j in range(0, len(columns)):       # write number of collumns 
            line += str(columns[j][i]) +'\t'
        filestream.writelines(line+'\n')       # write to the file 
#        print(line)                           # for debugging 
        
    filestream.close()    
    
def Poly_calc(data_x, data_y, polynomial_order):
    z = np.polyfit(data_x, data_y, polynomial_order)
    y_poply = []
    for i in range (len(data_x)):
        p = np.polyval(z, data_x[i])
        y_poply.append(p)
    return (y_poply)

def on_key_cut(event, fig_cut, ax_eta, file_matrix, eta_sel_line, Poly):
#    print('you pressed', event.key, event.xdata, event.ydata)
    # select range for analysis
    
    # check first file only
    # data = [f_info, data_x, eta, sig, ref]
    # file_matrix = [data[0], ... , data[n]]
    data_x = file_matrix[0][1]
    eta = file_matrix[0][2]
    
    global polynomial_order
    y_poply = []
    start_index = 0
    stop_index = 0
    if event.key == ' ' :
        sel_x.append(event.xdata)
        if len(sel_x)== 1:
            print('space is pressed first time')
        elif len(sel_x)== 2:
            print('space is pressed second time')
                       
            # search for index of sel_x in data_x - check first file only
            i=0
            while data_x[i] < min(sel_x[0], sel_x[1]):
                i=i+1
            start_index = i
            i = 0
            while data_x[i] < max(sel_x[0], sel_x[1]):
                i = i+1
            stop_index = i
#            print('start/stop_index =', start_index, '/', stop_index)
            index.append(start_index)
            index.append(stop_index)
            
            # plot selected eta data
            eta_sel_line.set_xdata(data_x[start_index:stop_index])
            eta_sel_line.set_ydata(eta[start_index:stop_index])
#            fig_cut.canvas.draw()
            
            # fit selected eta data
            global polynomial_order                
            y_poply = Poly_calc(data_x[start_index:stop_index], eta[start_index:stop_index], polynomial_order)
            Poly.set_xdata(data_x[start_index:stop_index])
            Poly.set_ydata(y_poply)
            fig_cut.canvas.draw()
        
        elif len(sel_x)== 3:
            print('space is pressed third time')
            sel_x.clear()
            index.clear()
            y_poply.clear()
            Poly.set_xdata(file_matrix[0][1])
            nul = zeros(len(file_matrix[0][1]))
            Poly.set_ydata(nul)
            plt.draw()
        
        else: print('else condition on_key_cut function - space counter >3')
        
    if event.key == 'up' :
#        print('if condition for key right')
        polynomial_order += 1
        print('polynomial_order =', polynomial_order)
        start_index = index[0]
        stop_index = index[1]
#        print('right arrow start/stop_index =', start_index, '/', stop_index)
        y_poply = Poly_calc(data_x[start_index:stop_index], 
                                eta[start_index:stop_index], polynomial_order)
        Poly.set_xdata(data_x[start_index:stop_index])
        Poly.set_ydata(y_poply)
            
        plt.draw()
#       fig_cut.canvas.draw()
        
    elif event.key == 'down' :
#        print('if condition for key left')
        polynomial_order -= 1
        print('polynomial_order =', polynomial_order)
        start_index = index[0]
        stop_index = index[1]
#        print('left arrow start/stop_index =', start_index, '/', stop_index)
        y_poply = Poly_calc(data_x[start_index:stop_index], 
                                eta[start_index:stop_index], polynomial_order)
        Poly.set_xdata(data_x[start_index:stop_index])
        Poly.set_ydata(y_poply)
        plt.draw()
   
    # ENTER # ENTER # ENTER # ENTER # ENTER # ENTER # ENTER # ENTER # ENTER #
    elif event.key == 'enter' : 
        c_matrix = []
        if len(index)==2:
            plt.close('cut Figure')
            start_index = index[0]
            stop_index = index[1]
            # data = [f_info, data_x, eta, sig, ref]
            # file_matrix = [data[0], ... , data[n]]
            for data in file_matrix : 
                f_info = data[0]
                x = data[1]
                eta = data[2]
                sig = data[3]
                ref = data[4]
                c_x = x[start_index:stop_index]
                c_eta = eta[start_index:stop_index]
                c_sig = sig [start_index:stop_index]
                c_ref = ref[start_index:stop_index]
                c_data = [f_info, c_x, c_eta, c_sig, c_ref] 
                c_matrix.append(c_data)
            step_data = [1, c_matrix] # the call of the "main" function
            cut_folder(step_data) # EXIT# EXIT# EXIT# EXIT# EXIT# EXIT# EXIT
            index.clear()
            sel_x.clear()
        # EXIT #  # EXIT # # EXIT # # EXIT # # EXIT # # EXIT # # EXIT # # EXIT        
        else:
            print('select region of interest')
        

def av_min(input_list):
    # the function to determine a minimum of the recorded signal 
    # to avoide occasional minimum due to noise increase n
    n = 0 # 2n+1 points will be averaged
    print('function av_min(input_list)')
    print('number of points acounted for signal minimum = ', 2*n+1)

    i = input_list.index(min(input_list))
#    print('i = ', i)
    # avoide an error if minimum at the end of the list
    j_list = []
    for ii in range(i-n, i+n+1):
        if ii < len(input_list):
            j_list.append(ii)
        elif i+3 >= len(input_list):
            j_list.append(ii-len(input_list))
        else:  print('else condition in av_min(), index selection')
#    print ('list j_list=', j_list)
    # averaging: 
    y_min_2 = 0
    for ii in j_list:
        y_min_2 += input_list[ii]/len(j_list)
    return y_min_2       

    
def open_LAS_data(dname):

    path = dname + '/'
    # print(f"{path=}")
    data_matrix = []
    listOffiles = [f for f in os.listdir(dname) if f.endswith("_signal.asc")]
    # print('listOffiles =', listOffiles)

    for i in range(len(listOffiles)):
#        print('listOffiles[i] = ', listOffiles[i])
        list_splt = listOffiles[i].split('_')
        meas_id = list_splt[0] + '_' + list_splt[1]
        
        fname = meas_id+ '_etalon.asc'
#        print('fname =', fname)
        etaData = read_QCL_file (path, fname)
        
        fname = meas_id+ '_reference.asc'
#        print('fname =', fname)
        refData = read_QCL_file (path, fname)

        fname = meas_id+ '_signal.asc'
#        print('fname =', fname)
        sigData = read_QCL_file (path, fname)
        
        f_info = [path, meas_id]
        data = [f_info, etaData[1], etaData[2], sigData[2] ,refData[2] ]
        data_matrix.append(data)
        
    print('files loaded')
    return data_matrix
     
def read_QCL_file (path, filename):
    fullname = path + filename
    # print(f"load file: {fullname}")
    h = 0 # number of header elements
    if os.path.isfile(fullname):
        f = open(fullname,'r')       # Open the file for read
        line = f.readlines()   # read all lines of the file
        f.close()
        l = int(len(line))  # number of lines elements
#        print('number of lines =', l)       # ptint number for debaging
        n = l - h
#         print('number of data lines =', n) # for debaging
        data_x = zeros(n)   # define the varible for wavenumbers and fill it with zeros 
        data_y = zeros(n) # define the varible for intensity and fill it with zeros
        for i in range(h, l):
            # change comma to point two times:
            line[i] =  line[i].replace(',','.',2)
            data_tmp = line[i].split('\t')
            data_x[i-h] = data_tmp[0]
            data_y[i-h] = data_tmp[1]

    else: 
        print('This file does not exist!')
        data_x = []
        data_y = []
    
    data = [filename, data_x, data_y]
#    print (data)
    return data

if __name__ == "__main__":
    cut_folder(0) # - it is "main function at the top"
