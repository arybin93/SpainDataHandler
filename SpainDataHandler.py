# -*- coding: utf-8 -*-
"""
Scripts for processing data from Spain shelf

Created on Tue Aug  9 22:21:58 2016

@author: Artem Rybin
@email: arybin93@gmail.com

inputs files formats:
During_storm_P1_salinity.dat
    LON    DEPTH    salinity  
-9.55848    -32.52     35.74

output files formats:

"""
import os
import numpy as np
import matplotlib.pyplot as plt

def read_file(fname):
    print("read file ", fname)
    data = np.loadtxt(fname)
    return data
    
def write_file(fname):
    print("write file %s", fname)

# Foffonoff state sea water 
def get_rho(temp, sal):
    print("get rho")
    R00 = 1000;
    A = 999.842594 + 6.793952e-2*temp - 9.09529e-3*temp**2 + 1.001685e-4*temp**3  - 1.120083e-6*temp**4 + 6.536332e-9*temp**5
    B = 8.24493e-1 - 4.0899e-3*temp + 7.6438e-5*temp**2 - 8.2467e-7*temp**3 + 5.3875e-9*temp**4
    C = -5.72466e-3  +  1.0227e-4*temp  -  1.6546e-6*temp**2
    D = 4.8314e-4
    E = 19652.21 + 148.4206*temp - 2.327105*temp**2 + 1.360477e-2*temp**3 - 5.155288e-5*temp**4
    F = 54.6746-.603459*temp  + 1.09987e-2*temp**2  -  6.167e-5*temp**3  
    G = 7.944e-2+1.6483e-2*temp - 5.3009e-4*temp**2
    H = 3.239908  +  1.43713e-3*temp  +  1.16092e-4*temp**2   -  5.77905e-7*temp**3
    I = 2.2838e-3 - 1.0981e-5*temp  -  1.6078e-6*temp**2
    J = 1.91075e-4
    M = 8.50935e-5  -  6.12293e-6*temp  +  5.2787e-8*temp**2
    N = -9.9348e-7  +  2.0816e-8*temp  +   9.1697e-10*temp**2
    R0 = A + B*sal + C*sal**1.5  +  D*sal**2
    P = 0
    RK = E+F*sal+G*sal**1.5+(H+I*sal+J*sal**1.5)*P+(M+N*sal)*P**2
    rho = R0/(1-P/RK)-R00
    return rho

if __name__ == '__main__':
    print("Start")
    
    os.chdir("D:\ScientificWork\Work\SpainData\source_data")
    filenamesTemperature = ["Pre_storm_P1_temperature.dat", "During_storm_P1_temperature.dat", "Post_storm_P1_temperature.dat",
                            "Pre_storm_S2_temperature.dat", "During_storm_S2_temperature.dat", "Post_storm_S2_temperature.dat"];
    
    filenamesSalinity = ["Pre_storm_P1_salinity.dat", "During_storm_P1_salinity.dat", "Post_storm_P1_salinity.dat",
                        "Pre_storm_S2_salinity.dat", "During_storm_S2_salinity.dat", "Post_storm_S2_salinity.dat"];
                        
    for i in range(1):
        dataTemp = read_file(filenamesTemperature[i]) 
        dataSal = read_file(filenamesTemperature[i])                    
        rho = get_rho(dataTemp[:,2], dataSal[:,2])                
                   
    print(rho)
    plt.plot(rho[1:15])