import shutil
import json, os
import numpy as np
import scipy
import ZachsModules as zm
import sys
import math
zm.zp.updateRCParams(**{'text.usetex':False})
plt = zm.plt

current_folder = os.getcwd()
all_txt_files = [f for f in os.listdir(current_folder) if f.endswith('.txt')]
# file_name = all_txt_files[1]

for file_name in all_txt_files:
    f = open(file_name, 'r+')
    data = f.readlines()
    
    if "8d9t" in file_name:
        num_lines = 30
    
    elif "4d5t" in file_name:
        num_lines = 24
    
    elif "0d5t" in file_name:
        num_lines = 18 
    
    elif "1d5t" in file_name:
        num_lines = 21
    
    elif "0d1t" in file_name:
        num_lines = 14
        
    cases = int(len(data)/num_lines)
    divide1 = "CLtarget: "
    divide2 = ",    h/b"
    coeff = "CL"
    points = []
    oge = []
    for i in range(1, cases):
        coeff_index = num_lines*i + 2
        CDi_index = num_lines*i + 6
        
        val = data[coeff_index].split(divide1)[1]
        if coeff!="h/b":            
            val = float(val.split(divide2)[0])
        elif coeff == "h/b":
            val = float(val)
        CDi = float(data[CDi_index].split("fun: ")[1])/10000
        points.append([val, CDi])
        oge.append([val, val**2/(math.pi*8)])
    points = np.array(points)
    oge = np.array(oge)
    #plt.plot(points[:, 0], points[:, 1], label=file_name, color="black")
    #plt.plot(oge[:, 0], oge[:, 1], label="Out of Ground Effect", color="black")
plot_linewidth = 2.5
plot_fontsize = 20
plot_font = "Times New Roman"
tick_size = 18
axis_thickness = 3
tick_length = 4
tick_width = 2
large_font = 20
plot_size = (6.7, 5.2)
dist_plot_size = (6.7, 5.2)
linestyle_list = ["-", "--", "-."]
grid_on = True
grid_width = 1.5
x_label = r"$\mathregular{C_L}$"
    
fig, cdi = plt.subplots(figsize=plot_size)
cdi.set_xlabel(x_label, fontsize=plot_fontsize, fontname=plot_font, fontstyle='italic')
cdi.set_ylabel(r"$\mathregular{C_{Di}}$", fontsize=plot_fontsize, fontname=plot_font, fontstyle='italic')
cdi.tick_params(axis='both', labelsize=tick_size, width=tick_width, length=tick_length)
cdi.spines['top'].set_linewidth(axis_thickness)  
cdi.spines['right'].set_linewidth(axis_thickness)  
cdi.spines['bottom'].set_linewidth(axis_thickness)  
cdi.spines['left'].set_linewidth(axis_thickness)
cdi.grid(grid_on, linewidth=grid_width)
cdi.minorticks_off()
cdi.set_xlim(min(points[:, 0]), max(points[:, 0]))
cdi.ticklabel_format(axis='y', style='plain', scilimits=(0, 0))

cdi.set_xlim(0.125, 1.0) 
cdi.set_xticks([0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0])
#xticks = np.linspace(0.125, 1.0, num=6)
#cdi.set_xticks(xticks)

cdi.set_ylim(0.00, 0.04)
y_ticks = np.linspace(0.00, 0.04, num=5)
cdi.set_yticks(y_ticks)

cdi.plot(points[:, 0], points[:, 1], color="black", linewidth=plot_linewidth, label=r"$\mathregular{h/b=0.25}$")
cdi.plot(oge[:, 0], oge[:, 1], label="Out of Ground Effect", color="black", linewidth=plot_linewidth, linestyle="--")
plt.legend(fontsize=14)
    