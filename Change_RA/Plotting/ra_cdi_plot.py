import shutil
import json, os
import numpy as np
import scipy
import ZachsModules as zm
import sys
from matplotlib.ticker import ScalarFormatter
zm.zp.updateRCParams(**{'text.usetex':False})
plt = zm.plt

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
linestyle_list = ["-.", "--", "-"]
grid_on = True
grid_width = 1.5
x_label = r"$\mathregular{R_{A}}$"

fig, cdi = plt.subplots(figsize=plot_size)
cdi.set_xlabel(x_label, fontsize=plot_fontsize, fontname=plot_font, fontstyle='italic')
cdi.set_ylabel(r"$\mathregular{\%}$ Drag Reduction", fontsize=plot_fontsize, fontname=plot_font, fontstyle='normal')
cdi.tick_params(axis='both', labelsize=tick_size, width=tick_width, length=tick_length)
cdi.spines['top'].set_linewidth(axis_thickness)  
cdi.spines['right'].set_linewidth(axis_thickness)  
cdi.spines['bottom'].set_linewidth(axis_thickness)  
cdi.spines['left'].set_linewidth(axis_thickness)
cdi.grid(grid_on, linewidth=grid_width)
cdi.minorticks_off()
cdi.ticklabel_format(axis='y', style='plain', scilimits=(0, 0))

# cdi.set_xlim(4, 14) 
# #cdi.set_xticks([0.125, 0.25, 0.5, 0.75, 1.0])
# xticks = np.linspace(4, 14, num=6)
# cdi.set_xticks(xticks)

# cdi.set_ylim(0.0, 70)
# ymin, ymax = cdi.get_ylim()
# y_ticks = np.linspace(0.0, 70, num=8)
# cdi.set_yticks(y_ticks)

current_folder = os.getcwd()
all_txt_files = [f for f in os.listdir(current_folder) if f.endswith('.txt')]
count = 0
# file_name = all_txt_files[0]
all_dict = {}
for file_name in all_txt_files:
    f = open(file_name, 'r+')
    data = f.readlines()
    
    if "8D9T" in file_name:
        num_lines = 30
    
    elif "4D5T" in file_name:
        num_lines = 24
    
    elif "0D5T" in file_name:
        num_lines = 18 
    
    elif "1D5T" in file_name:
        num_lines = 21
    
    elif "0D1T" in file_name:
        num_lines = 14
        
    cases = int(len(data)/num_lines)
    divide1 = "Ra: "
    divide2 = ",    Rt"
    coeff = "Ra"
    points = []
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
    points = np.array(points)

    label_name = file_name.split("cp")[1].split("_init")[0]
    all_dict[label_name] = points
    cdi.plot(points[:, 0], points[:, 1], color="black", linewidth=plot_linewidth, linestyle=linestyle_list[count], label=label_name)
    count+=1
#cdi.set_xlim(min(points[:, 0]), max(points[:, 0]))


# d0t1 = np.delete(np.delete(np.delete(np.delete(all_dict["0D1T"], (10), axis=0), (26), axis=0), (52), axis=0), 85, axis=0)
# d4t5 = all_dict["4D5T"]
# d0t5 = np.delete(np.delete(np.delete(np.delete(all_dict["0D5T"], (10), axis=0), (26), axis=0), (52), axis=0), 85, axis=0)


# reduc_45 = 100*abs(d4t5 - d0t1)/d0t1
# reduc_05 = 100*abs(d0t5 - d0t1)/d0t1
# cdi.plot(points[:, 0], reduc_45[:, 1], color="black", linewidth=plot_linewidth, linestyle=linestyle_list[2], label="4D5T")
# cdi.plot(points[:, 0], reduc_05[:, 1], color="black", linewidth=plot_linewidth, linestyle=linestyle_list[1], label="0D5T")
# cdi.legend(fontsize=14)