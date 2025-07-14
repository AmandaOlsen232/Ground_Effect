import shutil
import json, os
import numpy as np
import scipy
import ZachsModules as zm
import sys
zm.zp.updateRCParams(**{'text.usetex':False})
# plt = zm.plt
import matplotlib.pyplot as plt

current_folder = os.getcwd()
all_txt_files = [f for f in os.listdir(current_folder) if f.endswith('.txt')]

data = []
for file in all_txt_files:
    cp = (file.split("_cp")[1]).split("_LLsolver")[0]
    
    di_cp =int(cp.split("d")[0])
    t_cp = int((cp.split("d")[1]).split("t")[0])
    
    with open(file, 'r') as file:
        content = file.readlines()
    
    cdi = float(content[6].split("fun: ")[1])/10000
    
    data.append([di_cp, cdi])


data = np.array(data)
data = data[data[:, 0].argsort()]

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
    
# if coeff=="h/b":
#         divide1 = "h/b: "
#         x_label = r"$\mathregular{h/b}$"
        
# if coeff=="RT":
#     divide1 = "Rt: "
#     divide2 = ",    CLtarget: "
#x_label = r"$\mathregular{R_{T}}$"

# if coeff=="RA":
#     divide1 = "Ra: "
#     divide2 = ",    Rt:"
x_label = r"$\mathregular{N_{\Gamma}}$"
    
# if coeff=="CL":
#     divide1 = "CLtarget: "
#     divide2 = ",    h/b:"
#     x_label = r"$\mathregular{C_L}$"
        
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

cdi.set_xlim(min(data[:, 0]), max(data[:, 0]))
cdi.ticklabel_format(axis='y', style='plain', scilimits=(0, 0))
cdi.set_xticks([2, 4, 6,  8,  10,  12,  14, 16])

cdi.set_ylim(0.0025, 0.0075)
y_ticks = np.linspace(0.0025, 0.0075, num=6)
cdi.set_yticks(y_ticks)

# cdi.plot(data[:, 0], data[:, 1], color="black", linewidth=plot_linewidth)
# #cdi.scatter(data[:, 0], data[:, 1], color="black", s=100)
# cdi.scatter(data[:, 0], data[:, 1], s=100, c='black', linewidths=5)#, marker='o', edgecolors='black', )

cdi.plot(data[:, 0], data[:, 1], color="black", linewidth=plot_linewidth, marker='o', markersize=10, markeredgewidth=2)