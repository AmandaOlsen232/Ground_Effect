import shutil
import json, os
import numpy as np
import scipy
import ZachsModules as zm
import sys
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
# plot_size = (6.7, 5.2)
# dist_plot_size = (6.7, 5.2)
plot_size = (6.7, 5.2)
dist_plot_size = (6.7, 5.2)

linestyle_list = [":", "--", "-.", "-"]
grid_on = True
grid_width = 1.5
x_label = r"$\mathregular{h/b}$"

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

cdi.set_xlim(0.1, 1.0) 
#cdi.set_xticks([0.125, 0.25, 0.5, 0.75, 1.0])
xticks = np.linspace(0.1, 1.0, num=7)
cdi.set_xticks(xticks)

cdi.set_ylim(0.001, 0.011)
y_ticks = np.linspace(0.001, 0.011, num=6)
cdi.set_yticks(y_ticks)

cdi.ticklabel_format(axis='y', style='plain', scilimits=(0, 0))

#plt.gca().invert_yaxis()
#plt.axis('equal')

current_folder = os.getcwd()
all_txt_files = ["RA8.0_RT1.0_CL0.5_hbchange_cp0D1T_init0_LLsolvernonlinear_didistquadratic_grid100_optmethodSLSQP_dibounds(-90, 90)_opttol1e-12_optmaxiter1000_cases100.txt",
                 "RA8.0_RT1.0_CL0.5_hbchange_cp0D5T_init0_LLsolvernonlinear_didistquadratic_grid100_optmethodSLSQP_dibounds(-90, 90)_opttol1e-12_optmaxiter1000_cases100.txt",
                 "RA8.0_RT1.0_CL0.5_hbchange_cp1D5T_init0_LLsolvernonlinear_didistquadratic_grid100_optmethodSLSQP_dibounds(-90, 90)_opttol1e-12_optmaxiter1000_cases100.txt",
                 "RA8.0_RT1.0_CL0.5_hbchange_cp4D5T_initprevious_LLsolvernonlinear_didistquadratic_grid100_optmethodSLSQP_dibounds(-90, 90)_opttol1e-12_optmaxiter1000_cases100_full.txt"]

#[all_txt_files[3], all_txt_files[4], all_txt_files[6], all_txt_files[10]]#[all_txt_files[0], all_txt_files[1], all_txt_files[2], all_txt_files[6]]
#all_txt_files = [all_txt_files[3], all_txt_files[6]]#, all_txt_files[5], all_txt_files[6]]
count=0

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
    print(num_lines)  
    cases = int(len(data)/num_lines)
    divide1 = "h/b: "
    coeff = "h/b"
    points = []
    for i in range(1, cases):
        coeff_index = num_lines*i + 2
        CDi_index = num_lines*i + 6
        
        val = data[coeff_index].split(divide1)[1]
        val = float(val)
        
        CDi = float(data[CDi_index].split("fun: ")[1])/10000
        points.append([val, CDi])
    points = np.array(points)

      
    label_name = file_name.split("cp")[1].split("_init")[0]
    plt.plot(points[:, 0], points[:, 1], label=label_name, linewidth=plot_linewidth, linestyle=linestyle_list[count], color="black")
    count+=1
   

plt.legend(fontsize=14)