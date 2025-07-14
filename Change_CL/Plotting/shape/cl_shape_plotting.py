import shutil
import json, os
import numpy as np
import scipy
import ZachsModules as zm
import sys
zm.zp.updateRCParams(**{'text.usetex':False})
plt = zm.plt


def MachUp(x):
    ## x = [t0, t1, t2, t3, t4, t5, d1, d2, d3, d4, d5]
    tr = x[0]
    twi_points = x[1:t+1]
    di_points = x[t+1:]
    
    #n = int((len(x) - 1) / 2)
    n = int(len(twi_points))
    di_cp = np.linspace(0,1,len(di_points)+1)
    twi_cp = np.linspace(0,1,len(twi_points)+1)

    ## update mounting angle
    f = open('input.json', 'r+')
    data = f.readlines()
    f.seek(0)
    data[41] = ' '*12 + '"mounting_angle": {:.16f},\n'.format(tr)
    f.writelines(data)
    f.close()

    ## update washout file
    f = open('washout.txt', 'w')
    f.write('span,washout\n0.0,0.0\n')
    for i in range(len(twi_points)):
        f.write('{:.16f},{:.16f}\n'.format(twi_cp[i+1],twi_points[i]))
    f.close()
    
    if d != 0:
        di = di_dist(x, span_frac, di_cp, n)
    elif d == 0:
        di = np.zeros_like(span_frac)
        
    ## update dihedral file
    f = open('dihedral.txt', 'w')
    f.write('span,dihedral\n0.0,0.0\n')
    for i in range(len(span_frac)):
        f.write('{:.16f},{:.16f}\n'.format(span_frac[i], di[i]))
    
    f.close()
    
    ##################
    ## run MachUp
    os.system('MachUp.out input.json > out.txt')
MachUp.iter=0

def mb_list(x, s, n):
    mb = 0
    mb_list = []
    for b in s:
        '''determine bounds for equation'''
        b0, b1, g0, g1 = eq_bounds(b, x, s, n)
        
        '''solve for all mb value'''
        a = (g1 - g0 - mb*(b1 - b0))/((b1 - b0)**2)
        q = mb - 2*a*b0
        #c = g0 + a*b0**2 - mb*b0 
        mb = 2*a*b + q
        mb_list.append(mb)
    return np.array(mb_list)

def eq_bounds(b, x, s, n):
    for i in range(len(s) - 1):  
        if s[i] < b <= s[i+1]:
    
            b0 = s[i]
            b1 = s[i+1]
            g1 = x[i+n+1]
            
            if s[0] < b <= s[0+1]:
                g0 = 0.0
            else:
                g0 = x[i+n]
                
        elif s[i] == b and b == 0.0:
            b0 = 0.0 
            b1 = s[i+1]
            g0 = 0.0
            g1 = x[i+n+1]
            
    return b0, b1, g0, g1

def di_dist(x, span, s, n):
    slopes = mb_list(x, s, n)
    
    di_list = []
    for b in span:
        b0, b1, g0, g1 = eq_bounds(b, x, s, n)
        
        for i in range(len(s)):
            if b0 == s[i]:  
                mb = slopes[i]
    
    
        a = (g1 - g0 - mb*(b1 - b0))/((b1 - b0)**2)
        q = mb - 2*a*b0
        c = g0 + a*b0**2 - mb*b0 
        di = a*b**2 + q*b + c
        di_list.append(di)   
        
    return np.array(di_list)

def get_CD(x, runMU=True, jac=False):
    if runMU: MachUp(x)
    ##################
    ## get CD
    f = open('input_forces.json', 'r')
    data = json.load(f)
    f.close()
    if jac:
        return data['total']['MyAirplane']['stability_coordinates']['CD'] * 10000
    else:
        CD = data['total']['MyAirplane']['stability_coordinates']['CD']
        get_CD.iter += 1
        CL = get_CL_offset(x, runMU=False) + get_CL_offset.CLtarget
        hmin = get_min_z(x, runMU=False) + get_min_z.tol
        CLmax = get_sect_CL_max.CLmax - get_sect_CL_max(x, runMU=False)
        print('Iter: {:4d},    CD: {:23.16e},    CL: {:19.16f},  min h: {:19.16f},  sect CLmax: {:19.16f}'.format(get_CD.iter, CD, CL, hmin, CLmax))
        return CD * 10000
get_CD.iter = 0

def get_CL_offset(x, runMU=True, jac=False):
    if runMU: MachUp(x)
    ##################
    ## get CD
    f = open('input_forces.json', 'r')
    data = json.load(f)
    f.close()
    return (data['total']['MyAirplane']['stability_coordinates']['CL'] - get_CL_offset.CLtarget) * 1
get_CL_offset.CLtarget = 0.0

def get_min_z(x, runMU=True, jac=False):
    if runMU: MachUp(x)
    dist = MachUpDist()
    z        = np.array(dist['z'])
    dihedral = np.array(dist['dihedral'])
    twist    = np.array(dist['twist'])
    chord    = np.array(dist['chord'])
    te_z = max(z + 0.75*chord*np.sin(np.deg2rad(twist))*np.cos(np.deg2rad(dihedral)))
    le_z = max(z - 0.25*chord*np.sin(np.deg2rad(twist))*np.cos(np.deg2rad(dihedral)))
    return -(max(te_z, le_z) + get_min_z.tol)
get_min_z.tol = 0.0

def get_sect_CL_max(x, runMU=True, jac=False):
    if runMU: MachUp(x)
    dist = MachUpDist()
    return get_sect_CL_max.CLmax - max(dist['sect_CL'])
get_sect_CL_max.CLmax = 0.0

jac = {'jac':True}
def Jacobian(x, func):
    h = 1e-10
    l = x.size
    grad = np.zeros(l)
    grad[:] = None
    for i in range(l):
        xtemp = np.copy(x)
        xtemp[i] += h
        fpos = func(xtemp, **jac)
        xtemp[i] -= 2.*h
        fneg = func(xtemp, **jac)
        grad[i] = (fpos - fneg) / 2. / h
    return grad

def Jacobian_CD(x, ):
    return Jacobian(x, get_CD)

def Jacobian_CL(x):
    return Jacobian(x, get_CL_offset)

def Jacobian_z(x):
    return Jacobian(x, get_min_z)

def Jacobian_CLmax(x):
    return Jacobian(x, get_sect_CL_max)

def MachUpDist(first=True):
    f = open('input_distributions.txt', 'r')
    data = f.readlines()
    f.close()
    n = int((len(data)-1)/2)
    
    '''
    0 WingName
    1 ControlPoint(x)
    2 ControlPoint(y)
    3 ControlPoint(z)
    4 Chord
    5 Twist(deg)
    6 Sweep(deg)
    7 Dihedral(deg)
    8 Area
    9 Section_Alpha(deg)
    10 Section_CL
    11 Section_CD_parasitic
    12 Section_Cm
    13 CL(Ref)
    14 Section_alpha_L0(deg)
    '''
    
    x = [None]*n
    z = [None]*n
    y = [None]*n
    twist = [None]*n
    dihedral = [None]*n
    sect_CL = [None]*n
    CLref = [None]*n
    chord = [None]*n
    sweep = [None]*n
    area = [None]*n
    sect_alfa = [None]*n
    sect_CD = [None]*n
    sect_Cm = [None]*n
    sect_alfaL0 = [None]*n
    
    if first:
        s = slice(1,n+1)
    else:
        s = slice(n+1,2*n+2)
    
    for i,line in enumerate(data[s]):
        vals = [float(j) for j in line.split()[1:]]
        
        x[i]           = vals[0 ]
        y[i]           = vals[1 ]
        z[i]           = vals[2 ]
        chord[i]       = vals[3 ]
        twist[i]       = vals[4 ]
        sweep[i]       = vals[5 ]
        dihedral[i]    = vals[6 ]
        area[i]        = vals[7 ]
        sect_alfa[i]   = vals[8 ]
        sect_CL[i]     = vals[9 ]
        sect_CD[i]     = vals[10]
        sect_Cm[i]     = vals[11]
        CLref[i]       = vals[12]
        sect_alfaL0[i] = vals[13]
    
    zm.nm.zSort(y, x, z, chord, twist, sweep, dihedral, area, sect_alfa, sect_CL, sect_CD, sect_Cm, CLref, sect_alfaL0, verbose=False)
    
    return {
        'x': x,
        'y': y,
        'z': z,
        'chord': chord,
        'twist': twist,
        'sweep': sweep,
        'dihedral': dihedral,
        'area': area,
        'sect_alfa': sect_alfa,
        'sect_CL': sect_CL,
        'sect_CD': sect_CD,
        'sect_Cm': sect_Cm,
        'CLref': CLref,
        'sect_alfaL0': sect_alfaL0
        }

def get_file_name():
    if len(RA) > 1:
        save_RA = "change"
        cases = len(RA)
    else:
        save_RA = RA[0]
    if len(RT) > 1:
        save_RT = "change"
        cases = len(RT)
    else:
        save_RT = RT[0]
    if len(CL) > 1:
        save_CL = "change"
        cases = len(CL)
    else:
        save_CL = CL[0]   
    if len(hb) > 1:
        save_hb = "change"
        cases = len(hb)
    else:
        save_hb = hb[0]    
    if len(RA) == 1 and len(RT) == 1 and len(CL) == 1 and len(hb) == 1:
        cases = 1
    
    with open('input.json', 'r') as file:
        data = json.load(file)
    LLsolver = data["solver"]["type"]
    grid = data["wings"]["Main"]["grid"]

    file_name = "RA{}_RT{}_CL{}_hb{}_cp{}d{}t_init{}_LLsolver{}_didist{}_grid{}_optmethod{}_dibounds{}_opttol{}_optmaxiter{}_cases{}"
    file_name = file_name.format(save_RA, save_RT, save_CL, save_hb, d, t+1, init, LLsolver, didist, grid, opt_method, di_bounds, opt_tol, max_iter, cases)
    return file_name

def span_fraction():
    f = open('input.json', 'r+')
    data = f.readlines()
    f.close()
    ## get the span fraction locations
    grid = int(data[47][20:])
    theta = np.linspace(0, np.pi, grid+1)
    span_frac = (1 - np.cos(theta))/2
    count = 0
    mid_list = []
    for _ in range(len(theta)-1):
        new = (theta[count] + theta[count+1])/2
        mid_list.append(new)
        count+=1 
    theta_mid = np.array(mid_list)
    span_frac = (1 - np.cos(theta_mid))/2
    span_frac = np.append(span_frac, 1.0)
    return span_frac
    
def runCase(Ra, Rt, CLtarget, hb, CLmax=1.4, cw=1.):
    ## case setup
    #####################
    get_CL_offset.CLtarget = CLtarget
    get_sect_CL_max.CLmax = CLmax
    get_min_z.tol = 0.01*Ra
    
    cr, ct = 2.*cw/(1.+Rt), 2.*Rt*cw/(1.+Rt)
    
    h = hb * Ra
    
    get_CD.iter = 0
    
    f = open('input.json', 'r+')
    data = f.readlines()
    f.seek(0)
    data[14] = ' '*8  + '"area": {:.16f},\n'.format(Ra)
    data[16] = ' '*8  + '"lateral_length": {:.16f}\n'.format(Ra)
    data[20] = ' '*8  + '"ground": {:.16f}\n'.format(h)
    data[11] = ' '*8  + '"CGz": {:.16f}\n'.format(-h)
    data[34] = ' '*16 + '"dz": {:.16f},\n'.format(-h)
    data[37] = ' '*12 + '"span": {:.16f},\n'.format(Ra/2)
    data[44] = ' '*12 + '"root_chord": {:.16f},\n'.format(cr)
    data[45] = ' '*12 + '"tip_chord": {:.16f},\n'.format(ct)
    f.writelines(data)
    f.close()

if __name__ == '__main__':
    #baseline wing: RA=8.0, RT=1.0, CL=0.5, hb=0.25
    RA = [8.0]
    RT = [1.0]
    CL = [0.5]
    #hb = [0.25]#np.linspace(0.1, 1.0, 100) #if h/b not changing, make sure to change for loop
    d = 4
    t = 4 
    # ma = 6

    # opt_tol = 1e-12
    # opt_method = "SLSQP"
    # max_iter = 1000
    
    # didist = "quadratic" #linear, quadratic
    # di_bounds = (-90, 90)
    # init = "previous" #options: 0, "previous"
    #plt.plot(y_b, z_b)
    #plt.plot(y_b, dist["dihedral"])
    #plt.plot(y_b, dist["twist"])
    #plt.plot([-0.5, 0.5], [0, 0], linestyle="--")
    
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
    x_label = r"$\mathregular{z/b}$"
    
    fig, cdi = plt.subplots(figsize=plot_size)
    cdi.set_xlabel(x_label, fontsize=plot_fontsize, fontname=plot_font, fontstyle='italic')
    cdi.set_ylabel(r"$\mathregular{h/b}$", fontsize=plot_fontsize, fontname=plot_font, fontstyle='italic')
    cdi.tick_params(axis='both', labelsize=tick_size, width=tick_width, length=tick_length)
    cdi.spines['top'].set_linewidth(axis_thickness)  
    cdi.spines['right'].set_linewidth(axis_thickness)  
    cdi.spines['bottom'].set_linewidth(axis_thickness)  
    cdi.spines['left'].set_linewidth(axis_thickness)
    cdi.grid(grid_on, linewidth=grid_width)
    cdi.minorticks_off()
    # plt.gca().invert_yaxis()
    plt.axis('equal')
    # cdi.set_xlim(-0.55, 0.55)
    
    cdi.set_xlim(-0.5, 0.5)
    xticks = np.linspace(-0.5, 0.5, num=6)
    cdi.set_xticks(xticks)
    
    # cdi.set_ylim(-0.5, 8)
    # y_ticks = np.linspace(-0.5, 8, num=6)
    # cdi.set_yticks(y_ticks)

    #file_name = all_txt_files[0]
    
    file_name = "use_this.txt"#RA8.0_RT1.0_CL0.5_hbchange_cp4d5t_initprevious_LLsolvernonlinear_didistquadratic_grid100_optmethodSLSQP_dibounds(-90, 90)_opttol1e-12_optmaxiter1000_cases100_full.txt"
    #for file_name in all_txt_files:
    f = open(file_name, 'r+')
    opt_data = f.readlines()
    f.close() 
    
    linestyle_list = ["-", "--", "-."]
    num_lines = 24
    start = 15
    end = 23+1
    
    n = int(len(opt_data)/num_lines)
    cases = list(range(0, n))
    count=0
    for i in cases:
        coeff_index = num_lines*i + 2
        x = opt_data[num_lines*i+start:num_lines*i+end]
        
        coeff_string = opt_data[coeff_index]
        x = [float(val) for val in x]#this is the x that contains the twist and dihedral info
        
        RA = float(coeff_string.split("Ra: ")[1].split(",    Rt")[0])
        RT = float(coeff_string.split("Rt: ")[1].split(",    CLtarget")[0])
        CL = float(coeff_string.split("CLtarget: ")[1].split(",    h/b")[0])
        hb = float(coeff_string.split("h/b: ")[1].split("h/b: ")[0])
        
        legend = r"$\mathregular{C_{L}=" + f"{CL:.1f}" + "}$"
        span_frac = span_fraction()  
        #RA, RT, CL, h/b
        runCase(RA, RT, CL, hb)
        MachUp(x)
        dist = MachUpDist()

        
        #x = [ma] + [0]*t + [0]*d
        #file_name = get_file_name()
        span_frac = span_fraction()
        
        MachUp(x)
        runCase(RA, RT, CL, hb)
        
        dist = MachUpDist()
        f = open('input.json', 'r+')
        b_data = f.readlines()
        f.close()
        b = 2*float(b_data[37].split('"span": ')[1].split(",")[0])

        y_b = [y/b for y in dist["y"]]
        z_b = [-z/b for z in dist["z"]]
        

        cdi.plot(y_b, z_b, color="black", linewidth=plot_linewidth, linestyle=linestyle_list[count], label=legend)
        count+=1
    cdi.legend(fontsize=14)
