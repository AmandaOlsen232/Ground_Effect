import shutil
import json, os
import numpy as np
import scipy
import ZachsModules as zm
import sys
# zm.zp.updateRCParams(**{'text.usetex':False})
# plt = zm.plt

import matplotlib.pyplot as plt
import MyModules as my
plot_settings = my.MyPlot()

def MachUp():
    ##################
    ## run MachUp
    os.system('MachUp.out input.json > out.txt')
MachUp.iter=0

def get_CD(runMU=True, jac=True):
    if runMU: MachUp()
    ##################
    ## get CD
    f = open('input_forces.json', 'r')
    data = json.load(f)
    f.close()
    if jac:
        return data['total']['MyAirplane']['stability_coordinates']['CD'] * 10000
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
    f.truncate()
    f.close()





    


    
    
if __name__ == '__main__':

    #baseline wing: RA=8.0, RT=1.0, CL=0.5, hb=0.25
    RA = 8.0
    RT = 1.0
    CL = 0.5
    hb = 0.25 #if h/b not changing, make sure to change for loop
    d = 4
    t = 4 
    ma = 0
    nodes = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]

    CD_list = [0.00942497704042313, 0.0094395875856149, 0.00944301910251925, 0.00944426598994268, 0.00944485074568149, 0.00944517051405375, 0.00944536389325673, 0.00944548966741566, 0.00944557601693095, 0.00944563784143986, 0.00944568361626641, 0.00944571844977246, 0.00944574556916941, 0.00944576709425836, 0.00944578446390254, 0.00944579868253867, 0.00944581046854241, 0.00944582034668554, 0.0094458287075241, 0.00944583584657321]
    # runCase(RA, RT, CL, hb) #set h/b, RT, CL, and RA
    CD_500 = 0.00944589137741308
    # # CD_list = []
    # for _ in range(1):
    #     with open("input.json", "r+") as f:
    #         data = f.readlines()
    #         f.seek(0)
    #         data[47] = ' '*12 + '"grid": {:.16f}\n'.format(500)
    #         f.writelines(data)
    #         f.truncate()
    #         f.close()
    #     MachUp()
    #     CD = get_CD()/10000
    #     # CD_list.append(CD)
    #     print(500, CD)
        
    plt.plot(nodes, CD_list)
    plt.xlabel(r"$N_{\nu}$")
    plt.ylabel(r"$C_{Di}$")
    plt.xlim(10, 200) 
    plt.ylim(0.009425, 0.009447)
    y_ticks = np.linspace(0.009425, 0.009447, num=6)
    plt.yticks(y_ticks)
    #cdi.set_xticks([0.125, 0.25, 0.5, 0.75, 1.0])
    # xticks = np.linspace(4, 14, num=6)
    # plt.set_xticks(xticks)