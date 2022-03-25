import pandas as pd
import struct
import numpy as np


def read_pos(data):
    """ Loads an APT .pos file as a pandas dataframe.

    Columns:
        x: Reconstructed x position
        y: Reconstructed y position
        z: Reconstructed z position
        Da: mass/charge ratio of ion"""
    # read in the data
    n = int(len(data)/4)
    d = struct.unpack('>'+'f'*n,data)
                    # '>' denotes 'big-endian' byte order
    # unpack data
    pos = pd.DataFrame({'x': d[0::4],
                        'y': d[1::4],
                        'z': d[2::4],
                        'Da': d[3::4]})
    return pos


def read_epos(data):
    """Loads an APT .epos file as a pandas dataframe.

    Columns:
        x: Reconstructed x position
        y: Reconstructed y position
        z: Reconstructed z position
        Da: Mass/charge ratio of ion
        ns: Ion Time Of Flight
        DC_kV: Potential
        pulse_kV: Size of voltage pulse (voltage pulsing mode only)
        det_x: Detector x position
        det_y: Detector y position
        pslep: Pulses since last event pulse (i.e. ionisation rate)
        ipp: Ions per pulse (multihits)

     [x,y,z,Da,ns,DC_kV,pulse_kV,det_x,det_y,pslep,ipp].
        pslep = pulses since last event pulse
        ipp = ions per pulse

    When more than one ion is recorded for a given pulse, only the
    first event will have an entry in the "Pulses since last evenT
    pulse" column. Each subsequent event for that pulse will have
    an entry of zero because no additional pulser firings occurred
    before that event was recorded. Likewise, the "Ions Per Pulse"
    column will contain the total number of recorded ion events for
    a given pulse. This is normally one, but for a sequence of records
    a pulse with multiply recorded ions, the first ion record will
    have the total number of ions measured in that pulse, while the
    remaining records for that pulse will have 0 for the Ions Per
    Pulse value.
        ~ Appendix A of 'Atom Probe tomography: A Users Guide',
          notes on ePOS format."""
    # read in the data
    n = int(len(data)/4)
    rs = n / 11
    d = struct.unpack('>'+'fffffffffII'*int(rs),data)
                    # '>' denotes 'big-endian' byte order
    # unpack data
    pos = pd.DataFrame({'x': d[0::11],
                        'y': d[1::11],
                        'z': d[2::11],
                        'Da': d[3::11],
                        'ns': d[4::11],
                        'DC_kV': d[5::11],
                        'pulse_kV': d[6::11],
                        'det_x': d[7::11],
                        'det_y': d[8::11],
                        'pslep': d[9::11], # pulses since last event pulse
                        'ipp': d[10::11]}) # ions per pulse
    return pos


def read_rrng(f):
    """Loads a .rrng file produced by IVAS. Returns two dataframes of 'ions'
    and 'ranges'."""
    import re

    rf = open(f,'r').readlines()

    patterns = re.compile(r'Ion([0-9]+)=([A-Za-z0-9]+).*|Range([0-9]+)=(\d+.\d+) +(\d+.\d+) +Vol:(\d+.\d+) +([A-Za-z:0-9 ]+) +Color:([A-Z0-9]{6})')
    patterns2 = re.compile(r'ion([0-9]+)=([A-Za-z0-9]+).*|range([0-9]+)=(\d+.\d+) +(\d+.\d+) +vol:(\d+.\d+) +([A-Za-z:0-9 ]+) +color:([A-Z0-9]{6})')
    ions = []
    rrngs = []
    for line in rf:
        m = patterns.search(line)
        m2 = patterns2.search(line)
        if m:
            if m.groups()[0] is not None:
                ions.append(m.groups()[:2])
            else:
                rrngs.append(m.groups()[2:])
        if m2:
            if m2.groups()[0] is not None:
                ions.append(m2.groups()[:2])
            else:
                rrngs.append(m2.groups()[2:])
    ions = pd.DataFrame(ions, columns=['number','name'])
    ions.set_index('number',inplace=True)
    rrngs = pd.DataFrame(rrngs, columns=['number','lower','upper','vol','comp','colour'])
    rrngs.set_index('number',inplace=True)

    rrngs[['lower','upper','vol']] = rrngs[['lower','upper','vol']].astype(float)
    rrngs[['comp','colour']] = rrngs[['comp','colour']].astype(str)

    return ions,rrngs

def read_rng(f):
    import re
    rn = open(f,'r').readlines()
    start=rn[0]
    length=str(start).replace('\n','')
    length=length.split(' ')
    length=length[1]
    pattern_rng=re.compile('----------------- ')
    pattern_rng2=re.compile(' . ')
    ions = []
    rngs = []
    for line in rn:
        mr = pattern_rng.search(line)
        mr2=pattern_rng2.search(line)
        if mr is not None:
            ions.append(line)
        if mr2 is not None:
             rngs.append(line)
    ions2=str(ions).replace('\\n',"")
    ions2=ions2.replace("']","")
    ions2=str(ions2).split(' ')   
    
    kick=[]
    m=1E999
    for k in range (0,len(rn)):
        if [rn[k]]==ions:m=k
        if k>m and k<m+int(length)+1:kick.append(rn[k])

    rngs=kick
    colors = ['FF0000','0000FF','00FF00']
    color_index = 0 
    liste=[]
    for i in range(0,len(rngs)):
        comps=''
        number=i+1
        rngs1=str(rngs[i]).replace('\n',"")
        rngs2=rngs1.split(' ')
        for j in range(0, len(ions2)-1):
            if int(rngs2[4+2*j])!=0:
               comp=ions2[j+1]
               comps=comps+' '+comp+':'+rngs2[4+2*j]
        colour=colors[color_index]
        liste.append([number,float(rngs2[1]),float(rngs2[2]),0.0,comps,colour])
        color_index +=1
        if color_index>=7: color_index=0
    
    db_ions=[]
    for m in range(1,len(ions2)):
        db_ions.append([m,ions2[m]])

    ions_db = pd.DataFrame(db_ions, columns=['number','name'])
    ions_db.set_index('number',inplace=True)
    rngs_db = pd.DataFrame(liste, columns=['number','lower','upper','vol','comp','colour'])
    rngs_db.set_index('number',inplace=True)
    
    return ions_db,rngs_db
    
    
    

def label_ions(pos,rrngs):
    """labels ions in a .pos or .epos dataframe (anything with a 'Da' column)
    with composition and colour, based on an imported .rrng file."""
    rl=list(rrngs)
    rd=pd.DataFrame(rl.pop(1),columns=['lower','upper','vol','comp','colour'])
    pos['comp'] = ''
    pos['colour'] = '#FFFFFF'

    for n,r in rd.iterrows():
        pos.loc[(pos.Da >= r.lower) & (pos.Da <= r.upper),['comp','colour']] = [r['comp'],'#' + r['colour']]

    return pos


def deconvolve(lpos):
    """Takes a composition-labelled pos file, and deconvolves
    the complex ions. Produces a dataframe of the same input format
    with the extra columns:
       'element': element name
       'n': stoichiometry
    For complex ions, the location of the different components is not
    altered - i.e. xyz position will be the same for several elements."""

    import re

    out = []
    pattern = re.compile(r'([A-Za-z]+):([0-9]+)')

    for g,d in lpos.groupby('comp'):
        if g != '':
            for i in range(len(g.split(' '))):
                tmp = d.copy()
                cn = pattern.search(g.split(' ')[i]).groups()
                tmp['element'] = cn[0]
                tmp['n'] = cn[1]
                out.append(tmp.copy())
    return pd.concat(out)

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    
def unique(list1):
 
    # initialize a null list
    unique_list = []
     
    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    return unique_list
