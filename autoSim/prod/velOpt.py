import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize
from scipy.linalg import issymmetric
from scipy import signal
from scipy.ndimage import gaussian_filter
import sys

import psutil
import tracemalloc
import os

#####
def prep_window(length=50, sig=5):
    '''

      Apply window smoothing function

    '''
    winNorm = 1 / (2*np.pi*sig**2)**(1/2)
    window = signal.windows.gaussian(length, std=sig)
    return window, winNorm

#####
def enforce_sym(dots, trim):
    '''

      Enforces symmetry and real values 

    '''
    dim = dots[0].shape
    print(f'Dim: {dim}')
    for i in range(0,dim[0]):
        #for j in range(i, i+1):
        for j in range(0, dim[1]):
            #print(i,j)
            ct = dots[:,i,j]
            ct_w = np.fft.fft(ct)
            ct_w = ct_w#*window
            #plt.plot(ct_w)
            #plt.show()
            #plt.clf()
            ct_w = ct_w*np.conj(ct_w)
            ct_w = np.fft.ifft(ct_w)
            ct_w = ct_w.real #imaginary components gone (noise)
            dots[:,i,j] = ct_w

    #print(len(dots))
    for i, ct in enumerate(dots):
        dots[i] = (ct+ct.T)/2
        if issymmetric(dots[i]) == False:
            print(f'ERROR: Matrix at t({i}) not symmetric')
            sys.exit(1)

    #if len(dots) > 3*trim:
    #    dots = dots[0:trim]
    
    dots = np.array(list(np.flip(dots, axis=0)) + list(dots))
    #print(len(dots))
    #plt.plot(np.arange(0,len(dots[:,5,5])), dots[:,5,5])
    #plt.plot(np.arange(0,len(dots[:,4,4])), dots[:,4,4])
    #plt.show()
    return dots

#####
def calc_freq_mats(dots, dt=1, trim=75):
    '''

      Convert the velocity data to frequency domain
      Holds information for all cross correlations

    '''
    #print('Starting enforce_sym')
    #dots = enforce_sym(dots, -1*trim)
    
    dots = np.append(np.flip(dots, axis=0), dots, axis=0)

    #print('End enforce_sym')

    #print('Starting fft')
    
    ct_w = np.fft.fft(dots)
    dots = ct_w.real

    #print('End fft')
    
    freq = np.fft.fftfreq(len(dots), d=dt)

    return dots, freq

#####
def calc_spectra(w, freq_mat, ax=''):
    '''
      
      Calculate spectra from vacf
      Given by the trace of the produced freq mat

    '''
    
    #y = [np.trace(i[0]) for i in freq_mat]
    #y = freq_mat[0][0]

    #print(f'Spectra sum -inf to inf: {sum(y)}')

    y = np.sqrt(np.array(y)**2)
    
    #data, = ax.plot(w, gaussian_filter(y, sigma=2),lw=1,)
    #x = list(data.get_xdata())
    #y = list(data.get_ydata())
    
    y = gaussian_filter(y, sigma=2,)
    #y = list(y)
    #x = list(w)
    x = w

    #zero = x.index(0)
    #print(f'Spectra sum 0 to inf: {sum(y[zero:])}')
    
    return x, y

#####
def main(vels, max_dt=100, max_t=100, jump_t=100, dt=0.5e-15):
    raw_spectra = {}
    
    #print('Starting freq_mats')
    freq_mat, w = calc_freq_mats(vels, dt=dt, trim=75)
    #print('End freq_mats')

    #print('Start Calc Spectra')
    y = np.sqrt(freq_mat**2)
    y = gaussian_filter(y, sigma=2,) 
    #print('End Calc Spectra')

    return dict(zip(w,y))

#####
def main_vel(data='solute.data', max_dt=10000, max_t=10000, jump_t=1, dt=0.5e-15, root=''): 

    print(os.getcwd())

    with open(f'{root}/prod/ip.txt','r') as f:
        ip = f.readline().strip()
    
    with open(f'{root}/prod/port.txt','r') as f:
        port = f.readline().strip()

    print(ip, port)

    u = mda.Universe(data, f'imd://{ip}:{port}',)#buffer_size=100*1024**2)
   
    #ag = u.select_atoms('type 1 or type 2 or type 3 or type 4')
    ag = u.atoms

    spectra = {}

    print(len(ag))
    i = 0
    c = 0
    Ct = []

    #tracemalloc.start()
    #snap1 = tracemalloc.take_snapshot()

    n_atoms = len(u.atoms)
    dot = np.empty((n_atoms,n_atoms))
    upper_indices = np.triu_indices_from(dot)
    proc = psutil.Process(os.getpid())

    for ts in u.trajectory:
        if i%5 == 0: print(ts)
        for atom in ag.atoms:
            atom.velocity = atom.velocity * atom.mass**(1/2)
       
        #snap2 = tracemalloc.take_snapshot()
        #top_stats = snap2.compare_to(snap1, 'lineno')
        #print(f"Step {i}: Top memory growth:")
        #for stat in top_stats[:5]:  # top 5 lines causing growth
        #    print(stat)
        #snap1 = snap2  # reset snapshot for next step

        if i == 0: 
            v0 = ag.velocities

        #dot = np.matmul(v0, ag.velocities.T, out=dot)
        np.matmul(v0, ag.velocities.T, out=dot)
        
        #upper_tri = dot[upper_indices]
        diag = np.diagonal(dot)

        try:
            Ct[i] = (Ct[i] + diag) / 2
        except IndexError:
            Ct.append(diag.astype(np.float32))

        if i == max_dt:
            my_dot_prods = np.array(Ct) 
            print('Starting Spectra...')
            for atom_ind in range(0, n_atoms):
                vels = my_dot_prods[:,atom_ind]
    
                spectrum = main(vels, max_dt=max_dt, max_t=max_t, jump_t=jump_t, dt=dt)

                for k,v in spectrum.items():
                    try:
                        spectra[k] = (spectra[k]+spectrum[k]) / 2
                    except KeyError:
                        spectra[k] = spectrum[k]
            i = 0
            print('...End spectra')
            mem = proc.memory_info().rss / 1024**2
            print(f"{ts}: {mem:.2f} MB")
            c += 1

            with open('spectra.txt','w') as f:
                for k,v in spectra.items():
                    f.write(f'{k} {v}\n')

            #if c == 1: break
        #if i == max_dt: break

        i += 1

    with open('spectra.txt','w') as f:
        for k,v in spectra.items():
            f.write(f'{k} {v}\n')
        
#####

if __name__ == '__main__':
    import json
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--root",)
    args = parser.parse_args()

    settings = f'{args.root}/settings.json'

    with open(settings,'r') as j:
        settings = json.load(j)

    main_vel(
            max_dt=settings['spectra']['max_dt'], 
            max_t=settings['spectra']['max_t'], 
            jump_t=1, 
            dt=0.5e-15,
            root=args.root,
        )

