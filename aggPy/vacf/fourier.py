#fourier.py
"""

    Fourier transform 0 to +x values
    Return: Spectra in wavenumbers

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks

def ftransform(x, y, dt=1):
    
    #y = list(reversed(list(y))) + [0] + list(y)
    y = list(y)
    y = y[-2::-1] + y
    
    N = len(y)+1000000

    y = np.fft.fft(y, n=N)
    y = np.abs(y)
    y = np.fft.fftshift(y)
    
    x = np.fft.fftfreq(len(y), dt)
    x = np.fft.fftshift(x)
    x = x / 3e10

    print('Filtering...')
    for _ in range(3):
        y = gaussian_filter(y, 100)
    
    return x, y

