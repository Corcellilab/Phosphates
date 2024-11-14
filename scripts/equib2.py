import numpy as np
import matplotlib.pyplot as plt

class LAMMPS_Dump:
  def __init__(self, filename):
    self.title = filename
    #title = filename.split('/')
    #self.title = str(title[0] + ' ' + title[1].split('.')[0]) 
    with open(filename, mode='r') as f:
      lines = [line.rstrip('r\n').split() for line in f]
      lines = np.array(lines[1::])
    self.data = lines

  def get_data(self):   #Get variable types and data column - lammps specific
    variables = []
    for i in self.data[1]:
      if i == '=': break    #Lammps dummp deliminator
      variables.append(i)
    self.variables = variables

    length = len(variables)+1
    for i in range(0,len(variables)):
      setattr(self, variables[i], self.data[:, length+i].astype(float))

  def average(self, N):
    run_avgs = {}
    stdev = {}
    for i in self.variables:
      run_avgs[i] = np.convolve(getattr(self, i), np.ones(N)/N, mode='valid')
      stdev[i] = np.std(run_avgs[i])
    self.avgs = run_avgs
    self.stdevs = stdev

  def simple_avg(self, start):
    start = int(start)
    for i in range(1, len(self.variables)):
      name = self.variables[i]
      print(f'{name}: {np.average(getattr(self, name)[start::])} +/- {np.std(getattr(self, name)[start::])}')

def graphing(self, start, nrows, ncols):
    graphs = {}
    x = self.Step[start:]
    f, ax = plt.subplots(nrows, ncols, figsize=(15,15))
    c = 0
    r = 0
    for i in range(0, len(self.variables)):
      name = self.variables[i]
      if name == 'Step': continue           #Skip timestep v timestep plot
      y = getattr(self, name)[start:]
      min_y = np.min(y)*0.90
      max_y = np.max(y)*1.10
      ax[r, c].plot(x, y)
      ax[r, c].set_ylim(min_y, max_y)
      ax[r, c].set_title(f'{name}')
      c += 1
      if c >= ncols:
        c = 0
        r += 1
    f.suptitle(self.title, fontsize=25)
    plt.show()

########
def main(filename, start=0):
  start = int(start)
  data = LAMMPS_Dump(filename)
  data.get_data()
  data.Volume = data.Volume**(1/3)        #Convert volume to box side length

  data.average(1000)
  print(data.variables)
  #print(np.average(data.Pressure), np.std(data.Pressure))
  #print(data.avgs['Pressure'], data.stdevs['Pressure'])

  data.simple_avg(start)

  data.graphing(start, 2, 3)
