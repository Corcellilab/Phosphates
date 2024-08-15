import numpy as np

def timeCorr(self, min_dt=0, max_dt=20, skip_dt=1, t_prop=1):
  bonded = {}
  count = {}
  for i in range(0, self.length, t_prop):
    step_bonded = {}
    total = getattr(getattr(self, f'result{i}'), 'Network')
    step_bonded = [tuple(sorted(ID)) for ID in total] 
    for pair in step_bonded:
      try:
        bonded[pair] = bonded[pair] + 1
      except KeyError:
        bonded[pair] = 1
    
    unbound = set(bonded.keys()).difference(step_bonded)
    for u in unbound:
      try: 
        count[bonded[u]] += 1
        del bonded[u]
      except KeyError: 
        count[bonded[u]] = 1
        del bonded[u] 
  counts = dict(sorted(count.items()))

  time_corr = {}
  count = list(counts.values())
  for i in range(min_dt, max_dt, skip_dt):
    time_corr[i] = sum(count[i:])

  y = np.array(list(time_corr.values()))
  y = y / time_corr[0]
  x = np.array(list(time_corr.keys()))

  return x, y

