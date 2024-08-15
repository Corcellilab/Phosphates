import numpy as np

def timeCorr(self, min_dt=1, max_dt=20, skip_dt=1, t_prop=1):
  for i in range(0, self.length):
    total = getattr(getattr(self, f'result{i}'), 'Aggregate resids')

