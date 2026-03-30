import numpy as np
import json

#####

def check_convg(filename, col, start=1/3, t_value=None):
    with open(filename,'r') as f:
        data = json.load(f)

    roll_avg = np.mean(data[f'rollAvg_{col}'][int(len(data[f'rollAvg_{col}'])*start):])
    stdev = data[f'avg_{col}'][1]

    if t_value == None:
        return None, roll_avg, stdev

    if t_value == 'avg':
        t_value = float(roll_avg)

    maximum = roll_avg + 2*stdev
    minimum = roll_avg - 2*stdev

    if t_value <= maximum and t_value >= minimum:
        return True, roll_avg, stdev
    else:
        return False, roll_avg, stdev


#####

if __name__ == '__main__':
    convg, avg, stdev = check_convg('defrost.json', col='Temp', start=1/3, t_value=150)
    print(convg)

