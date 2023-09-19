import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from colour import Color
import matplotlib.colors
from itertools import islice

def graphing(datas, length, xlab='', ylab='', title='', mark='o',leg_title='',leg_labels=1,y_parms=(0,100)):
    data = dict(sorted(datas.items()))
    data = dict(islice(data.items(),length[0],length[1]))
    HSV_tuples = []
    leg = []
    counter = 0
    [leg.extend(list(v.keys())) for v in data.values()]
    leg = set(leg)
    
    color_code = []
    for k,v in data.items():
        ends = dict(sorted(v.items()))
        data[k] = ends
        color_code.extend(v.keys())
        
    red = Color('blue')
    color_code = set(color_code)
    colors = list(red.range_to(Color("red"),len(color_code)))
    color = [(matplotlib.colors.to_rgb(str(v))) for v in colors]
    colors = dict(zip((list(color_code)),color))
    
    for k,v in data.items():
        color_map = list(set(v.keys()) & set(colors.keys()))
        color_list = [colors[x] for x in color_map]
        k = [k] * len(v.values())
        plt.scatter(x=k,y=v.values(),c=color_list,marker=mark)
        if counter > len(data):
            break
        else:
            counter += 1
     
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.ylim(y_parms)
    handle = []
    for x in leg:
        if x % leg_labels == 0 or x == 0 or x == 1:
            try:
                handle.append(mlines.Line2D([], [], color=colors[x], marker=mark, ls='', label=str(x)))
            except:
                continue
    if list(leg)[-1] % leg_labels != 0:
        handle.append(mlines.Line2D([], [], color=colors[x], marker=mark, ls='', label=str(list(leg)[-1])))
    plt.legend(handles=handle,title=leg_title)
    plt.show()
    return None

with open('IonPair315.txt','r') as f:
    chain_len = eval(f.readline())
    cordination = eval(f.readline())

graph = graphing(
        cordination, [0,30], y_parms=(0,1),
        xlab='Interacton Number', ylab = '% of Na', title='P-Na interactions vs. P-P interactions',
        leg_title='Na interactions given P IN',
        leg_labels=1,
        mark='.'
        )

chain = graphing(
        chain_len, [0,80], y_parms=(0,1),
        xlab='# Phos in aggregate', ylab = '% of Na', title='Aggregate Length vs. P-Na interactions',
        leg_title='Na interactions given Agg length',
        leg_labels=10,
        mark='.'
        )





