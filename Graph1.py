import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from colour import Color
import matplotlib.colors
from itertools import islice

def graphing(datas, length, xlab='', ylab='', title='', mark='o',leg_title='',leg_labels=1,y_parms=(0,100), rang=[False]):
    data = dict(sorted(datas.items()))
    data = dict(islice(data.items(),length[0],length[1]))
    HSV_tuples = []
    leg = []
    counter = 0
    
    if rang[0] == True:
        rang_min = rang[1]
        rang_max = rang[2]
        rang_step = rang[3]
        rang = rang[0]
        lst = []

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
    plt.figure(figsize=(8,4.2)) 
    for k,v in data.items():
        color_map = list(set(v.keys()) & set(colors.keys()))
        color_list = [colors[x] for x in color_map]
        if rang == True:
            for i in range(rang_min,rang_max,rang_step):
                try:
                    plt.scatter(x=k,y=v[i],c=colors[i],marker=mark)
                    lst.append(i)
                except KeyError:
                    pass
        else:
            k = [k] * len(v.values())
            plt.scatter(x=k,y=v.values(),c=color_list,marker=mark)
        if counter > len(data):
            break
        else:
            counter += 1
    if rang == True:
        lst = sorted(set(lst))
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        #plt.title(title)
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
        #handle = [(mlines.Line2D([], [], color=colors[x], marker=mark, ls='', label=str(x))) for x in lst]
        plt.legend(handles=handle,title=leg_title)
        plt.savefig('aggSize.png',dpi=2000)
        plt.show()
    else:
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
        plt.legend(handles=handle,title=leg_title,loc='upper right')
        #plt.figure(figsize=(5,5))
        plt.savefig('aggSize.png',dpi=1200)
        plt.show()
        
    return None

with open('branches.txt','r') as f:
    #chain_len = eval(f.readline())
    #cordination = eval(f.readline())
    branches = eval(f.readline())

#graph = graphing(
#        cordination, [0,8], y_parms=(0,1),
#        rang=[True,2,80,3],
#        xlab='Interacton Number', ylab = '% of Na', title='P-Na interactions vs. P-P interactions at 0-3.15 angstrom',
#        leg_title='Na interactions given P IN',
#        leg_labels=1,
#        mark='o'
#        )

#chain = graphing(
#        chain_len, [0,80], y_parms=(0,0.4),
#        rang=[True,2,200,15],
#        xlab='# Phos in aggregate', ylab = '% of Na', title='Aggregate Length vs. P-Na interactions at 3.15-4.31 angstrom',
#        leg_title='Na interactions given Agg length',
#        leg_labels=15,
#        mark='.'
#        )

branch = graphing(
        branches, [0,80], y_parms=(0,100),
        rang=[True,2,80,3],
        xlab='Number of Molecules in Aggregate', ylab = '% of Ends', title='Ratio of Ends in Aggregate vs. Size',
        leg_title='# of Ends Given Size',
        leg_labels=5,
        mark='.'
        )




