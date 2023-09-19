import matplotlib.pyplot as plt
import scipy.optimize as opt
import numpy as np

def Graph(filename,Num_mols):
    data = []
    dic = {} 
    with open(filename, 'r') as f:
        for i in range(0,Num_mols-1):      #Number Molecules -1
            data = eval(f.readline())
            for k,v in data.items():
                try:
                    dic[k] += v
                except KeyError:
                    dic[k] = v
        norm = dic[1]
        for k,v in dic.items():
            dic[k] = v/norm
    x = np.array(list(dic.keys()))*0.1
    y = np.array(list(dic.values()))
    return x, y

def model_func(x,a,k1,b,k2,c,k3,k4):
    return a * np.exp(-1*x/k1) + b * np.exp(-1*x/k2) + c * np.exp(-1*x/k3) + (1-a-b-c) * np.exp(-1*x/k4) # * np.exp(-1*x/k5)

def Fitting(x,y):
    x = np.array(x)
    y = np.array(y)

    optimizedParameters, pcov = opt.curve_fit(model_func, x, y,p0=[1,1,1,1,1,1,0])
    residuals = y - model_func(x, *optimizedParameters)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y-np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    print('R-squared ' + str(r_squared))
  
    val = ['a', 'k1', 'b', 'k2', 'k3']
    lab = [(round(i,4)) for i in optimizedParameters]
    print(lab)
    fit = 'Fit: ' + '\n' + str(r'${c1}$'.replace('c1',str(lab[0]))
        + r'$e^{-x/{s1}}$'.replace('s1',str(lab[1])) 
        + r'$\cdot{c2}$'.replace('c2',str(lab[2]))
        + r'$e^{-x/{s2}}$'.replace('s2',str(lab[3]))
        + r'$\cdot{c3}$'.replace('c3',str(lab[4]))
        + r'$e^{-x/{s3}}$'.replace('s3',str(lab[5]))
        + r'$\cdot{c4}$'.replace('c4',str(round(1-lab[0]-lab[2]-lab[4],4))))
        #+ r'$e^{s4}$'.replace('s4',str(lab[6])))
    
    #r'*exp(-x/' + str(lab[1]) + ' + ' + str(lab[2]) + '*exp(-x/' + str(lab[3]) + ' + ' + str(lab[4]))
    plt.plot(x, model_func(x, *optimizedParameters),marker=',',label=fit) #+ '\n' + r'e^-kt'))
    plt.legend()
    plt.xlabel('Time, ps')
    plt.ylabel('Correlation')
    plt.title('Time Correlation of Hbonds')
    val = ['a', 'k1', 'b', 'k2', 'k3']
    print(val)
    [print(round(i,4)) for i in optimizedParameters]
    plt.show()

x,y = Graph('Corr.txt',192)
plt.plot(x,y, marker=',', label="Simulated Hbonds")
#plt.show()
Fitting(x,y)


