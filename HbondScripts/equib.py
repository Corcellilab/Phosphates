#Equib reader
import matplotlib.pyplot as plt

with open('equib_npt.txt', mode='r') as dat:
    lines = (line.rstrip('r\n') for line in dat)
    lst = []
    p_avg = 0
    t_avg = 0
    v_avg = 0
    s = []
    for line in lines:
        lst.append(line.split())
    lst.pop(0)
    p = []
    t = []
    v = []
    for i in lst:
        s.append(float(i[6]))
        p.append(float((i[-3])))
        t.append(float((i[-4])))
        v.append(float((i[-1])))
        p_avg += float((i[-3]))
        t_avg += float((i[-4]))
        v_avg += float((i[-1]))
    print('Avg Pressure ' + str((p_avg/len(s))))
    print('Avg Volume ' + str((v_avg/len(s))))
    print('Box side ' + str(((v_avg/len(s))**(1/3))))
    print('Avg Temp ' + str((t_avg/len(s))))


plt.scatter(s,p,marker='.',color='green')
plt.title('Pressure')
plt.show()
plt.clf()
plt.scatter(s,t,marker='.',color='red')
plt.title('Temperature')
plt.show()
plt.clf()
plt.scatter(s,v,marker='.',color='orange')
plt.title('Volume')
plt.show()
plt.clf()

