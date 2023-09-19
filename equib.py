#Equib reader
with open('prod.Conc1M.txt', mode='r') as dat:
    lines = (line.rstrip('r\n') for line in dat)
    lst = []
    p = 0
    t = 0
    s = []
    for line in lines:
        lst.append(line.split())
    lst.pop(0)
    for i in lst:
        s.append(i[0])
        p += float((i[-3]))
        t += float((i[-4]))
    print('Avg Pressure ' + str((p/len(s))))
    print('Avg Temp ' + str((t/len(s))))
