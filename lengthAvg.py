import matplotlib.pyplot as plt

with open('IonPair431.txt','r') as f:
    branches = eval(f.readline())

y = {}
for k,v in branches.items():
    y_temp = []
    y_temp = [(key*(value)) for key,value in v.items()]
    #for key,value in v.items():
    #    print(key,value)
    y_temp = sum(y_temp)/k
    y[k] = y_temp

x = list(y.keys())
y = list(y.values())
avg_y = round(sum(y)/len(y),4)

plt.scatter(x,y, label='Average Ratio = ' + str(avg_y), color='orange',marker='.')
plt.xlabel('# Phosphates in Aggregate')
plt.ylabel('Phosphates in Aggregate /' + '\n' + 'Average # Na Interactions')
plt.title('Ratio of Aggregate Size/Average amount of Na Interaction' + '\n' + 'vs. Phos in Aggregate')
plt.legend()
plt.show()

