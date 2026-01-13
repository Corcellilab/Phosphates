#plot_csv.py

import pandas as pd
import matplotlib.pyplot as plt
import sys

def plot_csv():
    filename = input('Filename: ')
    results = pd.read_csv(filename)
    print(results.columns)
    
    while True:
        col_name = input('Columnn to analyze: ')
        if col_name.lower() == 'exit': sys.exit(1)
        if col_name.lower() == 'quit': sys.exit(1)
        if col_name.lower() == 'new': 
            filename = input('Filename: ')
            results = pd.read_csv(filename)
            print(results.columns)
       
        print(results[col_name])

        x = list()
        y = list()
        for i, value in enumerate(results[col_name]):
            x.append(i)
            v = eval(value)
            y.extend(v)

        graph_type = input('Scatter plot or Histogram [s/h]?: ')

        if graph_type == 's':
            plt.scatter(x,y)
            plt.show()
        elif graph_type == 'h':
            bin_num = int(input('Number of bins: '))
            plt.hist(y, bins=bin_num)
            plt.show()

if __name__ == "__main__":
    plot_csv()

