import json
import matplotlib.pyplot as plt
import sys

def regraph(graph, plot_type, x, y, title='', xlabel='', ylabel='', num_bins=10, color='blue', label='series0'):
    fig = plt.figure()
    if plot_type == 'hist': fig = plt.hist(y, bins=num_bins, color=color, label=label)
    else: fig = graph(x, y, color=color, label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    return fig

####

def graphEdit(graph, plot_type, x, y, num_bins=20):
    if plot_type == 'hist': G = graph(y, bins=x)
    else: G = graph(x,y)
    title, xlabel, ylabel, label = ('',)*4
    color = 'blue'
    label = 'series0'

    while True:
        edit = input('Edit Property (list to show options): ')
        match edit:
            case 'list':
                print('title, xlabel, ylabel, bins, show, Exit')
            case 'title':
                title = (input('Title: '))
            case 'xlabel':
                xlabel = (input('xlabel: '))
            case 'ylabel':
                ylabel = (input('ylabel: '))
            case 'label':
                label = input('Series Label: ')
            case 'color':
                color = input('Color: ')
            case 'bins':
                num_bins = int(input('Num bins: '))
            case 'show':
                if plot_type == 'hist': 
                    G = regraph(graph, plot_type, x, y, title=title, xlabel=xlabel, ylabel=ylabel,
                        num_bins=num_bins, color=color, label=label,)
                else:
                    G = regraph(graph, plot_type, x, y, title=title, xlabel=xlabel, ylabel=ylabel,
                            color=color, label=label,)
                plt.show()
            case 'Exit':
                break
    return None
######
def plotting(plot_type, data):
    graph = getattr(plt, plot_type)      #plt.scatter/hist
    bin_num = 10

    if plot_type == 'scatter':
        keys = list(data.keys())
        keys.sort()
        sorted_data = {i: data[i] for i in keys}
        x, y = list(sorted_data.keys()),list(sorted_data.values())
        graph(x,y)
        plt.show() 
    elif plot_type == 'hist':
        bin_num = int(input('Number of bins: '))
        try:
            graph(data.values(), bins=bin_num)
            x, y = bin_num, data.values()
        except AttributeError:
            graph(data, bins=bin_num)
            x, y = bin_num, data
        plt.show()
    
    edit = input('Edit Graph? [y/n]: ')
    if edit == 'y':
        G = graphEdit(graph, plot_type, x, y, num_bins=bin_num)
        plt.show()

######
def jsonAnalyze():
    f = open(input('Filename: '))
    data = json.load(f)
    
    while True:
        print(f'Keys available: {list(data.keys())}')
        selection = input('Select data key (or Exit): ')
        if selection == 'Quit' or selection == 'Exit':
            break
        show = input('Display Data? [y/n]: ')
        if show == 'y': print(data[selection])
        
        plot_type = input('Plot type (scatter or hist): ')
        if plot_type == 'scatter': plotting(plot_type, data[selection])
        elif plot_type == 'hist': plotting(plot_type, data[selection])


#jsonAnalyze()

