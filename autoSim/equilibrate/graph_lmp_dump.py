import matplotlib.pyplot as plt
import numpy as np

#####

def graph_lmp_dump(data, window_size=5, outfile='lmp_dump_graphs.png'):
    n_keys = len(data.keys())
    if n_keys % 2 != 0 :
        n_keys += 1
        hide_final_plot = True

    fig, ax = plt.subplots(2, int(n_keys/2), figsize=(10,8))
    plt.tight_layout()

    avgs = {}
    keys = list(data.keys())
    values = list(data.values())
    for i, axis in enumerate(ax.flat):
        try:
            axis.plot(values[i], label='Instantaneous', color='cornflowerblue')
            kernel = np.ones(window_size) / window_size
            avg = np.convolve(values[i], kernel, mode='valid')
            avgs[f'rollAvg_{keys[i]}'] = list(avg)

            axis.plot(avg, ls='-', alpha=0.5, label='Running Avg.', color='indianred')
            
            t_avg = np.mean(values[i])
            plot_t_avg = [t_avg for _ in range(0,len(values[i]))]
            
            axis.plot(plot_t_avg, ls='-', alpha=0.5, label='Avg.', color='green')
            axis.set_title(keys[i])

            avgs[f'avg_{keys[i]}'] = (np.mean(values[i]), np.std(values[i]))

        except IndexError:
            ax.flat[i-1].legend(bbox_to_anchor=(1.05,1), loc='upper right')
            fig.delaxes(axis)
        
    #plt.show()

    plt.savefig(outfile,)

    return avgs

#####

if __name__ == '__main__':
    from read_lmp_dump import read_lmp_dump
    import json

    filename = 'defrost.txt'
    prefix = filename.split(".")[0]

    data = read_lmp_dump(filename)
    avgs = graph_lmp_dump(data, window_size=10, outfile=f'{prefix}_avg.png')

    with open(f'{prefix}.json','w') as j:
        json.dump(avgs, j, indent=4)

