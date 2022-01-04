# copied from `https://github.com/theislab/single-cell-tutorial/blob/master/latest_notebook/gprofiler_plotting.py`
# Plotting functions - 'GProfiler-official version'
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib import colors
from matplotlib import rcParams


def scale_data_5_75(data):
    mind = np.min(data)
    maxd = np.max(data)
    
    if maxd == mind:
        maxd=maxd+1
        mind=mind-1
        
    drange = maxd - mind
    return ((((data - mind)/drange*0.70)+0.05)*100)


def plot_enrich(data, n_terms=20, title=None, save=False, 
                fontsize=30, dpi=300, width=4, height=None): # 20 terms ~ default figsize
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sb
    from matplotlib import colors
    from matplotlib import rcParams

    # Test data input
    if not isinstance(data, pd.DataFrame):
        raise ValueError('Please input a Pandas Dataframe output by gprofiler.')
        
    if not np.all([term in data.columns for term in ['p_value', 'name', 'intersection_size']]):
        raise TypeError('The data frame {} does not contain enrichment results from gprofiler.'.format(data))
        
    data_to_plot = data.iloc[:n_terms,:].copy()    
    data_to_plot['go.id'] = data_to_plot.index
    

    min_pval = data_to_plot['p_value'].min()
    max_pval = data_to_plot['p_value'].max()
    
    # Scale intersection_size to be between 5 and 75 for plotting
    #Note: this is done as calibration was done for values between 5 and 75
    data_to_plot['scaled.overlap'] = scale_data_5_75(data_to_plot['intersection_size'])
    
    norm = colors.LogNorm(min_pval, max_pval)
    sm = plt.cm.ScalarMappable(cmap="coolwarm", norm=norm)
    sm.set_array([])

    sb.set(style="whitegrid")
    height = data_to_plot.shape[0]*0.4
    rcParams.update({'font.size': fontsize,'figure.figsize':[width,height]})

    path = plt.figure()
    plt.scatter(x='recall', y="name", c='p_value', cmap='coolwarm', 
                       norm=colors.LogNorm(min_pval, max_pval), 
                       data=data_to_plot, linewidth=1, edgecolor="grey", 
                       s=[(i+10)**1.5 for i in data_to_plot['scaled.overlap']])
    
    ax = plt.gca()
    
    ax.invert_yaxis()

    ax.set_ylabel('')
    ax.set_xlabel('Gene ratio', fontsize=14, fontweight='normal')
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)
    if title is not None:
        ax.set_title(title)

    # Shrink current axis by 20%
    box = ax.get_position()
    print('%d terms with height'%data_to_plot.shape[0], 'figsize', path.bbox_inches, 'box',box)
    
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Get tick marks for this plot
    #Note: 6 ticks maximum
    min_tick = np.floor(np.log10(min_pval)).astype(int)
    max_tick = np.ceil(np.log10(max_pval)).astype(int)
    tick_step = np.ceil((max_tick - min_tick)/6).astype(int)
    
    # Ensure no 0 values
    if tick_step == 0:
        tick_step = 1
        min_tick = max_tick-1
    
    ticks_vals = [10**i for i in range(max_tick, min_tick-1, -tick_step)]
    ticks_labs = ['$10^{'+str(i)+'}$' for i in range(max_tick, min_tick-1, -tick_step)]

    #Colorbar
    fig = plt.gcf()
    cbaxes = fig.add_axes([0.8, 0.1, 0.03, 0.3])
    cbar = ax.figure.colorbar(sm, ticks=ticks_vals, shrink=0.5, anchor=(0,0.1), cax=cbaxes)
    cbar.ax.set_yticklabels(ticks_labs)
    cbar.set_label("Adjusted p-value", fontsize=10, fontweight='normal')

    #Size legend
    min_olap = data_to_plot['intersection_size'].min()
    max_olap = data_to_plot['intersection_size'].max()
    olap_range = max_olap - min_olap
    
    #Note: approximate scaled 5, 25, 50, 75 values are calculated
    #      and then rounded to nearest number divisible by 5
    size_leg_vals = [np.round(i/5)*5 for i in 
                          [min_olap, min_olap+(20/70)*olap_range, min_olap+(45/70)*olap_range, max_olap]]
    size_leg_scaled_vals = scale_data_5_75(size_leg_vals)

    
    l1 = plt.scatter([],[], s=(size_leg_scaled_vals[0]+10)**1, edgecolors='none', color='black')
    l2 = plt.scatter([],[], s=(size_leg_scaled_vals[1]+10)**1, edgecolors='none', color='black')
    l3 = plt.scatter([],[], s=(size_leg_scaled_vals[2]+10)**1, edgecolors='none', color='black')
    l4 = plt.scatter([],[], s=(size_leg_scaled_vals[3]+10)**1, edgecolors='none', color='black')

    labels = [str(int(i)) for i in size_leg_vals]

    plt.legend([l1, l2, l3, l4], labels, ncol=1, frameon=False, title_fontsize=10, fontsize=8,
                     handlelength=1, loc = 'upper left', borderpad = 1, labelspacing = .6,
                     handletextpad=2, title='Gene overlap', scatterpoints = 1,  bbox_to_anchor=(-2, 8.2/3), 
                     facecolor='black')

    if save:
        plt.savefig(save, dpi=300, format='png', bbox='tight')
        
    else:
        plt.show()
        return plt.gca()


def patchplot(pal, width=7, y_ratio=0.8, labels=None):
    """Plot the values in a color palette as a horizontal array.

    Parameters
    ----------
    pal : sequence of matplotlib colors
        colors, i.e. as returned by seaborn.color_palette()
    size :
        scaling factor for size of plot

    """
    # set style
    # sns.set(palette='deep',style='white',font_scale=2)
    
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(width, y_ratio*width/n))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) + .05)
    ax.set_yticks([-.5, .5])
    
    ax.set_xticklabels(labels) if labels is not None else ax.set_xticklabels([])
    ax.set_yticklabels([])
    # set back
    # sns.set(context='notebook', style='whitegrid', palette='deep', font='sans-serif', font_scale=1.2, color_codes=True)
    return ax



