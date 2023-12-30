import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from atmopy.display import *

def draw_street_map(data, node, node_begin, node_end):

    width = 30 / 2.54
    #matplotlib.rcParams["figure.figsize"] = (0.9 * width, .5 * width)
    matplotlib.rcParams["figure.figsize"] = (width, width)
    matplotlib.rcParams["figure.subplot.left"] = 0.08 #0.05
    matplotlib.rcParams["figure.subplot.right"] = 0.75 # 0.8
    matplotlib.rcParams["figure.subplot.bottom"] = 0.08
    matplotlib.rcParams["figure.subplot.top"] = 0.95
    font_size = 11
    matplotlib.rcParams["font.size"] = font_size
    matplotlib.rcParams["axes.titlesize"] = font_size + 2
    matplotlib.rcParams["axes.labelsize"] = font_size +2
    matplotlib.rcParams["xtick.labelsize"] = font_size
    matplotlib.rcParams["ytick.labelsize"] = font_size
    matplotlib.rcParams["legend.fontsize"] = font_size
    matplotlib.rcParams["xtick.major.size"] = 2
    matplotlib.rcParams["xtick.minor.size"] = 1
    matplotlib.rcParams["ytick.major.size"] = 2
    matplotlib.rcParams["ytick.minor.size"] = 1

    
    
    vmax = max(data)

    # Display intervals
    ninterv = 15
    step = vmax / ninterv

    interv = np.arange(0,vmax,step)

    alpha = 1. # 0.5
    beta = 1.
    size = np.size(interv)
    interval = []
    interval_l = []

    # Interval for plotting
    for i in range(0,size):
        interval.append(alpha*interv[i])
        interval_l.append(alpha*interv[i]/beta)

    # Set the line thickness
    s = 800

    # PLOTTING
    # ---------

    w = 2
    w_0 = 1


    # Plot
    # ----------

    fig = plt.figure()
    fig.suptitle("Title")
    ax = plt.subplot(111)

    # Display the node numbers.
    # ------------------------
    disp_node_number = False  #True # False
    disp_street_number = False #True #False
    if (disp_node_number):
        for inode in node:
            ax.text(node[inode][0], node[inode][1], str(inode), size=10, color='b')



    lw_max = 2
    for i in range(len(node_begin)):
        xy_begin = node[node_begin[i]]
        xy_end = node[node_end[i]]


        # Display the street numbers.
        if (disp_street_number):
            x_ = (xy_begin[0] + xy_end[0]) / 2.0
            y_ = (xy_begin[1] + xy_end[1]) / 2.0
            ax.text(x_, y_, str(arc_id[i]), size=10, color='r')

    
        if data[i] != 1.0 and data[i] >= interval[0] and data[i] < interval[1]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'indigo',linestyle = '-',
                    lw = 0.5)

        elif data[i] >= interval[1] and data[i] < interval[2]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'midnightblue',linestyle = '-',
                    lw = 1.0) #data[i]/s)
        elif data[i] >= interval[2] and data[i] < interval[3]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'navy',linestyle = '-',
                    lw = 1.5) #data[i]/s)
        elif data[i] >= interval[3] and data[i] <interval[4] :
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'blue',linestyle = '-',
                    lw = 2.0) # data[i]/s)
        elif data[i] >= interval[4] and data[i] < interval[5]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'royalblue',linestyle = '-',
                    lw = 2.5) # data[i]/s)
        elif data[i] >= interval[5] and data[i] < interval[6]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'deepskyblue',linestyle = '-',
                    lw = 3.0) # data[i]/s)
        elif data[i] >= interval[6] and data[i] < interval[7]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'cyan',linestyle = '-',
                    lw = 3.5) #data[i]/s)
        elif data[i] >= interval[7] and data[i] < interval[8]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'lime',linestyle = '-',
                    lw = 4.0) #data[i]/s)
        elif data[i] >= interval[8] and data[i] < interval[9]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'yellow',linestyle = '-',
                    lw = 4.5) # data[i]/s)
        elif data[i] >= interval[9] and data[i] < interval[10]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'gold',linestyle = '-',
                    lw = 5.0) # data[i]/s)
        elif data[i] >= interval[10] and data[i] <interval[11]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'orange',linestyle = '-',
                    lw = 5.5) # data[i]/s)
        elif data[i] >= interval[11] and data[i] < interval[12]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'salmon',linestyle = '-',
                    lw = 6.0) #data[i]/s)
        elif data[i] >= interval[12] and data[i] < interval[13]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'red',linestyle = '-',
                    lw = 6.5) # data[i]/s)
        elif data[i] >= interval[13] and data[i] < interval[14]:
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color = 'firebrick',linestyle = '-',
                    lw = 7.0) #data[i]/s)
        else :
            ax.plot([xy_begin[0], xy_end[0]],
                    [xy_begin[1], xy_end[1]],
                    color='darkred',linestyle='-',
                    lw = 7.5) #data[i]/s)

    lw = 0.5
    intro = mlines.Line2D([],[],color = 'w',linestyle='-',lw = 0.5 )
    #intro = mlines.Line2D([],[],color = '#A4A4A4',linestyle='-',lw = 0.5 )
    in_0 = mlines.Line2D([],[],color = 'k',linestyle='-',lw = 0.5 )
    in_1 = mlines.Line2D([],[],color = 'indigo',linestyle='-',lw = lw * 1 )
    in_2 = mlines.Line2D([],[],color = 'midnightblue',linestyle='-',lw = lw * 2 )
    in_3 = mlines.Line2D([],[],color = 'navy',linestyle='-',lw = lw*3)
    in_4 = mlines.Line2D([],[],color = 'blue',linestyle='-',lw = lw *4)
    in_5 = mlines.Line2D([],[],color = 'royalblue',linestyle='-',lw = lw*5)
    in_6 = mlines.Line2D([],[],color = 'deepskyblue',linestyle='-',lw = lw*6)
    in_7 = mlines.Line2D([],[],color = 'cyan',linestyle='-',lw = lw*7)
    in_8 = mlines.Line2D([],[],color = 'lime',linestyle='-',lw = lw*8)
    in_9 = mlines.Line2D([],[],color = 'gold',linestyle='-',lw = lw*9)
    in_10 = mlines.Line2D([],[],color = 'orange',linestyle='-',lw = lw*10)
    in_11 = mlines.Line2D([],[],color = 'salmon',linestyle='-',lw = lw*11)
    in_12 = mlines.Line2D([],[],color = 'tomato',linestyle='-',lw = lw*12)
    in_13 = mlines.Line2D([],[],color = 'red',linestyle='-',lw = lw*13)
    in_14 = mlines.Line2D([],[],color = 'firebrick',linestyle='-',lw = lw *14)
    in_15 = mlines.Line2D([],[],color = 'darkred',linestyle='-',lw = lw*15)

    handles = [intro, in_0, in_1, in_2, in_3, in_4, in_5, in_6, in_7, in_8, in_9, in_10, in_11, in_12, in_13, in_14, in_15]


    interval_l = np.array(interval_l) 
    #interval_l /= 1000.
    # interval_l = interval_l_arr

    labels = ['in $\u03BC$g m$^{-1}$s$^{-1}$',
              '0.0',
              '[ '+str(interval_l[0])+', '+str(interval_l[1])+' [',
              '[ '+str(interval_l[1])+', '+str(interval_l[2])+' [',
              '[ '+str(interval_l[2])+', '+str(interval_l[3])+' [',
              '[ '+str(interval_l[3])+', '+str(interval_l[4])+' [',
              '[ '+str(interval_l[4])+', '+str(interval_l[5])+' [',
              '[ '+str(interval_l[5])+', '+str(interval_l[6])+' [',
              '[ '+str(interval_l[6])+', '+str(interval_l[7])+' [',
              '[ '+str(interval_l[7])+', '+str(interval_l[8])+' [',
              '[ '+str(interval_l[8])+', '+str(interval_l[9])+' [',
              '[ '+str(interval_l[9])+', '+str(interval_l[10])+' [',
              '[ '+str(interval_l[10])+', '+str(interval_l[11])+' [',
              '[ '+str(interval_l[11])+', '+str(interval_l[12])+' [',
              '[ '+str(interval_l[12])+', '+str(interval_l[13])+' [',
              '[ '+str(interval_l[13])+', '+str(interval_l[14])+' [',
              '[ '+str(interval_l[14])+', ~ [' ]

    # bbox_to_anchor(0.5,-0.1)
    plt.legend(handles,labels,bbox_to_anchor=(0.995,0.),loc = 'lower left', prop = {'size' : 20})
    # plt.subplots_adjust(right = 0.2)
    plt.subplots_adjust(bottom = 0.1)

    a = ax.get_yticks().tolist()
    ax.set_yticklabels(a)

