import csv
import numpy as np 
import pickle   
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import interpolate

import angelino_nozzle_design
import MOC

### NOTE: should return here for quick csv to np pickle code
def serialize_csv(filename):
    # read csv data
    with open(filename,'r') as dest_f:
        data_iter = csv.reader(dest_f, delimiter = '    ',  quotechar = '"')
        data = [data for data in data_iter]
    # serialize csv data
    with open(filename[:-4] + '.pickle','wb') as handle:
        pickle.dump(data, handle, protocol = pickle.HIGHEST_PROTOCOL)

def read_pickle(filename):
    with open(filename, 'rb') as handle:
        unserialized_data = pickle.load(handle)

    return unserialized_data


def strpickle_to_numpypickle(filename):
    data = read_pickle(filename)
    data = np.asarray(data)

    cols = data[0,:]


    cfd_dict = {}

    for i in range(len(cols)):
        data_column = [float(el) for el in data[1:,i]]
        cfd_dict[cols[i]] = data_column
    
    with open(filename[:-7]+'_np.pickle','wb') as handle:
        pickle.dump(cfd_dict,handle,protocol = pickle.HIGHEST_PROTOCOL)


def set_font_size(ax1):
    for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] + ax1.get_xticklabels() + ax1.get_yticklabels()):
        item.set_fontsize(40)

def contourf_plot(x,y,z,cax,plug,n=1000,method='linear'):

    x_grid = np.linspace(min(x),max(x),n)
    y_grid = np.linspace(min(y),0.06,n)
    x_grid,y_grid = np.meshgrid(x_grid,y_grid)

    z_grid = interpolate.griddata((x,y),z,(x_grid,y_grid),method = method)

    ### PLUG SPECIFIC
    def line(x,m,x1,y1): return m*(x-x1) + y1

    x_vals = np.linspace(plug.x[-1]+0.0001,x_grid.max(),10)
    y_vals = line(x_vals,-80,x_vals[0],plug.y[-1])
    x_interp = np.concatenate((plug.x,x_vals))
    y_interp = np.concatenate((plug.y,y_vals))


    tck = interpolate.splrep(x_interp,y_interp)
    #z_grid[y_grid < interpolate.splev(x_grid,tck)] = np.nan
    z_grid[y_grid < np.interp(x_grid,x_interp,y_interp)] = np.nan
    z_grid[x_grid < plug.x[0]] = np.nan
    # z_grid[z_grid < 0.1] = np.nan
    ###


    #z_grid[]
    levels = [el for el in np.linspace(0,3.4,15)]

    print(levels)
    fill = cax.contourf(x_grid,y_grid,z_grid,levels,extend = 'max', cmap = cm.jet)
    
    
    return fill
# data = read_pickle('asd0.pickle')
# cols = data[0];
# cfd_mesh = np.asarray(float(data[1:]))

if __name__ == '__main__':
    # ideal plug definition
    plug1 = angelino_nozzle_design.design_test_nozzle()

    # MOC for plug
    moc_mesh = MOC.chr_mesh(plug1,plug1.gamma,0,120,downstream_factor=3.0)
    moc_mesh.y = moc_mesh.y
    # read CFD mesh
    cfd_mesh = read_pickle('asd0_np.pickle')
    print(cfd_mesh.keys())

    fig, (ax1,ax2) = plt.subplots(2,1)
    plug1.plot_contour(ax1)
    fill1 = contourf_plot(cfd_mesh['Points:1'],cfd_mesh['Points:0'],cfd_mesh['Ma'],ax1,plug1)

    interp_x,idx = np.unique(moc_mesh.x[moc_mesh.ID_contour_chr],return_index = 1)
    interp_y = moc_mesh.y[moc_mesh.ID_contour_chr]
    interp_y = interp_y[idx]
    # ax2.plot(interp_x,interp_y*-1,'.')
    # interp_y,idx = np.unique(interp_y,return_index=1)
    # interp_x = interp_x[idx]

    idx = np.argsort(interp_x)

    interp_x = interp_x[idx]
    interp_y = interp_y[idx]

    tck = interpolate.splrep(interp_x,interp_y)
    
    plug1.plot_contour(ax2)

    X_plt = np.linspace(0,moc_mesh.x.max(),1000)
    Y_plt = np.linspace(0,moc_mesh.y.min(),1000)
    X_plt,Y_plt = np.meshgrid(X_plt,Y_plt)
    
    invalid_grid = interpolate.splev(X_plt.flatten(),tck)<Y_plt.flatten()
    invalid_grid = invalid_grid.reshape(X_plt.shape)
    levels = [el for el in np.linspace(0,3.4,15)]
    M_contour=interpolate.griddata((moc_mesh.x,moc_mesh.y),moc_mesh.M,(X_plt,Y_plt),method='linear')
   
    #ax1.plot(X_plt[valid_grid],Y_plt[invalid_grid],'.')

 
    M_contour[invalid_grid] = np.nan
    # ax1.plot(X_plt.flatten(),interpolate.splev(X_plt.flatten(),tck),'.')

    levels = [el for el in np.linspace(0,3.4,15)]
    M_fill = ax2.contourf(X_plt,-Y_plt,M_contour,levels,extend = 'max',cmap=cm.jet)
    # plt.colorbar(M_fill,ax=ax2,orientation = 'horizontal')

    ax1.set_xlim(ax2.get_xlim())
    ax1.set_ylim(ax2.get_ylim())
    # fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.1, 0.05, 0.8])
    set_font_size(cbar_ax)
    cbar = fig.colorbar(M_fill, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=20)
    # ax1.scatter(cfd_mesh['Points:1'],cfd_mesh['Points:0'],c=cfd_mesh['Ma'],cmap=cm.jet)
    # fig.colorbar(M_fill,ax=ax1)
    fig.set_size_inches(18.5,10.5)
    # Sfig.set_title('CFD vs. MOC')

    ax1.axes.get_xaxis().set_ticklabels([])
    ax1.axes.get_yaxis().set_ticklabels([])

    ax2.axes.get_xaxis().set_ticklabels([])
    ax2.axes.get_yaxis().set_ticklabels([])

    ax1.set_aspect('equal','box')
    ax2.set_aspect('equal','box')

    plt.savefig('CFDvMOC',dpi=100)
    plt.show() 
