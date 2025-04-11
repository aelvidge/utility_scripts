def regression_plot_and_statistics(x, y, xlabel="x", ylabel="y", fontsize=8, do_plot=1, annotation_location="NW", **kwargs):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.metrics import mean_absolute_error as MAE
    from sklearn.metrics import mean_squared_error as MSE

    if 'xlim' in kwargs:
        min_x = kwargs["xlim"][0]
        max_x = kwargs["xlim"][1]
    else:
        min_x = np.nanmin(x)
        max_x = np.nanmax(x)
    if 'ylim' in kwargs:
        min_y = kwargs["ylim"][0]
        max_y = kwargs["ylim"][1]
    else:
        min_y = np.nanmin(y)
        max_y = np.nanmax(y)
    mean_x = np.nanmean(x)
    mean_y = np.nanmean(y)
    std_x = np.nanstd(x)
    std_y = np.nanstd(y)
    
    idx_isfin = [i for i,_ in enumerate(x) if np.isfinite(x[i]) and np.isfinite(y[i])]
    x_isfin = x[idx_isfin]
    y_isfin = y[idx_isfin]
    
    slope, intercept = np.polyfit(x_isfin, y_isfin, 1)
    correlation_coefficient = np.corrcoef(x, y)[0][1]
    mean_bias = mean_y-mean_x
    mean_absolute_error = MAE(x_isfin, y_isfin)
    root_mean_squared_error = MSE(x_isfin, y_isfin, squared=False)

    stats_text = ("mean x, y = {:.2f}".format(mean_x) + ", {:.2f}".format(mean_y) + 
                  "\nstandard deviation x, y = {:.2f}".format(std_x) + ", {:.2f}".format(std_y) + 
                  "\nslope, intercept of regression = {:.2f}".format(slope) + ", {:.2f}".format(intercept) + 
                  "\ncorrelation coefficient = {:.2f}".format(correlation_coefficient) + 
                  "\nmean bias error of y = {:.2f}".format(mean_bias) + 
                  "\nmean absolute error of y = {:.2f}".format(mean_absolute_error) + 
                  "\nroot mean squared error of y = {:.2f}".format(root_mean_squared_error))
    
    if do_plot==1:
        fig, ax = plt.subplots()
        plt.scatter(x,y,s=0.1)
        x_lin = np.linspace(min_x,max_x,100)
        plt.plot(x_lin, slope*x_lin+intercept, color = 'r')
        xyl = [np.min([min_x,min_y]),np.max([max_x,max_y])]
        plt.plot(xyl,xyl, color = 'k', ls = '--')
        xl = [min_x,max_x]
        yl = [min_y,max_y]
        plt.xlim(xl)
        plt.ylim(yl)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if annotation_location=="NW":
            plt.annotate(stats_text, xy=(xl[0]+(xl[1]-xl[0])/100, yl[1]-(yl[1]-yl[0])/100), horizontalalignment='left',verticalalignment='top', fontsize=fontsize)
        elif annotation_location=="SE":
            plt.annotate(stats_text, xy=(xl[1]-(xl[1]-xl[0])/100, yl[0]+(yl[1]-yl[0])/100), horizontalalignment='right',verticalalignment='bottom', fontsize=fontsize)
        else:
            raise ValueError("annotation_location must be equal to either ""NW"" or ""SE""")
        ax.set_aspect('equal')
    else:
        print(stats_text)

    return  mean_x, mean_y, std_x, std_y, slope, intercept, correlation_coefficient, mean_bias, mean_absolute_error, root_mean_squared_error

def weibull(ws):
    import numpy as np
    import math
    import scipy.special as sc
    k = (math.sqrt(np.mean(abs(ws - np.mean(ws))**2))/np.mean(ws))**-1.089
    gamma_f = math.exp(sc.gammaln(1+(1/k)))
    A = (np.mean(ws)/gamma_f)
    ws_bins = np.linspace(np.min(ws),np.max(ws),1000)
    dist = (k / A) * (ws_bins / A)**(k - 1) * np.exp(-(ws_bins / A)**k)
    return [ws_bins, dist, A, k]