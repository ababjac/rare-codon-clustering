import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score



def make_scatterplot_data(X, Y, xlabel, ylabel, plt_title):
    plt.scatter(X, Y)
    plt.title(plt_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    idx = np.isfinite(X) & np.isfinite(Y)
    m,b = np.polyfit(X[idx], Y[idx], 1)
    Y_pred = m*X+b
    plt.plot(X, Y_pred, color='red')

    r2 = r2_score(Y[idx], Y_pred[idx])
    plt.annotate("r-squared = {:.3f}".format(r2), (min(X), max(Y)))
    plt.show()

df = pd.read_csv('../Data/Disorder/disorder_full_yeast.csv')

# make_scatterplot_data(df['phi_estimate'], df['disorder_size'], 'ML-Phi Estimate', 'Disorder', 'Phi Estimate vs Disorder for Yeast')
make_scatterplot_data(df['MM_avg'], df['disorder_size'], 'MM Average', 'Disorder', 'Average %MinMax vs Disorder for Yeast')
make_scatterplot_data(df['MM_min'], df['disorder_size'], 'MM Minimum', 'Disorder', 'Minimum %MinMax vs Disorder for Yeast')
make_scatterplot_data(df['MM_max'], df['disorder_size'], 'MM Maximum', 'Disorder', 'Maximum %MinMax vs Disorder for Yeast')
