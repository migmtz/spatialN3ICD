import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from main_functions import *

sns.set_theme(style="white")


if __name__ == "__main__":
    # Data reading
    df0 = pd.read_excel("aires_bac.xlsx", sheet_name=0)

    df = [df0]
    
    # Maximal radius to plot and proportion for creating groups
    r = 10.0
    q = 0.2
    
    # Reference mark for x-axis and marks to plot wrt reference
    mark_id = "Apical Area"
    mark_list = ["N3ICD mean", "N3ICD Sum"]

    aux_low = 0.0
    aux_mid = 0.0
    aux_high = 0.0
    voisinages = []
    n = 0
    r_list = [10.0]

    fig, ax = plt.subplots(2, 1, figsize=(7, 8), layout="tight")

    for group_id, df_1 in enumerate(df):
        marks = df_1[mark_id].to_numpy() 
        X, Y = df_1["X"].to_numpy(), df_1["Y"].to_numpy()

        # Determine window of study using average radius of nuclei
        areas = df_1['Apical Area']
        radius = np.sqrt(areas / np.pi)
        mean_radius = radius.mean()

        min_x, max_x = np.min(X), np.max(X)
        min_y, max_y = np.min(Y), np.max(Y)

        len_x, len_y = max_x - min_x + 2 * mean_radius, max_y - min_y + 2 * mean_radius
        for i, y_mark_id in enumerate(mark_list):
            # Marks for classification in groups.
            y_marks = df_1[y_mark_id].to_numpy()

            # Compute lower and upper bounds for marks
            lower_bound = np.quantile(y_marks, q)
            upper_bound = np.quantile(y_marks, 1 - q)

            # Compute signal inside circle of radius r for each individual cell
            signal, _ = individual_marked_K_func(X, Y, r_list, len_x, len_y, lower_bound, upper_bound, y_marks, y_marks, m_r2_func)

            # Colours according to class (y_marks)
            colors = [1.0 if mark > upper_bound else -1.0 if mark < lower_bound else 0.0 for mark in y_marks]

            # Indices to separate outliers (points with no surrounding signal)
            positive_id = signal > 0.0
            null_id = signal == 0.0

            # Plotting of points coloured by surrounding signal
            ax[i].scatter(marks[null_id], y_marks[null_id], c="k", edgecolors="w", linewidths=0.5, zorder=10)
            sc = ax[i].scatter(marks[positive_id], y_marks[positive_id], c=np.abs(signal[positive_id]), cmap="Reds_r", linewidths=0.5,
                                         vmax=np.sort(signal)[-1],
                                         edgecolors="k", zorder=11)
            cb = plt.colorbar(sc, ax=ax[i], label="Total surrounding signal")

            # Plot of groups according to lower and upper groups
            plot_with_confidence([np.min(marks), np.max(marks)], [np.min(y_marks), np.min(y_marks)],
                                 [lower_bound, lower_bound],
                                 ax=ax[i], c="#440154", alpha=0.4, coef_alpha=0.0)
            plot_with_confidence([np.min(marks), np.max(marks)], [upper_bound, upper_bound],
                                 [np.max(y_marks), np.max(y_marks)],
                                 ax=ax[i], c="#fde725", alpha=0.4, coef_alpha=0.0)

            plot_with_confidence([np.min(marks), np.max(marks)], [lower_bound, lower_bound],
                                 [upper_bound, upper_bound],
                                 ax=ax[i], c="#21918c", alpha=0.4, coef_alpha=0.0)

            ax[1].set_xlabel(mark_id)

    fig.suptitle("BAC")

    fig.savefig("images/bac_pairplot.pdf", format="pdf", bbox_inches="tight")

    plt.show()
