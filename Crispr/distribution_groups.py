import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.ndimage import gaussian_filter
from main_functions import mark_voisins#, fast_multi_periodogram_window_3d, fast_multi_periodogram_prewindow_3d

sns.set_theme()

group_name = ["Hemisphere 1", "Hemisphere 2"]


if __name__ == "__main__":
    # Data reading
    df0 = pd.read_excel("aires_crispr.xlsx", sheet_name=0)
    df1 = pd.read_excel("aires_crispr.xlsx", sheet_name=1)

    df = [df0, df1]

    # Parameters and axes initialisation
    r = 10.0
    q = 0.2
    mark_id = "N3ICD Sum"

    aux_low = 0.0
    aux_mid = 0.0
    aux_high = 0.0
    voisinages = []
    n = 0

    fig, ax = plt.subplots(1, 2, figsize=(4, 5.0), sharey=False)

    for group_id, df_1 in enumerate(df):
        marks = df_1[mark_id].to_numpy()

        # Compute lower and upper bounds for marks
        lower_bound = np.quantile(marks, q)
        upper_bound = np.quantile(marks, 1 - q)

        # Colours according to groups
        colors = [1.0 if mark > upper_bound else -1.0 if mark < lower_bound else 0.0 for mark in marks]

        # Plotting
        sns.swarmplot(marks, ax=ax[group_id], c=colors, cmap="viridis", edgecolor="k", linewidth=1.0)
        x0, x1 = ax[group_id].get_xlim()

        ax[group_id].plot([x0+0.1, x1-0.1], [upper_bound, upper_bound], c="k", alpha=0.5, linestyle="--")
        ax[group_id].plot([x0+0.1, x1-0.1], [lower_bound, lower_bound], c="k", alpha=0.5, linestyle="--")

        ax[group_id].set_ylabel("N3ICD Sum")

    fig.suptitle("Crispr")

    fig.savefig("images/crispr_groups_n3icd_sum.pdf", format="pdf", bbox_inches="tight")

    plt.show()