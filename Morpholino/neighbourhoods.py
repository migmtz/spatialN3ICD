import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.ndimage import gaussian_filter
from main_functions import mark_voisins#, fast_multi_periodogram_window_3d, fast_multi_periodogram_prewindow_3d

sns.set_theme()

colors = {1.0:"#026e0d", -1.0:"#7efc8b", 0.0:"#bcbcbc"}

group_name = ["MOctrl", "MOj1b"]

if __name__ == "__main__":
    # Data reading
    df = pd.read_excel("data_summary.xlsx")

    # Parameters and axes initialisation
    r = 10.0
    q = 0.2
    mark_id = "NOTCH sum"

    fig, ax = plt.subplots(2, 3, figsize=(8, 4), sharex=True, sharey=True)

    fig2, ax2 = plt.subplots(2, 1, figsize=(6, 11), sharex=True)

    for control_id in [0, 1]:
        # Initialisation of variables group-specific (Ctrl or J1b)
        aux_low = 0.0
        aux_mid = 0.0
        aux_high = 0.0
        voisinages = []
        n = 0

        for fish_id in range(1, 4):
            for hemisphere_id in range(1, 3):
                # Data reading to separate disjoint segmented windows.

                df_poisson = df.query(
                    '`Fish ID`==@fish_id and `Hemisphere ID`==@hemisphere_id and `J1b ?`==@control_id')

                if not np.all(df_poisson['Group ID'].isna()):
                    print("Careful")
                    df_poisson1 = df.query(
                        '`Fish ID`==@fish_id and `Hemisphere ID`==@hemisphere_id and `J1b ?`==@control_id and `Group ID` == 0')
                    df_poisson2 = df.query(
                        '`Fish ID`==@fish_id and `Hemisphere ID`==@hemisphere_id and `J1b ?`==@control_id and `Group ID` == 1')
                    df_list = [df_poisson1, df_poisson2]
                else:
                    df_list = [df_poisson]

                for df_1 in df_list:
                    n_1 = len(df_1)

                    X = df_1['X'].to_numpy()
                    Y = df_1['Y'].to_numpy()

                    marks = df_1[mark_id].to_numpy()

                    # Compute lower and upper bounds for marks
                    lower_bound = np.quantile(marks, q)
                    upper_bound = np.quantile(marks, 1-q)

                    # Compute relative neighbourhood for all points along with classified groups according to marks.
                    voisins_sample = mark_voisins(X, Y, r, marks, lower_bound, upper_bound)

                    aux_low += np.mean(np.array(voisins_sample[0])[:, 2] < 0.0)
                    aux_mid += np.mean(np.array(voisins_sample[0])[:, 2] == 0.0)
                    aux_high += np.mean(np.array(voisins_sample[0])[:,2] > 0.0)

                    try:
                        voisinages[0]
                        voisinages = [np.vstack((vois_prec, voisi)) for vois_prec, voisi in zip(voisinages, voisins_sample)]
                    except IndexError:
                        voisinages = [np.array(voisi) for voisi in voisins_sample]

                    n += 1

        # Plotting

        ax2[control_id].set_xlabel(r"$\mathrm{\mu}$m")
        ax2[control_id].set_ylabel(r"$\mathrm{\mu}$m")

        ax2[control_id].scatter(voisinages[3][:, 0], voisinages[3][:, 1], facecolor=None, edgecolors="k",
                                s=np.abs(voisinages[3][:, 2]) * 12,
                                alpha=0.75)
        ax2[control_id].scatter(voisinages[3][:, 0], voisinages[3][:, 1], c=voisinages[3][:, 2],
                                  s=(voisinages[3][:, 2] == 0) * 10,
                                  cmap="viridis", alpha=0.25)

        sc = ax2[control_id].scatter(voisinages[3][:, 0], voisinages[3][:, 1], c=voisinages[3][:, 2], s=np.abs(voisinages[3][:, 2])*10,
                                       cmap="viridis", alpha=0.75)

        ax2[control_id].set_title(group_name[control_id])

        for i in range(3):
            c = voisinages[i][:, 2]
            sc = ax[control_id, i].scatter(voisinages[i][:, 0], voisinages[i][:, 1], s=np.abs(voisinages[i][:, 2])*10, c=c, cmap="viridis", alpha=0.9)
            ax[control_id, i].scatter(voisinages[i][:, 0], voisinages[i][:, 1], s=(voisinages[i][:, 2] == 0)*10, c=c, cmap="viridis", alpha=0.25)
            ax[control_id, i].set_title(str(np.round(np.mean(voisinages[i][:, 2] < 0.0), 2)) + ",  " + str(np.round(np.mean(voisinages[i][:, 2] > 0.0), 2)) + ", n=" + str(len(c)))


    fig2.savefig("images/morpholino_neighbourhoods.pdf", format="pdf")#, bbox_inches="tight")
    plt.show()