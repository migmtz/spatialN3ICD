import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from main_functions import mark_voisins

sns.set_theme()

colors = {1.0:"#026e0d", -1.0:"#7efc8b", 0.0:"#bcbcbc"}

group_name = ["BAC"]

if __name__ == "__main__":
    # Data reading
    df = pd.read_excel("aires_bac.xlsx")

    # Parameters and axes initialisation
    r = 10.0
    q = 0.2
    mark_id = "N3ICD Sum"

    fig, ax = plt.subplots(1, 3, figsize=(8, 2.5), sharex=True, sharey=True)

    fig2, ax2 = plt.subplots(1, 1, figsize=(6, 5), sharex=True)

    aux_low = 0.0
    aux_mid = 0.0
    aux_high = 0.0
    voisinages = []
    n = 0

    X = df['X'].to_numpy()
    Y = df['Y'].to_numpy()

    marks = df[mark_id].to_numpy()

    # Compute lower and upper bounds for marks
    low_bound = np.quantile(marks, q)
    upper_bound = np.quantile(marks, 1 - q)

    # Compute relative neighbourhood for all points along with classified groups according to marks.
    voisins_sample = mark_voisins(X, Y, r, marks, low_bound, upper_bound)

    aux_low += np.mean(np.array(voisins_sample[0])[:, 2] < 0.0)
    aux_mid += np.mean(np.array(voisins_sample[0])[:, 2] == 0.0)
    aux_high += np.mean(np.array(voisins_sample[0])[:, 2] > 0.0)

    try:
        voisinages[0]
        voisinages_individuels = [np.array(voisi) for voisi in voisins_sample]
        voisinages = [np.vstack((vois_prec, voisi)) for vois_prec, voisi in zip(voisinages, voisins_sample)]
    except IndexError:
        voisinages_individuels = [np.array(voisi) for voisi in voisins_sample]
        voisinages = [np.array(voisi) for voisi in voisins_sample]


    # Plotting

    ax2.set_xlabel(r"$\mathrm{\mu}$m")
    ax2.set_ylabel(r"$\mathrm{\mu}$m")

    ax2.scatter(voisinages[3][:, 0], voisinages[3][:, 1], facecolor=None, edgecolors="k",
                            s=np.abs(voisinages[3][:, 2]) * 12,
                            alpha=0.75)
    ax2.scatter(voisinages[3][:, 0], voisinages[3][:, 1], c=voisinages[3][:, 2],
                              s=(voisinages[3][:, 2] == 0) * 10,
                              cmap="viridis", alpha=0.25)

    sc = ax2.scatter(voisinages[3][:, 0], voisinages[3][:, 1], c=voisinages[3][:, 2], s=np.abs(voisinages[3][:, 2])*10,
                                   cmap="viridis", alpha=0.75)

    ax2.set_title(group_name[0])

    for i in range(3):
        c = voisinages[i][:, 2]
        sc = ax[i].scatter(voisinages[i][:, 0], voisinages[i][:, 1], s=np.abs(voisinages[i][:, 2])*10, c=c, cmap="viridis", alpha=0.9)
        ax[i].scatter(voisinages[i][:, 0], voisinages[i][:, 1], s=(voisinages[i][:, 2] == 0)*10, c=c, cmap="viridis", alpha=0.25)
        ax[i].set_title(str(np.round(np.mean(voisinages[i][:, 2] < 0.0), 2)) + ",  " + str(np.round(np.mean(voisinages[i][:, 2] > 0.0), 2)) + ", n=" + str(len(c)))

    fig2.savefig("images/bac_neighbourhoods.pdf", format="pdf")#, bbox_inches="tight")
    plt.show()