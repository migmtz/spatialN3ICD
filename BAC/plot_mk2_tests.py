import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from main_functions import plot_with_confidence
import pickle


sns.set_theme()

# ['ID Imaris',                                     0
#  'N3ICD mean', 'N3ICD mediane', 'N3ICD Sum',      1-3
#  'PCNA', 'RFP',                                   4-5
#  'Sox2 mediane', 'nlsRFP mediane',                6-7
#  'Volume',                                        8
#  'Apical Area',                                   9
#  'X', 'Y', 'Z']                                   10-12

toSaveKeysSpatial = ['group_id',
                     'max_r', 'divisions',
                     'r_list', 'K_func', 'sq_mean_intensity']

toCheckKeysSpatial = ['group_id', 'max_r', 'divisions']

toSaveKeys = ['group_id', 'n_permut', 'q',
              'mark_id', 'max_r', 'divisions',
              'r_list', 'marked_K', 'mean_mark']

toCheckKeys = ['group_id', 'n_permut', 'q', 'mark_id', 'max_r', 'divisions']


if __name__ == "__main__":
    # Data reading (for averaging)
    df0 = pd.read_excel("aires_bac.xlsx", sheet_name=0)

    df = [df0]

    # Number of permutations
    n_permut = 199
    # List of marks to plot
    mark_id_list = ['N3ICD mean', 'N3ICD Sum', "Apical Area"]
    # Maximal radius and partitioning
    max_r = 30
    divisions = 201
    # Proportion for upper and lower groups
    q = 0.2

    # Mark-interaction function to consider
    func_name_list = ["func_inf"]
    mK_list = [[], [], []]

    colors = ['k', 'gray']

    fig, ax = plt.subplots(3, 1, figsize=(5, 9), sharex=True, sharey=False)

    for id_mark, mark_id in enumerate(mark_id_list):
        print(mark_id)
        K_func_mean = np.zeros((divisions))
        sq_mean_intensity_mean = np.zeros((1))

        M_func_mean = np.zeros((len(func_name_list), 1 + n_permut, divisions))
        mean_mark_mean = np.zeros((len(func_name_list), 1, 1))

        n = 0

        for group_id, df_1 in enumerate(df):
            print(" " * 5, "Hemisphere: ", group_id)
            n_1 = len(df_1)

            # Reading saved spatial estimations
            toCompare = [globals()[s] for s in toCheckKeysSpatial]

            fileSpatial = 'saved_estimations/spatial_K.pkl'

            with open(fileSpatial, 'rb') as f:
                saved_mK = pickle.load(f)

            saved_dic = {}

            for dic in saved_mK:
                dic_values = [dic[key] for key in toCheckKeysSpatial]
                if dic_values == toCompare:
                    saved_dic = dic

            r_list = saved_dic['r_list']
            K_func = saved_dic['K_func']
            sq_mean_intensity = saved_dic['sq_mean_intensity']

            # Averaging
            K_func_mean += K_func * n_1
            sq_mean_intensity_mean += sq_mean_intensity * n_1
            n += n_1

            for j, m_name in enumerate(func_name_list):
                # Reading saved marked estimations
                fileMarked = 'saved_estimations/permut_' + m_name + '_mK_pool.pkl'
                toCompare = [globals()[s] for s in toCheckKeys]

                with open(fileMarked, 'rb') as f:
                    saved_mK = pickle.load(f)

                saved_dic = {}

                for dic in saved_mK:
                    dic_values = [dic[key] for key in toCheckKeys]
                    if dic_values == toCompare:
                        saved_dic = dic
                M_func_list = saved_dic['marked_K']
                mean_mark = saved_dic['mean_mark']

                M_func_mean[j, :, :] += M_func_list
                mean_mark_mean[j, 0, 0] += mean_mark

        for j, m_name in enumerate(func_name_list):
            # Normalisation by points density and expected mark value
            M_aux = M_func_mean[j, :, :] / (mean_mark_mean[j, 0, 0] * (sq_mean_intensity_mean / n))
            mK_list[id_mark] += [M_aux[0, :]]

            # Compute variance-stabilised L-function
            L_aux = np.sqrt(M_aux / np.pi)
            # Extract observed L-function
            observed_Lfunc = L_aux[0, :]

            # Extract 100 permuted L-functions for estimated expected reference function.
            estimated_randomL = L_aux[1:101, :].mean(axis=0)

            # Compute deviations wrt estimated expected L-function using remaining permutations
            marked_L_envelope = np.abs(L_aux[101:, :] - estimated_randomL[np.newaxis, :])

            # Compute p-value
            u_0 = np.max(np.abs(observed_Lfunc - estimated_randomL))
            u_i = np.max(marked_L_envelope, axis=1)
            p_value = 1 - np.sum(u_i < u_0) / ((n_permut + 1) // 2)
            print(" "*10, "p-value = ", np.round(p_value, 3))

            # Compute maximal deviation and plot envelopes
            D_max = np.max(marked_L_envelope)

            min_plot, max_plot = estimated_randomL - D_max, estimated_randomL + D_max

            ax[id_mark].plot(r_list, observed_Lfunc - estimated_randomL, linestyle='-', c=colors[0], # c="#440154",
                                    alpha=1.0)
            ax[id_mark].plot(r_list, D_max * np.ones(r_list.shape), linestyle='--', c=colors[1],
                                    alpha=1.0)
            ax[id_mark].plot(r_list, - D_max * np.ones(r_list.shape), linestyle='--', c=colors[1],
                                    alpha=1.0)

            ##### Rejection zones
            upper_rejected =  (max_plot >= observed_Lfunc)
            under_rejected = (min_plot <= observed_Lfunc)

            upper_r_test = np.ones(r_list.shape) * (D_max)
            under_r_test = np.ones(r_list.shape) * (- D_max)
            upper_r_test[upper_rejected] = np.nan
            under_r_test[under_rejected] = np.nan

            up_curve_aux = (observed_Lfunc - estimated_randomL).copy()
            up_curve_aux[upper_rejected] = np.nan
            down_curve_aux = (observed_Lfunc - estimated_randomL).copy()
            down_curve_aux[under_rejected] = np.nan

            plot_with_confidence(r_list, up_curve_aux, upper_r_test , ax=ax[id_mark], c="r", alpha=0.5)
            plot_with_confidence(r_list, down_curve_aux, under_r_test, ax=ax[id_mark], c="r", alpha=0.5)


    ax[0].set_ylabel("N3ICD mean")
    ax[1].set_ylabel("N3ICD Sum")
    ax[2].set_ylabel("Apical Area")

    fig.suptitle("BAC | N3ICD Sum | Lower Expression group")

    fig.savefig("images/bac_test_of_randomness_lower.pdf", format="pdf", bbox_inches="tight")

    plt.show()
