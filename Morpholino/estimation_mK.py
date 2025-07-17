import pandas as pd
import seaborn as sns
from main_functions import *
from multiprocessing import Pool

sns.set_theme()

# ['ID Imaris',                                     0
#  'N3ICD mean', 'N3ICD mediane', 'N3ICD Sum',      1-3
#  'PCNA', 'RFP',                                   4-5
#  'Sox2 mediane', 'nlsRFP mediane',                6-7
#  'Volume',                                        8
#  'Apical Area',                                   9
#  'X', 'Y', 'Z']                                   10-12

toSaveKeysSpatial = ['control_id',
                     'max_r', 'divisions',
                     'r_list', 'K_func', 'sq_mean_intensity']

toCheckKeysSpatial = ['control_id', 'max_r', 'divisions']

toSaveKeys = ['control_id', 'n_permut', 'q',
              'mark_id', 'max_r', 'divisions',
              'r_list', 'marked_K', 'mean_mark']

toCheckKeys = ['control_id', 'n_permut', 'q', 'mark_id', 'max_r', 'divisions']


if __name__ == "__main__":
    # Data reading
    df = pd.read_excel("data_summary.xlsx", sheet_name=0)

    control_id_list = [0, 1]

    # Choice of number of random permutations
    n_permut = 199
    # Marks to use in estimation
    mark_id = 'NOTCH mean'#, 'NOTCH mean']

    # List of mark-interactions functions to use in estimationfunc_list = [m_func_inf]
    func_list = [m_func_inf]
    func_name_list = ["func_inf"]

    # Choice of radiuses to study
    divisions = 201
    max_r = 30
    r_list = np.linspace(0, max_r, divisions)

    # Proportion (quantile) to determine lower and upper bounds.
    q = 0.2

    for m_name, m_func in zip(func_name_list, func_list):

        for control_id in control_id_list:
            df_poisson = df.query('`J1b ?`==@control_id')
            n = len(df_poisson)

            K_func = np.zeros(divisions)
            sq_mean_intensity = np.zeros(1)
            marked_K = np.zeros((n_permut + 1, divisions))
            mean_mark = np.zeros(1)

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

                        # Determine window of study using average volume of nuclei
                        volume = df_1['Volume']
                        radius = (3 * volume / (4 * np.pi)) ** (1 / 3)
                        mean_radius = radius.mean()
                        print("Mean radius", mean_radius)

                        min_x, max_x = np.min(df_1['X']), np.max(df_1['X'])
                        min_y, max_y = np.min(df_1['Y']), np.max(df_1['Y'])

                        len_x, len_y = max_x - min_x + 2 * mean_radius, max_y - min_y + 2 * mean_radius

                        X = df_1['X'].to_numpy()
                        Y = df_1['Y'].to_numpy()

                        # Estimate spatial K-function and squared density.
                        K_func_aux, sq_mean_intensity_aux = spatial_K_func(X, Y, r_list, len_x, len_y)

                        K_func += K_func_aux * n_1 / n
                        sq_mean_intensity += sq_mean_intensity_aux * n_1 / n

                        # Marks for classification in groups.
                        marks_1 = df_1["NOTCH sum"].to_numpy()
                        # Marks to use in marked K-function
                        marks_2 = df_1[mark_id].to_numpy()

                        # Compute lower and upper bounds for marks
                        lower_bound, upper_bound = np.quantile(marks_1, q=q), np.quantile(marks_1, q=1 - q)

                        # Create list of mark permutations
                        np.random.seed(0)
                        marks_permut_1 = [marks_1] + [np.random.permutation(marks_1) for i in range(n_permut)]
                        np.random.seed(0)
                        marks_permut_2 = [marks_2] + [np.random.permutation(marks_2) for i in range(n_permut)]

                        # Estimation of marked K-function
                        with Pool(5) as p:
                            aux = zip([X] * (1 + n_permut), [Y] * (1 + n_permut), [r_list] * (1 + n_permut),
                                      [len_x] * (1 + n_permut), [len_y] * (1 + n_permut),
                                      [lower_bound] * (1 + n_permut), [upper_bound] * (1 + n_permut),
                                      marks_permut_1, marks_permut_2, [m_func] * (1 + n_permut))
                            res = p.starmap(marked_K_func, aux)

                        marked_K_aux = np.array([resi[0] for resi in res])
                        mean_mark_aux = res[0][1]

                        print(marked_K_aux.shape, n_1, n)

                        marked_K += marked_K_aux * n_1 / n
                        mean_mark += mean_mark_aux * n_1 / n

            # Save spatial and marked estimations
            fileSpatial = 'saved_estimations/spatial_K.pkl'

            toSaveSpatial = {s: globals()[s] for s in toSaveKeysSpatial}
            toCheckSpatial = [globals()[s] for s in toCheckKeysSpatial]

            update_pickle(fileSpatial, toSaveSpatial, toCheckSpatial, toCheckKeysSpatial)
            #print("*" * 150)

            fileMarked = 'saved_estimations/permut_' + m_name + '_mK_pool.pkl'

            toSave = {s: globals()[s] for s in toSaveKeys}
            toCheck = [globals()[s] for s in toCheckKeys]

            update_pickle(fileMarked, toSave, toCheck, toCheckKeys)
            print("*" * 150)
