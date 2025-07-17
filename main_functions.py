import numpy as np
from scipy.linalg import norm
import pickle


dict_for_columns = {(i, j):2*(i-1) + j-1 for i in range(1,4) for j in range(1,3)}

########################################
# Auxiliary functions
########################################

def mark_voisins(X, Y, r, marks, lower_bound, upper_bound):
    """
    Function computing the relative positions of all pair of points close at a distance smaller or equal to r.
    Point are classified in three categories along their marks, delimited by the lower and upper bounds of the middle
    class.
    Returns list of four lists, with relative distance and the mark of the second point in each pair.
    First list takes into account all points.
    Second one corresponds to centered neighbourhoods for the lower class.
    Third one corresponds to centered neighbourhoods for the middle class.
    Fourth one corresponds to centered neighbourhoods for the upper class.
    """
    voisinage = [[], [], [], []]
    for i in range(len(X)):
        x1, y1 = X[i], Y[i]
        for j in range(len(X)):
            x2, y2 = X[j], Y[j]
            flag = norm(np.array([x2 - x1, y2 - y1]))
            if 0.0 < flag <= r:
                if marks[j] < lower_bound:
                    aux_mark = -1.0
                elif marks[j] > upper_bound:
                    aux_mark = 1.0
                else:
                    aux_mark = 0.0

                voisinage[3] += [[x2 - x1, y2 - y1, aux_mark]]
                if marks[i] < lower_bound:
                    voisinage[0] += [[x2 - x1, y2 - y1, aux_mark]]
                elif marks[i] > upper_bound:
                    voisinage[2] += [[x2 - x1, y2 - y1, aux_mark]]
                else:
                    voisinage[1] += [[x2 - x1, y2 - y1, aux_mark]]

    return voisinage


def prepare_vertical_array(array):
    """
    Auxiliary function to get array in good format for other functions.
    """
    if len(array.shape) == 1:
        return array[:, np.newaxis]
    elif len(array.shape) == 2:
        if array.shape[1] == 1:
            return array
        elif array.shape[0] == 1:
            return array.T
        else:
             raise ValueError("Two dimensional array has to have at least one dimension of size 1.")
    else:
        raise ValueError("This is not a possible vertical array. Check dimensions. Has to be either one or two dimensional.")


def update_pickle(file, toUpdate, toCheck, toCheckKeys):
    """
    Auxiliary function to save and update saved estimations.
    """
    try:
        with open(file, 'rb') as f:
            old_data = pickle.load(f)

        idx_to_replace = []
        for i in range(len(old_data)):
            old_keys = [old_data[i][keys] for keys in toCheckKeys]
            if old_keys == toCheck:
                print('Already done before, replacing...')
                idx_to_replace += [i]

        for i in idx_to_replace:
            old_data.pop(i)

        with open(file, 'wb') as f:
            old_data += [toUpdate]
            pickle.dump(old_data, f)

    except FileNotFoundError:
        with open(file, 'wb') as f:
            toSave = [toUpdate]
            pickle.dump(toSave, f)

    return 0
########################################
# Mark-interaction functions
########################################


def m_multiplicative_func(x, y, m1, m2):
    return x * y

def m_r_func(x, y, m1, m2):
    return x

def m_r2_func(x, y, m1, m2):
    return y

def m_variogram_func(x, y, m1, m2):
    return 0.5 * ((x - y) ** 2)

def m_func_inf(x, y, m1, m2):
    return (x <= m1) * y

def m_func_mid(x, y, m1, m2):
    return (m1 < x) * (x < m2) * y

def m_func_sup(x, y, m1, m2):
    return (m2 <= x) * y


########################################
# Spatial and marked K-functions
########################################

def vW1W2_2d(x_diff, y_diff, len_x, len_y):
    """
    Edge-correction term using window translation.
    """
    x = np.positive(len_x - np.abs(x_diff))
    y = np.positive(len_y - np.abs(y_diff))
    return x * y


def spatial_K_func(X, Y, r_list, len_x, len_y):
    """
    Spatial K function for a set of points and list of radiuses and estimated squared mean intensity of the process.
    """
    X_0, Y_0 = prepare_vertical_array(X), prepare_vertical_array(Y)

    points = np.hstack((X_0, Y_0))

    distances = norm(points[np.newaxis, :, :] - points[:, np.newaxis, :], axis=-1)

    # id = distances > 0
    # dist = dist[id]
    vW1 = vW1W2_2d(X_0 - X_0.T, Y_0 - Y_0.T, len_x, len_y)  # .ravel()[id]
    # print(dist, vW1)
    K_func = np.array([np.sum((0 < distances) * (distances <= r) / vW1) for r in r_list])
    sq_mean_intensity = (X_0.shape[0] * (X_0.shape[0] - 1) / ((len_x * len_y) ** 2))

    return K_func, sq_mean_intensity


def individual_marked_K_func(X, Y, r_list, len_x, len_y, lower_bound, upper_bound, marks_1, marks_2, m_func):
    """
    Function computing the signal surrounding each point in a dataset of positions and marks for a mark-interaction
    function m_func.
    Lower and upper bounds are optional depending on the m_func.
    """
    X_0, Y_0 = prepare_vertical_array(X), prepare_vertical_array(Y)

    points = np.hstack((X_0, Y_0))

    distances = norm(points[np.newaxis, :, :] - points[:, np.newaxis, :], axis=-1)

    id = distances > 0
    # dist = dist[id]
    vW1 = vW1W2_2d(X_0 - X_0.T, Y_0 - Y_0.T, len_x, len_y)  # .ravel()[id]
    vW1 = np.ones(vW1.shape)
    # print(dist, vW1)

    m1, m2 = np.meshgrid(marks_1, marks_2)
    marks_2b2 = m_func(m1, m2, lower_bound, upper_bound)
    mean_mark = np.mean(marks_2b2)

    marked_K = np.sum(marks_2b2 * id * (distances <= r_list[0]) / vW1, axis=0)

    return marked_K, mean_mark


def marked_K_func(X, Y, r_list, len_x, len_y, lower_bound, upper_bound, marks_1, marks_2, m_func):
    """
    Marked K-function for a set of points and associated marks and list of radiuses and estimated squared mean
    intensity of the process.
    Lower and upper bounds are optional depending on the m_func.
    """
    X_0, Y_0 = prepare_vertical_array(X), prepare_vertical_array(Y)

    points = np.hstack((X_0, Y_0))

    distances = norm(points[np.newaxis, :, :] - points[:, np.newaxis, :], axis=-1)

    id = distances > 0
    # dist = dist[id]
    vW1 = vW1W2_2d(X_0 - X_0.T, Y_0 - Y_0.T, len_x, len_y)  # .ravel()[id]
    # print(dist, vW1)

    m1, m2 = np.meshgrid(marks_1, marks_2)
    marks_2b2 = m_func(m1, m2, lower_bound, upper_bound)
    mean_mark = np.mean(marks_2b2)

    marked_K = np.array([np.sum(marks_2b2 * id * (distances <= r) / vW1) for r in r_list])

    return marked_K, mean_mark


########################################
# Plotting related functions
########################################

def plot_with_confidence(x, y_down, y_up, ax, c='blue', alpha=0.2, coef_alpha=1.5):
    """
    Auxiliary function to plot filled region between two lines.
    """
    ax.plot(x, y_down, '--', color=c, alpha=alpha*coef_alpha)
    ax.plot(x, y_up, '--', color=c, alpha=alpha*coef_alpha)
    ax.fill_between(x, y_down, y_up, color=c, alpha=alpha)




