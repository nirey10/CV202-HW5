import numpy as np
from scipy.ndimage.interpolation import shift

# Computer Exercise 1
def G(row_s, Temp):
    return np.exp((1/Temp) * (row_s.T.dot(shift(row_s, -1, cval=0))))

# Computer Exercise 2
def F(row_s, row_t, Temp):
    return np.exp((1/Temp) * (row_s.T.dot(row_t)))

# Computer Exercise 3
def compute_2x2_Lattice_Ztemp(Temp):
    values = [-1, 1]
    res = 0
    for x_11 in values:
        for x_12 in values:
            for x_21 in values:
                for x_22 in values:
                    internal_sum = x_11*x_12 + x_11*x_21 + x_12*x_22 + x_21*x_22
                    res += np.exp((1/Temp) * internal_sum)
    return res

# Computer Exercise 4
def compute_3x3_Lattice_Ztemp(Temp):
    values = [-1, 1]
    res = 0
    for x_11 in values:
        for x_12 in values:
            for x_13 in values:
                for x_21 in values:
                    for x_22 in values:
                        for x_23 in values:
                            for x_31 in values:
                                for x_32 in values:
                                    for x_33 in values:
                                        internal_sum = x_11*x_12 + x_11*x_21 + x_12*x_22 + x_21*x_22\
                                                        + x_12*x_13 + x_13*x_23 + x_22*x_23 + x_21*x_31\
                                                        + x_22*x_32 + x_31*x_32 + x_32*x_33 + x_23* x_33

                                        res += np.exp((1/Temp) * internal_sum)
    return res


def y2row(y,width=8):
    """
    y: an integer in (0,...,(2**width)-1)
    """
    if not 0<=y<=(2**width)-1:
        raise ValueError(y)
    my_str=np.binary_repr(y,width=width)
    # my_list = map(int,my_str) # Python 2
    my_list = list(map(int,my_str)) # Python 3
    my_array = np.asarray(my_list)
    my_array[my_array==0]=-1
    row=my_array
    return row

# Computer Exercise 5
def compute_efficient_2x2_Lattice_Ztemp(Temp):
    y_s = [0, 1, 2, 3]
    res = 0
    for i in y_s:
        for j in y_s:
            y1 = y2row(i, width=2)
            y2 = y2row(j, width=2)
            res += G(y1, Temp)*G(y2, Temp)*F(y1, y2, Temp)
    return res

# Computer Exercise 6
def compute_efficient_3x3_Lattice_Ztemp(Temp):
    y_s = [0, 1, 2, 3, 4, 5, 6, 7]
    res = 0
    for i in y_s:
        for j in y_s:
            for k in y_s:
                y1 = y2row(i, width=3)
                y2 = y2row(j, width=3)
                y3 = y2row(k, width=3)
                res += G(y1, Temp)*G(y2, Temp)*G(y3, Temp)*F(y1, y2, Temp)*F(y2, y3, Temp)
    return res


# Computer Exercise 7
def get_T_arrays_and_Ztemp(lattice_size, Temp):
    T_arrays = [[] * pow(2, lattice_size)] * (lattice_size-1)
    final_res = 0
    T_arr = [1]*pow(2, lattice_size)
    for i in range(lattice_size-1):
        T_arr = calc_T(lattice_size, Temp, T_arr)
        T_arrays[i] = T_arr
        print('T' + str(i+1) + ': ' + str(T_arr))

    # Printing the last T, the normalizing factor Z_temp
    for i in range(pow(2, lattice_size)):
        y1_vector = y2row(i, width=lattice_size)
        final_res += T_arr[i] * G(y1_vector, Temp)
    print('T' + str(lattice_size) + ' = Z: ' + str(final_res))
    return final_res, T_arrays

def calc_T(lattice_size, Temp, prev_arr):
    T_vec = [0]*pow(2, lattice_size)
    res = 0
    for i in range(pow(2, lattice_size)):
        temp_sum = 0
        for j in range(pow(2, lattice_size)):
            y1_vector = y2row(j, width=lattice_size)
            y2_vector = y2row(i, width=lattice_size)
            temp_res = G(y1_vector, Temp) * F(y1_vector, y2_vector, Temp)*prev_arr[j]
            temp_sum += temp_res
            res += temp_res
        T_vec[i] = temp_sum
    return T_vec

def get_cond_prob_matrix(prob_k, Z_temp, T_arrays, lattice_size, Temp):
    k = prob_k-1
    if prob_k is lattice_size: # for p(y8) in lattice_size = 8 case
        res_matrix = [0]*pow(2, lattice_size)
        for i in range(pow(2, lattice_size)):
            y_vector = y2row(i, width=lattice_size)
            res_matrix[i] = T_arrays[k-1][i] * G(y_vector, Temp) / Z_temp
        return res_matrix
    else:
        if prob_k is 1: # for p(y1|y2)
            res_matrix = [[0] * pow(2, lattice_size)] * pow(2, lattice_size)
            for i in range(pow(2, lattice_size)):
                for j in range(pow(2, lattice_size)):
                    y1_vector = y2row(i, width=lattice_size)
                    y2_vector = y2row(j, width=lattice_size)
                    # notice that i are the rows (y1), j are the columns (y2)
                    res_matrix[i][j] = F(y1_vector, y2_vector, Temp) * G(y1_vector, Temp) / T_arrays[k][j]
            return res_matrix
        else:
            res_matrix = [[0] * pow(2, lattice_size)] * pow(2, lattice_size)
            for i in range(pow(2, lattice_size)):
                for j in range(pow(2, lattice_size)):
                    y1_vector = y2row(i, width=lattice_size)
                    y2_vector = y2row(j, width=lattice_size)
                    # notice that i are the rows (y1), j are the columns (y2)
                    res_matrix[i][j] = F(y1_vector, y2_vector, Temp) * G(y1_vector, Temp) * T_arrays[k-1][i]/ T_arrays[k][j]
            return res_matrix


# Testing
lattice_size = 8
Temp = 1
Z_temp, T_arrays = get_T_arrays_and_Ztemp(lattice_size, Temp)
p12 = get_cond_prob_matrix(1, Z_temp, T_arrays, lattice_size, Temp)
p23 = get_cond_prob_matrix(2, Z_temp, T_arrays, lattice_size, Temp)
p34 = get_cond_prob_matrix(3, Z_temp, T_arrays, lattice_size, Temp)
p45 = get_cond_prob_matrix(4, Z_temp, T_arrays, lattice_size, Temp)
p56 = get_cond_prob_matrix(5, Z_temp, T_arrays, lattice_size, Temp)
p67 = get_cond_prob_matrix(6, Z_temp, T_arrays, lattice_size, Temp)
p78 = get_cond_prob_matrix(7, Z_temp, T_arrays, lattice_size, Temp)
p8 = get_cond_prob_matrix(8, Z_temp, T_arrays, lattice_size, Temp)

