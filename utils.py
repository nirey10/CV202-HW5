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


# ------------------------------------------------------------------------------------------------------------------
def compute_efficient_generic_Lattice_Ztemp(lattice_size, Temp):
    y_s = range(pow(lattice_size, 2))
    y_arr = [0] * pow(lattice_size, 2)
    res = recursive(lattice_size, lattice_size, Temp, y_arr, 0)
    return res

def calc(lattice_size, index_arr, Temp):
    res = 1
    count = 1
    for i in range(lattice_size):
        print(i)
        res *= G(y2row(index_arr[i], width=lattice_size), Temp)
        if count < lattice_size:
            res *= F(y2row(index_arr[i], width=lattice_size), y2row(index_arr[i+1], width=lattice_size), Temp)
        count += 1
    return res

def recursive(counter, lattice_size, Temp, index_arr, res):
    if counter == -1:
        return calc(lattice_size, index_arr, Temp)
    for i in range(pow(lattice_size, 2)):
        index_arr[counter-1] = i
        res += recursive(counter - 1, lattice_size, Temp, index_arr, res)
        return res
# ------------------------------------------------------------------------------------------------------------------

def dynamic_programming(lattice_size, Temp):
    T_arr = [[] * lattice_size] * (lattice_size-1)
    res = 0
    res_vec = [0] * pow(2, lattice_size)
    for i in range(pow(2, lattice_size)):
        y1_vector = y2row(i, width=lattice_size)
        temp_res = recursive_dynamic(lattice_size - 1, lattice_size, Temp, i, T_arr)
        res += temp_res * G(y1_vector, Temp)
        res_vec[i] = temp_res

    for T_idx in range(lattice_size-2):
        print('T' + str(T_idx + 1) + ': ' + str(T_arr[T_idx + 1]))
    print('T' + str(lattice_size-1) + ': ' + str(res_vec))
    print('T' + str(lattice_size) + ' = Z: ' + str(res))
    return res

def recursive_dynamic(counter, lattice_size, Temp, y2, T_arr):
    y2_vector = y2row(y2, width=lattice_size)
    res = 0
    res_vec = [0] * pow(2, lattice_size)
    if counter is 1:
        for i in range(pow(2, lattice_size)):
            y1_vector = y2row(i, width=lattice_size)
            res += G(y1_vector, Temp) * F(y1_vector, y2_vector, Temp)
    else:
        for i in range(pow(2, lattice_size)):
            y1_vector = y2row(i, width=lattice_size)
            temp_res = recursive_dynamic(counter - 1, lattice_size, Temp, i, T_arr)
            res += temp_res * G(y1_vector, Temp) * F(y1_vector, y2_vector, Temp)
            res_vec[i] = temp_res
        T_arr[counter-1] = res_vec

    return res

# Testing
dynamic_programming(2, 1)