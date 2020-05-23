import numpy as np
from scipy.ndimage.interpolation import shift

print("------------------------------")


# Computer Exercise 1
def G(row_s, Temp):
    return np.exp((1/Temp) * (row_s.T.dot(shift(row_s, -1))))


print("Defined G.")
print("------------------------------")


# Computer Exercise 2
def F(row_s, row_t, Temp):
    return np.exp((1/Temp) * (row_s.T.dot(row_t)))


print("Defined F.")
print("------------------------------")


# Computer Exercise 3
def compute_2x2_Lattice_ZTemp(Temp):
    values = [-1, 1]
    res = 0
    for x_11 in values:
        for x_12 in values:
            for x_21 in values:
                for x_22 in values:
                    internal_sum = x_11*x_12 + x_11*x_21 + x_12*x_22 + x_21*x_22
                    res += np.exp((1/Temp) * internal_sum)
    return res


print("Computer Exercise 3 output:")
print(compute_2x2_Lattice_ZTemp(1.0))
print(compute_2x2_Lattice_ZTemp(1.5))
print(compute_2x2_Lattice_ZTemp(2.0))
print("------------------------------")


# Computer Exercise 4
def compute_3x3_Lattice_ZTemp(Temp):
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


print("Computer Exercise 4 output:")
print(compute_3x3_Lattice_ZTemp(1.0))
print(compute_3x3_Lattice_ZTemp(1.5))
print(compute_3x3_Lattice_ZTemp(2.0))
print("------------------------------")


def y2row(y, width=8):
    """
    y: an integer in (0,...,(2**width)-1)
    """
    if not 0 <= y <= (2**width)-1:
        raise ValueError(y)
    my_str = np.binary_repr(y, width=width)
    # my_list = map(int,my_str) # Python 2
    my_list = list(map(int, my_str))  # Python 3
    my_array = np.asarray(my_list)
    my_array[my_array == 0] = -1
    row = my_array
    return row


# Computer Exercise 5
def compute_efficient_2x2_Lattice_ZTemp(Temp):
    y_s = [0, 1, 2, 3]
    res = 0
    for i in y_s:
        for j in y_s:
            y1 = y2row(i, width=2)
            y2 = y2row(j, width=2)
            res += G(y1, Temp)*G(y2, Temp)*F(y1, y2, Temp)
    return res


print("Computer Exercise 5 output:")
print(compute_efficient_2x2_Lattice_ZTemp(1.0))
print(compute_efficient_2x2_Lattice_ZTemp(1.5))
print(compute_efficient_2x2_Lattice_ZTemp(2.0))
print("------------------------------")


# Computer Exercise 6
def compute_efficient_3x3_Lattice_ZTemp(Temp):
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


print("Computer Exercise 6 output:")
print(compute_efficient_3x3_Lattice_ZTemp(1.0))
print(compute_efficient_3x3_Lattice_ZTemp(1.5))
print(compute_efficient_3x3_Lattice_ZTemp(2.0))
print("------------------------------")


# Computer Exercise 7
def get_T_arrays(lattice_size, Temp):
    T_arrays = [[]] * lattice_size
    # final_res = 0
    prev_T = [1] * pow(2, lattice_size)  # T_0
    for k in range(lattice_size-1):
        curr_T = calc_T(lattice_size, Temp, prev_T)
        T_arrays[k] = curr_T
        # print('T' + str(k+1) + ': ' + str(curr_T))
        prev_T = curr_T

    # Printing the last T, the normalizing factor Z_temp
    T_arrays[lattice_size - 1] = 0
    for y_last in range(pow(2, lattice_size)):
        y_last_vector = y2row(y_last, width=lattice_size)
        T_arrays[lattice_size-1] += curr_T[y_last] * G(y_last_vector, Temp)
    # print('T' + str(lattice_size) + ' = Z: ' + str(T_arrays[lattice_size-1]))
    return T_arrays


def calc_T(lattice_size, Temp, prev_T):
    T_vec = [0]*pow(2, lattice_size)
    res = 0
    for y2 in range(pow(2, lattice_size)):
        temp_sum = 0
        for y1 in range(pow(2, lattice_size)):
            y1_vector = y2row(y1, width=lattice_size)
            y2_vector = y2row(y2, width=lattice_size)
            temp_sum += G(y1_vector, Temp) * F(y1_vector, y2_vector, Temp)*prev_T[y1]
            # res += temp_res
        T_vec[y2] = temp_sum
    return T_vec


def get_p_k(k, ZTemp, T, lattice_size, Temp):
    k -= 1  # converting to array indexes (was in range 1,..,lattice_size)
    if k == lattice_size-1:  # for p_last(y_last)
        res_matrix = np.ndarray((2**lattice_size , 1), dtype=float)
        # res_matrix = [0] * pow(2, lattice_size)
        for y_last in range(pow(2, lattice_size)):
            y_last_vector = y2row(y_last, width=lattice_size)
            res_matrix[y_last] = (T[lattice_size-2][y_last] * G(y_last_vector, Temp)) / ZTemp
        return res_matrix
    elif k == 0:  # p_(1|2)(y1|y2)
        # res_matrix = [[0] * pow(2, lattice_size)] * pow(2, lattice_size)
        res_matrix = np.ndarray((2**lattice_size, 2**lattice_size), dtype=float)
        for y1 in range(pow(2, lattice_size)):
            y1_vector = y2row(y1, width=lattice_size)
            for y2 in range(pow(2, lattice_size)):
                y2_vector = y2row(y2, width=lattice_size)
                res_matrix[y1][y2] = (F(y1_vector, y2_vector, Temp) * G(y1_vector, Temp)) / T[0][y2]
        return res_matrix
    else:
        # res_matrix = [[0] * pow(2, lattice_size)] * pow(2, lattice_size)
        res_matrix = np.ndarray((2 ** lattice_size, 2 ** lattice_size), dtype=float)
        for y_k in range(pow(2, lattice_size)):
            y_k_vector = y2row(y_k, width=lattice_size)
            for y_kplus1 in range(pow(2, lattice_size)):
                y_kplus1_vector = y2row(y_kplus1, width=lattice_size)
                res_matrix[y_k][y_kplus1] = (F(y_k_vector, y_kplus1_vector, Temp) * G(y_k_vector, Temp) * T[k-1][y_k]) / T[k][y_kplus1]
        return res_matrix


def calc_p(ZTemp, T, lattice_size, Temp):
    p = []
    for k in range(1, lattice_size+1):
        p.append(get_p_k(k, ZTemp, T, lattice_size, Temp))
    return p


def backward_sample(p, lattice_size):
    y = [0]*lattice_size
    y[lattice_size-1] = np.random.choice(pow(2, lattice_size), p=p[lattice_size-1][:, 0])
    k = lattice_size - 2
    while k >= 0:
        y[k] = np.random.choice(pow(2, lattice_size), p=p[k][:, y[k+1]])
        k -= 1
    # for k in range(lattice_size-1):
    #     prev_y = y_arr[lattice_size - k - 1]
    #     new_rand = np.random.choice(pow(2, lattice_size), p=p[lattice_size - k - 1][prev_y])
    #     y_arr[lattice_size-i - 2] = new_rand[0]
    return y


def convert_y_to_image(y, lattice_size):
    # image = [[]*lattice_size]*lattice_size
    image = np.ndarray((lattice_size, lattice_size))
    index = 0
    for i in range(len(y)):
        # row_values = [int(i) for i in list('{0:0b}'.format(row))]
        row = y2row(y[i], width=lattice_size)
        for j in range(len(row)):
            if row[j] == 0:
                row[j] = -1
        # missing_zeros = lattice_size - len(row_values)
        # padded_row_values = np.pad(row_values, (missing_zeros, 0), 'constant')
        # image[index] = padded_row_values
        # index += 1
        image[i] = row
    # np_image = np.array([arr for arr in image])
    return image


# Computer Exercise 8
def compute_empirical_expectation(i1, j1, i2, j2, samples):
    normalizing_factor = len(samples)
    sum = 0
    for sample in samples:
        sum += sample[i1][j1] * sample[i2][j2]
    return sum/normalizing_factor

print("Computer Exercise 7 output:")
# Testing
lattice_size = 3

for Temp in [1.0, 1.5, 2.0]:
    Ts = get_T_arrays(lattice_size, Temp)
    ZTemp = Ts[-1]
    print("Temp: " + str(Temp))
    # for i in range(lattice_size):
    #     print("T" + str(i+1) + ": " + str(Ts[i]))
    p = calc_p(ZTemp, Ts, lattice_size, Temp)
    # for i in range(lattice_size):
    #     print("p" + str(i+1) + ":")
    #     print(str(p[i]))
    y = backward_sample(p, lattice_size)
    print("y:")
    print(y)
    y_image = convert_y_to_image(y, lattice_size)
    print("image:")
    print(y_image)
    print("---------------")