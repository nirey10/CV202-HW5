import numpy as np


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


if __name__ == "__main__":
    print("Computer Exercise 4 output:")
    print(compute_3x3_Lattice_ZTemp(1.0))
    print(compute_3x3_Lattice_ZTemp(1.5))
    print(compute_3x3_Lattice_ZTemp(2.0))
    print("------------------------------")
