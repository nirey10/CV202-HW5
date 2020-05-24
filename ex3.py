import numpy as np


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


if __name__ == "__main__":
    print("Computer Exercise 3 output:")
    print(compute_2x2_Lattice_ZTemp(1.0))
    print(compute_2x2_Lattice_ZTemp(1.5))
    print(compute_2x2_Lattice_ZTemp(2.0))
    print("------------------------------")
