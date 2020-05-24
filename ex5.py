import numpy as np
from ex1 import G
from ex2 import F


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


if __name__ == "__main__":
    print("Computer Exercise 5 output:")
    print(compute_efficient_2x2_Lattice_ZTemp(1.0))
    print(compute_efficient_2x2_Lattice_ZTemp(1.5))
    print(compute_efficient_2x2_Lattice_ZTemp(2.0))
    print("------------------------------")
