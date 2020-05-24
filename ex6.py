import numpy as np
from ex1 import G
from ex2 import F
from ex5 import y2row


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


if __name__ == "__main__":
    print("Computer Exercise 6 output:")
    print(compute_efficient_3x3_Lattice_ZTemp(1.0))
    print(compute_efficient_3x3_Lattice_ZTemp(1.5))
    print(compute_efficient_3x3_Lattice_ZTemp(2.0))
    print("------------------------------")