import numpy as np
from scipy.ndimage.interpolation import shift


# Computer Exercise 1
def G(row_s, Temp):
    return np.exp((1/Temp) * (row_s.T.dot(shift(row_s, -1))))


if __name__ == "__main__":
    print("Defined G.")
    print("------------------------------")
