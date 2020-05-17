import numpy as np
from scipy.ndimage.interpolation import shift

# Computer Exercise 1
def G(row_s, Temp):
    return np.exp((1/Temp) * (row_s.T.dot(shift(row_s, -1, cval=0))))

# Computer Exercise 2
def F(row_s, row_t, Temp):
    return np.exp((1/Temp) * (row_s.T.dot(row_t)))





a = np.array([1,2,3])
b = np.array([1,2,3])

print(F(a, b, 1))