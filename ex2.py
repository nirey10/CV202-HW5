import numpy as np


# Computer Exercise 2
def F(row_s, row_t, Temp):
    return np.exp((1/Temp) * (row_s.T.dot(row_t)))


if __name__ == "__main__":
    print("Defined F.")
    print("------------------------------")
