import numpy as np
from ex7 import get_T_arrays, calc_p, backward_sample, convert_y_to_image
import time


# Computer Exercise 8
def compute_empirical_expectation(i1, j1, i2, j2, samples):
    normalizing_factor = len(samples)
    sum = 0
    for sample in samples:
        sum += sample[i1][j1] * sample[i2][j2]
    return sum/normalizing_factor


if __name__ == "__main__":
    print("Computer Exercise 7 output:")
    lattice_size = 8

    Temps = [1.0, 1.5, 2.0]
    NSAMPLES = 10000

    for Temp in Temps:
        print("Temp = " + str(Temp))
        Ts = get_T_arrays(lattice_size, Temp)
        ZTemp = Ts[-1]
        p = calc_p(ZTemp, Ts, lattice_size, Temp)

        print("Sampling:")
        s = time.time()
        x = [[]] * NSAMPLES
        for n in range(NSAMPLES):
            sample = backward_sample(p, lattice_size)
            x[n] = convert_y_to_image(sample, lattice_size)
        e = time.time()
        print("Done (" + str(round(e - s, 2)) + " secs).")

        print("Empirical Expectation 1:")
        E1 = compute_empirical_expectation(0, 0, 1, 1, x)
        print(E1)
        print("Empirical Expectation 2:")
        E2 = compute_empirical_expectation(0, 0, lattice_size-1, lattice_size-1, x)
        print(E2)
        print("---------------")