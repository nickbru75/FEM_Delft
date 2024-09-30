import numpy as np

def stiffness_matrix_element(x: np.array):
    # general data
    x_i = x[0]
    x_j = x[1]
    x_k = x[2]
    l = x_k - x_i
    # matrix terms
    K11 = (16/(3*l**4) * x_k**3 - 4/l**4 * (x_j+x_k) * x_k**2 + (x_j + x_k)**2/l**4 * x_k -
           (16/(3*l**4) * x_i**3 - 4/l**4 * (x_j+x_k) * x_i**2 + (x_j + x_k)**2/l**4 * x_i))
    K22 = (64/(3*l**4) * x_k**2 - 32/l**4 * (x_i+x_k) * x_k**2 + 16/l**4 * (x_i+x_k)**2 * x_k -
           (64/(3*l**4) * x_i**2 - 32/l**4 * (x_i+x_k) * x_i**2 + 16/l**4 * (x_i+x_k)**2 * x_i))
    K33 = (16/(3*l**4) * x_k**3 - 4/l**4 * (x_j+x_i) * x_k**2 + (x_j + x_i)**2/l**4 * x_k -
           (16/(3*l**4) * x_i**3 - 4/l**4 * (x_j+x_i) * x_i**2 + (x_j + x_i)**2/l**4 * x_i))
    K12 =