import numpy as np
from scipy.integrate import quad


# Define individual shape functions for the quadratic element
def shape_function_Ni(z, x_i, x_j, x_k):
    return 2 * (z - x_j) * (z - x_k) / (x_k - x_i)**2


def shape_function_Nj(z, x_i, x_j, x_k):
    return -4 * (z - x_i) * (z - x_k) / (x_k - x_i)**2


def shape_function_Nk(z, x_i, x_j, x_k):
    return 2 * (z - x_i) * (z - x_j) / (x_k - x_i)**2


# Define the load vector for a single element
def load_vector_element(a, b, x_i, x_j, x_k):
    # Define the distributed load function g(z) = a + bz
    g = lambda z: a + b * z

    # Initialize the load vector
    f = np.zeros(3)

    # Compute each component of the load vector by integrating N_i * g(z)
    f[0], _ = quad(lambda z: shape_function_Ni(z, x_i, x_j, x_k) * g(z), x_i, x_k)
    f[1], _ = quad(lambda z: shape_function_Nj(z, x_i, x_j, x_k) * g(z), x_i, x_k)
    f[2], _ = quad(lambda z: shape_function_Nk(z, x_i, x_j, x_k) * g(z), x_i, x_k)

    return f


# Example usage
a = 13
b = 0.13
x_i = 0
x_j = 250  # Midpoint
x_k = 500

f_element = load_vector_element(a, b, x_i, x_j, x_k)
print(f_element)
