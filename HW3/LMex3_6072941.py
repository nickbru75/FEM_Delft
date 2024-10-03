import numpy as np
import matplotlib.pyplot as plt


def stiffness_matrix_element(elastic_module: float, surface: float, l: float):
    A = np.array([[7, -8, 1], [-8, 16, -8], [1, -8, 7]])
    # A = np.array([[7, 1, -8], [1, 7, -8], [-8, -8, 16]])
    return elastic_module * surface / (3 * l) * A


def assemble_global_stiffness_matrix(nelem: int, ltot: float, elastic_module: float, surface: float):
    n = 3 + (nelem - 1) * 2
    l = ltot / nelem
    K_global = np.zeros((n, n))
    for elem in range(nelem):
        K_global[2 * elem:2 * elem + 3, 2 * elem:2 * elem + 3] += stiffness_matrix_element(elastic_module, surface, l)
    return K_global


def load_vector_simple(coeff_a: float, coeff_b: float, nelem: int, ltot: float):
    n = 3 + (nelem - 1) * 2
    nodes = np.linspace(0, ltot, n)
    load = np.zeros(n)
    dummy = coeff_a * nodes + coeff_b / 2 * nodes ** 2
    load[1:] = np.diff(dummy)
    return load


def f1(A, B, C, D, x, l_el):
    return 2/l_el**2 * (A*B*C*x + 0.5*A*B*D*x**2 - 0.5*A*C*x**2 - 1/3*A*D*x**3 - 0.5*B*C*x**2 - 1/3*B*D*x**3 + 1/3*C*x**3 + 1/4*D*x**4)


def f2(A, B, C, D, x, l_el):
    return -4/l_el**2 * (A*B*C*x + 0.5*A*B*D*x**2 - 0.5*A*C*x**2 - 1/3*A*D*x**3 - 0.5*B*C*x**2 - 1/3*B*D*x**3 + 1/3*C*x**3 + 1/4*D*x**4)


def f3(A, B, C, D, x, l_el):
    return 2/l_el**2 * (A*B*C*x + 0.5*A*B*D*x**2 - 0.5*A*C*x**2 - 1/3*A*D*x**3 - 0.5*B*C*x**2 - 1/3*B*D*x**3 + 1/3*C*x**3 + 1/4*D*x**4)


def load_vector_shape_functions(coeff_a: float, coeff_b: float, nelem: int, ltot: float):
    n = 3 + (nelem - 1) * 2
    l = ltot/nelem
    load = np.zeros(n)
    nodes = np.linspace(0, ltot, n)
    for ii in range(nelem):
        x_i = nodes[2 * ii]
        x_j = nodes[2 * ii + 1]
        x_k = nodes[2 * ii + 2]
        # load[2 * ii] += (2*coeff_a/l**2 * (
        #         1/3*(x_k**3 - x_i**3) - 1/2*(x_k**2 - x_i**2)*(x_j + x_k) + (x_k-x_i)*x_j*x_k) +
        #                 2*coeff_b/l**2 * (
        #                         1/4*(x_k**4 - x_i**4) - 1/3*(x_k**3 - x_i**3)*(x_j+x_k) + 1/2*(x_k**2 - x_i**2)*x_j*x_k))
        # load[2 * ii + 1] += -2*(
        #         2*coeff_a/l**2 * (1/3*(x_k**3 - x_i**3) - 1/2*(x_k**2 - x_i**2)*(x_i + x_k) + (x_k-x_i)*x_i*x_k) +
        #                 2*coeff_b/l**2 * (
        #                 1/4*(x_k**4 - x_i**4) - 1/3*(x_k**3 - x_i**3)*(x_i+x_k) + 1/2*(x_k**2 - x_i**2)*x_i*x_k))
        # load[2 * ii + 2] += (2*coeff_a/l**2 * (
        #         1/3*(x_k**3 - x_i**3) - 1/2*(x_k**2 - x_i**2)*(x_j + x_i) + (x_k-x_i)*x_j*x_i) +
        #                 2*coeff_b/l**2 * (
        #                             1/4*(x_k**4 - x_i**4) - 1/3*(x_k**3 - x_i**3)*(x_j+x_i) + 1/2*(x_k**2 - x_i**2)*x_j*x_i))
        load[2 * ii] += f1(x_j, x_k, coeff_a, coeff_b, x_k, l) - f1(x_j, x_k, coeff_a, coeff_b, x_i, l)
        load[2 * ii + 1] += f2(x_i, x_k, coeff_a, coeff_b, x_k, l) - f2(x_i, x_k, coeff_a, coeff_b, x_i, l)
        load[2 * ii + 2] += f3(x_i, x_j, coeff_a, coeff_b, x_k, l) - f3(x_i, x_j, coeff_a, coeff_b, x_i, l)
    return load


def apply_problem_speficic_bc(K_global, load):
    K_global = K_global[1:, 1:]
    load = load[1:]
    return K_global, load


def quadratic_elements_problem(ltot, surface, elastic_module, coeff_a, coeff_b, nelem):
    F = load_vector_shape_functions(coeff_a, coeff_b, nelem, ltot)
    K = assemble_global_stiffness_matrix(nelem, ltot, elastic_module, surface)
    K_reduced, F_reduced = apply_problem_speficic_bc(K, F)
    u_partial = np.linalg.solve(K_reduced, F_reduced)
    u_complete = np.zeros(3 + (nelem - 1) * 2)
    u_complete[1:] = u_partial
    return u_complete


def true_displacement_and_error(coeff_a, coeff_b, surface, elastic_module, sol, ltot, nelem):
    n = 3 + (nelem - 1) * 2
    nodes = np.linspace(0, ltot, n)
    true = 1/(surface*elastic_module)*(coeff_a/2 * nodes**2 + coeff_b/6 * nodes**3)
    error_ = np.sqrt(np.sum((true-sol)**2))
    return true, error_


# INPUT VALUES
L = 500
Area = 120
E = 70e3
a = 13
b = 0.
N = 100
# values with solution online
# L = 4
# A = 0.003
# E = 210e9
# F = np.array([0, 5e3, -10e3, -7e3, 10e3])
# N = 2

# COMPUTING FIXED NELEM PROBLEM
n = 3 + (N - 1) * 2
nodes = np.linspace(0, L, n)

u = quadratic_elements_problem(L, Area, E, a, b, N)
F = load_vector_shape_functions(a, b, N, L)
print(u[-1])
plt.figure(figsize=(8, 5))
F_dummy = a*nodes + b/2*nodes**2
F_true = np.zeros_like(nodes)
F_true[1:] = np.diff(F_dummy)
plt.plot(nodes, F, label='Computed')
plt.plot(nodes, F_true, label='Analytical')
plt.legend()
plt.grid()
plt.xlabel('Number of elements [-]')
plt.ylabel('Force [N]')
plt.savefig('force.png', dpi=200)

u_true, error = true_displacement_and_error(a, b, Area, E, u, L, N)
print(u_true[-1])
plt.figure(figsize=(8, 5))
plt.plot(nodes, u, label='Computed')
plt.plot(nodes, u_true, label='Analytical')
plt.legend()
plt.grid()
plt.xlabel('Number of elements [-]')
plt.ylabel('Displacement of last node [mm]')
plt.savefig('solution_vs_true.png', dpi=200)

# CONVERGENCE PROBLEM
# nmax = 100
# u_end = np.zeros(nmax)
# err = np.zeros(nmax)
# for i in range(1, nmax+1):
#     u = quadratic_elements_problem(L, Area, E, a, b, i)
#     u_end[i-1] = u[-1]
#     u_true, error = true_displacement_and_error(a, b, Area, E, u, L, i)
#     err[i-1] = error
# plt.figure(figsize=(8, 5))
# plt.plot(np.arange(1, nmax+1), u_end)
# plt.plot(np.arange(1, nmax+1), err)
# plt.grid()
# plt.xlabel('Number of elements [-]')
# plt.ylabel('Displacement of last node [mm]')
# plt.title('Convergence Study')
# plt.savefig('convergence.png', dpi=200)

