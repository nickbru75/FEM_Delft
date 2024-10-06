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


def true_displacement_and_error(coeff_a, coeff_b, surface, elastic_module, sol, ltot, nelem, number_points):
    n = 3 + (nelem - 1) * 2
    nodes = np.linspace(0, ltot, n)
    l = ltot/nelem
    nodes_dense = np.linspace(0, ltot, number_points*nelem - (nelem-1))
    sol_dense = np.zeros_like(nodes_dense)
    for ii in range(nelem):
        x_vals = nodes[2*ii:2*ii+3]
        y_vals = sol[2*ii:2*ii+3]

        A = np.array([
            [x_vals[0] ** 2, x_vals[0], 1],
            [x_vals[1] ** 2, x_vals[1], 1],
            [x_vals[2] ** 2, x_vals[2], 1]
        ])
        coeffs = np.linalg.solve(A, y_vals)
        a_, b_, c_ = coeffs

        def quadratic_polynomial(x):
            return a_ * x ** 2 + b_ * x + c_

        sol_dense[(number_points - 1) * ii:(number_points - 1) * ii + number_points] = quadratic_polynomial(
            nodes_dense[(number_points - 1) * ii:(number_points - 1) * ii + (number_points)])

    true = 1/(surface*elastic_module)*(-(coeff_a/2 * nodes_dense**2 + coeff_b/6 * nodes_dense**3)+nodes_dense*(coeff_a*ltot + coeff_b*ltot**2/2))
    error_ = np.sqrt(np.sum((true-sol_dense)**2)/len(true))

    return true, error_, sol_dense, nodes_dense


# INPUT VALUES
L = 500
Area = 120
E = 70e3
a = 13.
b = 0.13
N = 5


u = quadratic_elements_problem(L, Area, E, a, b, N)
u_true, err, u_dense, n_dense = true_displacement_and_error(a, b, Area, E, u, L, N, 500)
plt.figure(figsize=(8, 5))
plt.plot(n_dense, u_dense, label='Computed')
plt.plot(n_dense, u_true, 'r--', label='Analytical')
plt.legend()
plt.grid()
plt.xlabel('Number of elements [-]')
plt.ylabel('Displacement of last node [mm]')
plt.title('Computed solution with {} quadratic elements'.format(N))
plt.savefig('solution_vs_true.png', dpi=200)


# CONVERGENCE PROBLEM
nmax = 10
u_end = np.zeros(nmax)
err = np.zeros(nmax)
for i in range(1, nmax+1):
    u = quadratic_elements_problem(L, Area, E, a, b, i)
    u_end[i-1] = u[-1]
    u_true, error, u_dense, n_dense = true_displacement_and_error(a, b, Area, E, u, L, i, 500)
    err[i-1] = error
plt.figure(figsize=(8, 5))
plt.plot(np.arange(1, nmax+1), err)
plt.grid()
plt.xlabel('Number of elements [-]')
plt.ylabel('RMS difference from true solution [mm]')
plt.title('Convergence Study')
plt.savefig('convergence.png', dpi=200)

