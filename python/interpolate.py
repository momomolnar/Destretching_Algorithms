"""
A collection of 1D and 2D interpolation functions to be used in the code

Rewritten to be more easily parallelized.
"""

import numpy as np


def interp_1D_spline(grid, data):
    """
    Calculate the parameters for the 1D cubic spline interpolator 
    -- second derivatives on the grid. Assuming the second derivative to vanish
    at the boundary conditions.
    
    Following the ideas of implementation from the Numerical Recipes book.
    
    :param: grid: Grid on which the data is interpolated.
    :param data: Evaluations of the function to be interpolated.
    :return: interpolator_1D_spline: function to interpolate a location 
    _loc_ on the grid.
    evaluated on the grid. 
    """
    len_grid = len(grid)

    A = np.zeros((len_grid-2, len_grid-2))

    for j in range(1, len_grid-1):
        A[j-1][j-1] = (grid[j+1] - grid[j-1]) /3 #Diagonal elements
    for j in range(1, len_grid-2):
        A[j][j-1] = (grid[j] - grid[j-1])/6
    for j in range(1, len_grid - 2):
        A[j-1][j] = (grid[j] - grid[j - 1]) / 6

    b1 = np.array([(data[j+1]-data[j])/(grid[j+1] - grid[j])
                   for j in range(1, len_grid-1)])
    b2 = np.array([(data[j]-data[j-1])/(grid[j] - grid[j-1])
                   for j in range(1, len_grid-1)])
    b = b1 - b2

    second_der_t = np.linalg.solve(A, b)
    print(second_der_t)

    second_der = np.zeros(len_grid)
    second_der[1:-1] = second_der_t

    def interpolator_1D_spline(loc, data=data, second_der=second_der):
        '''
        Function to interpolate a function using a cubic spline,
        provided the derivatives needed.

        :param: loc: where the function to be estimated
        :param data: evaluations of the functions
        :param second_der: second derivatives needed for the
            spline

        :return: y : evaluation of the function
        '''

        #find the neighbor
        j1_ind = int(np.argmax(1/(grid-loc))) #1 over the smalles positive number
                                         # gives the biggest positive number
                                         # and the right neighbor
        j_ind = int(j1_ind - 1)

        xj  = grid[j_ind]
        xj1 = grid[j1_ind]
        yj  = data[j_ind]
        yj1 = data[j1_ind]
        ddyj = second_der[j_ind]
        ddyj1 = second_der[j1_ind]

        A = (xj1 - loc) / (xj1 - xj)
        B = 1 - A
        C = (A**3 - A) * (xj1 - xj)**2 / 6
        D = (B**3 - B) * (xj1 - xj)**2 / 6
        y = A*yj + B*yj1 + C*ddyj + D*ddyj1

        return y


    return interpolator_1D_spline

def cubic_fn(x, C):
    return C[0] * x**3 + C[1] * x**2 + C[2] * x**1 + C[3]

if __name__ == "__main__":
    x = np.linspace(-10, 10, num=9)
    y = cubic_fn(x, [0, -1, 0 , 1])
    f = interp_1D_spline(x, y)

    x_tight = np.linspace(10, -10, num=1000)
    y_true = cubic_fn(x_tight, [0, -1, 0 , 1])
    y_interp = np.array([f(el) for el in x_tight])

    import matplotlib.pyplot as pl

    pl.plot(x, y, 'ro')
    pl.plot(x_tight, y_true, 'b--', label="True")
    pl.plot(x_tight, y_interp, 'g.--', label="Interp")
    pl.grid(alpha=0.5)
    pl.legend()
    pl.show()