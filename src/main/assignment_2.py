import numpy as np


def nevilles_method(x, y, number):
    # create a matrix
    matrix1 = np.zeros((len(x), len(x)))

    # initialize matrix1 with the y values
    for index, row in enumerate(matrix1):
        row[0] = y[index]

    # find the values for matrix1
    i: int = 0
    j: int = 0
    for i in range(1, len(x)):
        for j in range(1, i+1):
            product1 = (number - x[i-j]) * matrix1[i][j-1]
            product2 = (number - x[i]) * matrix1[i-1][j-1]
            denominator = x[i] - x[i - j]

            matrix1[i][j] = (product1 - product2)/denominator

    return matrix1


def divided_difference_table(x, y):
    # create a matrix
    matrix2 = np.zeros((len(x), len(x)))

    # initialize matrix2 with the y values
    for index, row in enumerate(matrix2):
        row[0] = y[index]

    # find the values for matrix2
    i: int = 0
    j: int = 0
    for i in range(1, len(x)):
        for j in range(1, i+1):
            # numerator is result of left - diagonal left matrix values
            numerator = matrix2[i][j-1] - matrix2[i-1][j-1]
            denominator = x[i] - x[i-j]
            matrix2[i][j] = numerator/denominator

    return matrix2


def find_approximation(table, x, number):
    # initialize value for product of (number - x) terms
    product = 1
    # initialize approximation result
    result = table[0][0]

    i: int = 0
    for i in range(1, len(x)):
        polynomial_coefficient = table[i][i]
        product *= (number - x[i-1])
        multiplication = polynomial_coefficient * product
        result += multiplication

    return result


def hermite(x, y, slopes):
    num_points = len(x)
    matrix4 = np.zeros((2*len(x), 2*len(x)))

    i: int = 0
    j: int = 0
    for i in range(0, (2*len(x))-1, 2):
        matrix4[i][0] = x[j]
        matrix4[i][1] = y[j]
        matrix4[i+1][0] = x[j]
        matrix4[i+1][1] = y[j]
        j += 1
    i = 0
    j = 0
    for i in range(1, 2*len(x), 2):
        matrix4[i][2] = slopes[j]
        if i+1 < 2*len(x):
            matrix4[i+1][2] = ((matrix4[i+1][1] - matrix4[i][1])/(matrix4[i+1][0] - matrix4[i][0]))
        j += 1

    matrix4 = divided_difference(matrix4)
    return matrix4


def divided_difference(matrix4_result):
    i: int = 0
    for i in range(2, len(matrix4_result)):
        for j in range(3, i+2):
            # if a value already exists, skip
            if (j >= len(matrix4_result[i])) or matrix4_result[i][j] != 0:
                continue
            left: float = matrix4_result[i][j-1]
            diagonal_left: float = matrix4_result[i-1][j-1]
            denominator = matrix4_result[i][0] - matrix4_result[i-j+1][0]

            matrix4_result[i][j] = ((left - diagonal_left)/denominator)

    return matrix4_result


def cubic_spline_interpolation(x, y):
    size = len(x)
    matrix5 = np.zeros((size, size))

    # matrix A
    i: int = 0
    j: int = 0
    for i in range(0, size):
        if i == 0:
            matrix5[i][0] = 1
        elif i < len(x)-1:
            matrix5[i][i-1] = x[i-1 + 1] - x[i-1]
            matrix5[i][i] = 2 * ((x[i-1 + 1] - x[i-1]) + (x[i + 1] - x[i]))
            matrix5[i][i+1] = x[i + 1] - x[i]
        else:
            matrix5[i][i] = 1
    print(matrix5)
    print("")

    # vector b
    b = np.zeros(len(x))
    i = 0
    for i in range(1, len(x)-1):
        b[i] = ((3/(x[i+1] - x[i]))*(y[i+1] - y[i])) - ((3/(x[i-1+1] - x[i-1]))*(y[i] - y[i-1]))
    print(b)
    print("")

    # vector x
    # Ax = b
    x = np.zeros((len(x)))
    x = np.linalg.inv(matrix5).dot(b)
    print(x)
    print("")


if __name__ == "__main__":
    # 1.) Using Neville’s method, find the 2nd degree interpolating value for f(3.7)
    # for the following set of data
    #     x      f(x)
    #    3.6     1.675
    #    3.8     1.436
    #    3.9     1.318
    x_values1 = [3.6, 3.8, 3.9]
    y_values1 = [1.675, 1.436, 1.318]
    number_for_approximation1 = 3.7
    approximation1 = nevilles_method(x_values1, y_values1, number_for_approximation1)
    print(approximation1[len(x_values1)-1][len(x_values1)-1])
    print("")

    # 2.) Using Newton’s forward method, print out the polynomial approximations
    # for degrees 1, 2, and 3 using the following set of data
    #     x      f(x)
    #    7.2     23.5492
    #    7.4     25.3913
    #    7.5     26.8224
    #    7.6     27.4589
    x_values2 = [7.2, 7.4, 7.5, 7.6]
    y_values2 = [23.5492, 25.3913, 26.8224, 27.4589]
    table2 = divided_difference_table(x_values2, y_values2)
    problem2_answer_length = len(x_values2)-1
    problem2 = np.array(np.zeros(problem2_answer_length))
    for i in range(0, len(problem2)):
        problem2[i] = table2[i+1][i+1]
    np.set_printoptions(precision=15)
    print(np.array2string(problem2, separator=", "))
    print("")

    # 3.) Using the results from 3, approximate f(7.3)?
    number_for_approximation3 = 7.3
    approximation3 = find_approximation(table2, x_values2, number_for_approximation3)
    print('{0:.7f}'.format(approximation3).rstrip('0'))
    print("")

    # 4.)  Using the divided difference method, print out the
    # Hermite polynomial approximation matrix
    #     x      f(x)        f'(x)
    #    3.6     1.675      -1.195
    #    3.8     1.436      -1.188
    #    3.9     1.318      -1.182
    x_values4 = [3.6, 3.8, 3.9]
    y_values4 = [1.675, 1.436, 1.318]
    slopes4 = [-1.195, -1.188, -1.182]
    table4 = hermite(x_values4, y_values4, slopes4)
    np.set_printoptions(suppress=True, formatter={'float_kind': '{:f}'.format})
    print(table4)
    print("")

    # 5.) Using cubic spline interpolation, solve for the following using this set of data:
    #     x      f(x)
    #     2      3
    #     5      5
    #     8      7
    #    10      9
    # Find matrix A
    # Vector b
    # Vector x
    x_values5 = [2, 5, 8, 10]
    y_values5 = [3, 5, 7, 9]
    cubic_spline_interpolation(x_values5, y_values5)