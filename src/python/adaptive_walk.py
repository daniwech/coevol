"""
Visualization of fitness landscape (internal fitness) and simulation
of an adaptive walk on the landscape.

python -m science.miscellaneous.experiments.continuous_coevolution.adaptive_walk
"""


import numpy
import scipy.integrate
import scipy.misc
import matplotlib.pyplot as plt
import pandas
from ....plotting import plotting_helpers
from ....plotting.styles import basic_style
from mpl_toolkits.mplot3d import Axes3D



MAX = 1.0


###################################################
# FITNESS FUNCTION (Using absolute function)
###################################################

def my_i(A, X, i, c_i):

    u_i = c_i[i] + (A[i,:] * X)[0,0]
    return u_i

def gamma_i(A, X, i, c_i):
    x_i = X[i][0]
    gamma_i = x_i - my_i(A, X, i, c_i)
    return gamma_i

def sigma_i(A, X, i, c_i):
    """
    Returns the fitness contribution of gene i.
    :param A:
    :param X:
    :param i:
    :return:
    """
    sigma_i = MAX / (1.0 + numpy.fabs(gamma_i(A, X, i, c_i)))
    return sigma_i

def f(A, X, c_i):
    """
    Fitness function. For each X returns a value in range [0, MAX]
    :param A:
    :param X:
    :return:
    """
    sum = 0.0
    N = numpy.shape(X)[0]

    for i in range(0, N):
        sum += sigma_i(A,X,i, c_i)

    return sum / float(N)

###################################################
# FITNESS FUNCTION (Square)
###################################################

FALLOFF = 0.1

def sqr_sigma_i(A, X, i, c_i):

    sigma_i = MAX / (1.0 + FALLOFF * numpy.power(gamma_i(A, X, i, c_i), 2.0))
    return sigma_i


def sqr_f (A, X, c_i):

    """
    Fitness function (using square). For each X returns a value in range [0, MAX]
    :param A:
    :param X:
    :return:
    """
    sum = 0.0
    N = numpy.shape(X)[0]

    for i in range(0, N):
        sum += sqr_sigma_i(A,X,i, c_i)

    return sum / float(N)



###################################################
# FITNESS FUNCTION GRADIENT (Absolute)
###################################################

def d_sigma_i_i(A, X, i, c_i):
    """
    The derivative of sigma_i with respect to i.
    :param A:
    :param X:
    :param i:
    :return:
    """
    gamma_i_ = gamma_i(A, X, i, c_i)
    dx_i_sigma = -MAX / numpy.power(1.0 + numpy.fabs(gamma_i_), 2.0) * numpy.sign(gamma_i_) * (1.0 - A[i,i])
    return dx_i_sigma


def d_sigma_j_i(A, X, i, j, c_i):
    """
    The derivative of sigma_j, with respect to i.
    :param A:
    :param X:
    :param i:
    :param j:
    :return:
    """
    gamma_i_ = gamma_i(A, X, j, c_i)
    dx_i_sigma = -MAX / numpy.power(1.0 + numpy.fabs(gamma_i_), 2.0) * numpy.sign(gamma_i_) * (-A[j,i])
    return dx_i_sigma


def df_x_i(A, X, i, c_i):
    """
    Change of fitness with respect to i.
    :param A:
    :param X:
    :param i:
    :return:
    """
    sum = 0.0
    N = numpy.shape(X)[0]
    for j in range(0, N):

        if i == j:
            sum += d_sigma_i_i(A,X,i, c_i)
        else:
            sum += d_sigma_j_i(A,X,i,j, c_i)


    return sum / float(N)


###################################################
# FITNESS FUNCTION GRADIENT (Absolute)
###################################################



def sqr_d_sigma_i_i(A, X, i, c_i):
    """
    The derivative of sigma_i with respect to i.
    :param A:
    :param X:
    :param i:
    :return:
    """
    gamma_i_ = gamma_i(A, X, i, c_i)
    x_i = X[i][0]
    tmp = 2.0*x_i - 2.0*c_i[i] - 2.0 * (A[i,:] * X)[0,0]
    dx_i_sigma = -MAX / numpy.power(1.0 + FALLOFF * numpy.power(gamma_i_, 2.0), 2.0) * FALLOFF * tmp
    return dx_i_sigma


def sqr_d_sigma_j_i(A, X, i, j, c_i):
    """
    The derivative of sigma_j, with respect to i.
    :param A:
    :param X:
    :param i:
    :param j:
    :return:
    """

    gamma_i_ = gamma_i(A, X, j, c_i)
    x_j = X[j][0]
    tmp = -2.0*x_j*A[j,i] + 2.0*c_i[j]*A[j,i] + 2.0*A[j,i] * (A[j,:] * X)[0,0]
    dx_i_sigma = -MAX / numpy.power(1.0 + FALLOFF * numpy.power(gamma_i_, 2.0), 2.0) * FALLOFF * tmp
    return dx_i_sigma


def sqr_df_x_i(A, X, i, c_i):
    """
    Change of fitness with respect to i.

    :param A:
    :param X:
    :param i:
    :return:
    """
    sum = 0.0
    N = numpy.shape(X)[0]
    for j in range(0, N):

        if i == j:
            sum += sqr_d_sigma_i_i(A,X,i, c_i)
        else:
            sum += sqr_d_sigma_j_i(A,X,i,j, c_i)

    return sum / float(N)



###################################################
# GRADIENT ASCENT
###################################################

def dX (X_, t, A, beta, c_i):

    N = len(X_)
    X = numpy.ndarray(shape=(N, 1))
    for i in range(0, N):
        X[i] = [X_[i]]

    dX_ = list()
    for i in range(0, N):
        dX_.append(df_x_i(A, X, i, c_i)[0] * beta)

    return dX_

def sqr_dX (X_, t, A, beta, c_i):

    N = len(X_)
    X = numpy.ndarray(shape=(N, 1))
    for i in range(0, N):
        X[i] = [X_[i]]

    dX_ = list()
    for i in range(0, N):
        dX_.append(sqr_df_x_i(A, X, i, c_i) * beta)

    return dX_





###################################################
# START TESTS
###################################################

if __name__ == "__main__":

    #
    # Plot fitness functions and derivatives.
    # Using absolute fitness function
    #
    if True:

        c_i = [4.0, -4.0, -4.0]
        A = [[0.0, 0.5, 0.6],
             [0.1, 0.0, 0.7],
             [0.1, 0.3, 0.0]]
        A = numpy.matrix(A)

        X = [[2.0],
             [1.0],
             [4.0]]

        i = 0
        j = 1
        x_i_range = numpy.arange(-60, 60, 0.1)

        # Absolute fitness function
        sigma_xi = []       # Sigma_i for different xi
        sigma_xj = []       # Sigma_j for different xi
        f_xi = []           # Fitness for different xi
        d_f_xi = []         # Change of fitness for different xi
        d_sigma_i_xi = []   # Change of sigma_i for different xi
        d_sigma_j_xi = []   # Change of sigma_j for different xi

        # Square fitness function
        sqr_sigma_xi = []       # Sigma_i for different xi
        sqr_sigma_xj = []       # Sigma_j for different xi
        sqr_f_xi = []           # Fitness for different xi
        sqr_d_f_xi = []         # Change of fitness for different xi
        sqr_d_sigma_i_xi = []   # Change of sigma_i for different xi
        sqr_d_sigma_j_xi = []   # Change of sigma_j for different xi

        for x_i in x_i_range:

            X[i] = [x_i]
            # Absolute fitness function
            sigma_xi.append(sigma_i(A, X, i, c_i))
            sigma_xj.append(sigma_i(A, X, j, c_i))
            f_xi.append(f(A, X, c_i))
            d_f_xi.append(df_x_i(A,X,i, c_i))
            d_sigma_i_xi.append(d_sigma_i_i(A, X, i, c_i))
            d_sigma_j_xi.append(d_sigma_j_i(A, X, i, j, c_i))

            # Square fitness function
            sqr_sigma_xi.append(sqr_sigma_i(A, X, i, c_i))
            sqr_sigma_xj.append(sqr_sigma_i(A, X, j, c_i))
            sqr_f_xi.append(sqr_f(A, X, c_i))
            sqr_d_f_xi.append(sqr_df_x_i(A,X,i, c_i))
            sqr_d_sigma_i_xi.append(sqr_d_sigma_i_i(A, X, i, c_i))
            sqr_d_sigma_j_xi.append(sqr_d_sigma_j_i(A, X, i, j, c_i))



        fig, ax = plt.subplots(nrows=3, ncols=2, sharex=False, sharey=False,squeeze=False, dpi=120, figsize=(8, 8))

        # sigma_i and d_sigma_i
        ax[0,0].plot(x_i_range, sigma_xi, "r-", label="2")
        ax[0,0].plot(x_i_range, d_sigma_i_xi, "b-", label="2")
        ax[0,0].set_xlabel("$x_i$", labelpad=15, fontsize=16)
        plotting_helpers.setLegend(ax[0,0], ["$ \sigma_i(x_i) $",
                                             "$\\frac{\partial \sigma_i(x)}{\partial x_i} $"],
                                   fontsize=10,
                                   location="upper right")

        # sigma_i and d_sigma_i
        ax[1,0].plot(x_i_range, d_sigma_j_xi, "r-", label="2")
        ax[1,0].plot(x_i_range, sigma_xj, "b-", label="2")
        ax[1,0].set_xlabel("$x_i$", labelpad=15, fontsize=16)
        plotting_helpers.setLegend(ax[1,0], ["$ \sigma_j(x_i) $",
                                            "$\\frac{\partial \sigma_j(x)}{\partial x_i} $"],
                                   fontsize=10,
                                   location="upper right")

        # f and d_f_i
        ax[2,0].plot(x_i_range, f_xi, "r-")
        ax[2,0].plot(x_i_range, d_f_xi, "g-")
        ax[2,0].set_xlabel("$x_i$", labelpad=15, fontsize=16)
        plotting_helpers.setLegend(ax[2,0], ["$ f(x_i) $",
                                             "$\\frac{\partial f(x)}{\partial x_i} $"],
                                   fontsize=10,
                                   location="upper left")




        # sigma_i and d_sigma_i
        ax[0,1].plot(x_i_range, sqr_sigma_xi, "r-", label="2")
        ax[0,1].plot(x_i_range, sqr_d_sigma_i_xi, "b-", label="2")
        ax[0,1].set_xlabel("$x_i$", labelpad=15, fontsize=16)
        plotting_helpers.setLegend(ax[0,0], ["$ \sigma_i(x_i) $",
                                             "$\\frac{\partial \sigma_i(x)}{\partial x_i} $"],
                                   fontsize=10,
                                   location="upper right")


         # sigma_i and d_sigma_i
        ax[1,1].plot(x_i_range, sqr_d_sigma_j_xi, "r-", label="2")
        ax[1,1].plot(x_i_range, sqr_sigma_xj, "b-", label="2")
        ax[1,1].set_xlabel("$x_i$", labelpad=15, fontsize=16)
        plotting_helpers.setLegend(ax[1,0], ["$ \sigma_j(x_i) $",
                                            "$\\frac{\partial \sigma_j(x)}{\partial x_i} $"],
                                   fontsize=10,
                                   location="upper right")

        # f and d_f_i
        ax[2,1].plot(x_i_range, sqr_f_xi, "r-")
        ax[2,1].plot(x_i_range, sqr_d_f_xi, "g-")
        ax[2,1].set_xlabel("$x_i$", labelpad=15, fontsize=16)
        plotting_helpers.setLegend(ax[2,0], ["$ f(x_i) $",
                                             "$\\frac{\partial f(x)}{\partial x_i} $"],
                                   fontsize=10,
                                   location="upper left")

        fig.tight_layout()
        plt.show()


    #
    # Plot sqr fitness landscape
    #
    if True:

        c_i = [4.0, -0.0, -10.0, 0.0]
        A = [[0.0, 0.1],
             [0.0, 0.0]]
        A = numpy.matrix(A)

        X = [[2.0],
             [1.0]]

        A = [[0.0, 0.1, -0.2, 0.0],
             [0.0, 0.0, -0.1, 0.0],
             [0.2, 0.1, 0.0, 0.0],
             [0.3, 0.1, -2.0, 0.0]]
        A = numpy.matrix(A)

        X = [[2.0],
             [1.0],
             [1.0],
             [1.0]]

        X_0_range = numpy.arange(-50, 50, 5.0)
        X_1_range = numpy.arange(-50, 50, 5.0)
        X_0, X_1 = numpy.meshgrid(X_0_range, X_1_range)

        f_ = numpy.ndarray(shape=(len(X_0), len(X_1)))
        f_i_ = numpy.ndarray(shape=(len(X_0), len(X_1)))
        df_x0 = numpy.ndarray(shape=(len(X_0), len(X_1)))

        dfFitness = pandas.DataFrame(index=X_0_range, columns=X_1_range)
        i = 0
        for x_0 in X_0_range:
            j = 0
            for x_1 in X_1_range:

                X[0] = [x_0]
                X[1] = [x_1]
                f_ij = sqr_f(A, X, c_i)
                f_[i,j] = f_ij
                df_x0[i,j] = sqr_df_x_i(A, X, 0, c_i)
                f_i_[i,j] = sqr_sigma_i(A, X, 0, c_i)

                dfFitness.ix[x_0][x_1] = sqr_sigma_i(A, X, 0, c_i)
                j += 1

            i += 1

        """
        fig = plt.figure(figsize=(12, 12))
        ax = fig.gca(projection='3d')
        ax.plot_surface(X_0, X_1, f_, cmap=plt.get_cmap('rainbow'),rstride=1, cstride=1)
        ax.set_xlabel("$X_1$", fontsize=16)
        ax.set_ylabel("$X_0$", fontsize=16)
        plt.show()


        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot_surface(X_0, X_1, df_x0, cmap=plt.get_cmap('rainbow'),rstride=1, cstride=1)
        plt.show()
        """

        beta = 2.2 # Mutation Rate
        X0 = [5.0, 5.0, 5.0, 5.0]

        t = numpy.arange(0, 200, 0.1)
        sol = scipy.integrate.odeint(sqr_dX, X0, t, args=(A, beta, c_i))
        df = pandas.DataFrame(sol, index=t)
        fig, ax = plt.subplots(nrows=3, ncols=1, sharex=False, sharey=False,squeeze=False, dpi=120, figsize=(6, 8))
        df.plot(ax=ax[0,0])

        x0 = list()
        x1 = list()
        f = list()
        for t in t:
            X_t = numpy.array(df.ix[t].tolist()).reshape(4, 1)
            print X_t
            f_t = sqr_f(A, X_t, c_i)
            f.append(f_t)

        plt.show()
        print df


    #
    # Gradient ascent
    #
    if True:

        c_i = [4.0, -5.0, -10.0]

        A = [[0.0, -2.0, 0.5],
             [1.0, 0.0, 0.7],
             [0.0, 0.3, 0.0]]
        A = numpy.matrix(A)

        beta = 2.2 # Mutation Rate
        X0 = [2.0, 1.0, 4.0]



        c_i = [4.0, -40.0, -10.0, 3.0]
        A = [[2.0, 0.1, -0.2, 0.0],
             [3.0, 0.0, -0.1, 0.0],
             [0.2, 0.1, 0.0, 0.0],
             [0.3, 0.1, -2.0, 0.0]]
        A = numpy.matrix(A)

        X0 = [200.0, 1.0, 4.0, 5.0]


        t = numpy.arange(0, 10000, 1.0)
        sol = scipy.integrate.odeint(sqr_dX, X0, t, args=(A, beta, c_i))

        df = pandas.DataFrame(sol, index=t)

        f_t = list()
        for t in t:
            X_t = numpy.array(df.ix[t].tolist()).reshape(len(X0), 1)
            f_t.append(sqr_f(A, X_t, c_i))

        df['F'] = f_t

        fig, ax = plt.subplots(nrows=2, ncols=1, sharex=False, sharey=False,squeeze=False, dpi=120, figsize=(6, 8))
        df[range(0, len(X0))].plot(ax=ax[0,0])

        df['F'].plot(ax=ax[1,0])
        plt.show()





