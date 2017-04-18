"""
python -m science.miscellaneous.coevol.src.python.coevolution
"""

import numpy
import pandas
import scipy.integrate
import matplotlib.pyplot as plt


class EcoCoEvol:


    def __init__(self, I, n, rho, lambd, eta, gamma, M, beta, r, c, h):

        # Paramters
        self.n = n
        self.I = I

        if numpy.trace(I) != 0:
            print "Error: No self connections between traits allowd."
            raise Exception()

        I1, I2 = numpy.shape(I)
        if I1 != I2 or I1 % n != 0:
            raise Exception()
        self.k = I1 / n


        self.rho = rho
        self.lambd = lambd
        self.eta = eta
        self.gamma = gamma
        self.M = M
        self.beta = beta
        self.r = r
        self.c = c
        self.h = h

        # State
        self.N_ = None
        self.V_ = None

    ##############################################################################
    # INTEGRATION
    ##############################################################################

    def dX (self, X, t):
        """
        X :=    The state vector which is a concatenation of N and V
                X = [N_1, N_2, N_n, T_1_1, .. T_1_k, .. T_n_1, .. T_n_k]
        """
        self.N_ = X[0:self.n]
        self.V_ = X[self.n:len(X)]

        dN_dt_ = self.dN_dt()

        dV_dt_ = self.dV_dt()
        # Extract abundance and trait vector

        dX_dt = dN_dt_ + dV_dt_
        return dX_dt



    def integrate(self, N0, V0, T):
        """
        Integrates the system using the given initial conditions.
        :param N0:
        :param V0:
        :param T:
        :return:
        """

        X0 = N0 + V0
        sol = scipy.integrate.odeint(self.dX, X0, T, mxstep=50000000)

        return pandas.DataFrame(sol, index=T, columns=self.getNColNames() + self.getVColNames())



    ##############################################################################
    # HELPER FUNCTIONS TO ACCESS INTERACTION MATRICES
    ##############################################################################

    def A(self, i, j, p, q):
        """
        Returns the effect of trait q of species j on trait p of species i
        :param i:
        :param j:
        :param p:
        :param q:
        :return:
        """

        return self.I[i*self.k+p,j*self.k+q]

    def B(self, i, p, q):

        return self.I[i*self.k+p,i*self.k+q]

    def N(self, i):
        return self.N_[i]

    def V(self, i, p = None):
        if p is None:
            return self.V_[i*self.k:i*self.k+self.k]

        return self.V_[i*self.k+p]


    ##############################################################################
    # ECOLOGICAL DYNAMICS
    ##############################################################################


    def g(self, i, j, p):
        """
        The net effect of the traits of species j on the trait p of
        species i.

        :param V_i_p:
        :param V_j:
        :return:
        """

        V_i_p = self.V(i,p)
        V_j = self.V(j)

        sum = 0.0
        for q in range(0, self.k):
            sum += self.F_match(i, j, p, q) * self.A(i,j,p,q)

        return sum

    def G(self, i, j):
        """
        The net effects of the traits of species j on the traits of
        species i.

        :param V_i:
        :param V_j:
        :param k:
        :return:
        """

        sum = 0.0
        for p in range(0, self.k):
            sum += self.g(i, j, p)

        return sum

    def dN_dt(self):

        # Extract abundance and trait vector

        dN_dt_ = []

        # Compute growth of each species (1)
        for i in range(0, self.n):

            # Compute effect of interactions on growth (2)
            sum = 0.0
            N_i = self.N(i)
            for j in range(0, self.n):
                if i == j:
                    continue

                N_j = self.N(j)
                sum += self.G(i, j)* N_j / (1.0 + self.h * N_j)

            dN_i_dt = N_i * (self.r[i] - self.beta*N_i + self.gamma*self.Psi(i) + sum)
            dN_dt_.append(dN_i_dt)

        return dN_dt_

    ##############################################################################
    # INTERNAL FITNESS FUNCTION
    ##############################################################################

    def mu(self, i, p):
        """

        :param i:
        :param p:
        :return:
        """
        sum = 0.0
        for q in range(0, self.k):
            B_i_pq = self.B(i, p, q)
            V_i_q = self.V(i, q)
            sum += B_i_pq * V_i_q


        return self.c[i][p] + sum

    def phi(self, i, p):
        """
        The fitness contribution (internal) of trait p on of species i.
        :param i:
        :param p:
        :return:
        """

        V_i_p = self.V(i,p)
        return self.M / (1.0 + self.rho * numpy.power(V_i_p - self.mu(i, p), 2.0))

    def Psi(self, i):
        """
        The internal fitness of species i.
        :param i:
        :return:
        """
        sum = 0.0
        for p in range (0, self.k):
            sum += self.phi(i, p)

        return 1.0 / float(self.k) * sum

    ##############################################################################
    # TRAIT EVOLUTION
    ##############################################################################


    def dV_dt (self):

        dV_dt_ = []
        for i in range(0, self.n):
            for p in range(0, self.k):
                dV_dt_.append(self.dV_i_p_dt(i, p))

        return dV_dt_

    def dV_i_p_dt (self, i, p):
        """
        Temporal change of trait p of species i due to internal and external forces.
        :param i:
        :param p:
        :return:
        """

        N_i = self.N(i)
        sum = 0.0
        for j in range(0, self.n):
            if i == j:
                continue
            N_j = self.N(j)
            sum += N_j * self.Pg_PV_p(i, j, p)

        dV_i_p_dt_ = self.lambd * N_i * (sum + self.eta * self.PPsi_PV_p(i, p))
        return dV_i_p_dt_

    def PPsi_PV_p (self, i, p):
        """
        Partial derivative of Psi with respect to trait p of species i (internal)
        :param p:
        :return:
        """

        V_i_p = self.V(i, p)

        denominator = numpy.power(1.0 + self.rho * numpy.power(V_i_p - self.mu(i, p), 2.0), 2.0)
        firstTerm = (-V_i_p + self.mu(i, p)) / denominator

        sum = 0.0

        for q in range(0, self.k):
            if p == q:
                continue
            B_i_qp = self.B(i, q, p)
            V_i_q = self.V(i, q)

            denominator = numpy.power(1.0 + self.rho * numpy.power(V_i_q - self.mu(i, q), 2.0), 2.0)
            sum += (B_i_qp * (V_i_q - self.mu(i, q)) / denominator)

        return 2.0 * self.rho * self.M * (1.0 / self.k) * (firstTerm + sum)

    def Pg_PV_p (self, i, j, p):
        """

        :param i:
        :param j:
        :param p:
        :return:
        """

        sum = 0.0
        for q in range(0, self.k):
            A_ij_pq = self.A(i, j, p, q)

            sum += A_ij_pq * self.PF_match_PV_p(i, j, p, q)

        return -sum



    ##############################################################################
    # TRAIT MATCH FUNCTIONS
    ##############################################################################

    def F_match (self, i, j, p, q):
        """
        Trait matching function returning the strength of the effect [0,1]
        of trait q of species j on trait p of specie i
        :return:
        """
        V_i_p = self.V(i,p)
        V_j_q = self.V(j,q)



        return numpy.exp(-numpy.power(V_i_p-V_j_q, 2.0))

    def PF_match_PV_p (self, i, j, p, q):
        """
        Partial derivative of trait match function with respect to trait p of species i.
        :param i:
        :param j:
        :param p:
        :param q:
        :return:
        """

        V_i_p = self.V(i,p)
        V_j_q = self.V(j,q)

        return (V_i_p-V_j_q) * self.F_match(i, j, p, q)


    ##############################################################################
    # OTHER HELPER FUNCTIONS
    ##############################################################################


    def getNColNames(self):
        N_cols = []
        for i in range(0, self.n):
            N_cols.append("N-" + str(i))
        return N_cols

    def getVColNames(self, i=None):

        V_cols = []
        if i is None:
            for i in range(0, self.n):
                for p in range(0, self.k):
                    V_cols.append("V-" + str(i) + "_" + str(p))
        else:
            for p in range(0, self.k):
                V_cols.append("V-" + str(i) + "_" + str(p))

        return V_cols

    def printC (self, i, j, p, q):

        print "i=" + str(i) + " j=" + str(j) + " p=" + str(p) + " q=" + str(q)



if __name__ == "__main__":



    numpy.set_printoptions(linewidth=200)

    n_ = 2
    k_ = 2
    N0 = [0.5, 0.6]
    V0 = [0.4, 0.5,0.6,0.7]
    rho_ = 1.0
    lambda_ = 1.0
    eta_ = 1.0
    gamma_ = 0.0
    M_ = 1.0
    r_ = [0.5, 0.5]
    beta_ = 1.0
    c_ = [[1.0, 1.0],
          [1.0, 1.0],
        ]
    h_ = 0.5

    AB = numpy.matrix([[0.0, 0.0, 0.0, 1.0],
                       [1.0, 0.0, 0.0, 1.0],
                       [1.0, 0.0, 0.0, 0.0],
                       [1.0, 1.0, 1.0, 0.0]
                       ])


    ecv = EcoCoEvol(AB, n_, rho_, lambda_, eta_, gamma_, M_, beta_, r_, c_, h_)
    T = numpy.arange(0, 20, 0.1)


    df = ecv.integrate(N0, V0, T)

    fig, ax = plt.subplots(nrows=n_+1, ncols=1, sharex=False, sharey=False,squeeze=False, dpi=120, figsize=(8, 8))
    df[ecv.getNColNames()].plot(ax=ax[0,0], legend=None)
    for i in range(0, n_):
        df[ecv.getVColNames(i)].plot(ax=ax[i+1,0], legend=None)

    plt.show()

