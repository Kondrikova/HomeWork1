import sys
import math
import numpy
import numpy.linalg
import pylab

h = math.pi / 5
tau = 0.25
error = 0
n = int(math.pi/h)
m = int(10/tau)

u1 = numpy.zeros((m+1,n+1))
f = numpy.zeros((m+1,n+1))
err = numpy.zeros((m+1,n+1))

def exactSolution(x, t):
    return numpy.sin(x) + numpy.log(t ** 2 + 1)

def g(x, t):
    return math.sin(x)+(2 * t) / (t ** 2 + 1)

def explicit(u):
    for j in range(0, m):
        for i in range(1, n):
            u[j+1][i] = (tau/h**2)*(u[j][i+1]-2*u[j][i]+u[j][i-1])+u[j][i]+tau*g(i*h,j*tau)
    return u

def implicit(u):
    a = b = -tau
    c = h ** 2 + 2 * tau
    for j in range(0, m):
        alfa = [0] * (n+1)
        beta = [0] * (n+1)

        alfa[1] = 0
        beta[1] = u[j+1][0]

        for i in range(1, n):
            f = (h ** 2) * u[j][i] + tau * (h ** 2) * g(i * h, (j + 1) * tau)
            alfa[i + 1] = (-b) / (a * alfa[i] + c)
            beta[i + 1] = (f - a * beta[i]) / (a * alfa[i] + c)
        for i in reversed(range(1, n)):
            u[j + 1][i] = u[j + 1][i + 1] * alfa[i + 1] + beta[i + 1]
    return u

def KNmethod(u):
    a = b = -tau
    c = 2*(h ** 2) + 2 * tau
    for j in range(0, m):
        alfa = [0] * (n+1)
        beta = [0] * (n+1)

        alfa[1] = 0
        beta[1] = u[j+1][0]

        for i in range(1, n):
            f = tau*u[j][i-1] + (2* (h**2) - 2*tau)*u[j][i] + tau* u[j][i+1] +2*tau * (h ** 2) * g(i * h, j * tau)
            alfa[i + 1] = (-b) / (a * alfa[i] + c)
            beta[i + 1] = (f - a * beta[i]) / (a * alfa[i] + c)
        for i in reversed(range(1, n)):
            u[j + 1][i] = u[j + 1][i + 1] * alfa[i + 1] + beta[i + 1]
    return u

def calcErr(u,f):
    for i in range(0, n + 1):
        for j in range(0, m + 1):
            err[j][i] = f[j][i]-u[j][i]
    error = numpy.absolute(err).max()
    return error

def makeData():
    x = numpy.arange(0, math.pi+h, h)
    t = numpy.arange(0, 10+tau, tau)
    xgrid, tgrid = numpy.meshgrid(x, t)
    f = exactSolution(xgrid, tgrid)
    return f

def main():
    u = numpy.zeros((m + 1, n + 1))

    for i in range(0, n + 1):
        u[0][i] = math.sin(i * h)

    for j in range(1, m + 1):
        u[j][0] = math.log((j * tau) ** 2 + 1)

    for j in range(1, m + 1):
        u[j][n] = math.log((j * tau) ** 2 + 1)

    uu = u.copy()
    uuu = u.copy()
    u1 = explicit(u)
    f = makeData()
    print calcErr(u1,f)
    u2 = implicit(uu)
    print calcErr(u2, f)
    u3 = KNmethod(uuu)
    print calcErr(u3,f)

if __name__ == '__main__':
    main()

