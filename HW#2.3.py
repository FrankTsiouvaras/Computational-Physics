
# smf_cosmos.dat

# x_i+1 = xi - gamma del(f)

# a_i+1 = a_i - gamma[[f(a_i,b_i-1)-f(a_i-1,b_i-1)]/[a_i - a_i-1]

# b_i+1 = b_i - gamma[[f(a_i-1,b_i)-f(a_i-1,b_i-1)]/[b_i - b_i-1]

# chi^2 = sum_j=1^N = [f(x_j)-y_j]^2/sigma_jj^2

# aim: chi^2/N \approx 1



import matplotlib.pyplot as plt
import numpy as np
from numpy import empty, linspace
from pylab import contour, plot, xlabel, ylabel, minorticks_on, show



def f(a, b):
    return (a - 2)**2 + (b - 2)**2

gamma = 0.1

a = [0.0, 0.1]     # a_0, a_1
b = [0.0, 0.1]     # b_0, b_1
fs = [f(a[0], b[0]), f(a[1], b[1])]

for i in range(1, 30):
    dfda = (f(a[i], b[i-1]) - f(a[i-1], b[i-1])) / (a[i] - a[i-1])
    dfdb = (f(a[i-1], b[i]) - f(a[i-1], b[i-1])) / (b[i] - b[i-1])
    a.append(a[i] - gamma * dfda)
    b.append(b[i] - gamma * dfdb)
    fs.append(f(a[i+1], b[i+1]))

print(a[-1], b[-1])

plt.plot(fs)
plt.xlabel("step i")
plt.ylabel("f(a, b)")
plt.yscale("log")
plt.minorticks_on()
plt.show()

N = 100
avals = linspace(-1, 4, N)
bvals = linspace(-1, 4, N)
Z = empty([N, N])

for i in range(N):
    for j in range(N):
        Z[i, j] = f(avals[j], bvals[i])

contour(avals, bvals, Z)
plot(a, b, "ro-")
xlabel("a")
ylabel("b")
minorticks_on()
show()


#---------------------------------------------------------------------------------

data = np.loadtxt("smf_cosmos.dat")
logM = data[:, 0]
n    = data[:, 1]
err  = data[:, 2]

def schechter(logM, phi, logMstar, alpha):
    x = 10**(logM - logMstar)          # M / M*
    return phi * x**(alpha + 1) * np.exp(-x) * np.log(10)

def chi2(phi, logMstar, alpha):
    model = schechter(logM, phi, logMstar, alpha)
    return np.sum(((n - model) / err)**2)


gamma_phi   = 1e-9
gamma_logM  = 1e-5
gamma_alpha = 1e-5

eps = 1e-10

def fit(phi0, logM0, alpha0):
    phis   = [phi0, phi0 * 1.01]
    logMs  = [logM0, logM0 + 0.01]
    alphas = [alpha0, alpha0 + 0.01]
    chis = [chi2(phis[0], logMs[0], alphas[0]), chi2(phis[1], logMs[1], alphas[1])]

    i = 1
    converged = False
    while not converged and i < 50000:
        old = chi2(phis[i-1], logMs[i-1], alphas[i-1])
        dphi   = (chi2(phis[i], logMs[i-1], alphas[i-1]) - old) / (phis[i] - phis[i-1])
        dlogM  = (chi2(phis[i-1], logMs[i], alphas[i-1]) - old) / (logMs[i] - logMs[i-1])
        dalpha = (chi2(phis[i-1], logMs[i-1], alphas[i]) - old) / (alphas[i] - alphas[i-1])

        phis.append(phis[i] - gamma_phi * dphi)
        logMs.append(logMs[i] - gamma_logM * dlogM)
        alphas.append(alphas[i] - gamma_alpha * dalpha)
        chis.append(chi2(phis[i+1], logMs[i+1], alphas[i+1]))

        converged = (abs(phis[i+1] - phis[i])     < eps * abs(phis[i]) or
                     abs(logMs[i+1] - logMs[i])   < eps * abs(logMs[i]) or
                     abs(alphas[i+1] - alphas[i]) < eps * abs(alphas[i]))
        i = i + 1

    return phis[-1], logMs[-1], alphas[-1], chis

starts = [(0.005, 10.5, -1.2),
          (0.001, 10.0, -1.5),
          (0.004, 11.2, -0.9)]

for phi0, logM0, alpha0 in starts:
    phi, logMstar, alpha, chis = fit(phi0, logM0, alpha0)
    print(phi0, logM0, alpha0, "->", phi, logMstar, alpha, chis[-1])
    plt.plot(chis, label=f"start {phi0}, {logM0}, {alpha0}")

plt.xlabel("step i")
plt.ylabel(r"$\chi^2$")
plt.yscale("log")
plt.minorticks_on()
plt.legend()
plt.show()



logM_fine = linspace(logM[0], logM[-1], 200)
model = schechter(logM_fine, phi, logMstar, alpha)

plt.errorbar(10**logM, n, yerr=err, fmt="o", label="COSMOS data")
plt.plot(10**logM_fine, model, label="best-fit Schechter")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$M_{gal}$")
plt.ylabel(r"$n(M_{gal})$ [1/dex/(Mpc/h)$^3$]")
plt.minorticks_on()
plt.legend()
plt.show()





print("chi^2/N =", chis[-1] / len(n))








