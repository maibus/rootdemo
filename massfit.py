import numpy as np
from scipy.optimize import curve_fit
from rootget import loadhist
import matplotlib.pyplot as plt

filename = "Run2012BC_DoubleMuParked_Muons.root"
treename = "Events"
filters = ["nMuon", "/Muon_(eta|phi|pt)/"]
label = "mass [GeV]"

background_hist = np.array(loadhist(filename, treename, filters, 40, 60, 300, label=label, display=True, sample_size=200000))

def kexp(E, a, k, d):
    return a * np.exp(-k * (E - d))

E = np.linspace(40, 60, 300)
pbackground, pcov = curve_fit(kexp, E, background_hist, p0 = [50, 0.1, 40])
plt.scatter(E, background_hist, s=2)
plt.plot(E, kexp(E, *pbackground))
plt.show()

masshist = np.array(loadhist(filename, treename, filters, 75, 105, 1000, label=label, display=True, sample_size=2000000))

def breit_wigner(E, a, gamma, m):
    return a * gamma ** 2 / ((E - m) ** 2 + gamma ** 2 / 4)

E = np.linspace(75, 105, 1000)
masshist -= kexp(E, *pbackground)

hbar = 1.055 * 10 ** -34
e = 1.6 * 10 ** -19

p, pcov = curve_fit(breit_wigner, E, masshist, p0=[200, 2, 90])
print(p[-1], "Z0 mass in GeV")
print(hbar / (p[-2] * e * 10 ** 9), "rough estimate mean particle lifetime in seconds")

plt.scatter(E, masshist, s=2)
plt.plot(E, breit_wigner(E, *p))
plt.show()
