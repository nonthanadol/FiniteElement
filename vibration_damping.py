import numpy as np
import matplotlib.pyplot as plt

m = 3
k = 40000
c = 100
x0 = 0
x0dot = 12
Cc = 2 * np.sqrt(m * k)
Z = c / Cc
wn = np.sqrt(k / m)
wd = wn * (np.sqrt(1 - Z ** 2))
t = np.arange(0.0001, 1, 0.001)


def Underdamped(x0, x0dot, Z, wn, wd):
    X = np.sqrt((x0) ** 2 + ((x0dot + Z * wn * x0) / wd) ** 2)
    Angle = np.arctan((x0 * wd) / (x0dot + Z * wn * x0))
    x = X * np.e ** (-Z * wn * t) * np.sin(wd * t + Angle)
    return x


def Overdamped(x0, x0dot, Z, wn, wd):
    s1 = -Z * wn + wn * np.sqrt(Z ** 2 - 1)
    s2 = -Z * wn - wn * np.sqrt(Z ** 2 - 1)
    c1 = (x0dot - s2 * x0) / (s1 - s2)
    c2 = (x0dot - s1 * x0) / (s2 - s1)
    x = c1 * np.e ** (s1 * t) + c2 * np.e ** (s2 * t)
    return x


def Criticallydamped(x0, x0dot, Z, wn, wd):
    x = (x0 + (x0dot + wn * x0) * t) * np.e ** (-wn * t)
    return x


if (0 < Z < 1):
    print("damping ratio is {0:.4f}".format(Z))
    print("System is Underdamped ")
    x = Underdamped(x0, x0dot, Z, wn, wd)  # เรียกใช้function
elif (Z == 1):
    print("damping ratio is {0:.4f}".format(Z))
    print("System is Criticallydamped ")
    x = Criticallydamped(x0, x0dot, Z, wn, wd)
else:
    print("damping ratio is {0:.4f}".format(Z))
    print("System is Overdamped ")
    x = Overdamped(x0, x0dot, Z, wn, wd)

plt.plot(t, x)
plt.show()