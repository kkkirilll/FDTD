
import requests
import pathlib
import json
import matplotlib.pyplot as plt
import numpy as np

from scipy.constants import pi, speed_of_light
from scipy.special import spherical_jn as jn
from scipy.special import spherical_yn as yn


def hn(n, x):
    return jn(n, x) + 1j * yn(n, x)


def bn(n, x):
    return (x * jn(n - 1, x) - n * jn(n, x)) / (x * hn(n - 1, x) - n * hn(n, x))


def an(n, x):
    return jn(n, x) / hn(n, x)


link="https://jenyay.net/uploads/Student/Modelling/task_02_01.txt"
variant=1

content_link=requests.get(link)
content_string=content_link.text
date=content_string.split("\n")[variant-1].split()

D=float(date[1][date[1].find("=")+1:-1])
fmin=float(date[2][date[2].find("=")+1:-1])
fmax=float(date[3][date[3].find("=")+1:])
fstep = 1e6

r = D / 2
freq = np.arange(fmin, fmax, fstep)
lambd = speed_of_light / freq
k = 2 * pi / lambd

arr_sum = [((-1) ** n) * (n + 0.5) * (bn(n, k * r) - an(n, k * r)) for n in range(1, 20)]
summa = np.sum(arr_sum, axis=0)
rcs = (lambd ** 2) / pi * (np.abs(summa) ** 2)

result = {
    "freq": freq.tolist(),
    "lambda": lambd.tolist(),
    "rcs": rcs.tolist(),
}

dir_ = pathlib.Path("results")
dir_.mkdir(exist_ok=True)
file = dir_ / "task2.json"
with file.open("w") as f:
    json.dump(result, f, indent=4)

plt.plot(freq / 10e6, rcs)
plt.xlabel("$f, МГц$")
plt.ylabel(r"$\sigma, м^2$")
plt.grid()
plt.savefig("results/task2.png")
plt.show()
