import numpy as np
import matplotlib.pyplot as plt  # type: ignore

data = [
    (50, 0.23874107),
    (100, 0.11912461),
    (500, 0.049564756),
    (1000, 0.06429757),
    (5000, 0.021661203),
]

N = np.array([N_ for N_, _ in data])
sigma = np.array([sigma_ for _, sigma_ in data])

p, cov = np.polyfit(1 / np.sqrt(N), sigma, deg=1, cov=True)

k_best_fit = p[0]
k_error = np.sqrt(cov[0][0])
k_true = np.sqrt(np.pi * (4 - np.pi))
print(f"Best fit {k_best_fit:.3f} +/- {k_error:.3f}")
print(f"True value: {k_true:.3f}")

fig, ax = plt.subplots()

ax.scatter(N, sigma, label="Data")
N_linspace = np.linspace(30, 6000)
ax.plot(
    N_linspace,
    k_best_fit / np.sqrt(N_linspace),
    label="Best-fit",
    color="red",
)
ax.fill_between(
    N_linspace,
    (k_best_fit - k_error) / np.sqrt(N_linspace),
    (k_best_fit + k_error) / np.sqrt(N_linspace),
    color="red",
    alpha=0.1,
)
ax.plot(
    N_linspace,
    k_true / np.sqrt(N_linspace),
    label="True",
    linestyle="--",
    color="green",
)

ax.set_yscale("log")
ax.set_xscale("log")
ax.legend()
ax.set_xlabel("Number of throws")
ax.set_ylabel("$ \\sigma $, $ \\pi $ value estimation error")


fig.savefig("out/ex3/pi-error-fit.png")
