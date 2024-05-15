import math
import time

import random as rng


TWO_PI = 2.0 * math.pi


def sample_circle_analytic() -> tuple[float, float]:
    ampl = math.sqrt(rng.uniform(0, 1))
    phase = rng.uniform(0, TWO_PI)
    return (ampl * math.cos(phase), ampl * math.sin(phase))


def sample_circle_rejection() -> tuple[float, float]:
    while True:
        x = rng.uniform(-1, 1)
        y = rng.uniform(-1, 1)
        if x**2 + y**2 <= 1.0:
            return (x, y)


def ex2():
    sample_size = 10_000_000
    print(f"\nSampling {sample_size} points inside the unit circle")
    print("Analytic method...")
    start = time.time()
    for _ in range(sample_size):
        sample_circle_analytic()
    print(f"... done in {time.time() - start:.4} sec")
    print("Rejection method...")
    start = time.time()
    for _ in range(sample_size):
        sample_circle_rejection()
    print(f"... done in {time.time() - start:.4} sec")


if __name__ == "__main__":
    ex2()
