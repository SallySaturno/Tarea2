# Tutorial: Cuadratura de Gauss-Legendre para la integración numérica

Este código utiliza la cuadratura de Gauss-Legendre para aproximar la integral de una función en un intervalo arbitrario. La función a integrar se define y luego se aproxima usando los puntos y pesos de Gauss-Legendre.

## Ejemplo de uso

A continuación se muestra un ejemplo de uso de las funciones definidas para calcular la integral de la función `f(x) = x^6 - x^2 sin(2x)` en el intervalo [1, 3].

```python

import matplotlib.pyplot as plt
import numpy as np

def func(x):
    return x**6 - x**2 * np.sin(2*x)

def gaussxw(N):
    a = np.linspace(3, 4 * (N - 1), N) / ((4 * N) + 2)
    x = np.cos(np.pi * a + 1 / (8 * N * N * np.tan(a)))

    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = np.ones(N, dtype=float)
        p1 = np.copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
        dp = (N + 1) * (p0 - x * p1) / (1 - x * x)
        dx = p1 / dp
        x -= dx
        delta = np.max(np.abs(dx))

    w = 2 * (N + 1) * (N + 1) / (N * N * (1 - x * x) * dp * dp)
    return x, w

def gaussxwab(a, b, x, w):
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def integrar(function, wk, xk):
    return np.sum(wk * function(xk))

val = gaussxw(4)
val_escal = gaussxwab(1, 3, val[0], val[1])
val_int = integrar(func, val_escal[1], val_escal[0])
```
