import matplotlib.pyplot as plt
import numpy as np

def gaussxw(N):
    """
    Calcula los puntos y pesos de Gauss-Legendre para la integración numérica.

    Parameters:
    N (int): Número de puntos de muestreo.

    Returns:
    tuple: Arrays de puntos de muestreo (x) y pesos (w).
    """
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
    """
    Escala los puntos y pesos de Gauss-Legendre a un intervalo arbitrario [a, b].

    Parameters:
    a (float): Límite inferior del intervalo.
    b (float): Límite superior del intervalo.
    x (array): Puntos de muestreo de Gauss-Legendre.
    w (array): Pesos de Gauss-Legendre.

    Returns:
    tuple: Arrays de puntos escalados y pesos escalados.
    """
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def func(x):
    """
    Define la función a integrar.

    Parameters:
    x (float or array): Variable independiente.

    Returns:
    float or array: Valor de la función evaluada en x.
    """
    return x**6 - x**2 * np.sin(2*x)

def integrar(function, wk, xk):
    """
    Calcula la integral numérica usando la cuadratura de Gauss-Legendre.

    Parameters:
    function (callable): Función a integrar.
    wk (array): Pesos de Gauss-Legendre.
    xk (array): Puntos de muestreo de Gauss-Legendre.

    Returns:
    float: Valor aproximado de la integral.
    """
    return np.sum(wk * function(xk))

val_real = 317.34424667382635
i = 1
error = 1

while error > 10e-12:
    val = gaussxw(i)
    val_escal = gaussxwab(1, 3, val[0], val[1])
    val_int = integrar(func, val_escal[1], val_escal[0])
    error = np.abs(val_int - val_real) / val_real
    print(f'El valor de la integral para N = {i} es: {val_int}')
    i += 1
