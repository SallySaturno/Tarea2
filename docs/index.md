# Introducción

Este script implementa la cuadratura de Gauss-Legendre para la integración numérica.

El objetivo es aproximar la integral de una función en un intervalo dado mediante
el método de cuadratura de Gauss. Se busca determinar la cantidad mínima de puntos
de muestreo necesarios para que el error relativo sea menor a 10⁻¹²(en mi caso).

Se realiza lo siguiente:

- Se calculan los puntos y pesos de Gauss-Legendre.
- Se escalan los valores al intervalo de integración deseado.
- Se evalúa la función en los puntos escalados y se calcula la integral aproximada.
- Se repite el proceso aumentando el número de puntos hasta alcanzar el error deseado.

La integral a calcular utilizada para este código es:

$$
I = \int_{1}^{3} \left( x^6 - x^2 \sin(2x) \right) \,dx
$$

El resultado es una aproximación precisa del valor de la integral con un error controlado.
