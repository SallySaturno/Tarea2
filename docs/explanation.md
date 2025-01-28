## Descripción del Método Numérico: Cuadratura Gaussiana

Uno de los métodos más poderosos para evaluar integrales de forma numérica es la **cuadratura Gaussiana**.

### Idea Principal

La cuadratura Gaussiana busca aproximar una integral definida de la siguiente manera:

\[
\int_a^b {\rm{d}}x f(x) \approx \sum_{k=1}^{N+1} w_k f(x_k).
\]

donde:

  - $ w_k $ son los "pesos".
  - $ x_k $ son los puntos de muestreo. Se usan $ N+1 $ puntos (es decir, $ N $ subregiones o subintervalos).

### Comparación con el Método de Newton-Cotes

En el método de Newton-Cotes:

  - Los puntos de muestreo son **equidistantes**.
  - Una ecuación de Newton-Cotes de orden $ N $ es *exacta* para un polinomio de grado $ N $.
  - Un polinomio de orden $ N $ aproxima mejor una función bien comportada que un polinomio de orden $ N-1 $.

Por el contrario, en la **cuadratura Gaussiana**:

  - Los puntos de muestreo **no son equidistantes**, lo que introduce más grados de libertad para la misma discretización en $ N $ subregiones.
  - Es exacta para un polinomio de orden $ 2N - 1 $, lo que significa que proporciona la misma precisión que un polinomio de orden $ 2N - 1 $ en Newton-Cotes.

### Elección de Pesos y Puntos de Muestreo

Existe una **regla universal** para elegir los valores de $ x_k $ y $ w_k $:

  - $ x_k $ corresponden a las $ N $ raíces (ceros) de los polinomios de Legendre $ P_N(x) $ de orden $ N $.
  - Los pesos $ w_k $ se calculan como:

\[
w_k = \left[\frac{2}{1-x^2}\left(\frac{dP_N}{dx}\right)^{-2}\right]_{x={x_k}},
\]

donde $ x_k $ son los ceros de $ P_N(x) $.

### Pros y Contras de la Cuadratura Gaussiana

#### Pros:
  - Aunque la ecuación para evaluar los errores es compleja, la aproximación mejora con un error que decrece por un factor \( \text{const} / N^2 \) cuando se incrementa el número de subregiones en uno.
  - Ejemplo: Pasar de \( N=10 \) a \( N=11 \) mejora la estimación por un factor de \( \approx 100 \), lo que indica una convergencia rápida con pocos puntos de muestreo.

#### Contras:
  - Funciona bien solo si la función a integrar es relativamente bien comportada. Si no lo es, se requieren más puntos de muestreo cerca de las regiones problemáticas.
  - Es complicado evaluar el error de manera precisa si se necesita con exactitud.

Este método es una herramienta fundamental en la integración numérica debido a su precisión y eficiencia, aunque su aplicación requiere un análisis cuidadoso de la función a integrar.
