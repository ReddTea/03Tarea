# Tarea nro. 3
## FI3104B - Metodos Numericos Para la Ciencia y la Ingenieria
#### Prof. Valentino González

## P1

El oscilador de van der Pool fue propuesto para describir la dinámica de algunos
circuitos eléctricos. La ecuación es la siguiente:

<img src='eqs/van_der_pool_1.png' alt='van de Pool eqn.' height='40'>

> Latex:
    `\dfrac{d^2x}{dt^2} = - k x - \mu (x^2 - a^2) \dfrac{dx}{dt}`

donde k es la constante elástica y &mu; es un coeficiente de roce. Si |x| > a el
roce amortigua el movimiento, pero si |x| < a el roce inyecta energía. Se puede
hacer un cambio de variable para convertir la ecuación a:

<img src='eqs/van_der_pool_2.png' alt='van de Pool eqn. transformada' height='40'>

> Latex:
    `\dfrac{d^2y}{ds^2} = - y - \mu^* (y^2 - 1) \dfrac{dy}{ds}`

con lo cual ahora la ecuación sólo depende de un parámetro, &mu;<sup>\*</sup>.
Indique cuál es el cambio de variable realizado.

Integre la ecuación de movimiento usando el método de Runge-Kutta de orden 3
visto en clase. Se pide que Ud. implemente su propia versión del algoritmo,
describa la discretización usada y el paso de tiempo. Use
&mu;<sup>\*</sup>=1.RRR, donde RRR son los 3 últimos dígitos de su RUT antes
del guión.

> Nota: No olvide incluir su Nombre y RUT en el informe.

Integre la solución hasta T=20&pi; (aproximadamente 10 períodos) para las
siguientes condiciones iniciales:

<img src='eqs/iniciales.png' alt='condiciones iniciales' height='80'>

>Latex:
    `1) \dfrac{dy}{ds} = 0; y = 0.1\\
     2) \dfrac{dy}{ds} = 0; y = 4.0`

Grafique y(s) y la trayectoria en el espacio (y, dy/ds).

## P2

El sistema de Lorenz es un set de ecuaciones diferenciales ordinarias conocido
por tener algunas soluciones caóticas, la más famosa el llamado atractor de
Lorenz. El sistema de ecuaciones es el siguiente:

<img src='eqs/lorenz.png' alt='Lorenz system' height='120'>

> Latex:
    ```
    \begin{flalign*}
    \dfrac{dx}{ds} &= \sigma (y - x)\\
    \dfrac{dy}{ds} &= x (\rho - z) - y\\
    \dfrac{dz}{ds} &= xy - \beta z
    \end{flalign*}
    ```

La solución más famosa se obtiene con los parámetros &sigma;=10, &beta;=8/3 y
&rho;=28. Utilice esos parámetros, elija un set de condiciones iniciales
(x<sub>0</sub>, y<sub>0</sub>, z<sub>0</sub>) e integre la ecuación por un
tiempo que estime conveniente. Esta vez se pide que utilice un algoritmo RK4
pero no necesita implementarlo, puede usar los algoritmos disponibles en
`scipy.integrate` o cualquier otro que encuentre y que sea de uso libre.

Plotee en 3D la solución (x(t), y(t), z(t)).

> Nota: Este repositorio incluye un archivo 3D.py que demuestra como hacer un
  plot en 3D usando `matplotlib`. Puede usar este archivo como guía de ayuda.

__Otras Notas.__

- Utilice `git` durante el desarrollo de la tarea para mantener un historial de
  los cambios realizados. La siguiente [*cheat
  sheet*](https://education.github.com/git-cheat-sheet-education.pdf) le puede
  ser útil. Esta vez empezaremos a evaluar el uso efectivo de `git`. Recuerde
  hacer cambios significativos pero relativamente pequeños y guardar seguido.
  Evite hacer `commits` de código que no compila y deje mensajes que permitan
  entender los cambios realizados.

- Tambien comenzaremos a revisar su uso correcto de python. Si define una
  función relativamente larga o con muchos parámetros, recuerde escribir el
  *doctsring* que describa los parametros y que es lo que hace la funcion.
  También recuerde usar nombres explicativos para las variables y las funciones.
  El mejor nombre es aquel que permite entender que hace la funcion sin tener
  que leer su implementacion.

- La tarea se entrega como un *pull request* en github. El *pull request* debe
  incluir todos los códigos usados además de su informe.

- El informe debe ser entregado en formato *pdf*, este debe ser claro sin
  información ni de más ni de menos.
