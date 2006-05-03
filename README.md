# Computational Physics
## Numerical recipes applied to thermodynamics

This is a bunch of old C programs based on the good ol' **Numerical Recipes** [http://www.nr.com/](http://www.nr.com/). Vectors' components start in **1** and not in **0** as we are used to.

**Example:**
```c
	for (i=1; i<=N; i++) {
		printf("Enter initial value of parameter %d: ", i);
		scanf("%lf", &a[i]);
	}
```

## Folder structure

There are three computational physics problems proposed:

* **Folder `breitwigner`:** Solving the Breit-Wigner formula, actually a Nuclear Physics problem.
* **Folder `pvdiagrams`:** Pressure-Volume (PV) diagrams and the *van der Waals* equation.
* **Folder `gibbs`:** Vapor-liquid equilibrium and calculation of the excess Gibbs potential.

The **Numerical Recipes** libraries used for this calculations are stored in the folder `nr`, but you can also install them in your system and load them normally, if you prefer. Check the NR documentation to learn more about it.

The libraries in the `nr` folder used for these computational physics problems are NOT included due to [license restrictions](http://www.nr.com/licenses/redistribute.html). The snippets might not work without them.

## Compiling and running

The `-lm` is needed to link the math libraries.

### Compile:

```bash
$ gcc filename.c -o filename -lm
```

### Run:

```bash
$ ./filename
```

## Dedicated site

Learn more at http://nrthermo.tk

## Licenses

[![License](https://img.shields.io/badge/gnu-license-green.svg?style=flat)](https://opensource.org/licenses/GPL-2.0)
All code is [GNU licensed](https://opensource.org/licenses/GPL-2.0).

[![License](https://img.shields.io/badge/NR3-license-green.svg?style=flat)](http://www.nr.com/licenses/redistribute.html)
Numerical Recipes libraries are licensed according to [these terms](http://www.nr.com/licenses/redistribute.html).

