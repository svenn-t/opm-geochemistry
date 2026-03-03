# Equilibrium with exchange sites - `IEXCHANGE` keyword
<!-- Table of contents: Run pandoc with --toc option -->

To model ion exchange one need to use the `IEXCHANGE` keyword. For each ion exchange basis specie in the database, `GEOCHEMX` reads two items which must be on the same line: name of the exchanger followed by the exchange capacity in units of eq/L. A default exchanger $X=X^{-}$ is hard-coded into the database, and may be entered either with or without the charge included. For custom exchangers (see examples further below), the charge should be included, e.g., `Z+ 0.01`, to model positive exchange sites.

### Input file

```
GEOCHEM

SOLUTION 0
pH 7 charge
K                0.2e-3
NO3             1.2e-3
Na               1.0e-3
/end

iexchange
X 1.1e-3
/ end

/end

```

*Explanation of input file:*
Note that the solution keyword is described in the [solution example](../../solution/markdown/main_solution.md). To model exchange capacity of e.g. a clay mineral or organic matter, we need to invoke the `IEXCHANGE` keyword.
`IEXCHANGE .. /end`:
  :    
  Defines an ion exchange calculation block.
`X`:
  :    
  is the name of the ion exchanger
`1.1e-3`:
  :    
  is the total exchange capacity (in mol/kg water equivalent).

When run, `GeoChemX` will:

1. Create a solution at pH 7.
2. Add the dissolved ions $\mathrm{K^{+}}$, $\mathrm{NO_3^{-}}$, and $\mathrm{Na^{+}}$.
3. Perform a charge-balance adjustment using $\mathrm{H^{+}}$.
4. Add an ion-exchange phase with $1.1\cdot 10^ {-3}$ mol of exchange sites.
5. Compute how positive ions ($\mathrm{K^{+}}$, $\mathrm{Na^{+}}$, $\mathrm{Ca^{2+}}$, $\mathrm{H^{+}}$, etc.) partition between:
  * aqueous phase
  * exchange sites


### Description of output
The $pH$ of the solution is close to 3, this is because positive ions are pulled from the solution onto the ion exchanger.
In the file `<root name>_solution.out`, the following table is given (for description of the other output see [solution example](../../solution/html/solution-readable.md))

<table class="dotable" border="1">
<thead>
<tr><th align="center">﻿Species</th> <th align="center">Concentration</th> <th align="center">Equivalents</th> <th align="center">Equivalent fraction</th> <th align="center">Log gamma </th> <th align="center">Log K</th> </tr>
</thead>
<tbody>
<tr><td align="center">   NaX         </td> <td align="center">   9.07E-04         </td> <td align="center">   9.07E-04       </td> <td align="center">   8.24E-01               </td> <td align="center">   -0.0165962    </td> <td align="center">   0        </td> </tr>
<tr><td align="center">   KX          </td> <td align="center">   1.96E-04         </td> <td align="center">   1.96E-04       </td> <td align="center">   1.79E-01               </td> <td align="center">   -0.0167785    </td> <td align="center">   -0.7     </td> </tr>
</tbody>
</table>
This table describes the composition of exchange sites (`X`) in an ion-exchange model.
Each row represents an exchange complex, such as `NaX` or `KX`, where a cation ($\mathrm{Na^{+}}$, $\mathrm{K^{+}}$) is bound to an exchange site (`X⁻`).

`Species`:
  :    
  name of the exchange complex
`Concentration`:
  :    
  Total amount (mol/kg water) of this surface species
`Equivalents`:
  :    
  Because $\mathrm{Na^{+}}$ and $\mathrm{K^{+}}$ both have charge +1, equivalents = concentration. If you had CaX$_2$ , equivalents would be 2 $\times$ concentration.
`Equivalent fraction`:
  :    
  Fraction of total exchange capacity occupied by that cation. They should sum to one, decreasing `chemtol` would make the sum closer to one.  
  * NaX = 0.824 $\to$ 82.4 $\%$ of exchange sites are Na-bound
  * KX = 0.179 $\to$  17.9 $\%$ are K-bound

`Log gamma`:
  :    
  Logarithm of the activity coefficient of the exchange species. Slightly negative values indicate small non-ideality on the exchange surface. Very small magnitude, behavior is close to ideal.
`Log K`:
  :    
  The log-equilibrium constant for the exchange reaction. For KX, $\log K=−0.7$. This is often defined relative to NaX, which typically has $\log K = 0$ (reference state). Interpretation: $\mathrm{K^{+}}$ binds slightly less strongly to exchange sites than $\mathrm{Na^{+}}$ under the model assumptions.

### Equilibrium without changing solution composition
In many cases one wants to calculate the ion exchange composition *without* changing the solution cations. By using the `equilibrate` keyword `GeoChemX` will adjust the exchange composition, and at the same time keep bulk concentration constant.    


```
GEOCHEM

equilibrate 1

SOLUTION 0
pH 7 charge
K            0.2e-3
NO3          1.2e-3
Na           1.0e-3
/end

iexchange
X 1.1e-3
/ end

chemtol
1e-8  1e-5
/end

/end
```

Note that we have increases the accuracy on the mass balance, by the `chemtol` keyword. Running `GeoChemX` gives a final pH close to 7 and the following exchange composition

<table class="dotable" border="1">
<thead>
<tr><th align="center">﻿Species</th> <th align="center">Concentration</th> <th align="center">Equivalents</th> <th align="center">Equivalent fraction</th> <th align="center">Log gamma </th> <th align="center">Log K</th> </tr>
</thead>
<tbody>
<tr><td align="center">   NaX         </td> <td align="center">   5.49E-04         </td> <td align="center">   5.49E-04       </td> <td align="center">   4.99E-01               </td> <td align="center">   -0.0165936    </td> <td align="center">   0        </td> </tr>
<tr><td align="center">   KX          </td> <td align="center">   5.51E-04         </td> <td align="center">   5.51E-04       </td> <td align="center">   5.01E-01               </td> <td align="center">   -0.0167758    </td> <td align="center">   -0.7     </td> </tr>
</tbody>
</table>
Note that the sum of equivalent fractions are now 1.

## Ion exchange with positive exchange sites

It is possible to add new species to model e.g. positive exchange sites, which might be important for chalk [@nodland2024new]. To illustrate we introduce a positive exchange site `Z+`


```
GEOCHEM

equilibrate 1
include "ions.txt"

SOLUTION 0
pH 7 charge
K               0.2e-3
NO3             1.2e-3
SO4		0.4e-3
Na              1.8e-3
/end

iexchange
Z+ 1.1e-3
/ end

chemtol
1e-8  1e-5
/end

/end

```

New specie definition has to be included in a separate file `ions.txt`, using the `include` keyword.

*Explanation of input file `ions_txt`:*


```
EXCHANGE_SPECIES
#name 
Z+  
/end

SECONDARY_SPECIES
SO4Z2 = 2Z+ + SO4-2 / ANA  -0.8 /
ZNO3   = Z+ + NO3- /ANA -0.7/
/end
```

`EXCHANGE_SPECIES`:
  :    
  Defines the exchange master species name(s) used by the exchanger, note the charge needs to be included.
`SECONDARY_SPECIES`:
  :    
  Secondary species that involve the exchange site `Z+` and aqueous anions ($\mathrm{SO_{4}^{2-}}$, $\mathrm{NO_3^{-}}$). Each line is a mass-action reaction with an associated $\log K$ (formation constant).
`ANA`:
  :    
  Analytical formula for $\log K$, defined in the equation below , where $T$ is in Kelvin. 

$$
\begin{equation}
\log_{10} K = a_1 + a_2 T + \frac{a_3}{T} + a_4\log_{10}{T} + \frac{a_5}{T^2} + a_6 T^2\,
\end{equation}
$$

Running `GeoChemX` gives a final pH close to 7 and the following exchange composition

<table class="dotable" border="1">
<thead>
<tr><th align="center">﻿Species</th> <th align="center">Concentration</th> <th align="center">Equivalents</th> <th align="center">Equivalent fraction</th> <th align="center">Log gamma </th> <th align="center">Log K</th> </tr>
</thead>
<tbody>
<tr><td align="center">   SO4Z2       </td> <td align="center">   4.88E-04         </td> <td align="center">   9.75E-04       </td> <td align="center">   8.87E-01               </td> <td align="center">   -0.0920841    </td> <td align="center">   -0.8     </td> </tr>
<tr><td align="center">   ZNO3        </td> <td align="center">   1.25E-04         </td> <td align="center">   1.25E-04       </td> <td align="center">   1.13E-01               </td> <td align="center">   -0.0233732    </td> <td align="center">   -0.7     </td> </tr>
</tbody>
</table>
<br />
## Bibliography

 1. <div id="nodland2024new"></div> **O. Nodland and A. Hiorth**.  A New Formulation of the Surface Charge/surface Potential Relationship in Electrolytes With Valence Less Than Three, *Computational Geosciences*, 28(2), pp. 289-304, 2024.


