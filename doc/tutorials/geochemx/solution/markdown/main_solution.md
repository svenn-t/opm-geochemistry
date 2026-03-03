# Seawater speciation and the `SOLUTION` keyword
<!-- Table of contents: Run pandoc with --toc option -->

## Seawater speciation
This tutorial showcase the `SOLUTION` keyword, and we use a typical seawater composition to demonstrate how the calculations are done. The distribution of the major species in seawater is shown in the tables below, see [@garrels1965solutions] (p. 105) at ambient conditions.

<table class="dotable" border="1">
<thead>
<tr><th align="center">Metal (Me)</th> <th align="center">Molality</th> <th align="center">Free ion (%)</th> <th align="center">Me-SO4 (%)</th> <th align="center">Me-HCO3 (%)</th> <th align="center">Me-CO3 (%)</th> </tr>
</thead>
<tbody>
<tr><td align="center">   Na            </td> <td align="center">   0.48        </td> <td align="center">   99              </td> <td align="center">   1             </td> <td align="center">   0.04           </td> <td align="center">   0.006         </td> </tr>
<tr><td align="center">   Mg            </td> <td align="center">   0.054       </td> <td align="center">   86              </td> <td align="center">   13            </td> <td align="center">   0.6            </td> <td align="center">   0.4           </td> </tr>
<tr><td align="center">   Ca            </td> <td align="center">   0.01        </td> <td align="center">   90              </td> <td align="center">   9             </td> <td align="center">   0.6            </td> <td align="center">   0.1           </td> </tr>
<tr><td align="center">   K             </td> <td align="center">   0.01        </td> <td align="center">   98              </td> <td align="center">   12            </td> <td align="center">   -              </td> <td align="center">   -             </td> </tr>
</tbody>
</table>
<table class="dotable" border="1">
<thead>
<tr><th align="center">Anion (An)</th> <th align="center">Molality</th> <th align="center">Free ion (%)</th> <th align="center">Ca-An (%)</th> <th align="center">Mg-An (%)</th> <th align="center">Na-An (%)</th> <th align="center">K-An  (%)</th> </tr>
</thead>
<tbody>
<tr><td align="center">   SO4           </td> <td align="center">   0.028       </td> <td align="center">   52              </td> <td align="center">   3            </td> <td align="center">   26           </td> <td align="center">   19           </td> <td align="center">   0.6          </td> </tr>
<tr><td align="center">   HCO3          </td> <td align="center">   0.0024      </td> <td align="center">   76              </td> <td align="center">   3            </td> <td align="center">   14           </td> <td align="center">   8            </td> <td align="center">   -            </td> </tr>
<tr><td align="center">   CO3           </td> <td align="center">   0.00027     </td> <td align="center">   9               </td> <td align="center">   7            </td> <td align="center">   72           </td> <td align="center">   12           </td> <td align="center">   -            </td> </tr>
<tr><td align="center">   Cl            </td> <td align="center">   0.56        </td> <td align="center">   100             </td> <td align="center">   -            </td> <td align="center">   -            </td> <td align="center">   -            </td> <td align="center">   -            </td> </tr>
</tbody>
</table>
The geochemical calculations are performed by  giving the *total concentration* of  *basis species*. The default basis species is a choice, and in `GeoChemX` we have chosen to use the basis species listed in  [Appendix: List of basis species](#sec:basis_species).

*changing basis species.* 
In `GeoChemX` it is possible to change the basis, see the example towards the end.



### Input file


```
GEOCHEM
Temp 25
Pres 1.01325e5
debug 0

SOLUTION 0
pH 7 charge
Na 0.48
K 0.01
Mg 0.054
Ca 0.01
SO4 0.028
HCO3 0.00267
Cl 0.56
/end

chemtol
1e-8  1e-5
/end
/end

```

*Explanation of input file:*

`GEOCHEM .. /end`:
  :    
  Defines a GEOCHEM calculation block.
`Temp`:
  :    
  temperature of solution, units: Celsius, default 25 $^\circ$C
`Pres`:
  :    
  pressure of solution, units: Pascal, default 10$^5$ Pa
`debug`:
  :    
  controls level of debugging
  * 0 No additional debug information
  * 1 Prints debug data to screen and some extra files

`SOLUTION <name>`:
  :    
  Defines the aqueous solution to be modeled. `<name>` must be a unique identifier (in this case 0):
  * Specie names are listed (without charge), units: mol/kg water (Molality).
  * The solution block is terminated by `\end`.

`pH 7 charge`:
  :    
  the activity of H$^+$ is determined by charge balance and 7 tells the solver to use $10^{-7}$ as initial guess. Alternatively one can enter the total concentration of $\mathrm{H^{+}}$. However as the total concentration of *free hydrogen ions* can be negative in the calculations, the solver adds 1 to the hydrogen concentration. To be consistent one should then write `H 1.0000001`.
`chemtol`:
  :    
  Sets numerical convergence criteria, default is `1e-5 1e-5`.
  * Solver uses one loop for mass balance, and a second outer loop for pH or charge balance. The first number controls the convergence of mass balance and the second, charge balance.
  * `chemtol` block should be put at the end of the file. 


The input file described above saved in `seawater.dat`, and the code can be run

```
Terminal> GeoChemX EQSOLVER seawater.dat
```

`GeoChemX` produce four files with data, described below and one file named `converged.txt` containg 1 if converged and zero otherwise. 

### Description of output
If debug flag is zero there are four output files. The naming convention is root name and then extension `_solution.out`, `_aq.out`, `_buffer.out`, and `_species.out`. Each of the files are described below. Very briefly, the `_solution.out` file contains the most comprehensive description of the calculation. The other files contains subset of information, but in a tab separated format that can be easily read in pandas.

### `seawater_solution.out` Brief Explanation

*1. Header: Surface Charge and Electrostatic Info.*


```
    pH   Sch[C/m^2]   Sch_DL[C/m^2]   psi[mV]   F[C/mol]   S[m^2/L]   cbal[eq/L]
    6.48826	 0	 0	 0	 96484.6	 1	3.6198588848623814e-09
```

`pH`:
  :    
  6.48826
`Sch, Sch_DL`:
  :    
  Surface charge densities (C/m$^2$), here zero (no
    surface).
`psi`:
  :    
  Surface potential (mV).
`F`:
  :    
  Faraday constant.
`S`:
  :    
  Specific surface area (m$^2$/L).
`cbal`:
  :    
  Charge balance error (eq/L).

*2. Physical Properties.*

```
    Temp   Pres[bar]   Io_   permittivity   density_water[kg/m^3]
    25     1.01325     0.649459   78.2439   997.048
```

`Temp`:
  :    
  25 C
`Pres`:
  :    
  1.01325 bar
`Io_`:
  :    
  Ionic strength
`Permittivity`:
  :    
  Dielectric constant
`Water density`:
  :    
  density of water (not brine)

*3. Double-Layer / Ion-Group Parameters.*

```
    n(1) n(-1) n(2) n(-2) g(1) g(-1) g(2) g(-2) frac_DL
    0    0     0    0     0    0     0    0     1
```

* Defines charge groups for electrostatic models.
* Here all are zero, no surface-complexation model.
* *frac_DL = 1*,  everything is in the bulk solution (volume of diffusive layer is zero)

* 4. Species List and Type Flags.*

```
    Ca+2  Cl-  H+  H2O  HCO3-  K+  Mg+2  Na+  SO4-2  CaOH+ ...
    type 0    0    0    0    0     0    0     0     0       ...
```

* Lists all *aqueous species* considered, basis and secondary species.
* `type = 0` for all, i.e.  all are dissolved aqueous species.

*5. Thermodynamic Data.*

`log_a`:
  :    
  Logarithm of species activity.
`log_m`:
  :    
  Logarithm of molality (mol/kg).
`log_g`:
  :    
  Logarithm of the activity coefficient $\gamma$
`logK`:
  :    
  Equilibrium constant for the formation reaction of each species.
`charge / scharge`:
  :    
  Electrical charge; surface charge contribution (0 here).

*6. Species Concentrations.*

Example:

```
    Species   ch   c_tot   c_tot_aq   c_tot_SC   c_tot_IO   c_dl
    Ca+2       2   1.0e-02   1.0e-02    0            0           0
    Cl-       -1   5.6e-01   5.6e-01    0            0           0
    H+         1   1.00000010e+00	 9.99327309e-01   0   0   0
    H2O        0   0                 0.999941391  0  0  0
```

Field meanings:

`c_tot`:
  :    
  total input concentration
`c_tot_aq`:
  :    
  dissolved (aqueous) concentration
`c_tot_SC / IO / dl`:
  :    
  surface, inner/outer layer, diffuse layer (here all zero)

*Notice.* 

GEOCHEM offsets total hydrogen by +1 mol/kg,  hence `1.0000001`
    for H+ total



*7. Mineral Saturation Indices (SI).*

Example:

```
    Name      SI        logK        reaction
    CALCITE   0.248389  1.84866     1Ca+2 -1H+ 1HCO3-
    DOLOMITE  2.4845    2.5135      1Ca+2 -2H+ 2HCO3- 1Mg+2
    HALITE   -2.61068   1.5855      1Cl- 1Na+
    ANHYDRITE -0.9547  -4.30642	   +1Ca+2+1SO4-2
```

* `SI` Saturation Index
  * $SI > 0$,  supersaturated
  * $SI = 0$, equilibrium
  * $SI < 0$, undersaturated

* `logK` Equilibrium constant for mineral.
* reaction: Dissolution/precipitation reaction for the mineral, note that the

*Notice.* 
All reactions are written in terms of basis species, the `logK` values are only valid for this specific basis combination.



### `seawater_aq.out` Brief Explanation
The file contains a tab separated table that lists aqueous species and ion pairs along with their:

`log_m`:
  :    
  log molality (mol/kg water)
`log_a`:
  :    
  log activity
`log_g`:
  :    
  log activity coefficient ($\gamma$)
`logK`:
  :    
  equilibrium constant of formation (where applicable)
`mol_volume`:
  :    
  partial molal volume
`reaction`:
  :    
  formation reaction from basis species

The species are listed in descendent order of concentration.

*$\log K$ values.* 
The $\log K$ values are *only valid for the specific reaction*, e.g. if one would like the $\log K$ of the reaction $\mathrm{MgCO_{3}^0} {\rightleftharpoons} \mathrm{Mg^{2+}} + \mathrm{CO_2(aq)}$, one need subtract $\log K$ for the reaction $\mathrm{MgCO_{3}^0}$ and $\mathrm{CO_2(aq)}$ in the datafile.


### `seawater_buffer.out` Brief Explanation

The file contains a tab separated table that lists minerals and gasses along with their:

`SI`:
  :    
  saturation index, for gases this would be equal to the partial pressure
`logK`:
  :    
  equilibrium constant of formation (where applicable)
`reaction`:
  :    
  formation reaction from basis species

The minerals are listed in descendent order of saturation, i.e. the most supersaturated mineral at the top.

*$\log K$ values.* 
The $\log K$ values are *only valid for the specific reaction*, e.g. if one would like the $\log K$ of the reaction $\mathrm{MgCO_{3}^0} {\rightleftharpoons} \mathrm{Mg^{2+}} + \mathrm{CO_3^{2-}}$ (magnesite), one need subtract $\log K$ for the reaction $\mathrm{MgCO_{3}^0}$ in *this file* and $\mathrm{CO_3^{2-}}$ in the `_aq.out` datafile.


### `seawater_species.out` Brief Explanation

This file lists the total and partitioned concentrations of dissolved *basis* species in a geochemical solution. 

The file contains a tab separated table that lists basis species along with their:

`ch`:
  :    
  electrical charge 
`c_tot`:
  :    
  total concentration at the start of simulation, units: mol/kg water
`c_tot_aq`:
  :    
  total aqueous concentration after the simulation, units: mol/kg water
`c_tot_SC`:
  :    
  total surface complex concentration after the simulation, units: mol/kg water
`c_tot_IO`:
  :    
  total ion exchange concentration after the simulation, units: mol/kg
`c_dl`:
  :    
  total concentration in the diffusive layer after the simulation, units: mol/kg 

## Changing basis species `as` keyword

The choice of basis species is flexible and can be modified by the user. This flexibility can be particularly useful when modeling brines or other systems in which the redox potential differs significantly from standard conditions. Under sufficiently reducing conditions, sulfur is preferentially present in reduced forms such as $\text{HS}^⁻$ and $\text{H}_2\text{S}$, whereas under oxidizing conditions sulfate, $\mathrm{SO_{4}^{2-}}$, is the dominant sulfur species. Using the same example as above we can now do the calculation using $\text{HS}^⁻$

```
SOLUTION 0
pH 7 charge
Na 0.48
K 0.01
Mg 0.054
Ca 0.01
SO4 0.028 as HS
HCO3 0.00267
Cl 0.56
/end
```

This will a transform all reactions to use $\text{HS}^⁻$ as basis specie and update $\log K$ correspondingly. Note that in this case some minerals will not form, such as anhydrite. This is because now anhydrite cannot be written in terms of the basis species used.

## Specifying `pe`
pe is a dimensionless parameter that represents the electron activity (or redox potential) of a solution. It is analogous to pH, but instead of describing proton activity, it describes the tendency of a solution to accept or donate electrons.

In geochemistry:
$$
\begin{equation}
\text{pe} = -\log a_e,
\end{equation}
$$
where $a_e$ is the electron activity in moles per liter. By changing the pe one can control sulfate reduction reactions
$$
\begin{equation}
\mathrm{SO_{4}^{2-}} + 8e^- + 10\mathrm{H^{+}} \text{HS}^⁻{\rightleftharpoons} \text{HS}^- + 4 \mathrm{H_2O}.
\end{equation}
$$
A low pe will move the reaction to the right. In `GeoChemX` we can specify the pe as shown below, mathematically the electron is treated as a basis specie and the by defining the pe we just keep the activity of the electron constant as $10^-\text{pe}$.

```
SOLUTION 0
pH 7 charge
pe 2
Na 0.48
K 0.01
Mg 0.054
Ca 0.01
SO4 0.028 as HS
HCO3 0.00267
Cl 0.5
/end
```

## Bibliography

 1. <div id="garrels1965solutions"></div> **R. Garrels and C. Christ**.  *Solutions, Minerals, and Equilibria*, *Harper international student reprint*, Harper  Row, 1965, <https://books.google.no/books?id=OT5RAAAAMAAJ>.

## Appendix: List of basis species
<div id="sec:basis_species"></div>

<!-- The table below shows all the possible aqueous species that can appear as -->
<!-- basis species. However, as discussed in the main text, minerals are sometimes -->
<!-- also swapped into the basis. -->

<table class="dotable" border="1">
<thead>
<tr><th align="center">  Name </th> <th align="center">Nick Name</th> <th align="center">Mol Weight</th> <th align="center"> Name </th> <th align="center">Nick Name</th> <th align="center">Mol Weight</th> </tr>
</thead>
<tbody>
<tr><td align="center">   Ag+        </td> <td align="center">   Ag           </td> <td align="center">   107.868       </td> <td align="center">   MoO4-2    </td> <td align="center">   MoO4         </td> <td align="center">   159.9376      </td> </tr>
<tr><td align="center">   Al+3       </td> <td align="center">   Al           </td> <td align="center">   26.98154      </td> <td align="center">   NO3-      </td> <td align="center">   NO3          </td> <td align="center">   62.0049       </td> </tr>
<tr><td align="center">   Am+3       </td> <td align="center">   Am           </td> <td align="center">   243           </td> <td align="center">   Na+       </td> <td align="center">   Na           </td> <td align="center">   22.98977      </td> </tr>
<tr><td align="center">   Ar         </td> <td align="center">   Ar(aq)       </td> <td align="center">   39.948        </td> <td align="center">   NbO3-     </td> <td align="center">   NbO3         </td> <td align="center">   140.9046      </td> </tr>
<tr><td align="center">   Au+        </td> <td align="center">   Au           </td> <td align="center">   196.9665      </td> <td align="center">   Nd+3      </td> <td align="center">   Nd           </td> <td align="center">   144.24        </td> </tr>
<tr><td align="center">   BO2-       </td> <td align="center">   BO2          </td> <td align="center">   42.8088       </td> <td align="center">   Ne        </td> <td align="center">   Ne(AQ)       </td> <td align="center">   20.179        </td> </tr>
<tr><td align="center">   Ba+2       </td> <td align="center">   Ba           </td> <td align="center">   137.33        </td> <td align="center">   Ni+2      </td> <td align="center">   Ni           </td> <td align="center">   58.7          </td> </tr>
<tr><td align="center">   Be+2       </td> <td align="center">   Be           </td> <td align="center">   9.01218       </td> <td align="center">   Pb+2      </td> <td align="center">   Pb           </td> <td align="center">   207.2         </td> </tr>
<tr><td align="center">   Bi+3       </td> <td align="center">   Bi           </td> <td align="center">   208.9804      </td> <td align="center">   Pd+2      </td> <td align="center">   Pd           </td> <td align="center">   106.4         </td> </tr>
<tr><td align="center">   Br-        </td> <td align="center">   Br           </td> <td align="center">   79.904        </td> <td align="center">   Pm+3      </td> <td align="center">   Pm           </td> <td align="center">   145           </td> </tr>
<tr><td align="center">   Ca+2       </td> <td align="center">   Ca           </td> <td align="center">   40.08         </td> <td align="center">   Pr+3      </td> <td align="center">   Pr           </td> <td align="center">   140.90765     </td> </tr>
<tr><td align="center">   Cd+2       </td> <td align="center">   Cd           </td> <td align="center">   112.41        </td> <td align="center">   Pt+2      </td> <td align="center">   Pt           </td> <td align="center">   195.078       </td> </tr>
<tr><td align="center">   Ce+3       </td> <td align="center">   Ce           </td> <td align="center">   140.116       </td> <td align="center">   Ra+2      </td> <td align="center">   Ra           </td> <td align="center">   222           </td> </tr>
<tr><td align="center">   Cl-        </td> <td align="center">   Cl           </td> <td align="center">   35.453        </td> <td align="center">   Rb+       </td> <td align="center">   Rb           </td> <td align="center">   85.4678       </td> </tr>
<tr><td align="center">   Co+3       </td> <td align="center">   Co           </td> <td align="center">   58.9332       </td> <td align="center">   Rh+3      </td> <td align="center">   Rh           </td> <td align="center">   102.9055      </td> </tr>
<tr><td align="center">   Cr+3       </td> <td align="center">   Cr           </td> <td align="center">   51.9961       </td> <td align="center">   Rn        </td> <td align="center">   Rn(AQ)       </td> <td align="center">   222           </td> </tr>
<tr><td align="center">   Cs+        </td> <td align="center">   Cs           </td> <td align="center">   132.9054      </td> <td align="center">   RuO4-2    </td> <td align="center">   RuO4         </td> <td align="center">   165.1676      </td> </tr>
<tr><td align="center">   Cu+        </td> <td align="center">   Cu           </td> <td align="center">   63.546        </td> <td align="center">   SO4-2     </td> <td align="center">   SO4          </td> <td align="center">   96.0626       </td> </tr>
<tr><td align="center">   Dy+3       </td> <td align="center">   Dy           </td> <td align="center">   162.5         </td> <td align="center">   Sc+3      </td> <td align="center">   Sc           </td> <td align="center">   44.9559       </td> </tr>
<tr><td align="center">   Er+3       </td> <td align="center">   Er           </td> <td align="center">   167.259       </td> <td align="center">   SeO4-2    </td> <td align="center">   SeO4         </td> <td align="center">   142.9576      </td> </tr>
<tr><td align="center">   Eu+3       </td> <td align="center">   Eu           </td> <td align="center">   151.964       </td> <td align="center">   SiO2      </td> <td align="center">   Si           </td> <td align="center">   56.171        </td> </tr>
<tr><td align="center">   F-         </td> <td align="center">   F            </td> <td align="center">   18.998403     </td> <td align="center">   Sm+3      </td> <td align="center">   Sm           </td> <td align="center">   150.36        </td> </tr>
<tr><td align="center">   Fe+3       </td> <td align="center">   Fe           </td> <td align="center">   55.847        </td> <td align="center">   Sn+2      </td> <td align="center">   Sn           </td> <td align="center">   118.69        </td> </tr>
<tr><td align="center">   Fr+        </td> <td align="center">   Fr           </td> <td align="center">   223           </td> <td align="center">   Sr+2      </td> <td align="center">   Sr           </td> <td align="center">   87.62         </td> </tr>
<tr><td align="center">   Ga+3       </td> <td align="center">   Ga           </td> <td align="center">   69.735        </td> <td align="center">   Tb+3      </td> <td align="center">   Tb           </td> <td align="center">   158.92534     </td> </tr>
<tr><td align="center">   GdO+       </td> <td align="center">   GdO          </td> <td align="center">   173.2494      </td> <td align="center">   TcO4-     </td> <td align="center">   TcO4         </td> <td align="center">   162.9038      </td> </tr>
<tr><td align="center">   H+         </td> <td align="center">   H            </td> <td align="center">   1.0079        </td> <td align="center">   Tl+3      </td> <td align="center">   Tl           </td> <td align="center">   204.37        </td> </tr>
<tr><td align="center">   H2O        </td> <td align="center">   H2O          </td> <td align="center">   18.0152       </td> <td align="center">   Tm+3      </td> <td align="center">   Tm           </td> <td align="center">   168.93421     </td> </tr>
<tr><td align="center">   HAsO4-2    </td> <td align="center">   HAsO4        </td> <td align="center">   139.92714     </td> <td align="center">   U+3       </td> <td align="center">   U            </td> <td align="center">   238.029       </td> </tr>
<tr><td align="center">   HCO3-      </td> <td align="center">   HCO3         </td> <td align="center">   61.0171       </td> <td align="center">   VO4-3     </td> <td align="center">   VO4          </td> <td align="center">   114.9391      </td> </tr>
<tr><td align="center">   HPO4-2     </td> <td align="center">   HPO4         </td> <td align="center">   95.979301     </td> <td align="center">   WO4-2     </td> <td align="center">   WO4          </td> <td align="center">   247.8476      </td> </tr>
<tr><td align="center">   He         </td> <td align="center">   He(AQ)       </td> <td align="center">   4.0026        </td> <td align="center">   Xe        </td> <td align="center">   Xe(aq)       </td> <td align="center">   131.3         </td> </tr>
<tr><td align="center">   Hf+4       </td> <td align="center">   Hf           </td> <td align="center">   178.49        </td> <td align="center">   Y+3       </td> <td align="center">   Y            </td> <td align="center">   88.9059       </td> </tr>
<tr><td align="center">   Hg+2       </td> <td align="center">   Hg           </td> <td align="center">   200.59        </td> <td align="center">   Yb+3      </td> <td align="center">   Yb           </td> <td align="center">   173.04        </td> </tr>
<tr><td align="center">   HoO+       </td> <td align="center">   HoO          </td> <td align="center">   180.92972     </td> <td align="center">   Zn+2      </td> <td align="center">   Zn           </td> <td align="center">   65.38         </td> </tr>
<tr><td align="center">   I-         </td> <td align="center">   I            </td> <td align="center">   126.9045      </td> <td align="center">   Zr+4      </td> <td align="center">   Zr           </td> <td align="center">   91.22         </td> </tr>
<tr><td align="center">   In+3       </td> <td align="center">   In           </td> <td align="center">   114.82        </td> <td align="center">   e-        </td> <td align="center">   e-           </td> <td align="center">   0             </td> </tr>
<tr><td align="center">   K+         </td> <td align="center">   K            </td> <td align="center">   39.0983       </td> <td align="center">   X-        </td> <td align="center">   X            </td> <td align="center">   0             </td> </tr>
<tr><td align="center">   Kr         </td> <td align="center">   Kr(AQ)       </td> <td align="center">   83.8          </td> <td align="center">   GCO3-     </td> <td align="center">   GCO3         </td> <td align="center">   0             </td> </tr>
<tr><td align="center">   La+3       </td> <td align="center">   La           </td> <td align="center">   138.9055      </td> <td align="center">   GCa+      </td> <td align="center">   GCa          </td> <td align="center">   0             </td> </tr>
<tr><td align="center">   Li+        </td> <td align="center">   Li           </td> <td align="center">   6.941         </td> <td align="center">   E         </td> <td align="center">   E            </td> <td align="center">   0             </td> </tr>
<tr><td align="center">   Lu+3       </td> <td align="center">   Lu           </td> <td align="center">   174.967       </td> <td align="center">   GSiOH     </td> <td align="center">   GSiOH        </td> <td align="center">   0             </td> </tr>
<tr><td align="center">   Mg+2       </td> <td align="center">   Mg           </td> <td align="center">   24.305        </td> <td align="center">   HSd-      </td> <td align="center">   HSd-         </td> <td align="center">   33.0729       </td> </tr>
<tr><td align="center">   Mn+2       </td> <td align="center">   Mn           </td> <td align="center">   54.938        </td> <td align="center">   Hdg       </td> <td align="center">   Hdg          </td> <td align="center">   2.0158        </td> </tr>
<tr><td align="center">              </td> <td align="center">                </td> <td align="center">                 </td> <td align="center">   Ndg       </td> <td align="center">   Ndg          </td> <td align="center">   28.0134       </td> </tr>
</tbody>
</table>

