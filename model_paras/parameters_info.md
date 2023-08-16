### Model parameters info



Unit used in Model:

Time: $h$; Length: $\mu m$; Concentration: $mM$; Mass: $pg$





| Note                | Description                                                  | Value                                             | Unit                              |
| ------------------- | ------------------------------------------------------------ | ------------------------------------------------- | --------------------------------- |
| $\lambda_s$         | `maxGrowthRate` Batch culture growth rate.                   | 0.9$^{*}$                                         | $h^{-1}$                          |
| $Y$                 | `yield` yield factor                                         | 0.5                                               | $\mu g /\mu g$                    |
| $\rho_0$            | `Rho0` local cell mass density                               | 0.102                                             | $p g \cdot \mu m ^{-3}$           |
| $l_{div}$           | `L_divide`                                                   |                                                   | $\mu m$                           |
| $r_{cell}$          | `Radius` $w_0 = 2 r_{cell}$                                  |                                                   | $\mu m$                           |
| $\sigma$            | `SigmaMoL` the convert factor for cell elongation rate to   mass growth rate. | $\frac{\ln{(2 l_{div} / (l_{div} - w_0))}}{\ln2}$ | -                                 |
| $K_c$               | `KC` capacity constant (the Monod constant) for sugar        | 0.02                                              | $mM$                              |
| $C_{\mathrm{max}}$  | `maxCarbon` Maximum concentration of carbon source           | 11.1                                              | $mM$                              |
| $\Delta t$          | `dt` step length for cell movements                          | 0.00005                                           | $h^{-1}$                          |
| $\Delta t_c$        | `Cdt` step length for update carbon source concentration     | 0.05                                              | $h^{-1}$                          |
| $N_x$               | `Box_x` number of grids in the $x$ direction                 | 200                                               | -                                 |
| $N_y$               | `Box_y ` number of grids in the $y$ direction                | 200                                               | -                                 |
| $N_{zc}$            | `Box_z` number of grids in the $z$ direction in colony       | 50                                                | -                                 |
| $N_{za}$            | `Box_z_agar`number of grids in the $z$ direction in agar     | 20                                                | -                                 |
| $D_{+}$             | `DiffColony`                                                 | 0.025                                             | $\mu m ^{2} \cdot h^{-1}$         |
| $D_{-}$             | `DiffAgar`                                                   | 0.167                                             | $\mu m ^{2} \cdot h^{-1}$         |
|                     |                                                              |                                                   |                                   |
| $C_{\mathrm{rate}}$ | `C_rate` carbon consumption rate                             | $\frac{\lambda_s \cdot \rho_0}{Y \cdot M_r}$      | $10^{9} \mathrm{mM \cdot h^{-1}}$ |







|                 |                                       |                                     |       | Unit                              |                                                              |
| --------------- | ------------------------------------- | ----------------------------------- | ----- | --------------------------------- | ------------------------------------------------------------ |
| $h_{grid}$      | `BoxLength`                           | physical size of each box           | NA    | $\mu m$                           | We set a dynamic spatial size of grids, cell length + cell width. |
| $\phi_0^{max}$  | `Density_Threshold`                   | maximum volume fractions of cells   | 0.70  | NA                                | see sup.2 in paper, 0.68 + 0.03                              |
| $\gamma_{surf}$ | `Surface_Tension`, `tension ` in code | surface tension constant            | 150   | $pg \cdot h^{-1}$                 | see Sup.2                                                    |
| $k_cc$          | `k_cc`                                | cell-cell Hertzian elastic constant | 30000 | $pg \cdot\mu m^{-1} \cdot h^{-2}$ | see Sup.2                                                    |
|                 |                                       |                                     |       |                                   |                                                              |

$^{*}$:  This is a representative value.



* carbon consumption rate:
  $$
  C_{\mathrm{rate}} = \frac{\lambda_s \cdot \rho_0}{Y \cdot M_r}
  $$

​	$M_r = 180$  $ \mathrm{g \cdot mol^{-1}}$ Unit: $10^{-6} \mathrm{mmol^{-1} \cdot h^{-1}} \cdot \mu m ^{-3}$ = $10^{9} \mathrm{mM \cdot h^{-1}}$

