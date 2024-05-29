# cahn-hilliard-2D-spectral
Solution of coupled Cahn-Hilliard equations in multiple dimensions.

The dynamics of the state of the two mixed systems under consideration, $u$ and $v$ evolves to minimize the value of
an energy functional in the following expression: 



$$
F_{\epsilon_u, \epsilon_v ,\sigma }(u,v) = \int_{\Omega} \left( \frac{\epsilon_u^2}{2} |\nabla u|^2 + \frac{\epsilon_v^2}{2} |\nabla v|^2 + W(u,v) + \frac{\sigma}{2} \left| (-\Delta)^{-1/2} (v - \overline{v}) \right|^2 \right) \, dr
$$



where

$$W\left( u,v\right) =\frac{\left( u^{2}-1\right) ^{2}}{4}+\frac{\left(
v^{2}-1\right) ^{2}}{4}+b_1 uv-b_{2}\frac{uv^{2}}{2}$$

