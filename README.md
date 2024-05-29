# cahn-hilliard-2D-spectral
Solution of coupled Cahn-Hilliard equations in multiple dimensions.

The dynamics of the state of the two mixed systems under consideration, $u$ and $v$ evolves to minimize the value of
an energy functional in the following expression: 



$$
F_{\epsilon_u, \epsilon_v ,\sigma }\left( u,v\right) = \int\limits_{\Omega }\left\{ \frac{\epsilon _{u}^{2}}{2}\left\vert \nabla u\right\vert ^{2}+\frac{\epsilon _{v}^{2}}{2}\left\vert \nabla v\right\vert ^{2}+W\left(u,v\right) +\frac{\sigma }{2}\left\vert \left( -\Delta \right)^{-1/2}\left( v-\overline{v}\right) \right\vert ^{2}\right\} \, dr,
$$


where

$$W\left( u,v\right) =\frac{\left( u^{2}-1\right) ^{2}}{4}+\frac{\left(
v^{2}-1\right) ^{2}}{4}+b_1 uv-b_{2}\frac{uv^{2}}{2}$$

