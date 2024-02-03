## Classical Laminate Plate Theory (or) Kirchhoff–Love theory

The Kirchhoff–Love theory of plates is a two-dimensional mathematical model that is used to determine the stresses and deformations in thin plates subjected to forces and moments. This theory is an extension of Euler-Bernoulli beam theory and was developed in 1888 by Love using assumptions proposed by Kirchhoff. The theory assumes that a mid-surface plane can be used to represent a three-dimensional plate in two-dimensional form.

The following kinematic assumptions that are made in this theory:

1) straight lines normal to the mid-surface remain straight after deformation
2) straight lines normal to the mid-surface remain normal to the mid-surface after deformation
3) the thickness of the plate does not change during a deformation.

| Euler-Bernoulli Theory      | Kirchoff-Love Theory |
| ----------- | ----------- |
| Header      | Title       |
| Paragraph   | Text        |


<!-- Kirchoff Plate Theory Schematic --> 
<p align="center">
  <img width="240" src="https://github.com/MekaSaiKrishna/CodesComposites/assets/93347557/c5e84442-fe5e-448d-bbd8-cf9c842b2668">
</p>

___
$$\begin{equation}
\begin{bmatrix}
\sigma_{11} \\
\sigma_{22} \\
\sigma_{12} \\
\end{bmatrix} = 
\begin{bmatrix}
C_{11} & C_{12} & C_{13} \\
C_{12} & C_{22} & C_{23} \\
C_{13} & C_{23} & C_{33} \\
\end{bmatrix}
\begin{bmatrix}
\varepsilon_{11} \\
\varepsilon_{22} \\
\varepsilon_{12} \\
\end{bmatrix}
\end{equation}$$

___
<iframe src="interactive_plot.html" width="800" height="600"></iframe>

