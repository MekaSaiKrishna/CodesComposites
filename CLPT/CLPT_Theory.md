## Classical Laminate Plate Theory (or) Kirchhoff–Love theory

The Kirchhoff–Love theory of plates is a two-dimensional mathematical model that is used to determine the stresses and deformations in thin plates subjected to forces and moments. This theory is an extension of Euler-Bernoulli beam theory and was developed in 1888 by Love using assumptions proposed by Kirchhoff. The theory assumes that a mid-surface plane can be used to represent a three-dimensional plate in two-dimensional form.

 Checkout the [Wikipedia Link](https://en.wikipedia.org/wiki/Kirchhoff%E2%80%93Love_plate_theory) about Kirchoff-Love Theory for more details and derivations

| Euler-Bernoulli Theory      | Kirchoff-Love Theory |
| ----------- | ----------- |
| analyze the behavior of slender beams subjected to bending loads      | analyze the behavior of thin plates subjected to in-plane loads, such as bending, stretching, and twisting       |
| **Assumptions:** <br> * The plane sections remain plane assumption i.e. the cross-section remains plane and perpendicular to the neutral axis during deformation. [link](https://learnaboutstructures.com/Bernoulli-Euler-Beam-Theory) <br> ![image](https://github.com/MekaSaiKrishna/CodesComposites/assets/93347557/32e6c2b2-4112-4a42-b86a-3c863d24be75) <br> *hello| **Assumptions:** <br> * straight lines normal to the mid-surface remain straight after deformation <br> * straight lines normal to the mid-surface remain normal to the mid-surface after deformation <br> *  the thickness of the plate does not change during a deformation <br> <p align="center"> <img width="240" src="https://github.com/MekaSaiKrishna/CodesComposites/assets/93347557/76329864-d536-4fc4-8e9c-e410fffb2b3f"> </p>|


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

<p align="center">
  <img width="240" src="https://github.com/MekaSaiKrishna/CodesComposites/blob/main/CLPT/interactive_plot.html">
</p>

<iframe src="https://github.com/MekaSaiKrishna/CodesComposites/blob/main/CLPT/interactive_plot.html" width="800" height="600"></iframe>

