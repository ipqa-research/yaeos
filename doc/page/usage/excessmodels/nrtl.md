---
title: NRTL
---

# NRTL

The Non-Random-Two-Liquid model presented by Renon and Prausnitz.
This model uses local compositions to represent the excess Gibbs energy.

\[
   G^E = nRT \cdot \sum_i x_i \frac{\sum_j x_j \tau_{ji} G_{ji}}{\sum_j x_j G_{ji}}
\]

with:

\[\tau_{ij} = A_{ij} + \frac{B_{ij}}{T}\]

\[G_{ij} = exp(-C\tau_{ij})\]

Where \(A_{ij}, A_{ji}, B_{ij}, B_{ji} \text{ and } C\) are fittable parameters.

1. H. Renon and J. M. Prausnitz, “Local compositions in thermodynamic excess
functions for liquid mixtures,” AIChE Journal, vol. 14, no. 1, pp. 135–144,
1968, doi: ![10.1002/aic.690140124.](https://onlinelibrary.wiley.com/doi/abs/10.1002/aic.690140124)
