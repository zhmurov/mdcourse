Harmonic potential
==================

    .. math::

        V_{bond}=\sum_{(i,j)\in bonds}\frac{k^s_{ij}}{2}(r_{ij}-r_{ij}^0)^2=\sum_{(i,j)\in bonds}V_{bond}^{ij}

    .. math::

        V_{bond}^{ij}=\frac{k^s_{ij}}{2}(r_{ij}-r_{ij}^0)^2

    .. math::

        \nabla_iV_{bond}^{ij}=\frac{\partial V_{bond}}{\partial r_{ij}}\nabla r_{ij} = -k^s_{ij}(r_{ij}-r_{ij}^0)\frac{\mathbf{r}_{ij}}{r_{ij}}