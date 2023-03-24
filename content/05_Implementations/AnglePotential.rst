Three-body potential
====================

    .. figure:: img/angle_eng.png
        :name: fig:ff-angle-1
  
        Schematic representation of the angle potential (left panel) and the shape of the potential function as a function of angle :math:`\theta` (right panel). The minimum of the harmonic potential corresponds to the angle :math:`\theta_0`.

Introducing a triplet of atoms :math:`(i,j,k)`, forming a particular angle :numref:`fig:ff-angle-1` and a set of angles in a system, :math:`A`, the angle potential can be written as:

    .. math::
        :label: eq:ff-angles2
        
        V_{angles} = \sum_{(i,j,k) \in A}\frac{k_{ijk}^{\theta}}{2}(\theta_{ijk}-\theta^0_{ijk})^2

The force, acting on a particular atom :math:`l` due to the angle potential is then

    .. math::
        :label: eq:ff-angles3

        \begin{split}
        \mathbf{f}_{l}=-\nabla_{l}V_{angles} &= -\nabla_{l}\left(\sum_{(i,j,k) \in A}\frac{k_{ijk}^{\theta}}{2}(\theta_{ijk}-\theta^0_{ijk})^2\right)\\
        &= -\sum_{(i,j,k) \in A}k_{ijk}^{\theta}(\theta_{ijk}-\theta^0_{ijk})\nabla_{l}\theta_{ijk}
        \end{split}

Since, from the atomic coordinates, it is easier to compute :math:`\cos\theta_{ijk}` and :math:`\sin\theta_{ijk}`, let us switch to a :math:`\cos` and :math:`\sin` representation as follows:

    .. math::
        :label: eq:ff-angles4
        
        \theta_{ijk}=\arccos\left(\cos\theta_{ijk}\right),

and
    .. math::
        :label: eq:ff-angles5

        \nabla_{l}\theta_{ijk} = \nabla_{l}\left(\arccos\left(\cos\theta_{ijk}\right)\right) = -\frac{1}{\sqrt{1-\cos^{2}\theta_{ijk}}}\nabla_{l}\cos\theta_{ijk} = -\frac{\nabla_{l}\cos\theta_{ijk}}{\sin\theta_{ijk}}.

Let us denote vectors that connect particles :math:`i`, :math:`j` and :math:`k` in a triplet as :math:`\mathbf{r}_{ji}=\mathbf{r}_{i}-\mathbf{r}_{j}` and :math:`\mathbf{r}_{jk}=\mathbf{r}_{k}-\mathbf{r}_{j}`; and the respective distances as :math:`r_{ji}=|\mathbf{r}_{ji}|` and :math:`r_{jk}=|\mathbf{r}_{jk}|` Then

    .. math::
        :label: eq:ff-angles6

        \cos\theta_{ijk}=\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{|\mathbf{r}_{ji}||\mathbf{r}_{jk}|}=\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}r_{jk}}

Here, :math:`\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}` is a dot product of vectors :math:`\mathbf{r}_{ji}` and :math:`\mathbf{r}_{jk}`. It is clear, that :math:`\nabla_{l}\cos\theta_{ijk}` will vanish if :math:`l` does not belong to the triplet :math:`(i,j,k)`. Considering all three cases (:math:`l=i`, :math:`l=j` and :math:`l=k`) separately and using Eq. :eq:`eq:ff-angles6` we obtain:

    .. math::
        :label: eq:ff-angles7


        \begin{split}
        \nabla_{l}\cos\theta_{ijk}\bigg|_{l=i}&=\nabla_{l}\left(\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}r_{jk}}\right)\bigg|_{l=i}=\left(\frac{1}{r_{ji}r_{jk}}\nabla_{l}\left(\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}\right)-\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}^{2}r_{jk}}\nabla_{l}r_{ji}\right)\bigg|_{l=i}=\\
        &=\frac{1}{r_{ji}r_{jk}}\mathbf{r}_{jk}-\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}^{2}r_{jk}}\frac{\mathbf{r}_{ji}}{r_{ji}}=\frac{1}{r_{ji}}\left[\frac{\mathbf{r}_{jk}}{r_{jk}}-\cos\theta_{ijk}\frac{\mathbf{r}_{ji}}{r_{ji}}\right]\text{,}\\
        \end{split}

    .. math::
        :label: eq:ff-angles8
        
        \begin{split}
        \nabla_{l}\cos\theta_{ijk}\bigg|_{l=k}&=\nabla_{l}\left(\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}r_{jk}}\right)\bigg|_{l=k}=\left(\frac{1}{r_{ji}r_{jk}}\nabla_{l}\left(\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}\right)-\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}r_{jk}^{2}}\nabla_{l}r_{jk}\right)\bigg|_{l=k}=\\
        &=\frac{1}{r_{ji}r_{jk}}\mathbf{r}_{ji}-\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}r_{jk}^{2}}\frac{\mathbf{r}_{jk}}{r_{jk}}=\frac{1}{r_{jk}}\left[\frac{\mathbf{r}_{ji}}{r_{ji}}-\cos\theta_{ijk}\frac{\mathbf{r}_{jk}}{r_{jk}}\right]\text{, and}\\
        \end{split}

    .. math::
        :label: eq:ff-angles9
        
        \begin{split}
        \nabla_{l}\cos\theta_{ijk}\bigg|_{l=j}&=\nabla_{l}\left(\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}r_{jk}}\right)\bigg|_{l=j}=\\
        &=\left(\frac{1}{r_{ji}r_{jk}}\nabla_{l}\left(\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}\right)-\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}^{2}r_{jk}}\nabla_{l}r_{ji}-\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}r_{jk}^{2}}\nabla_{l}r_{jk}\right)\bigg|_{l=j}=\\
        &=-\frac{1}{r_{ji}r_{jk}}\mathbf{r}_{jk}-\frac{1}{r_{ji}r_{jk}}\mathbf{r}_{ji}+\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}^{2}r_{jk}}\frac{\mathbf{r}_{ji}}{r_{ji}}+\frac{\mathbf{r}_{ji}\cdot\mathbf{r}_{jk}}{r_{ji}r_{jk}^{2}}\frac{\mathbf{r}_{jk}}{r_{jk}}=\\
        &=\frac{1}{r_{ji}}\left[\cos\theta_{ijk}\frac{\mathbf{r}_{ji}}{r_{ji}}-\frac{\mathbf{r}_{jk}}{r_{jk}}\right] + \frac{1}{r_{jk}}\left[\cos\theta_{ijk}\frac{\mathbf{r}_{jk}}{r_{jk}}-\frac{\mathbf{r}_{ji}}{r_{ji}}\right]=\\
        &=-\nabla_{l}\cos\theta_{ijk}\bigg|_{l=i}-\nabla_{l}\cos\theta_{ijk}\bigg|_{l=k}
        \end{split}

Summirizing Eqs. :eq:`eq:ff-angles3`, :eq:`eq:ff-angles5` and :eq:`eq:ff-angles7` to :eq:`eq:ff-angles9` we see that, one can compute three components of force :math:`\mathbf{f}_{i}`, :math:`\mathbf{f}_{j}` and :math:`\mathbf{f}_{k}`, acting on each atom in the triplet :math:`(i,j,k)` due to the angle potential between atoms :math:`i`, :math:`j` and :math:`k` using the following relations:

    .. math::
        :label: eq:ff-angles10
        
        \begin{split}
        \mathbf{f}_{i}&=k_{ijk}^{\theta}(\theta_{ijk}-\theta^0_{ijk})\left(-\frac{1}{\sin\theta_{ijk}}\right)\frac{1}{r_{ji}}\left[\cos\theta_{ijk}\frac{\mathbf{r}_{ji}}{r_{ji}}-\frac{\mathbf{r}_{jk}}{r_{jk}}\right]\\
        \mathbf{f}_{k}&=k_{ijk}^{\theta}(\theta_{ijk}-\theta^0_{ijk})\left(-\frac{1}{\sin\theta_{ijk}}\right)\frac{1}{r_{jk}}\left[\cos\theta_{ijk}\frac{\mathbf{r}_{jk}}{r_{jk}}-\frac{\mathbf{r}_{ji}}{r_{ji}}\right]\\
        \mathbf{f}_{j}&=-\mathbf{f}_{i}-\mathbf{f}_{k}
        \end{split}


Hence, the forces acting on three atoms connected by the angle potential are related. Consequently it is more efficient to use the modification of the potential pair based parallelization approach.
