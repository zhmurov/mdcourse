Numerical integration of equations of motion 
============================================

    .. math::

        m\frac{d^{2}\mathbf{r}_{i}}{dt^{2}}=\mathbf{f}_{i}
        
        \frac{d\mathbf{r}_i}{dt}=\mathbf{v}_{i}

Using Taylor expansion:

    .. math::

        \mathbf{r}_{i}(t+\Delta t)=\mathbf{r}_{i}(t)+\mathbf{v}_{i}(t)\Delta t+\frac{\mathbf{f}_{i}(t)}{2m}\Delta t^{2}+\frac{1}{6}\frac{d^{3}\mathbf{r}_{i}(t)}{dt^{3}}\Delta t^{3}+\mathbf{o}(\Delta t^{4})
        
        \mathbf{r}_{i}(t-\Delta t)=\mathbf{r}_{i}(t)-\mathbf{v}_{i}(t)\Delta t+\frac{\mathbf{f}_{i}(t)}{2m}\Delta t^{2}-\frac{1}{6}\frac{d^{3}\mathbf{r}_{i}(t)}{dt^{3}}\Delta t^{3}+\mathbf{o}(\Delta t^{4})

Adding these up:

    .. math::

        \mathbf{r}_{i}(t+\Delta t)=2\mathbf{r}_{i}(t)-\mathbf{r}_{i}(t-\Delta t)+\frac{\mathbf{f}_{i}(t)}{m}\Delta t^{2}+\mathbf{o}(\Delta t^{4}).

To compute velocities, one can use:

    .. math::

        \mathbf{v}_{i}(t)=\frac{\mathbf{r}_{i}(t+\Delta t)-\mathbf{r}_{i}(t-\Delta t)}{2\Delta t}+\mathbf{o}(\Delta t^{2}).
