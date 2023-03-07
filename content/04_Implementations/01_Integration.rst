Numerical integration of equations of motion 
============================================

    .. math::

        m\frac{d^{2}\vec{r}_{i}}{dt^{2}}=\vec{f}_{i}
        
        \frac{d\vec{r}_i}{dt}=\upsilon_i

Using Taylor expansion:

    .. math::

        \vec{r}_{i}(t+\Delta t)=\vec{r}_{i}(t)+\vec{\upsilon}-{i}(t)\Delta t+\frac{\vec{f}_{i}(t)}{2m}\Delta t^{2}+\frac{1}{6}\frac{d^{3}\vec{r}_{i}(t)}{dt^{3}}\Delta t^{3}+\vec{o}(\Delta t^{4})
        
        \vec{r}_{i}(t-\Delta t)=\vec{r}_{i}(t)-\vec{\upsilon}_{i}(t)\Delta t+\frac{\vec{f}_{i}(t)}{2m}\Delta t^{2}-\frac{1}{6}\frac{d^{3}\vec{r}_{i}(t)}{dt^{3}}{6}\Delta t^{3}+\vec{o}(\Delta t^{4})

Adding these up:

    .. math::

        \vec{r}_{i}(t+\Delta t)=2\vec{r}_{i}(t)-\vec{r}_{i}(t-\Delta t)+\frac{\vec{f}_{i}(t)}{m}\Delta t^{2}+\vec{o}(\Delta t^{4}).

To compute velocities, one can use:

    .. math::

        \vec{\upsilon}_{i}(t)=\frac{\vec{r}_{i}(t+\Delta t)-\vec{r}_{i}(t-\Delta t)}{2\Delta t}+\vec{o}(\Delta t^{2}).
