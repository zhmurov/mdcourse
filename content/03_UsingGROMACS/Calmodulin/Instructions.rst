MD simulations of calmodulin 
============================

Set up gmx if not in the path:

    .. code-block:: shell

        GMX=/usr/local/gromacs/bin/gmx 

Get the CHARMM36 forcefield:

    .. code-block:: shell

        wget http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz -O charmm36-jul2022.ff.tgz
        tar -xvf charmm36-jul2022.ff.tgz
        rm charmm36-jul2022.ff.tgz

Get the PDB file from https://www.rcsb.org/structure/1cfd :

    .. code-block:: shell

        wget https://files.rcsb.org/download/1CFD.pdb

Pre-process the file, select forcefield:

    .. code-block:: shell
        
        $GMX pdb2gmx -f 1CFD.pdb -o processed.gro -ignh

Chose CHARMM from local folder (usually option 1), select TIP3P for water.

Define the solvation box:

    .. code-block:: shell

        $GMX editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic

Solvate:

    .. code-block:: shell
        
        $GMX solvate -cp newbox.gro -o solv.gro -p topol.top

Add ions:

    .. code-block:: shell
        
        touch ions.mdp
        $GMX grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
        $GMX genion -s ions.tpr -o solv_CAL.gro -p topol.top -np 4 -pname CAL -pq 2 -nn 0
        $GMX grompp -f ions.mdp -c solv_CAL.gro -p topol.top -o ions.tpr
        $GMX genion -s ions.tpr -o solv_ions.gro -conc 0.15 -p topol.top -pname NA -nname CL -neutral

Select group SOL

Get the mdp files for energy minimization, nvt and npt equilibration and production md.

Run energy minimization:

    .. code-block:: shell
        
        $GMX grompp -f emin-charmm.mdp -c solv_ions.gro -o em.tpr
        $GMX mdrun -v -deffnm em

Get the energy plot (to get plain .dat file, use `-xvg none` option):

    .. code-block:: shell

        $GMX energy -f em.edr -o potential.xvg

Run temperature equilibration:

    .. code-block:: shell

        $GMX grompp -f nvt-charmm.mdp -c em.gro -r em.gro -o nvt.tpr
        $GMX mdrun -v -deffnm nvt

Check the temperature:

    .. code-block:: shell

        $GMX energy -f nvt.edr -o temperature.xvg

Equilibrate pressure:

    .. code-block:: shell

        $GMX grompp -f npt-charmm.mdp -c nvt.gro -r nvt.gro -o npt.tpr
        $GMX mdrun -v -deffnm npt

Check the pressure:

    .. code-block:: shell

        $GMX energy -f npt.edr -o pressure.xvg

Production run:

    .. code-block:: shell

        $GMX grompp -f md-charmm.mdp -c npt.gro -t npt.cpt -o md.tpr
        $GMX mdrun -v -deffnm md




