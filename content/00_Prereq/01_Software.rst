Software
========

.. _gromacs-installation:

Build and install GROMACS
-------------------------

1. Get the code from GROMACS web page:

    .. code-block:: shell

        wget https://ftp.gromacs.org/gromacs/gromacs-2022.5.tar.gz
        tar -xvf gromacs-2022.5.tar.gz
        cd gromacs-2022.5

Alternatively, you can checkout required version using Git:

    .. code-block:: shell

        git clone https://gitlab.com/gromacs/gromacs.git
        cd gromacs
        git checkout v2022.5

2. Build and install

    .. code-block:: shell
        
        mkdir build
        cmake -S. -Bbuild
        cd build
        make -j4
        sudo make install

Optionally, one can set the installation path by adding ``-DCMAKE_INSTALL_PREFIX=/your/desired/path`` to the ``cmake`` command.

3. You can optionally check the build before installing by running the tests:

    .. code-block:: shell
        
        make check -j4

4. It is convenient to set a variable for GROMACS executable. If you install into the default location:

    .. code-block:: shell

        GMX=/usr/local/gromacs/bin/gmx


Force-fields
------------

.. _trappeua-installation:

TraPPE-UA
^^^^^^^^^

Get the forcefield files:

    .. code-block:: shell

        git clone https://github.com/wesbarnett/trappeua.git

``gmx pdb2gmx`` utility uses enumeration for the water model selection.
Hence, one need the water mode to have one of the pre-defined names, which is not the case in the ``trappeua`` implementation we just downloaded.
In order to use pdb2gmx later on, we now need to rename the water model file for tip4p model and change this name in the watermodels.dat file before we copy forcefield to the GROMACS installation folder:

    .. code-block:: shell

        mv trappeua/trappeua.ff/tip4p2005.itp trappeua/trappeua.ff/tip4p.itp
        sed -i 's/tip4p2005/tip4p/' trappeua/trappeua.ff/watermodels.dat

To install the forcefield, copy the result to GROMACS installation folder.
Assuming that GROMACS in installed at ``/usr/local/gromacs``:

    .. code-block:: shell
        
        sudo cp -pr trappeua/trappeua.ff /usr/local/gromacs/share/gromacs/top/

Note that you can also keep force-field files in your local folder without installing.

.. _charmm36-installation:

CHARMM36
^^^^^^^^

Download and extract the force-field files:

    .. code-block:: shell

        wget https://www.charmm.org/archive/charmm/resources/charmm-force-fields/download.php?filename=CHARMM_ff_params_files/archive/charmm36-mar2019.ff.tgz
        tar -xvf charmm36-mar2019.ff.tgz

Or clone `this repository <https://gitlab.com/artemzhmurov/charmm36>`_:

    .. code-block:: shell

        git clone git@gitlab.com:artemzhmurov/charmm36.git

To install the forcefield, copy it to your GROMACS installation folder.
Assuming that GROMACS in installed at ``/usr/local/gromacs``:

    .. code-block:: shell
        
        sudo cp -pr trappeua/charmm36.ff /usr/local/gromacs/share/gromacs/top/

Note that you can also keep force-field files in your local folder without installing.

.. _packmol-installation:

PackMol
-------

To create the coordinates for a box of molecules, we can use Packmol software.
You will need ``gfortran``, which you can install by running ``sudo apt install gfortran``.
To get and install Packmol:

    .. code-block:: shell

        git clone https://github.com/m3g/packmol.git
        cd packmol
        git checkout v20.3.5
        ./configure
        make
        PACKMOL=$(pwd)/packmol

.. _vmd-installation:

VMD
---

Visual Molecular Dynamics (VMD) is a visualization program, that should be installed locally.
You can get a copy `here <https://www.ks.uiuc.edu/Research/vmd/>`_.