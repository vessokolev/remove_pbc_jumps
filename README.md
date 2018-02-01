## ``REMOVE_PBC_JUMPS`` - PBC jumps remover

### This is a "bonzai" style Python 3 demo code which goal is to remove the PBC jumps from trajectory frames saved as GRO files.

#### Author: Veselin Kolev <vesso.kolev@gmail.com>
#### License: GPLv2
#### Version: 2018013100

#### Content:

#### 1. Introduction.
#### 2. How to download the source code of the project.
#### 3. Requirements.
#### 4. Executing the code.


_1. Introduction_

Imposing periodic boundary conditions (PBC) in one or more dimensions, is a popular technique to create a virtual system of infinite size by periodically repeating one sample (simulation box) of a finite size. That idea is widely implemented in almost any modern molecular dynamics simulation software, like GROMACS. When a molecule is passing through the virtual wall of the simulation box it gets its chain split to be able to "re-enter" the simulation box through the opposite wall. That split does not represent any real physical process and does not affect the energies or the bond lengts. But when saving the simulation frames that PBC jump might complicate or distort the analysis. To avoid that the PBC jumps need to be removed and the molecule - assembled back to a single object. Of course, one can always resolve fast and easy the PBC jumps by means of the GROMACS command line tools. That is why the goal of providing this simple Python 3 code is only to help everyone interested in removing of the PBC jumps to understand the idea in details.

Feel free to modify and adapt the code to fit your requirements.


_2. How to download the source code of the project._

The preferable method for obtaining the code is to use the tool ``git`` and clone the source tree locally to your file system:

```
git clone https:////github.com/vessokolev/remove_pbc_jumps.git
```

You may download the source as a ZIP-archive by pressing the button "Clone or download" in the web-interface on GitHub, or by using wget:

```
wget https://github.com/vessokolev/remove_pbc_jumps/archive/master.zip
```


_3. Requirements._

The code requires recent versions of both Python 3 and NumPy.


_4. Executing the code._

Sumply make the file ``remove_pbc_jumps.py`` executable:

$ chmod 755

and execute it:

```
$ ./remove_pbc_jumps.py
```

The script can be always executed by passing the path to the file as an argument to the ``python3`` executable.

```
$ python3 remove_pbc_jumps.py
```

Use VMD to compare the initial and created GRO files as frames of a trajectory:

```
$ vmd split.gro assembled.gro
```

