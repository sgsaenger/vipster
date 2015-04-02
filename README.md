#PW Toolbox

Visualization of PWScf, Cube and xyz files.

Depends on Python2.7, Numpy, f2py, gfortran, PyQT4 and the Python OpenGL bindings.
(https://pypi.python.org/pypi/PyOpenGL/3.0.2#downloads)

##Installation:

Compile the fortran part and you're good to go:
```
f2py -c faster.f90 -m mol_f
```

##Usage:

For startup:
```
./ptb_main.py -h
```

The GUI should be pretty self-explanatory (or will be at some point...)

###Mouse:

**Normal mode** (R):
Left-click: Rotate
Middle-click: Move
Right-click: Set back view

**Select mode** (S):
Left-click: Add atom to selection
Middle-click: Remove last atom from selection
Right-click: Empty selection

###Tools:

**Pick:**

Displays information for selected atoms


**Script:**

Interface for more complex actions:

- **def**ine (list) name: define a named group of atoms
- **rep**eat num: repeat the following operation num times
- **shi**ft (list) (vec): shift the given list of atoms along a given direction
- **rot**ate (list) angle (vec) (shift): rotates around vector vec with offset shift
- **mir**ror (list) (vec1) (vec2) (shift): mirrors at a plane given by the vectors 1 and 2, offset by shift

Arguments can be given as follows:

(list) of atoms can be either:
- defined name
- single index: idx/[idx]
- range: i1-i2/[i1-i2]
- explicit list: [i1,i2,i3,i4-i5]

mandatory vector (vec) or optional vector (shift) can be one of:
- position vector of atom: idx
- path between atoms: i1-i2
- explicit vector in bohr: (x,y,z)
- explicit vector in other formats: (x,y,z,'format')
  with format being one of bohr,angstrom,crystal,alat


**Mult. Cell:**

Permanently multiplies the unit cell by the given multiplicators.
