# TODO for first release V0.1:

- [x] Parse shell arguments for file-loading
- [x] Bonds:
  - [x] reset upon modification
- [x] Edit window:
  - [x] reenable PW Parameter tab
    - [x] reenable KPoint handling
  - [x] make coordinate editor not waste resources if not shown
  - [x] delete atoms in the editor
- [x] Fix PWI Write function

# V0.2:

- [x] Speed!
  - [x] remove coordination calculation with offset from paintGL routine
- [x] PBC Bonds
- [x] Make GUI optional
- [x] Fix C&P and Delete functions for multiple selection ranges

# V0.3:

- [x] Reduce bond calculations
- [x] Recalculate bonds upon edit when necessary
- [x] Checkboxes that applies editing to all molecules in mol-list
- [x] Actually handle multiple KPoint formats
- [x] Add more editing capabilities:
  - [x] Coord editor right click menu
  - [x] PW Editor editing capabilities
  - [x] Better commandline parsing

# V0.4:

- [x] Improve rotation
- [x] Improved Keyboard controls:
  - [x] GLWidget can retain focus
  - [x] Rotate with arrow keys
- [x] Screenshots
- [x] Vastly improved speed of bond-setting routine

# V0.5:

- [x] Improved startup by removing unnecessary bond calculations
- [x] Fixed bug asking for PW parameters when saving to xyz
- [x] Select Atoms graphically
- [x] parse different bravais lattice specifications from PWInput
- [x] privatised functions in backend
- [x] added function to get all atoms at once
- [x] exported model of sphere to external file
- [x] exported model of bonds to external file
- [ ] added doxygen documentation
- [ ] moved classes to seperate files
- [ ] moved k-points to molecule (including molecule edit area)
- [ ] added VBOs and instanced rendering:
  - [x] for atoms
  - [x] for bonds
  - [ ] for everything else
- [ ] some kind of visual representation of volume data:
  - [ ] 2D
  - [ ] 3D point cloud
  - [ ] 3D isosurface



# In the long run, wishlist:

- [ ] Calculate distances and angles
- [ ] Small scripting language with point&click support
- [ ] Better picture generation
  - [ ] Better screenshot dialog
  - [ ] Export to Povray
- [ ] Bond setting has to be faster
  - [ ] Move to fortran code
  - [ ] Threads to improve responsiveness
- [ ] Improve orthogonal projection
  - [x] display without zooming first
  - [ ] make it zoom on the same scale as perspective projection
- [ ] Commandline creation of new Files
- [ ] Default PW Parameter sets
- [ ] PW Parameters from PW Output
- [ ] Input validators
- [ ] Reduce Model polygons with increasing model count
