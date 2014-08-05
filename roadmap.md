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
- [ ] Recalculate bonds upon edit when necessary (always?)
- [x] Checkboxes that applies editing to all molecules in mol-list
- [x] Actually handle multiple KPoint formats
- [ ] Add more editing capabilities:
  - [ ] multiply cells automatically
  - [x] Coord editor right click menu
  - [x] PW Editor editing capabilities
  - [ ] Commandline creation of new Files

# In the long run, wishlist:

- [ ] Improve orthogonal projection
  - [ ] display without zooming first
  - [ ] make it zoom on the same scale as perspective projection
- [ ] Improve rotation
- [ ] Default PW Parameter sets
- [ ] PW Parameters from PW Output
- [ ] Input validators
- [ ] Select Atoms graphically
- [ ] Reduce Model polygons with increasing model count
