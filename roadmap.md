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

# In the long run, wishlist:

- [x] Make GUI optional
- [ ] Improve orthogonal projection
  - [ ] display without zooming first
  - [ ] make it zoom on the same scale as perspective projection
- [ ] Add more editing capabilities:
  - [ ] multiply cells automatically
  - [ ] Coord editor right click menu
  - [ ] PW Editor editing capabilities
- [ ] Speed!
  - [ ] remove coordination calculation with offset from paintGL routine
- [ ] Improve rotation
- [ ] Bonds:
  - [ ] between periodic images
    - [ ] orientation! check if +/-
  - [ ] speed up
    - maybe initalize upon reading the file, not loading the configuration?
- [ ] Checkboxes that applies editing to all molecules in mol-list
- [ ] Default PW Parameter sets
- [ ] PW Parameters from PW Output
- [ ] Actually handle multiple KPoint formats
- [ ] Input validators
