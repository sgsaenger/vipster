FSources = fortran/iso_surf_lut.f90 fortran/make_iso_surf.f90 fortran/set_bonds_f.f90

mol_f.so: $(FSources)
	f2py -c $(FSources) -m mol_f
