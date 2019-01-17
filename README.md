## Particle Reducer

### About

This software provides reducing the number of particles in the ndf file using the Voronoi algorithm[1]

The software is written in python 3, python 2 is not supported. For the rest of the document by python we assume some python3 binary.
### Installing

Clone this repository by running
```
git clone https://github.com/KseniaBastrakova/particle_reduction
```

or download an archive from the github page and unpack it.

The only dependency is h5py, which can be installed e.g., by
pip install h5py


### Particle reducing 
```
run voronoi_reduction(hdf_file, hdf_file_reduction, tolerance_momentum, tolerance_position)
```

where
* hdf_file is the path to an input HDF file;
* hdf_file_reduction is the path to an output reduction hdf5 file, by default result.h5;
* tolerance_momentum is the tolerance for momentum (We consider all values of different particles for momentum to be one particle);
* tolerance_position is the tolerance for position (We consider all values of different particles for position to be one particle);
	
[1] Voronoi particle merging algorithm for PIC codes Phuc T.Luu T.Tückmantel A.Pukhov
