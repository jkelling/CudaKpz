# CudaKPZ

The CudaKPZ package includes program to simulate [Kardar-Parisi-Zhang][KPZ] Surface
growth using the [Octahedron][oct] model using random-sequential (RS) or
stochastic-cellular-automaton (SCA) updates. The programs calculate the
roughness, moments of the height distribution and autocorrelation functions.

## License

When using results generated using programs from this package, please cite

- Kelling, Jeffrey, Ódor, Géza, Gemming, Sibylle: Suppressing correlations in massively
  parallel simulations of lattice models, <https://arxiv.org/abs/1705.01022>, 2017

for RS simulations and

- Kelling, Jeffrey, Ódor, Géza, Gemming, Sibylle: Bit-Vectorized GPU
  Implementation of a Stochastic Cellular Automaton Model for Surface Growth,
  2016 IEEE International Conference on Intelligent Engineering Systems, 2016.
  INES '16, IEEE, <https://doi.org/10.1109/INES.2016.7555127>, 2016 

for SCA simulations.

Copyright 2011 - 2017 Jeffrey Kelling <j.kelling@hzdr.de>
               Helmholtz-Zentrum Dresden-Rossendorf
               Institute of Ion Beam Physics and Materials Research

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Foobar is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

## Compiling

To compile and run programs from this package you need the following libraries
installed with paths set appropriately:
- CUDA
- HDF5
- libSplash https://github.com/ComputationalRadiationPhysics/libSplash
- zlib
It is recommended to use at least CUDA 7.x and GCC 5.x

The top-level Makefile will build static libraries, located in, inc/ and tools,
located in tools/ . The programs are located in kpz/ and the Makefile will will
the target kpzScaling, which can run on both CPU and GPU. The GPU part of the
resulting program will be optimized for NVIDIA Kepler and later architectures.

## Usage of kpzScaling

The following command runs an RS simulation with DTrDB TC=1,1 for 100kMCS and a
system size of 2^16 x 2^16 . It will compute AC with waiting times 30, 100 and
500MCS . It will run for at most 6.6h and save the system state to the file
"myState.h5" . The simulation will be in the KPZ case (p=1 , q=0).

        kpzScaling 100000 16 --noSave --scheduler dt -l dis --sampleRate .01 \
            --noCountTagsAsMeasurement --corrStart 30 --corrStart 100 --corrStart 500 \
            --saveH5 myState.h5 --maxH 6.6

Please refer to the file kpz/kpzFrontend.cpp to learn about the available
schedulers (defining the type of simulation) and their specific options. Note,
that the parallel CPU implementations were only created for performance
comparisons and not for productive use. Most of them generate parallel random
numbers in an unsave way (using random seeds for each thread instead of
independent sequences).

## References

[KPZ]: http://link.aps.org/doi/10.1103/PhysRevLett.56.889 "Kardar, Mehran,
Parisi, Giorgio, Zhang, Yi-Cheng: Dynamic Scaling of Growing Interfaces, Phys.
Rev. Lett. 56, American Physical Society, 889–892, Mar 1986"

[oct]: https://journals.aps.org/pre/abstract/10.1103/PhysRevE.79.021125 "Ódor,
Géza, Liedke, B., Heinig, K.-H.: Mapping of 2+1 dimensional KPZ growth onto
driven lattice gas model of dimers, Phys. Rev. E 79(021125), 021125, 2009"
