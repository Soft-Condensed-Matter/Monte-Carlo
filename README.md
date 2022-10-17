
<h1 align="center">
  <br>
  <a href="https://github.com/Soft-Condensed-Matter/Monte-Carlo/blob/master/MC.png"><img src="https://github.com/Soft-Condensed-Matter/Monte-Carlo/blob/master/MC.png" alt="Markdownify" width="200"></a>
  <br>
  Monte Carlo
  <br>
</h1>

<h3 align="center">Monte Carlo homogeneous fluid simulation  </h3>
<h4 align="center">FORTRAN code to simulate Hard-spheres in three dimensions  </h4>

<p align="center">
  <a href="#code">Code</a> •
  <a href="#compile">Compile</a> •
  <a href="#execute">Execute</a> •
    <a href="#visualization">Visualization</a> •
  <a href="#observables">Observables</a>
</p>

## Code
* MC_HS.f90
  - Hard-sphere fluid in three dimensions simulated with Monte Carlo using the Metropolis alghorithm. Within the code the radial distribution function, pressure and compressibility factor averages are computed. The code also makes a file with frames that are used to build a the simulation movie with vmd.

* MC_SW.f90
  - Square-well fluid in three dimensions simulated with Monte Carlo using the Metropolis alghorithm. Within the code the density profile is computed. The code also makes a file with frames that are used to build a the simulation movie with vmd.

## Compile
```bash
Serial:
   # GNU compiler
   $ gfortran -O3 MC_HS.f90 -o mc
   
   # Intel oneAPI
   $ ifort -O3 MC_HS.f90 -o mc
   
   # Nvidia HPC SDK
   $ nvfortran -O3 MC_HS.f90 -o mc

Parallel (share memory):
   # GNU compiler

```   

## Execute
* Simulation parameters
The simulation conditions as the number of particles, the volume fraction, temperature and number of Monte Carlo cycles are specified in the MC.inp file. Only for the Square-Well fluid code the attractive range must be specified at the end of file.
<i>Each time the code is executed files with results are rewritten</i>

* Run simulation
```bash
# Run the code
$./mc
```

## Visualization
* Simulation movie file
By default the simulation file creation is not active, to active uncomment line 401 in file MC_HS.f90. Larger simulations and/or great number of particles will create larger files (>MB) that will be hard to manipulate in personal computers

* Simulation movie
The file is created in <i> *.xyz</i> format that could be visualizated in xmakemol or vmd
```bash
# Load movie file
$ vmd -f MCMovie.xyz
```

* Initial configuration
The code also creates a snapshot of the initial configuration, to see it run
```bash
# Load snapshot file
$ vmd -f MCPic.xyz
```

## Observables
* Thermodynamic properties
Averages of thermodynamic properties are summarized in file MCAvr.dat
	
* Microscopic properties
Radial distribution function is saved in file MCGr.dat and can be plotted with any graphical sortware like gnuplot
```bash
gnuplot> p "MCGr.dat" title '{/Symbol f}=0.4' with line line style 1 line width 2
```

Density profile is saved in file MCGz.dat and can be plotted with any graphical sortware like gnuplot. The second column is the profile direct computation and the third one is the profile center at the center of mass
```bash
gnuplot> p "MCGz.dat" using 1:3 title '{/Symbol T}=1.0' with line line style 1 line width 2
```


## License

GNU LGPL-v3


## Cite this code
[![DOI](https://zenodo.org/badge/542237113.svg)](https://zenodo.org/badge/latestdoi/542237113)

* BibTex
```
@software{MC_HS2022,
  author = {Alexis Torres-Carbajal},
  doi = {10.5281/zenodo.7117628,
  month = {10},
  title = {{Hard-Sphere Monte Carlo}},
  url = {https://github.com/Soft-Condensed-Matter/Monte-Carlo},
  version = {1.0},
  year = {2022}
}
```

---

> Twitter [@Alpixels](https://twitter.com/Alpixels)
