# Disk Anisotropy


## Authors
Chong-Chong He (che1234 @ umd.edu)


## Description

_DiskAnisotropy_ is a program to compute anisotropy factors of radiations from accreting 
neutron stars, taking into account geometries effects, realistic surface radiative transfer models, and GR effects.

The underlining physics is explained in [C.-C. He & L. Keek, 2016, ApJ, 819, 47](https://ui.adsabs.harvard.edu/abs/2016ApJ...819...47H/abstract) with the following additional features: 
(1) allowing a gap between the inner edge of the disk and the NS surface, (2) allowing GR effects

## How to run

Complie the code with a simple `make`, or by running

```
gfortran -O2 -o main vectors.f95 main.f95
```

Then, run the executable (main). The input parameters are applied by changing the namelist file "input.nml". The program will read parameters from "input.nml", do the calculation, and print out the results on the screen. You may want to use `./main > output.out` to store the outputs into a file. The code is tested with gfortran 7.2.0 on macOS 10.13.2. 


### Namelist

- GEOMETRY

  | Parameters | Type    | Default  | Description                                                                          |
  |------------|---------|----------|--------------------------------------------------------------------------------------|
  | shape      | integer | 1        | 1, 2, 3 corresponds to disk a, b, c in Fig. 1 of He & Keek 2016                      |
  | para1      | real    | 0        | If shape = 2 or 3, this is the slope of the disk surface                             |
  | diskRadius | real    | 4000     | Disk radius in units of NS radius                                                    |
  | r0         | real    | 1        | Inner radius of the disk in units of NS radius. Must be in the range [1, diskRadius) |
  | h0         | real    | 0        | The height of the inner edge of the disk.                                            |
  | los        | string  | 'degree' | List of observational angles as 'degree' or 'cos' (the cos of the angle)             |


- PHYSICS

  | Parameters | Type    | Default | Description                                                                                                                                                                                                                 |
  |------------|---------|---------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
  | lambertian | logical | .true.  | true: Lambertian surface (cos function); false: H-function surface                                                                                                                                                          |
  | GR         | logical | .false. | whether or not to take into account GR effects, i.e. light bending and gravitational redshift. Current version only supports shape = 1. Note that you have to set rNS to a reasonable number before you can see GR effects. |
  | rNS        | real    | 3.0     | NS radius in units of gravitational radius r\_g when GR is enabled. For a 1.4 solar mass NS with radius = 12 km, rNS ~ 3 rg.                                                                                                |

## TODO

1. Add GR effects to the calculation of radiation onto the disk.
2. Add GR effects to inclined disks.
3. Make h0 work in the function radiation

