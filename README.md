# Disk Anisotropy


## Authors
Chong-Chong He


## Description
A program to compute anisotropy factors of radiations from accreting 
neutron stars, considering idealized geometry and GR effects.

## How to run

Complie the code with 

```
gfortran -O2 -o main vectors.f95 main.f95
```

There is also a makefile to do this. Then, run the executable (main). No parameter is required. The parameters are set by the namelist file "input.nml". The code is tested with gfortran 7.2.0 on macOS 10.13.2. The program will read parameters from namelist "input.nml", do the calculation, and print out the results on the screen. You may want to use `./main > output.out` to store the outputs into file.


### Namelist

- GEOMETRY

  | Parameters | Type    | Default  | Description                                                  |
  | ---------- | ------- | -------- | ------------------------------------------------------------ |
  | shape      | integer | 1        | 1, 2, 3 corresponds to disk a, b, c in Fig. 1 of He & Keek 2016 |
  | para1      | real    | 0        | If shape = 2 or 3, this is the slope of the disk surface     |
  | diskRadius | real    | 4000     | Disk radius in units of NS radius                            |
  | r0         | real    | 1        | Inner radius of the disk in units of NS radius. Must in the range [1, diskRadius) |
  | h0         | real    | 0        | The height of the inner edge of the disk. Has no effect for now. |
  | los        | string  | 'degree' | List of observational angles as 'degree' or 'cos(angle)'     |


- PHYSICS

  | Parameters | Type    | Default | Description                                                  |
  | ---------- | ------- | ------- | ------------------------------------------------------------ |
  | lambertian | logical | .true.  | false for H-function radiation                               |
  | GR         | logical | .false. | whether or not to take into account GR effects, i.e. light bending and gravitational redshift. Current version only supports shape = 1. Note that you have to set rNS to a reasonable number before you can see GR effects. |
  | rNS        | real    | huge    | NS radius in units of gravitational radius r\_g when GR is enabled. For a 1.4 solar mass NS with radius = 12 km, rNS ~ 3 rg. The default is set to infinity so that it converges to non-GR case, which should give the same results as when GR is turned off. |

## TODO

1. Make GR disk anisotropy work for inclined disks.
2. Add persistant anisotropy
3. Add GR to the calculation of radiation onto the disk.
