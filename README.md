# TRPL-PVK
Matlab implementation for the determination of carrier (electron or hole) diffusion coefficient (diffusion length) using transient photoluminescence decay (TRPL) experimental data and sum squared error minimization between the theoretical one-dimensional finite element diffusion model and the experimental decay. 

Code and data from paper "[Effect of Mesostructured Layer upon Crystalline Properties and Device Performance on Perovskite Solar Cells](https://doi.org/10.1021/acs.jpclett.5b00483)" *J. Phys. Chem. Lett.* **2015**, 6, 1628-1637.




##Dependencies:
* [Matlab](https://www.mathworks.com/)
* used version R2012b

##Usage

1. Determination of &beta;<sub>s</sub> and &tau;<sub>s</sub> parameters guessing a stretched exponential decay for TRPL data of blocked device (natural recombination without extraction layer). Run `stretched_exponential_fitting.m` routine in the `1_recombination_blocked_device` folder using a TRPL normalized data in the same format as `block_device_data.csv` data file. Grab the output &beta;<sub>s</sub> and &tau;<sub>s</sub> values for the next step.

2. Determination of diffusion coefficient. Open `minimization_proc.m` in folder `2_diffusion_coeficient` and set beta;<sub>s</sub> and &tau;<sub>s</sub> parameters as obtained in step 1. TRPL data obtained from device containing quenched layer (hole or electron extraction layer) is included normalized as `quenched_device_data.csv` file. Run `minimization_proc.m`. See comments in main routine `main_diffusion_funk.m` to study the implementation done here and another options as illumination and quenched side position.

3. Diffusion Length is obtained using &beta;<sub>s</sub> and &tau;<sub>s</sub> parameters from step 1 and D coefficient from step 2. Run `Diffusion_Length.m`

##Citation
If you find this work useful for your research, please cite:
```
@Article{Listorti2015Effect1628,
  Title                    = {Effect of Mesostructured Layer upon Crystalline Properties and Device Performance on Perovskite Solar Cells},
  Author                   = {Listorti, Andrea and Juarez-Perez, Emilio J. and Frontera, Carlos and Roiati, Vittoria and Garcia-Andrade, Laura and Colella, Silvia and Rizzo, Aurora and Ortiz, Pablo and Mora-Sero, Ivan},
  Journal                  = {J. Phys. Chem. Lett.},
  Pages                    = {1628--1637},
  Volume                   = {6},
  Year                     = {2015},
  Doi                      = {10.1021/acs.jpclett.5b00483},
}

```

![Figure 1](https://octodex.github.com/images/yaktocat.png)



## Contact
Feel free to contact me if there is any question (remove white spaces in email address)
```
e j juarez perez (at) gmail.co m 
```
