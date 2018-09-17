# TRPL-PVK
Matlab implementation for the determination of carrier (electron or hole) diffusion coefficient (diffusion length) using transient photoluminescence decay or time resolved photoluminescence (TRPL) experimental data.  Diffusion coefficient is obtained by the sum squared error minimization procedure between the theoretical curve obtained from the one-dimensional finite element diffusion model and the experimental decay. 


![Figure 1](https://github.com/ej-jp/TRPL-PVK/blob/master/img/Graph1.png)


Code and data from paper "[Effect of Mesostructured Layer upon Crystalline Properties and Device Performance on Perovskite Solar Cells](https://doi.org/10.1021/acs.jpclett.5b00483)" *J. Phys. Chem. Lett.* **2015**, 6, 1628-1637.
See [Support information file](https://pubs.acs.org/doi/suppl/10.1021/acs.jpclett.5b00483/suppl_file/jz5b00483_si_001.pdf) of the paper for details in this implementation. 



## Dependencies:
* [Matlab](https://www.mathworks.com/) (used version R2012b)

## Usage

The procedure consist of three steps:

1. Determination of &beta;<sub>s</sub> and &tau;<sub>s</sub> parameters guessing a stretched exponential decay for TRPL data of blocked device (natural recombination without extraction layer). Run `stretched_exponential_fitting.m` routine in the `1_recombination_blocked_device` folder using a TRPL normalized experimental data in the same format as `block_device_data.csv` data file. Grab the output &beta;<sub>s</sub> and &tau;<sub>s</sub> values for the next step.

2. Determination of diffusion coefficient. Open `minimization_proc.m` in folder `2_diffusion_coeficient` and set beta;<sub>s</sub> and &tau;<sub>s</sub> parameters as obtained in step 1. TRPL experimental data obtained from device containing quenched layer (hole or electron extraction layer) is included normalized in the `quenched_device_data.csv` file. Run `minimization_proc.m`. See comments in main routine `main_diffusion_funk.m` to study the implementation done here and another options as illumination and quenched side position.

3. Diffusion Length is obtained using &beta;<sub>s</sub> and &tau;<sub>s</sub> parameters from step 1 and D coefficient from step 2. Run `Diffusion_Length.m` to obtain the diffusion length (nm)


## Citation

If you find this work useful for your research, please cite:

A. Listorti, E. J. Juarez-Perez, C. Frontera, V. Roiati, L. Garcia-Andrade, S. Colella, A. Rizzo, P. Ortiz, I. Mora-Sero, [Effect of Mesostructured Layer upon Crystalline Properties and Device Performance on Perovskite Solar Cells](https://doi.org/10.1021/acs.jpclett.5b00483)" *J. Phys. Chem. Lett.* **2015**, 6, 1628-1637.

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





## Contact
Feel free to contact me here in github if there is any question or by email (remove white spaces in email address)
```
e j juarez perez (at) gmail.co m 
```
