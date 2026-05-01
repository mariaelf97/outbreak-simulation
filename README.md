



## Covasim 

### Installation

`pip install covasim`

### Simulation
Code : `src/run_sim.py`
The goal here is to simulate an Covid-19 outbreak in a defined region and calculate the viral shedding rates into wastewater. Here we describe a case where the region is a rectangle with 4 sub-regions. 2 rows and 2 columns. Let's break this down into pieaces.

1. We start with defining parameters for the simulation. These include `pop_size`,`pop_type`,`location`, `n_days`, `pop_infected` and `n_imports`. There are other parameters available for simulation through covasim that can be directly added to `define_sim_parameters()` function as well. For more explanation on what each of these parameters mean, please visit https://docs.idmod.org/projects/covasim/en/latest/tutorials/tut_intro.html




