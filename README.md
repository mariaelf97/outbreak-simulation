



## Covasim 

### Installation

`pip install -r requirements.txt`

### Simulation
Code : `src/run_sim.py`
The goal here is to simulate an Covid-19 outbreak in a defined region and calculate the viral shedding rates into wastewater. Here we describe a case where the region is a rectangle with 4 sub-regions. 2 rows and 2 columns. Let's break this down into pieaces.

1. We start with defining parameters for the simulation. These include `pop_size`,`pop_type`,`location`, `n_days`, `pop_infected` and `n_imports`. There are other parameters available for simulation through covasim that can be directly added to `define_sim_parameters()` function as well. For more explanation on what each of these parameters mean, please visit https://docs.idmod.org/projects/covasim/en/latest/tutorials/tut_intro.html

The following are the list of all available parameters for the simulation:

```
  --pop_size POP_SIZE   Population size
  --pop_type POP_TYPE   Population type (e.g. random, hybrid, etc.)
  --n_days N_DAYS       Number of simulation days
  --location LOCATION   Location name
  --pop_infected POP_INFECTED
                        Initial infected population (default: 0)
  --n_imports N_IMPORTS
                        Number of imported infections (default: 0)
  --n_rows N_ROWS       Number of rows in the region (default: 2)
  --n_cols N_COLS       Number of columns in the region (default: 2)
  --n_init_inf N_INIT_INF
                        Number of initial infections in the region (default: 20)
  --r_init_inf R_INIT_INF
                        Row index for where the infection initiated (default: 0)
  --c_init_inf C_INIT_INF
                        column index for where the infection initiated (default: 0)
```

Please note that covasim by default include each country's age distribution and household size. That is the reason for including the name of the country as a parameter. Function `check_age_household_dist()` checks whether the user input country has available data so if you get an error from this function, it indicates that you either have a typo in the name or the data is not available for that geolocation.

2. Once the parameters are set, the code will try to assign each individual to a random region.
This is done by function `assign_people()`, we use numpy's random shuffle function to ensure the assignment is random. Depending on the number of rows and columns you pass for your region, you will have differing number of regions. This function returns an array of length N where N is number of people and the values will be 0-X where X is Total number of defined regions.
3. We then initiate number of N number of infections in one of the regions using `--n_init_inf N_INIT_INF`,`--r_init_inf R_INIT_INF` and `--r_init_inf R_INIT_INF`. You define the region by its row and column indices. Please note that we started with 0 `--pop_infected` infected individual and 0 `n_imports`(average daily number of imported cases).
4. We then run the simulation, and calculate number of people infected in each region in each time point based on the parameters we passed earlier. Briefly, the function `calculate_new_infections()` first start with finding when each individual becomes infected and then use that information plus where each people is assigned to to create an array of total number of infected people in each region and each time point.
4. The code will calculate the viral shedding for each region and time point using two different methods, simple and using the viral loads calculated internally in covasim.
### viral shedding using "Eclipse phase" model
This model is a naive model that assumes all individuals, regardless of their age and other model parameters. The parameters for the viral shedding function is assigned as previously described(here).



