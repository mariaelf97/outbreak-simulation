import covasim as cv
import covasim.data as cvdata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import epyestim.covid19 as covid19
import networkx as nx
from networkx.drawing.nx_pydot import graphviz_layout
import argparse

def define_sim_parameters(pop_size, pop_type,
                          n_days, location,
                          pop_infected, n_imports):
    pars = dict(
    pop_size = pop_size,
    pop_type = pop_type,
    n_days = n_days,
    location = location,
    pop_infected = pop_infected,   # no automatic initial infections
    n_imports = n_imports,      # no imported infections
    )
    return pars


def check_age_household_dist(location):
    try:
        age_data = cvdata.get_age_distribution(location)
        print(f"Age distribution for {location}:")
        print(age_data)
    except ValueError as e:
        print(f"No age data for {location}, using default. ({e})")

    try: 
        household_size = cvdata.get_household_size(location) 
        print(f"Household size distribution for {location}:") 
        print(household_size) 
    except ValueError as e: 
        print(f"No household data for {location}, using default. ({e})")


def assign_people(sim, n_regions, pop_size):
    # assign each person to a region
    region_ids = np.repeat(np.arange(n_regions), pop_size // n_regions)
    np.random.shuffle(region_ids)
    sim.people.region = region_ids
    return sim


def initiate_infection(sim, row_index, col_index, n_infected):
    sim.people.infect(inds=np.where(sim.people.region == row_index)[col_index][:n_infected])
    return sim


def calculate_new_infections(sim, n_regions):
    new_cases = np.zeros((len(sim.results['new_infections']), n_regions))
    for t in range(len(sim.results['new_infections'])):
        newly_infected = np.where(sim.people.date_infectious == t)[0]
        for r in range(n_regions):
            new_cases[t, r] = np.sum(sim.people.region[newly_infected] == r)
    return new_cases


def eclipse_model(y, t, b, k, delta, p, mu, c):
    T, I1, I2, Vi, Vni = y
    dydt = [-b*Vi*T, b*Vi*T - k*I1, k*I1 - delta*I2, p*mu*I2 - c*Vi, p*(1.-mu)*I2 - c*Vni]
    return dydt


def solve_eclipse():
    b = 5e-5
    c_clear = 10
    k = 6
    mu = 1e-4
    p = 1e5
    delta = 0.5
    t = np.linspace(0, 40, 41)
    y0 = [1.33e5, 0, 1/30, 0, 0]
    # solve ode
    sol = odeint(eclipse_model, y0, t, args=(b, k, delta, p, mu, c_clear))
    c = sol[:,3:].sum(axis=1)
    c = c/c.sum()
    return c


def viral_shedding_simple(new_cases):
    ans = solve_eclipse()
    # create an array of zeros the same simze of new cases
    ww_signal = np.zeros_like(new_cases)
    # this calculates viral shedding in each day per region
    for i in range(new_cases.shape[1]):
        ww_signal[:, i] = np.convolve(new_cases[:, i], ans[::-1], mode='same')
    return ww_signal


def calculate_viral_load_matrix(sim):
    viral_load_matrix = np.zeros((sim['n_days'], sim.n))
    for t in range(sim['n_days']):
        vl = cv.utils.compute_viral_load(
            t,
            sim.people.date_infectious,
            sim.people.date_recovered,
            sim.people.date_dead,
            sim['viral_dist']['frac_time'],
            sim['viral_dist']['load_ratio'],
            sim['viral_dist']['high_cap'],
        )
        viral_load_matrix[t, :] = vl  # store full vector
    return viral_load_matrix

        
def viral_shedding_covasim(sim):
    viral_load_matrix = calculate_viral_load_matrix(sim)
    n_days, n_people = viral_load_matrix.shape
    n_regions = sim.people.region.max() + 1
    regional_viral_load = np.zeros((n_days, n_regions))
    for t in range(n_days):
        for p in range(n_people):
            r = sim.people.region[p]  # Get the region for this specific person
            regional_viral_load[t, r] += viral_load_matrix[t, p]
    return regional_viral_load


def plot_shedding(data,name):
    fig, ax = plt.subplots(figsize=(6,10))
    im = ax.imshow(data, aspect='auto', cmap='viridis')

    # Discrete x-axis labels
    x_labels = ['1', '2', '3', '4']   # replace with your actual category names
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels)

    ax.set_xlabel("regions")
    ax.set_ylabel("days")

    fig.colorbar(im, ax=ax, label="viral load")
    plt.savefig(f'{name}.pdf')


def main():
    parser = argparse.ArgumentParser(
        description="Define simulation parameters"
    )

    parser.add_argument("--pop_size", type=int, required=True,
                        help="Population size")

    parser.add_argument("--pop_type", type=str, required=True,
                        help="Population type (e.g. random, hybrid, etc.)")

    parser.add_argument("--n_days", type=int, required=True,
                        help="Number of simulation days")

    parser.add_argument("--location", type=str, required=True,
                        help="Location name")

    parser.add_argument("--pop_infected", type=int, default=0,
                        help="Initial infected population (default: 0)")
    parser.add_argument("--n_imports", type=int, default=0,
                        help="Number of imported infections (default: 0)")
    parser.add_argument("--n_rows", type=int, default=2,
                        help="Number of rows in the region (default: 2)")
    parser.add_argument("--n_cols", type=int, default=2,
                        help="Number of columns in the region (default: 2)")
    parser.add_argument("--n_init_inf", type=int, default=20,
                        help="Number of initial infections in the region (default: 20)")
    parser.add_argument("--r_init_inf", type=int, default=0,
                        help="Row index for where the infection initiated (default: 0)")
    parser.add_argument("--c_init_inf", type=int, default=0,
                        help="column index for where the infection initiated (default: 0)")

    args = parser.parse_args()
    ### STEP 1 ####
    pars = define_sim_parameters(
        pop_size=args.pop_size,
        pop_type=args.pop_type,
        n_days=args.n_days,
        location=args.location,
        pop_infected=args.pop_infected,
        n_imports=args.n_imports,
    )
    check_age_household_dist(args.location)
    # initialize the simulation with the parameters
    sim = cv.Sim(pars)
    sim.initialize()
    # Define regions for the simulation
    n_rows, n_cols = args.n_rows, args.n_cols
    n_regions = n_rows * n_cols
    # Assign people to each region
    sim = assign_people(sim, n_regions, args.n_days)
    # initiate infection
    # 20 initial infections in Region_0_0
    sim = initiate_infection(sim, args.r_init_inf, args.c_init_inf, args.n_init_inf)
    sim.run()
    #### STEP 2 ####
    # calculate number of new cases in each region
    new_cases = calculate_new_infections(sim, n_regions)
    # calculate wastewater shedding per region per time point using
    # basic shedding model
    shedding_simple = viral_shedding_simple(new_cases)
    # covasim viral load model
    shedding_covasim = viral_shedding_covasim(sim)
    # plot values
    plot_shedding(shedding_simple,"simple")
    plot_shedding(shedding_covasim,"covasim")


if __name__ == "__main__":
    main()