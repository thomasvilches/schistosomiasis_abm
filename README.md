# Agent-based model for schistosomiasis transmission and treatment

The model is used in the manuscript called 

> **Add a name** - *Vilches et al.*

### Brief description

The model is implemented in [Julia Language](https://julialang.org/project/), and a genetic algorithm (GA) is used to find the unknown parameters

* infection probabilities, for humans and snails;
* reduction on infection probability due to acquired immunity;
* period until the development of immunity;
* the worm mortality;

of the model based on the data of age-stratified reported prevalence of schistosomiasis infections in the district of Candeal, municipality of Estância-SE, Brazil. Please, refeer to [*Lidholtz et al*](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0006274) for more information about the data.

### Running a single scenario

Both the files `run_local.jl` and `run_cluster.jl` have the same function. The difference is that the first one can be used to run the simulations in a local computer, while the second one will distribute the simulations in a cluster. Please refeer to [ClusterManagers Julia Package](https://github.com/JuliaParallel/ClusterManagers.jl).

> :warning: **We do not recommend these simulations to be done in a local machine. A high capacity computational source is strongly recommended**


&rarr; **Firstly**, make sure that the path `folder0`, inside functions `file_names`, is pointing to an existing path in your system. You can use the command line

```julia
include("run_local.jl")
```

or


```julia
include("run_cluster.jl")
```

to compile the code. 

Change the parameters inside functions `run1`, `run2` or `run3` to run the simulations for 500, 1000 or 2000 snail-population size, respectively. These functions they will receive eight argumentsd in the following order:

<!-- idx,wld=1,treat=false,interv=[0],roundss=[0],strat=:dg,eff=0.0,et=0 -->
- **idx** &rarr; *Int64*: the index that will be used in the name of folder that will store the results.
- **wld** &rarr; *Int64*: the minimum number of worm pairs that the diagnostic method can find.
- **treat** &rarr; *Boolean*: if the treatment will be applied or not.
- **interv** &rarr; *Vector{Float64}*: the time between treatment, in years. If more than one element, the function runs the simulation for each element of this vector.
- **rounds** &rarr; *Vector{Int64}*: number of treatment rounds. If more than one element, the function runs the simulation for each element of this vector.
- **strat** &rarr; *Symbol*: if `:dg` the treatment is applied only to diagnosed population. if `:total` the treatment is applied to the entire population.
- **eff** &rarr; *Float64*: the efficacy of the treatment.

For instance, if one wants to run the simulation for 1000 snails, with the tratment of the diagnosed population, applied every 0.5 year, during 4 rounds. Considering that the treatment is able to diagnose one pair of worms in the human body and the drug efficacy is 90%. One should run:

```julia
run2(1,1,true,[0.5],[4],:dg,0.9)
```

### Output of the simulations

After running the functions `run*`, the ABM will create a folder, in the chosen path, named:

> result\_*nn*\_method\_2\_*rr*\_*ii*

if the argument *treat* is true, in which *nn* is snail population size, *rr* is the number of treatment rounds, and *ii* is the time between two rounds.

If the argument *treatment* is false, then the folder has the name

> result\_*nn*\_method_2

Inside those folders, there will be several files:

- **age\_data\_r\__idx_.dat**: Matrix containing the age of each individual at the end of the simulation. Each Monte-Carlo simulation corresponds to a column.
- **group_data_r\__idx_.dat**: Matrix containing the age group of each individual at the end of the simulation. Each Monte-Carlo simulation corresponds to a column.
- **inf_data_r\__idx_.dat**: Matrix containing 0 (if the individual is not infected at the end of the simulation) and 1 (if the individual is infected at the end of the simulation). Each Monte-Carlo simulation corresponds to a column.
- **cercaria_data_r\__idx_.dat**: Temporal evolution of the cercaria reservoir.
- **inf_time_series_r\__idx_.dat**: Prevalence of infections. 
- **inf_time_series_r_found\__idx_.dat**: Prevalence of infections that the chosen diagnosis could find. Each Monte-Carlo simulation corresponds to a column.
- **worms_c_data_r\__idx_.dat**: Matrix containing the number of worm couples per individual per simulation. Each Monte-Carlo simulation corresponds to a column.
- **worms_data_r\__idx_.dat**: Matrix containing the number of worms per individual per simulation. Each Monte-Carlo simulation corresponds to a column.

The parameter **_idx_** stands for the argument passed to the function `run*` and can be used to easily differ between scenarios, such as different drug efficacies.
