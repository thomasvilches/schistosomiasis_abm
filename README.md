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

### Running the genetic algorithm

Both the files `AG.jl` and `AG_Cluster.jl` have the same function. The difference is that the first one can be used to run the simulations in a local computer, while the second one will distribute the simulations in a cluster. Please refeer to [ClusterManagers Julia Package](https://github.com/JuliaParallel/ClusterManagers.jl).

> :warning: **We do not recommend these simulations to be done in a local machine. A high capacity computational source is strongly recommended**


&rarr; **Firstly**, make sure that the path inside functions `running_generations` and `run_ga` are pointing to an existing path in your system. The reader can use the command line

```julia
include("AG.jl")
```

or


```julia
include("AG_Cluster.jl")
```

to compile the code. The function `run_ga` receives three arguments:

1. The diagnostic method that is taken into account when reading the data from `Age.dat`, being 
    1. Kato-Katz (KK).
    1. Helmintex (HTX).
2. The number of pairs of worms that the diagnostic is supposed to detect the infection.
3. The snail population size.

For instance, one can run the GA for the HTX, which is able to detect the infection with one pair of worms, and a snail-population size of 2000 snails using.


```julia
run_ga(2,1,2000)
```

### Output of the genetic algorithm

For every generation of the GA[^1], two files called `GA_result_fitness_ii_xx_zz.dat` and `GA_result_pop_ii_xx_zz.dat` are saved. The name *pop* stands for the population of *chromossomes*, corresponding to the parameters of interest, in that generation, while the name *fitness* stands for the score that the given chromossome achieved. The indeces *ii* represents the generation, *xx* represents the method that was chosen in the simulation, and *size* represents the snail-population size.

After finding the best fit, one can use the code in the branch [singlesims](https://github.com/thomasvilches/schistosomiasis_abm/tree/singlesims) to run the simulations for the desired parameters.


[^1]: The number of generations can be changed in the header of `AG*.jl`.