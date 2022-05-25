####This file builds an Genetic Algorithm in order to estimate the parameters of a function
# based on the Least Square Method
using Distributed
addprocs(4)

using Distributions
using LatinHypercubeSampling
@everywhere include("main.jl")
@everywhere const n_sim = 8

include("prevalence.jl")

const n_parameters = 5
const n_crom = 10
const n_generations = 10


function generating_parameters(n_parameters::Int64,n_crom::Int64)
    #This function must generate the parameters distribuition
    
    population = zeros(Float64,n_parameters,n_crom)
    A = LHCoptim(n_crom,n_parameters,5)
    B = A[1]/n_crom
    ###In this case, we have 5 parameters (infection_human,infection_snail,beta,L,mu_w)
    population[1,:] = map(x-> (0.001-0.0001)*x+0.0001,B[:,1])#rand(Uniform(0,10),n_crom) 
    population[2,:] = map(x-> (0.0001-0.00001)*x+0.00001,B[:,2])#rand(Uniform(0,500),n_crom)
    population[3,:] = map(x-> x,B[:,3])#rand(Uniform(0,1),n_crom)
    population[4,:] = map(x-> (15-1)*x+1,B[:,4])#rand(Uniform(0,15),n_crom)
    population[5,:] = map(x-> (1-0.05)*x+0.05,B[:,5])#rand(Uniform(0.1,1),n_crom)

    return population
end

function least_square(estimated_prevalence::Array{Float64,1},prevalence_::Array{Float64,1})

    F = estimated_prevalence-prevalence_
    F = sum(map(x->x^2,F))
    F = 1/sqrt(F)
    return F

end

function fitness_crom(cromoss::Array{Float64,1},prevalence_diag::Array{Float64,1})
    P = SCHparameters(infection_human = cromoss[1],
        infection_snail = cromoss[2],
        beta = cromoss[3],
        immunity_parameter = cromoss[4],
        mu_w = cromoss[5]
    )

    run(P,n_sim)

    real_prevalence,estimated_prevalence,n_sim_zero = calc_prevalence_sim(P)
    
    f::Float64 = least_square(estimated_prevalence,prevalence_diag)
    if sum(cromoss) < 0.0005
        f = f/10.0
    end

    return estimated_prevalence,f

end

function finding_fitness_pop(population::Array{Float64,2},P::SCHparameters)

    n_parameters,n_crom = size(population)
    F = zeros(Float64,n_crom)
    groups,prevalence_diag = assimetric_prevalence() #first column = Kato-Katz, second = HTX

    for i = 1:n_crom
        println("Finding fitness pop = $i")
        estimated_prevalence,F[i] =  fitness_crom(population[:,i],prevalence_diag[:,P.method])
    end

    return F,prevalence_diag[:,P.method]

end

function sorting_crom_roulette(F::Array{Float64,1})

    F1 = zeros(Float64,length(F))
    for i = 1:length(F)
        F1[i] = F[i]
    end

    F2 = cumsum(F1)/sum(F1)
    r = rand()
    first_one = findfirst(x-> x > r,F2)

    return first_one

end

function crossover(crom1::Array{Float64,1},crom2::Array{Float64,1},prevalence_diag::Array{Float64,1},P::SCHparameters)

    r = rand(1:(length(crom1)-1))

    crom3 = [crom1[1:r];crom2[r+1:length(crom1)]]
    crom4 = [crom2[1:r];crom1[r+1:length(crom1)]]

    estimated_prevalence,f3 =  fitness_crom(crom3,prevalence_diag)
    estimated_prevalence,f4 =  fitness_crom(crom4,prevalence_diag)

    aux_vec = [f3;f4]
    aux_matrix = [crom3 crom4]

    p = sortperm(aux_vec,rev = true)

    return aux_matrix[:,p[1]],aux_vec[p[1]]
end

function mutation(crom1::Array{Float64,1},F::Float64,P::SCHparameters,prevalence_diag::Array{Float64,1})
    
    mut::Int64 = 0
    if rand() < P.mutation_prob
        crom1[1] = rand(Uniform(0,10))
        mut = 1
    end
    if rand() < P.mutation_prob
        mut = 1
        crom1[2] = rand(Uniform(0,500))
    end
    if rand() < P.mutation_prob
        mut = 1
        crom1[3] = rand(Uniform(0,1))
    end
    if rand() < P.mutation_prob
        mut = 1
        crom1[4] = rand(Uniform(0,15))
    end
    if rand() < P.mutation_prob
        mut = 1
        crom1[5] = rand(Uniform(0.1,1))
    end

    if mut == 1
        println("Mutation--")
        estimated_prevalence,F =  fitness_crom(crom1,prevalence_diag)
    end

    return crom1,F

end

function running_generations(population::Array{Float64,2},F::Array{Float64,1},P::SCHparameters,n_generations::Int64,prevalence_diag::Array{Float64,1})
    n_parameters,n_crom = size(population)
    population_aux = zeros(Float64,size(population))
    F_aux = zeros(Float64,length(F))
    c1::Int64 = -1
    c2::Int64 = -1
    for i = 1:n_generations
        println("generation = $i")
       #= for kk = 1:n_crom
            F_aux[kk] = F[kk]
            for k = 1:n_parameters
                population_aux[k,kk] = population[k,kk]
            end
        end=#
        for j = 1:Int((n_crom/2))
            c1 = sorting_crom_roulette(F)
            population_aux[:,j] = population[:,c1]
            F_aux[j] = F[c1]
            F[c1] = 0.0
        end

        for j = 1:Int((n_crom/2))
            println("pop_fill = $i $j")
            c1 = -1
            c2 = -1
            while c1 == c2
                c1 = Int(rand(1:(n_crom/2)))
                c2 = Int(rand(1:(n_crom/2)))
            end

            population_aux[:,Int((n_crom/2)+j)],F_aux[Int((n_crom/2)+j)] = crossover(population_aux[:,c1],population_aux[:,c2],prevalence_diag,P)

        end

        for k = 1:n_crom
            population_aux[:,k],F_aux[k] = mutation(population_aux[:,k],F_aux[k],P,prevalence_diag)
        end

        for kk = 1:n_crom
            F[kk] = F_aux[kk]
            for k = 1:n_parameters
                population[k,kk] = population_aux[k,kk]
            end
        end
        writedlm("/data/thomas/GA_results/GA_result_pop_$(i)_$(P.method).dat",population)
        writedlm("/data/thomas/GA_results/GA_result_fitness_$(i)_$(P.method).dat",F)
    end

    return population,F
end


function run_ga(met::Int64,lim::Int64)
    P = SCHparameters(method=$(met),worms_lim_diag = $(lim))

    population = generating_parameters(n_parameters,n_crom)

    Fit_Pop,prevalence_diag = finding_fitness_pop(population,P)

    population,Fit_Pop = running_generations(population,Fit_Pop,P,n_generations,prevalence_diag)

   # writedlm("GA_result_pop.dat",population)
    #writedlm("GA_result_fitness.dat",Fit_Pop)
end