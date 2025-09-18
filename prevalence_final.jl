
using Distributed
addprocs(10)
@everywhere using DelimitedFiles
@everywhere using Parameters
@everywhere using Statistics
@everywhere include("parameters.jl")

method = 2
snail_pop = 2000
file = 0
n_pop_ga = 500
n_gen = 20
n_boots = 150
limite = 1

rnd = 1
intv = 0.0


P = SCHparameters(file_index = file,method = method, grid_size_snail = snail_pop, rounds = rnd,Interval = intv,treatment = true)

folder = string("result_",snail_pop,"_method_",method)

data_worm = readdlm(string(folder,"/worms_c_data_r_",file,".dat"),Int64,header = false)
data_group = readdlm(string(folder,"/group_data_r_",file,".dat"),Int64,header = false)

mg = maximum(data_group)
m = zeros(Float64,n_boots,mg)
l = 1:size(data_group,2)
n_inf = Array{Float64,1}(undef,mg)
n_ind = Array{Float64,1}(undef,mg)

@everywhere function prevalence_matrix(data_worm,data_group,l,mg,i,P)
    println(i)
    n_inf = [0 for i = 1:mg]#length(n_inf)]
    n_ind = [0 for i = 1:mg]#length(n_ind)]
    pos = rand(l,1000)
    m = Array{Float64,1}(undef,mg)
    for w = 1:length(pos)
        j = pos[w]
        for k = 1:size(data_worm,1)
            #println("$i,$w,$k")
            n_ind[data_group[k,j]] +=  1
            if data_worm[k,j] >=  P.worms_lim_diag
                n_inf[data_group[k,j]] += 1
            end
        end
    end
    m = n_inf./n_ind
    return m
end

m = pmap(1:n_boots) do i
    prevalence_matrix(data_worm,data_group,l,mg,i,P)
end
m = hcat(m...)

writedlm(string(folder,"/matrix_$(file).dat"),m)
