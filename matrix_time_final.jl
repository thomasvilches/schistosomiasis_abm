
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

folder0 = "./"

println("$intv $rnd")
P = SCHparameters(file_index = file,method = method, grid_size_snail = snail_pop, rounds = rnd,Interval = intv,treatment = false)

    
if !P.treatment
    folder = string(folder0,"result_$(P.grid_size_snail)_method_$(P.method)/")#"Cluster/fixed_seed/size_500_method_2/"#
else
    folder = string(folder0,"result_$(P.grid_size_snail)_method_$(P.method)_$(P.rounds)_$(P.Interval)/")
end
data_time = readdlm(string(folder,"inf_time_series_r_",file,".dat"),Int64,header = false)
#data_group = readdlm(string(folder,"/group_data_r_",file,".dat"),Int64,header = false)

m = Array{Float64,2}(undef,size(data_time,1),n_boots)
l = 1:size(data_time,2)

m = pmap(1:n_boots) do i
    #println(i)
    mean(data_time[:,rand(l,1000)],dims=2)
end

m = hcat(m...)
writedlm(string(folder,"/matrix_time_$(file).dat"),m)

