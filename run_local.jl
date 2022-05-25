using DelimitedFiles
using DataFrames
using Parameters
#using Plots
#using Colors
#using PyPlot
#using Plotly
#using GR
#using Cairo,Compose,Fontconfig
using Distributions
using Distributed
addprocs(4)
@everywhere using Random
@everywhere using Parameters
@everywhere using DelimitedFiles
@everywhere include("parameters.jl")
@everywhere include("schisto_abm.jl")
#@everywhere include("main.jl")
@everywhere const n_sim = 1000

include("prevalence.jl")


function file_names(P::SCHparameters)

    age_prevalence_data = "Age_r_$(P.file_index).dat"

    folder0 = ""
    
    if !P.treatment
        folder = string(folder0,"result_$(P.grid_size_snail)_method_$(P.method)/")#"Cluster/fixed_seed/size_500_method_2/"#
    else
        folder = string(folder0,"result_$(P.grid_size_snail)_method_$(P.method)_$(P.rounds)_$(P.Interval)/")
    end
    
    time_data = "inf_time_series_r_$(P.file_index).dat"
    age_data = "age_data_r_$(P.file_index).dat"
    inf_data = "inf_data_r_$(P.file_index).dat"
    group_data = "group_data_r_$(P.file_index).dat"
    worms_data = "worms_data_r_$(P.file_index).dat"
    cercaria_data = "cercaria_data_r_$(P.file_index).dat"
    worms_c_data = "worms_c_data_r_$(P.file_index).dat"
    time_data_found = "inf_time_series_r_found_$(P.file_index).dat"

    return age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data,time_data_found
end

function data_process(results,n_sim,P::SCHparameters)
    
    age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data,time_data_found = file_names(P)

    if !Base.Filesystem.isdir(folder)
        Base.Filesystem.mkpath(folder)
    end

    time_series = zeros(Int64,P.total_sim_time,n_sim)
    cercaria_series = zeros(Int64,P.total_sim_time,n_sim)
    age_pop = zeros(Int64,P.grid_size_human,n_sim)
    group_pop = zeros(Int64,P.grid_size_human,n_sim)
    inf_pop = zeros(Int64,P.grid_size_human,n_sim)
    worms_pop = zeros(Int64,P.grid_size_human,n_sim)
    worms_c_pop = zeros(Int64,P.grid_size_human,n_sim)
    time_series_found = zeros(Int64,P.total_sim_time,n_sim)
    for i = 1:n_sim
        time_series[:,i] = results[i][1]
        age_pop[:,i] = results[i][2]
        group_pop[:,i] = results[i][3]
        inf_pop[:,i] = results[i][4]
        worms_pop[:,i] = results[i][5]
        cercaria_series[:,i] = results[i][6]
        worms_c_pop[:,i] = results[i][7]
        time_series_found[:,i] = results[i][8]
    end

    writedlm(string(folder,time_data),time_series)
    writedlm(string(folder,age_data),age_pop)
    writedlm(string(folder,inf_data),inf_pop)
    writedlm(string(folder,group_data),group_pop)
    writedlm(string(folder,worms_data),worms_pop)
    writedlm(string(folder,cercaria_data),cercaria_series)
    writedlm(string(folder,worms_c_data),worms_c_pop)
    writedlm(string(folder,time_data_found),time_series_found)
end


function run_s(P,n_sim)

    age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data,time_data_found = file_names(P)
    results = pmap(x->main(x,P),1:n_sim)
    data_process(results,n_sim,P)

end


P = SCHparameters(infection_human = 0.00493120,#6.688e-04,
        infection_snail = 0.00017434,
        beta = 0.01065208,
        immunity_parameter = 7.49600000,
        mu_w = 0.29320000,
        grid_size_snail = 500,
        method=2,
        worms_lim_diag = 1,
        file_index = 1
    )

    
age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data,time_data_found = file_names(P)

#P = SCHparameters(beta = 0.01,file_index = 2,grid_size_human = 1000,infection_human = 0.000003,infection_snail = 0.0000001)
n_sim = 1000
run_s(P,n_sim)



P = SCHparameters(infection_human = 2.6542e-03,#6.688e-04,
        infection_snail = 3.0898e-04,
        beta = 3.2000e-02,
        immunity_parameter = 1.1836e+01,
        mu_w = 6.7130e-01,
        grid_size_snail = 1000,
        method=2,
        worms_lim_diag = 1,
        file_index = 1
    )


age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)

#P = SCHparameters(beta = 0.01,file_index = 2,grid_size_human = 1000,infection_human = 0.000003,infection_snail = 0.0000001)
n_sim = 1000
run(P,n_sim)



P = SCHparameters(infection_human = 4.96e-04,#6.688e-04,
        infection_snail = 9.3862e-04,
        beta = 1.2000e-01,
        immunity_parameter = 1.4076e+01,
        mu_w = 8.955e-01,
        grid_size_snail = 2000,
        method=2,
        worms_lim_diag = 1,
        file_index = 1
    )


age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)

#P = SCHparameters(beta = 0.01,file_index = 2,grid_size_human = 1000,infection_human = 0.000003,infection_snail = 0.0000001)
n_sim = 1000
run_s(P,n_sim)



#=---------------------------------------------=#


P = SCHparameters(infection_human = 0.00493120,#6.688e-04,
        infection_snail = 0.00017434,
        beta = 0.01065208,
        immunity_parameter = 7.49600000,
        mu_w = 0.29320000,
        grid_size_snail = 500,
        method=2,
        worms_lim_diag = 1,
        file_index = 1,
        rounds = 10,
        Interval = 1,
        treatment = true
    )

    
age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data,time_data_found = file_names(P)

#P = SCHparameters(beta = 0.01,file_index = 2,grid_size_human = 1000,infection_human = 0.000003,infection_snail = 0.0000001)
n_sim = 1000
run_s(P,n_sim)



P = SCHparameters(infection_human = 2.6542e-03,#6.688e-04,
        infection_snail = 3.0898e-04,
        beta = 3.2000e-02,
        immunity_parameter = 1.1836e+01,
        mu_w = 6.7130e-01,
        grid_size_snail = 1000,
        method=2,
        worms_lim_diag = 1,
        file_index = 1
    )


age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)

#P = SCHparameters(beta = 0.01,file_index = 2,grid_size_human = 1000,infection_human = 0.000003,infection_snail = 0.0000001)
n_sim = 1000
run(P,n_sim)



P = SCHparameters(infection_human = 4.96e-04,#6.688e-04,
        infection_snail = 9.3862e-04,
        beta = 1.2000e-01,
        immunity_parameter = 1.4076e+01,
        mu_w = 8.955e-01,
        grid_size_snail = 2000,
        method=2,
        worms_lim_diag = 1,
        file_index = 1
    )


age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)

#P = SCHparameters(beta = 0.01,file_index = 2,grid_size_human = 1000,infection_human = 0.000003,infection_snail = 0.0000001)
n_sim = 1000
run_s(P,n_sim)












real_prevalence,estimated_prevalence,n_sim_zero = calc_prevalence_sim(P)

Plots.plot([real_prevalence estimated_prevalence],label = ["real" "estimated"])


cercaria = readdlm(string(folder,cercaria_data))
cerc = map(x->mean(cercaria[x,:]),1:size(cercaria)[1])
Plots.plot(cerc)


#=
par = readdlm("Cluster/GA/Old_datas/GA_result_pop_19_2_1000.dat",header=false)
f = readdlm("Cluster/GA/Old_datas/GA_result_fitness_19_2_1000.dat",header=false)[:,1]
pos=sortperm(f,rev=true)
sort(f,rev=true)

NN = 4

for i = 2:2
    println(i)
    v = par[:,pos[i]]
    P = SCHparameters(infection_human = v[1],
            infection_snail = v[2],
            beta = v[3],
            immunity_parameter = v[4],
            mu_w = v[5],
            grid_size_snail = 1000,
            method=2,
            worms_lim_diag = 1,
            file_index = i
        )

    age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)
    #P = SCHparameters(beta = 0.01,file_index = 2,grid_size_human = 1000,infection_human = 0.000003,infection_snail = 0.0000001)
    n_sim = 40
    run(P,n_sim)
end
=#

#@everywhere someething
results = pmap(x-> main(x,data),1:1000)