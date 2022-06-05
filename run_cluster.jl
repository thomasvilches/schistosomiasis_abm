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
using ClusterManagers
#addprocs(4)
addprocs(SlurmManager(250), N=8, topology=:master_worker, exeflags="--project=.")
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

    folder0 = "/data/thomas/schisto_treat/"
    #folder0 = "."
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


#=---------------------------------------------=#


function run1(idx,wld=1,treat=false,interv=[0],roundss=[0],strat=:dg,eff=0.0,et=0)
    
    for inter = interv,roun = roundss
        P = SCHparameters(infection_human = 0.00493120,#6.688e-04,
            infection_snail = 0.00017434,
            beta = 0.01065208,
            immunity_parameter = 7.49600000,
            mu_w = 0.29320000,
            grid_size_snail = 500,
            method=2,
            worms_lim_diag = wld,
            file_index = idx,
            med_efficacy = eff,
            rounds = roun,
            treat_strat = strat,
            Interval = inter,
            treatment = treat,
            eff_type = et
        )


        age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data,time_data_found = file_names(P)

        #P = SCHparameters(beta = 0.01,file_index = 2,grid_size_human = 1000,infection_human = 0.000003,infection_snail = 0.0000001)
        n_sim = 1000
        run_s(P,n_sim)
    end
end




function run2(idx,wld=1,treat=false,interv=[0],roundss=[0],strat=:dg,eff=0.0,et=0)
        
    for inter = interv,roun = roundss
        P = SCHparameters(infection_human = 2.6542e-03,#6.688e-04,
        infection_snail = 3.0898e-04,
        beta = 3.2000e-02,
        immunity_parameter = 1.1836e+01,
        mu_w = 6.7130e-01,
        grid_size_snail = 1000,
        method=2,
        worms_lim_diag = wld,
        file_index = idx,
        med_efficacy = eff,
        rounds = roun,
        treat_strat = strat,
        Interval = inter,
        treatment = treat,
        eff_type = et
        )


        age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)

        #P = SCHparameters(beta = 0.01,file_index = 2,grid_size_human = 1000,infection_human = 0.000003,infection_snail = 0.0000001)
        n_sim = 1000
        run_s(P,n_sim)
    end

end


function run3(idx,wld=1,treat=false,interv=[0],roundss=[0],strat=:dg,eff=0.0,et=0)
    for inter = interv,roun = roundss

        P = SCHparameters(infection_human = 4.96e-04,#6.688e-04,
        infection_snail = 9.3862e-04,
        beta = 1.2000e-01,
        immunity_parameter = 1.4076e+01,
        mu_w = 8.955e-01,
        grid_size_snail = 2000,
        method=2,
        worms_lim_diag = wld,
        file_index = idx,
        med_efficacy = eff,
        rounds = roun,
        treat_strat = strat,
        Interval = inter,
        treatment = treat,
        eff_type = et
        )


        age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)

        #P = SCHparameters(beta = 0.01,file_index = 2,grid_size_human = 1000,infection_human = 0.000003,infection_snail = 0.0000001)
        n_sim = 1000
        run_s(P,n_sim)
    end
end

#run1(0)
#run2(0)
#run3(0)
#= 
run1(4,2,true,[0.0],[1],:dg,0.8,0)
run2(4,2,true,[0.0],[1],:dg,0.8,0)
run3(4,2,true,[0.0],[1],:dg,0.8,0)
run1(5,2,true,[0.0],[1],:dg,1.0,1)
run2(5,2,true,[0.0],[1],:dg,1.0,1)
run3(5,2,true,[0.0],[1],:dg,1.0,1) =#

#= 
run1(4,2,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0)
run2(4,2,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0)
run3(4,2,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0)

run1(5,2,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1)
run2(5,2,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1)
run3(5,2,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1) =#

#= 
run1([0],[1],:dg,0.8,1,0)
run1([0.5;1;2],[2;4;6;8;10],:dg,0.8,1,0)

run2([0],[1],:dg,0.8,1,0)
run2([0.5;1;2],[2;4;6;8;10],:dg,0.8,1,0)

run3([0],[1],:dg,0.8,1,0)
run3([0.5;1;2],[2;4;6;8;10],:dg,0.8,1,0)

#run1([0],[1],:dg,1.0,2,1)
#run1([0.5;1;2],[2;4;6;8;10],:dg,1.0,2,1)


#run2([0],[1],:dg,1.0,2,1)
#run2([0.5;1;2],[2;4;6;8;10],:dg,1.0,2,1)

#run3([0],[1],:dg,1.0,2,1)
#run3([0.5;1;2],[2;4;6;8;10],:dg,1.0,2,1)


run1([0],[1],:total,0.8,3,0)
run1([0.5;1;2],[2;4;6;8;10],:total,0.8,3,0)


run2([0],[1],:total,0.8,3,0)
run2([0.5;1;2],[2;4;6;8;10],:total,0.8,3,0)

run3([0],[1],:total,0.8,3,0)
run3([0.5;1;2],[2;4;6;8;10],:total,0.8,3,0)


#run1([0],[1],:total,1.0,4,1)
#run1([0.5;1;2],[2;4;6;8;10],:total,1.0,4,1)


#run2([0],[1],:total,1.0,4,1)
#run2([0.5;1;2],[2;4;6;8;10],:total,1.0,4,1)

#run3([0],[1],:total,1.0,4,1)
#run3([0.5;1;2],[2;4;6;8;10],:total,1.0,4,1) =#