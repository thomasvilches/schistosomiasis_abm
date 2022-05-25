using DelimitedFiles
#using DataFrames
using Parameters
#using Plots
#using Colors
#using PyPlot
#using Plotly
#using GR
#using Cairo,Compose,Fontconfig
using Distributions

#@everywhere include("parameters.jl")
#@everywhere include("schisto_abm.jl")

function file_names(P::SCHparameters)

    age_prevalence_data = "Age_r_$(P.file_index).dat"
    folder = "result_$(P.grid_size_snail)_method_$(P.method)/"#"Cluster/fixed_seed/size_500_method_2/"#
    time_data = "inf_time_series_r_$(P.file_index).dat"
    age_data = "age_data_r_$(P.file_index).dat"
    inf_data = "inf_data_r_$(P.file_index).dat"
    group_data = "group_data_r_$(P.file_index).dat"
    worms_data = "worms_data_r_$(P.file_index).dat"
    cercaria_data = "cercaria_data_r_$(P.file_index).dat"
    worms_c_data = "worms_c_data_r_$(P.file_index).dat"
    time_data_found = "inf_time_series_r_found_$(P.file_index).dat"

    return age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data
end

function data_process(results,n_sim,P::SCHparameters)
    
    age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)

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

function run(P,n_sim)

    age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)
    results = pmap(x->main(x,P),1:n_sim)
    data_process(results,n_sim,P)

end