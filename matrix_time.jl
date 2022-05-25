
using Distributed
addprocs(4)
@everywhere using DelimitedFiles
@everywhere using Parameters
@everywhere using Statistics
@everywhere include("parameters.jl")

method = 2
snail_pop = 500
file = 1
n_pop_ga = 500
n_gen = 20
n_boots = 150
limite = 1

rnd = 1
intv = 0.0

folder0 = "Cluster/treatment/"

println("$intv $rnd")
P = SCHparameters(file_index = file,method = method, grid_size_snail = snail_pop, rounds = rnd,Interval = intv,treatment = true)

    
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


data_time = readdlm(string(folder,"inf_time_series_r_found_",file,".dat"),Int64,header = false)
#data_group = readdlm(string(folder,"/group_data_r_",file,".dat"),Int64,header = false)

m = Array{Float64,2}(undef,size(data_time,1),n_boots)
l = 1:size(data_time,2)

m = pmap(1:n_boots) do i
    #println(i)
    mean(data_time[:,rand(l,1000)],dims=2)
end

m = hcat(m...)
writedlm(string(folder,"/matrix_time_found_$(file).dat"),m)




r = [2;4;6;8;10]
i = [0.5;1.0;2.0]


for intv = i,rnd = r

        println("$intv $rnd")
    P = SCHparameters(file_index = file,method = method, grid_size_snail = snail_pop, rounds = rnd,Interval = intv,treatment = true)

        
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


    data_time = readdlm(string(folder,"inf_time_series_r_found_",file,".dat"),Int64,header = false)
    #data_group = readdlm(string(folder,"/group_data_r_",file,".dat"),Int64,header = false)

    m = Array{Float64,2}(undef,size(data_time,1),n_boots)
    l = 1:size(data_time,2)

    m = pmap(1:n_boots) do i
        #println(i)
        mean(data_time[:,rand(l,1000)],dims=2)
    end

    m = hcat(m...)
    writedlm(string(folder,"/matrix_time_found_$(file).dat"),m)

end



##########################

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

writedlm(string(folder,"/matrix_$(file)_$snail_pop.dat"),m)




####################################### zeros
P = SCHparameters(file_index = file,method = method, grid_size_snail = snail_pop, rounds = 10,Interval = 0.5,treatment = true)

if !P.treatment
    folder = string(folder0,"result_$(P.grid_size_snail)_method_$(P.method)/")#"Cluster/fixed_seed/size_500_method_2/"#
else
    folder = string(folder0,"result_$(P.grid_size_snail)_method_$(P.method)_$(P.rounds)_$(P.Interval)/")
end

data_worm = readdlm(string(folder,"/worms_data_r_",file,".dat"),Int64,header = false)
data_group = readdlm(string(folder,"/group_data_r_",file,".dat"),Int64,header = false)

n = 0

for i = 1:size(data_worm,2)
    if sum(data_worm[:,i]) == 0
        global n +=1
    end
end
n


################################################################################################33



for snail_pop in [500;1000;2000]
    for file = 1:4
        method = 2
        #snail_pop = 500
        #file = 1
        n_pop_ga = 500
        n_gen = 20
        n_boots = 150
        limite = 1

        rnd = 1
        intv = 0.0

        folder0 = "Cluster/treatment/"

        println("$intv $rnd $file $snail_pop")
        P = SCHparameters(file_index = file,method = method, grid_size_snail = snail_pop, rounds = rnd,Interval = intv,treatment = true)

            
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


        data_time = readdlm(string(folder,"inf_time_series_r_found_",file,".dat"),Int64,header = false)
        #data_group = readdlm(string(folder,"/group_data_r_",file,".dat"),Int64,header = false)

        m = Array{Float64,2}(undef,size(data_time,1),n_boots)
        l = 1:size(data_time,2)

        m = pmap(1:n_boots) do i
            #println(i)
            mean(data_time[:,rand(l,1000)],dims=2)
        end

        m = hcat(m...)
        writedlm(string(folder,"/matrix_time_found_$(file).dat"),m)




        r = [2;4;6;8;10]
        i = [0.5;1.0;2.0]


        for intv = i,rnd = r

                println("$intv $rnd")
            P = SCHparameters(file_index = file,method = method, grid_size_snail = snail_pop, rounds = rnd,Interval = intv,treatment = true)

                
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


            data_time = readdlm(string(folder,"inf_time_series_r_found_",file,".dat"),Int64,header = false)
            #data_group = readdlm(string(folder,"/group_data_r_",file,".dat"),Int64,header = false)

            m = Array{Float64,2}(undef,size(data_time,1),n_boots)
            l = 1:size(data_time,2)

            m = pmap(1:n_boots) do i
                #println(i)
                mean(data_time[:,rand(l,1000)],dims=2)
            end

            m = hcat(m...)
            writedlm(string(folder,"/matrix_time_found_$(file).dat"),m)

        end
    end
end

#= 

function assimetric_prevalence()

    Data = readdlm("Age.dat",header = false)
    
    group_lim = limit_of_ages()
    group_lim[size(group_lim)[1],2] = Int(maximum(Data[:,1]))
    groups = Array{String,1}(undef,size(group_lim)[1])
    for i = 1:(size(group_lim)[1]-1)
        groups[i] = string("$(group_lim[i,1])","-","$(group_lim[i,2])")
    end
    groups[end] = string("$(group_lim[end,1])","+")

    n_neg = zeros(Float64,length(groups),size(Data)[2]-1)
    n_pos = zeros(Float64,length(groups),size(Data)[2]-1)

    for i=1:size(Data)[1]
        pos = findfirst(x->x>=Data[i,1],group_lim[:,2])
        for j = 1:(size(Data)[2]-1)
            
            if Data[i,j+1] == "Positivo"
                n_pos[pos,j]+=1
            else
                n_neg[pos,j]+=1
            end
        end
    end
    prev = n_pos./(n_neg+n_pos)
    return groups,prev
end



function limit_of_ages() 
    groups_ages = [
        0 4;
        5 9;
        10 14;
        15 19;
        20 24;
        25 29;
        30 34;
        35 39;
        40 49;
        50 64;
        65 1000
    ]
    return groups_ages
end

groups,prevalence = assimetric_prevalence()
groups_label = ["0-1";"2-4";"5-9";"10-14";"15-19";"20-24";"25-29";"30-34";"35-39";"40-44";"45-49";"50-54";"55-59";"60-64";"65+"]

writedlm("prevalence_field.dat",[groups prevalence])




function calc_prevalence_sim(P::SCHparameters)

    age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)
    age_ = readdlm(string(folder,age_data),Int64)
    group_ = readdlm(string(folder,group_data),Int64)
    worms_ = readdlm(string(folder,worms_data),Int64)
    worms_c = readdlm(string(folder,worms_c_data),Int64)
    inf_ = readdlm(string(folder,inf_data),Int64)

    n_groups = maximum(group_)
    mean_worms_group = zeros(Float64,n_groups)
    #variance_worms_group = zeros(Float64,n_groups)
    n_people_group = zeros(Float64,n_groups)
    estimated_n_people_inf_group = zeros(Float64,n_groups)
    real_n_people_inf_group = zeros(Float64,n_groups)
    n_sim_zero::Int64 = 0

    estimated_prevalence = zeros(Float64,n_groups)
    real_prevalence = zeros(Float64,n_groups)
    for j = 1:size(worms_)[2]
        if sum(inf_[:,j]) > 0
            for i = 1:size(worms_)[1]
                n_people_group[group_[i,j]] += 1
                real_n_people_inf_group[group_[i,j]] += inf_[i,j]
                if worms_c[i,j] >= P.worms_lim_diag
                    estimated_n_people_inf_group[group_[i,j]] += inf_[i,j]
                end
                mean_worms_group[group_[i,j]] += worms_[i,j]
            end
        else
            n_sim_zero += 1

            for i = 1:size(worms_)[1]
                n_people_group[group_[i,j]] += 1
            end 
        end
    end

    for i = 1:n_groups
        if n_people_group[i] > 0
            real_prevalence[i] = real_n_people_inf_group[i]/n_people_group[i]
            estimated_prevalence[i] = estimated_n_people_inf_group[i]/n_people_group[i]
        else 
            real_prevalence[i] = 0.0
            estimated_prevalence[i] = 0.0
        end
    end
    #writedlm(string(folder,"Summary_data_","$(P.file_index)",".dat"),[real_prevalence estimated_prevalence mean_worms_group])
    
    return real_prevalence,estimated_prevalence,n_sim_zero

end


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


    return age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data
end

real_prevalence,estimated_prevalence,n_sim_zero = calc_prevalence_sim(P)

estimated_prevalence =#