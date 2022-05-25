
function Plotting_sims(P::SCHparameters)
    #P = SCHparameters(max_worms = 21)
    age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)

    ###reading time series for infection prevalence
    
    #time_series = readdlm(string(folder,time_data))
    ##in teste.jl

    ###the next analysis are done considering the last day of simulation
    #make sure to run  it long enough in order to reach de equilibrium point
    age_prevalence_data,folder,time_data,age_data,inf_data,group_data,worms_data,cercaria_data,worms_c_data = file_names(P)

    age_ = readdlm(string(folder,age_data),Int64)
    group_ = readdlm(string(folder,group_data),Int64)
    worms_ = readdlm(string(folder,worms_data),Int64)
    inf_ = readdlm(string(folder,inf_data),Int64)

    n_groups = maximum(group_)
    mean_worms_group = zeros(Float64,n_groups)
    variance_worms_group = zeros(Float64,n_groups)
    n_people_group = zeros(Float64,n_groups)
    n_people_inf_group = zeros(Float64,n_groups)
    n_sim_zero = 0
    for j = 1:size(worms_)[2]
        if sum(inf_[:,j]) > 0
            for i = 1:size(worms_)[1]
                n_people_group[group_[i,j]] += 1
                n_people_inf_group[group_[i,j]] += inf_[i,j]
                mean_worms_group[group_[i,j]] += worms_[i,j]
            end
        else
            global n_sim_zero+=1
        end
    end
    mean_worms_group = mean_worms_group./n_people_group

    for j = 1:size(worms_)[2]
        for i = 1:size(worms_)[1]
            if inf_[i,j]==1
                variance_worms_group[group_[i,j]] += (worms_[i,j]-mean_worms_group[group_[i,j]])^2
            end
        end
    end

    for i = 1:length(variance_worms_group)
        variance_worms_group[i] = variance_worms_group[i]/(n_people_inf_group[i]-1)
    end
    variance_worms_group = sqrt.(variance_worms_group)
    prevalence = n_people_inf_group./n_people_group

    writedlm(string(folder,"Summary_data_","$(P.file_index)",".dat"),[prevalence mean_worms_group variance_worms_group])

    Plots.plot(prevalence,label = "", xlabel="age group",ylabel="Prevalence",lw = 3, colour = :red)
   # Plots.plot(mean_worms_group,label="",xlabel="age group",ylabel="mean worm burden",lw = 3, colour = :red)
    #Plots.histogram(worms_,label="",xlabel="worm burden",ylabel="frequency")
   # Plots.plot(variance_worms_group,label="",xlabel="age group",ylabel="worm burden distribution",lw = 3, colour = :red)
    return n_sim_zero
    
end

Plotting_sims(P)




###Plotting barplot

using DataFrames, VegaLite, VegaDatasets
group_lim = limit_of_ages()
Data = readdlm("Age.dat",header = false)
group_lim[size(group_lim)[1],size(group_lim)[2]] = maximum(Data[:,1])
group_lim_2 = zeros(Float64,size(group_lim)[1],2)
for i=1:size(group_lim)[1]
    for j=1:size(group_lim)[2]
        group_lim_2[i,j] = group_lim[i,j]
    end
end
groups = [group_lim_2; group_lim_2]
method = Array{String,1}(undef,length(groups[:,1]))
for i = 1:Int(length(method)/2)
    method[i] = "Kato-Katz"
    method[length(method)-i+1] = "HTX"
end

groups[:,2]=groups[:,2].+1
prevalence2 = [prevalence[:,1]; prevalence[:,2]-prevalence[:,1]]

df = DataFrame(bin_start=groups[:,1], bin_end=groups[:,2], count=prevalence2,method = method)

p=df |> @vlplot(
:bar,width = 450,height = 450,
x={:bin_start, bin=:binned, title = "Idade (anos)",axis={tickStep=1,titleFontSize=30,labelFontSize=15,},type = :quantitative},
x2=:bin_end, 
y={"sum(count)",title = "Prevalência",axis={titleFontSize=30,labelFontSize=15,}},
color={
    :method,
    scale = {
        domain = ["Kato-Katz","HTX"],
        range = [:orange,:blue]
    },
    legend = {
        title = "Método de Diagnóstico",
        labelFontSize = 20
    }
},
background=:white)|>FileIO.save("DiagMethod.png")

using Distributions,Bootstrap,DelimitedFiles,Plots
using VegaLite
using Colors
using Images
using FileIO
#time_series = readdlm("results_teste/inf_time_series_r_1.dat")
time_series = readdlm(string(folder,time_data))
n = size(time_series)[1]
data1 = zeros(Float64,n,3)

for i = 1:n
    vec = time_series[i,:]
    mbst = bootstrap(mean,vec,BasicSampling(100))
    sci = confint(mbst,BasicConfInt(0.95))
    data1[i,1] = sci[1][1]
    data1[i,2] = sci[1][2]
    data1[i,3] = sci[1][3]
    #data1[i,2] = minimum(vec)
    #data1[i,3] = maximum(vec)
end
dt = DataFrame(Time = (1:n)/365,Mean = data1[:,1],CIL = data1[:,2],CIH = data1[:,3])
p_prev = dt |>@vlplot( 
    width = 450,
    height = 450,
    x = {
        field="Time",
        axis={title = "Tempo (anos)",
        titleFontSize = 32,
        labelFontSize=15,
        }
    },
    background = :white
)+@vlplot(mark={
    type=:errorband,
    opacity = 0.5
    },
    encoding={ 
        y = {field = "CIL",
        type = :quantitative,
            axis={
                title = "Número de indivíduos infectados",
                titleFontSize = 24,
                labelFontSize=15,
            },   
            
           # scale = {"zero",false}
        },
        y2 = {
            field = "CIH",
            type = :quantitative
        },
        color = {value=:lightgrey}
    },
)+@vlplot(mark=:line,
    encoding ={
        y={
            field="Mean",
            type = :quantitative,
            lw = 1
        },
        color = {value=:black}
    },
    background = :white
)|>FileIO.save("time_series.png")

####Prevalence nice plot

prevalence = readdlm(string(folder,"Summary_data_","$(P.file_index)",".dat"))[:,1]
position = 1:length(prevalence)
groups_label = ["0-1";"2-4";"5-9";"10-14";"15-19";"20-24";"25-29";"30-34";"35-39";"40-44";"45-49";"50-54";"55-59";"60-64";"65+"]

dt = DataFrame(Position = position,Group = groups_label,Prevalence = prevalence)
dt|> @vlplot(
    width = 450,
    height = 450,
    x = {
        field = "Position",
        type = "nominal",
        axis={
            title = "Grupo de idade",
            titleFontSize = 32,
            #labelAngle = -45,
            labelFontSize = 15,
        },
    },
    background = :white
)+@vlplot(mark=:line,
    encoding ={
        y={
            field="Prevalence",
            type = "quantitative",
            lw = 1
        },
        color = {value=:lightgrey}
    },
    background = :white
)+@vlplot(
    
    mark={
        type=:point,
        size=50,
        filled = true,
    },
    encoding={ 
        y = {
            field = "Prevalence",
            type = "quantitative",
            axis={
                title = "Prevalência",
                titleFontSize = 24,
                labelFontSize = 15,
            },   
            
           # scale = {"zero",false}
        },
        color = {value=:red}
    },
    background = :white
)|>FileIO.save("prevalence.png")

#################################
worms = [1; 11; 21]

prevalence = zeros(Float64,15,length(worms))
MWB = zeros(Float64,15,length(worms))
WBD = zeros(Float64,15,length(worms))

for i = 1:length(worms)
    data = readdlm(string(folder,"Summary_data_","$(worms[i])",".dat"))
    prevalence[:,i] = data[:,1]
    MWB[:,i] = data[:,2]
    WBD[:,i] = data[:,3]
end

labels = ["MW=1";"MW=11";"MW=21"]
cls = [:red :blue :green]
Plots.plot(prevalence,label = labels,xlabel="age group",ylabel="Prevalence",lw = 3,color = cls)
Plots.plot(MWB,label=labels,xlabel="age group",ylabel="mean worm burden",lw = 3, colour = cls)
Plots.plot(WBD,label=labels,xlabel="age group",ylabel="worm burden distribution",lw = 3, colour = cls)


age_ = zeros(Int64,P.grid_size_human)

for i = 1:P.grid_size_human
    age_[i] = humans[i].age
end

Plots.histogram(age_)



d = Exponential(1/0.015)

r = rand(d,1000)

Plots.histogram(r)

age_death_v = zeros(Int64,P.grid_size_human)

for i = 1:P.grid_size_human
    age_death_v[i] = humans[i].death_age
end

Plots.histogram(age_death_v)

r = rand(10000)

Plots.histogram(age_,labels ="")



age_ = zeros(Int64,P.grid_size_human)

for i = 1:P.grid_size_human
    age_[i],trash = age_death(P,0)
end


n = size(age_)[1]*size(age_)[2]
reshape(age_,n,1)
age_ = age_[:,1]
Plots.histogram(age_)




real_prevalence,estimated_prevalence,n_sim_zero = calc_prevalence_sim(P)
groups,prevalence = assimetric_prevalence()
groups_label = ["0-1";"2-4";"5-9";"10-14";"15-19";"20-24";"25-29";"30-34";"35-39";"40-44";"45-49";"50-54";"55-59";"60-64";"65+"]

Plots.plot([prevalence[:,P.method] estimated_prevalence],label = ["Campo" "in-silico"],linewidth=4,xlabel = "Idade",ylabel = "Prevalência",ylim=(0,1),color = [:red :blue])
Plots.xticks!(1:15,groups)|>FileIO.save("prevalence-poster.png")

Plots.plot([prevalence[:,P.method] estimated_prevalence real_prevalence],label = ["Campo" "Estimada in-silico" "\"Real\" in-silico"],linewidth=4,xlabel = "Idade",ylabel = "Prevalência",ylim=(0,1),color = [:red :blue :green])