
function Prevalence_Calc(age_prevalence_data::String,folder::String,P::SCHparameters)

    Data = readdlm("Age.dat",header = false)

    max_value = maximum(Data[:,1])
    n_bins = Int(ceil(max_value/P.age_interval))
    
    if Int(max_value%P.age_interval) > 0
        n_bins += 1
        groups = Vector{String}(undef,n_bins)

        for i = 1:n_bins
            groups[i] = "$((i-1)*5)-$(i*5-1)"
        end
    else
        groups = Vector{String}(undef,n_bins)
        for i = 1:(n_bins-1)
            groups[i] = "$((i-1)*5)-$(i*5-1)"
        end
        groups[n_bins] = "$((n_bins-1)*5)-$(n_bins*5)"
    end

    n_neg = zeros(Float64,n_bins,size(Data)[2]-1)
    n_pos = zeros(Float64,n_bins,size(Data)[2]-1)
    for i=1:size(Data)[1]
        for j = 1:(size(Data)[2]-1)
            if Data[i,j+1] == "Positivo"
                n_pos[Int(ceil(Data[i,1]/P.age_interval)),j]+=1
            else
                n_neg[Int(ceil(Data[i,1]/P.age_interval)),j]+=1
            end
        end
    end

    prev = n_pos./(n_neg+n_pos)
   

    writedlm(string(folder,"prevalence.dat"),[groups prev])

    return groups,prev
end

#groups,prevalence = Prevalence_Calc(age_prevalence_data,folder,P)
#=
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
    

end=#

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

#groups,prevalence = assimetric_prevalence(age_prevalence_data,folder)

function Plotting_prev(folder,groups,prevalence)
    color_ =  [:blue :red]
    gr()
    Plots.plot(prevalence,xticks = (1:length(groups), groups),linetype=:line,linestyle = :solid,linecolor = color_,lw=3,label = ["Kato-Katz" "HTX"],xrotation = rad2deg(pi/3))
    Plots.plot(prevalence,seriestype=:bar,layout="grouped-bar",linecolor = color_,lw=3,label = ["Kato-Katz" "HTX"],xrotation = rad2deg(pi/3))
    Plots.savefig("prevalence.png")
end
