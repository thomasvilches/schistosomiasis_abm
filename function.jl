function random_initial_inf(h1::Array{Human{Int64},1},P::SCHparameters)
    r = rand(1:P.grid_size_human)
    h1[r].n_worms_f = Int(P.init_n_worms/2)
    h1[r].n_worms_m = Int(P.init_n_worms/2)
    h1[r].health = INF
    h1[r].health_immunity = INF
    h1[r].time_infected += 1
    return r
end

function event_test(h::Human,P::SCHparameters,cercaria_reservoir::Int64)

    lambda::Float64 = 0.0
    event::Int64 = 0
    inf_::Int64 = 0
    die_::Int64 = 0
    
    if h.time_infected > P.immunity_period
        lambda = P.beta*P.infection_human
    else
        lambda = P.infection_human
    end
    prob_inf::Float64 = 1-exp(-lambda/365)
    
    r = rand()
    if (h.n_worms_f+h.n_worms_m) > 0 && r <= (1-exp(-(h.n_worms_f+h.n_worms_m)*P.mu_w/365))
        if rand() < h.n_worms_f/(h.n_worms_f+h.n_worms_m)
            h.n_worms_f-= 1
        else
            h.n_worms_m-= 1
        end
    end

    if rand() < 1-(1-prob_inf)^cercaria_reservoir
        #prob_inf = 1-exp(-lambda/365)
        r = Binomial(max(0,Int(cercaria_reservoir-1)),prob_inf)
        r1 = rand(r)+1
        r2 = Int(rand(Binomial(max(0,Int(r1)),0.5)))
        h.n_worms_f += r2
        h.n_worms_m += Int(r1-r2)
        cercaria_reservoir -= Int(r1) 
    end
  
    return cercaria_reservoir

end

function update_population(h1::Array{Human{Int64},1},P::SCHparameters,t::Int64)
    n_infected::Int64 = 0
    for i = 1:P.grid_size_human
        if (h1[i].n_worms_f+h1[i].n_worms_m) > 0
            h1[i].health = INF
            n_infected += 1
            h1[i].health_immunity = INF
        else 
            h1[i].health = SUSC
        end

        if h1[i].health_immunity == INF
            h1[i].time_infected += 1
        end
        
        if h1[i].age_days == 365
            h1[i].age += 1
            h1[i].age_days = 0

            ages_lim = limit_of_ages()
            g = findfirst(x-> x == h1[i].age,ages_lim[:,1])
            if g != nothing
                h1[i].group = g
            end

        else 
            h1[i].age_days += 1
        end

    end
    return n_infected
end

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


function update_cercaria!(P::SCHparameters,cercaria_inc_vec::Array{Int64,1},t::Int64)
    
    mortes::Int64 = 0

    for i = 1:16
        j = t-i
        if j > 0
            pt = exp(P.a/P.b*(1-exp(P.b*i*1000/16)))
            pt = 1-pt
            mortes += Int(round(pt*cercaria_inc_vec[j]))
            cercaria_inc_vec[j] -= Int(round(pt*cercaria_inc_vec[j]))
        else
            break
        end
    end

    i = 17
    j = t-i
    if j > 0
        mortes += Int(round(cercaria_inc_vec[j]))
        cercaria_inc_vec[j] = 0
    end

    return mortes
end