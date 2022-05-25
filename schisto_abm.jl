include("population.jl")
include("function.jl")
function init_human(P::SCHparameters)
    [Human(i,P) for i = 1:P.grid_size_human]
end

function init_snail(P::SCHparameters)
    [Snail(i,P) for i = 1:P.grid_size_snail]
end

function main(sim_index::Int64,P::SCHparameters)
    #println("$sim_index")

    Random.seed!(276*sim_index)

    humans = init_human(P)
    snails = init_snail(P)
    setup_demographic_human(humans,P)
    f_i = random_initial_inf(humans,P)

    ###running time
    infected_n = zeros(Int64,P.sim_time)
    cercaria_vec = zeros(Int64,P.sim_time)
    cercaria_inc_vec = zeros(Int64,P.sim_time)
    miracidium_reservoir::Int64 = 0
    cercaria_reservoir::Int64 = 0
    number_of_pairs::Float64 = 0.0
    number_of_inf_snails::Float64 = 0.0
    for t = 1:P.sim_time
        
        number_of_pairs = 0.0
        number_of_inf_snails = 0.0
        

        for i=1:P.grid_size_human #calculating the MWB over time
            number_of_pairs += min(humans[i].n_worms_f,humans[i].n_worms_m)
        end

        for i=1:P.grid_size_snail
            #=if snails[i].health == INF
                number_of_inf_snails += 1
            end=#
            y = Int(floor(Int(snails[i].health)/2))
            number_of_inf_snails += y
        end

        r = Poisson(P.oviposition_rate/365*P.development_probability*number_of_pairs)
        miracidium_reservoir = rand(r)
        
        r = Poisson(P.shedding_rate/365*number_of_inf_snails)
        cercaria_reservoir = rand(r)
        cercaria_vec[t] = cercaria_reservoir
        #snail dynamic
        for i = 1:P.grid_size_snail
            if snails[i].health == SUSC
                lambda_s = 1-exp(-P.infection_snail/365)
                if rand() <= 1-(1-lambda_s)^miracidium_reservoir#(1-exp(-P.infection_snail/365*miracidium_reservoir))
                    snails[i].health = LAT
                    miracidium_reservoir -= 1
                end

            elseif snails[i].health == LAT
                if rand() <= (1-exp(-P.mu_line_s/365))
                    snails[i].health = SUSC
                    snails[i].time_latent = 0
                elseif snails[i].time_latent >= P.latent_period
                    snails[i].health = INF
                else 
                    snails[i].time_latent += 1
                end
            else 
                if rand() <= (1-exp(-P.mu_2line_s/365))
                    snails[i].health = SUSC
                    snails[i].time_latent = 0
                end
            end
        end
        #human dynamic
        for i=1:P.grid_size_human ##since there is no contact structure, it doesn't matter the update sequence
            ##First, let's test if the human will die
            r = rand(1:P.grid_size_human)
            if (humans[r].age >= humans[r].death_age && humans[r].age_days >= humans[r].death_days)
                #if so, they are replaced by a newborn
                human_return(humans[r],P)
            else
                cercaria_reservoir = event_test(humans[r],P,cercaria_reservoir)  
            end
        end

        infected_n[t] = update_population(humans,P,t)
        #=cercaria_vec[t] = cercaria_reservoir

        p = 1-exp(-P.miracidium_mortality/365)
        miracidium_reservoir -= min(miracidium_reservoir,rand(Binomial(max(0,miracidium_reservoir),p)))
        

        #p = 1-exp(-P.cercaria_mortality/365)
        #cercaria_reservoir -= min(cercaria_reservoir,rand(Binomial(max(0,cercaria_reservoir),p)))
        mortes = update_cercaria!(P,cercaria_inc_vec,t)
        cercaria_reservoir =- min(mortes,cercaria_reservoir)

        r = Poisson(P.oviposition_rate/365*P.development_probability*number_of_pairs)
        miracidium_reservoir += rand(r)

        r = Poisson(P.shedding_rate/365*number_of_inf_snails)
        y = rand(r)
        cercaria_reservoir += y
        cercaria_inc_vec[t] = y=#
    end

    #=idade_p=zeros(Float64,15)
    idade_n=zeros(Float64,15)
    for i = 1:P.grid_size_human
        g::Int64 = humans[i].group
        if humans[i].n_worms>0
            idade_p[g]+=1
        else
            idade_n[g]+=1
        end
    end=#
    age_sim = zeros(Int64,P.grid_size_human)
    inf_sim = zeros(Int64,P.grid_size_human)
    worms_sim = zeros(Int64,P.grid_size_human)
    worms_c_sim = zeros(Int64,P.grid_size_human)
    group_sim = zeros(Int64,P.grid_size_human)
    for i = 1:P.grid_size_human
        age_sim[i] = humans[i].age
        group_sim[i] = humans[i].group

        if (humans[i].n_worms_f+humans[i].n_worms_m) > 0
            worms_sim[i] = (humans[i].n_worms_f+humans[i].n_worms_m)
            inf_sim[i] = 1
            worms_c_sim[i] = min(humans[i].n_worms_f,humans[i].n_worms_m)
        end
    end


    return infected_n,age_sim,group_sim,inf_sim,worms_sim,cercaria_vec,worms_c_sim #idade_p,idade_n,
end
