mutable struct Snail{Int64} ## mutable structs are stored on the heap 
    index::Int64
    #age::Union{Int64, Nothing}
    #age_days::Int64   
    #group::Union{Int64, Nothing}  
    time_latent::Int64
    #n_worms::Int64
    health::HEALTH
    #health_immunity::HEALTH
    #death_age::Int64
    #death_days::Int64
    #Here I added the strain matrix, a vector for the time the strain showed up, number of strains in the body, and
    #the efficacy of the vaccine against the strain that was transmitted (necessary for the asymp-symp trial)
    Snail{Int64}(P::SCHparameters) = new(0,-1,SUSC)

    function Snail(idx,P::SCHparameters)
        s = Snail{Int64}(P)
        s.index = idx
        return s
    end
end


mutable struct Human{Int64} ## mutable structs are stored on the heap 
    index::Int64
    age::Union{Int64, Nothing}
    age_days::Int64   
    group::Union{Int64, Nothing}  
    time_infected::Int64
    n_worms_f::Int64
    n_worms_m::Int64
    health::HEALTH
    health_immunity::HEALTH
    death_age::Int64
    death_days::Int64
    #Here I added the strain matrix, a vector for the time the strain showed up, number of strains in the body, and
    #the efficacy of the vaccine against the strain that was transmitted (necessary for the asymp-symp trial)
    Human{Int64}(P::SCHparameters) = new(0,nothing,0,0,-1,0,0,SUSC,SUSC,-1,-1)

    function Human(idx,P::SCHparameters)
        h = Human{Int64}(P)
        h.index = idx
        return h
    end
end



function human_return(h::Human,P::SCHparameters)
    
    h.age = 0
    h.age_days = 0
    h.group = 1
    h.time_infected = -1
    h.health = SUSC
    h.health_immunity = SUSC
    h.n_worms_f = 0
    h.n_worms_m = 0
    age_1 = -(1/P.mu_h)*log(1-rand())
    age_2 = Int(floor(age_1))
    age_1 = age_1-age_2
    h.death_age = min(age_2,100)
    h.death_days = Int(floor(365*age_1))
end

function age_death(P::SCHparameters,age::Int64,day::Int64)

    #d = Exponential(1/P.mu_h)
    age_d::Int64 = -1
    day_d::Int64 = -1
    aux::Bool = false

  #=  for i = 1:(P.max_trials)
        r = rand(d)
        if r > age
            if r > P.age_max
                age_d = P.age_max
                days_d = Int(floor(365*rand()))
                aux = true
            else
                age_d = Int(floor(r))
                days_d = Int(floor(365*(r-age_d)))
                aux = true
            end
            break
        end
    end
    
    if aux == false
        age_d = P.age_max
        days_d = rand(1:365)
    
    end=#

    if rand() <= (1-exp(-P.mu_h*age))
        age_d = age
        days_d = rand(day:365)
    else
        age_d = P.max_age
        days_d = rand(1:365)
    end


    return age_d,days_d
end

function setup_demographic_human(h1,P::SCHparameters)
    dist, AgeMin, AgeMax = distribution_age()
    for i = 1:length(h1)
        if h1[i].age == nothing
           # rn = rand()
           # g = findfirst(x -> rn <= x, dist) ## THERE IS A CLOSURE BUG HERE. DO NOT USE INBOUNDS. SEE FTEST() testing. This was fixed. See my issue on github
            #if g !== nothing
                age_1 = -(1/P.mu_h)*log(1-rand())
                age_2 = Int(floor(age_1))
                age_1 = age_1-age_2
                h1[i].death_age = min(age_2,100)
                h1[i].death_days = Int(floor(365*age_1))
                h1[i].age = rand(0:h1[i].death_age)
                h1[i].age_days = rand(1:365)

                ages_lim = limit_of_ages()
                h1[i].group = findfirst(x-> x >= h1[i].age,ages_lim[:,2])
                #h1[i].death_age = P.max_age
                #h1[i].death_days = Int(floor(rand()*365))
           # end
        end

    end
end

function distribution_age()
    ProbBirthAge = Vector{Float64}(undef, 15)
    SumProbBirthAge = Vector{Float64}(undef, 15)
    AgeMin = Vector{Int64}(undef, 15)
    AgeMax = Vector{Int64}(undef, 15)
     ProbBirthAge[1] = 0.006211180124223602
     ProbBirthAge[2] = 0.13442768411712513
     ProbBirthAge[3] = 0.30789707187222715
     ProbBirthAge[4] = 0.4622892635314996
     ProbBirthAge[5] = 0.5496894409937888
     ProbBirthAge[6] = 0.611357586512866
     ProbBirthAge[7] = 0.673469387755102
     ProbBirthAge[8] = 0.7213842058562555
     ProbBirthAge[9] = 0.7724046140195208
    ProbBirthAge[10] = 0.8176574977817213
    ProbBirthAge[11] = 0.8487133984028392
    ProbBirthAge[12] = 0.8877551020408162
    ProbBirthAge[13] = 0.9139307897071871
    ProbBirthAge[14] = 0.9418811002661933
    ProbBirthAge[15] = 1.0
    SumProbBirthAge = ProbBirthAge #cumsum(ProbBirthAge)
    #SumProbBirthAge[end] = 1.0

    AgeMin[1] = 0;
    AgeMax[1] = 1;
    
    AgeMin[2] = 2;
    AgeMax[2] = 4;

    AgeMin[3] = 5;
    AgeMax[3] = 9;

    AgeMin[4] = 10;
    AgeMax[4] = 14;

    AgeMin[5] = 15;
    AgeMax[5] = 19;

    AgeMin[6] = 20;
    AgeMax[6] = 24;

    AgeMin[7] = 25;
    AgeMax[7] = 29;

    AgeMin[8] = 30;
    AgeMax[8] = 34;

    AgeMin[9] = 35;
    AgeMax[9] = 39;

    AgeMin[10] = 40;
    AgeMax[10] = 44;

    AgeMin[11] = 45;
    AgeMax[11] = 49;

    AgeMin[12] = 50;
    AgeMax[12] = 54;

    AgeMin[13] = 55;
    AgeMax[13] = 59;

    AgeMin[14] = 60;
    AgeMax[14] = 64;

    AgeMin[15] = 65;
    AgeMax[15] = 90;

    return SumProbBirthAge,AgeMin,AgeMax
end
