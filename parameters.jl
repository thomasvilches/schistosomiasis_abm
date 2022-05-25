@enum HEALTH SUSC=0 INF=2 LAT=1
#@enum SNAIL_HEALTH SUSC=1 LAT=2 INF=3

@with_kw mutable struct SCHparameters @deftype Int64

    file_index = 1
    init_n_worms = 2
    age_interval = 5
    grid_size_human = 500
    grid_size_snail = 1000
    years_sim = 60
    sim_time = years_sim*365
    method = 2 #1 - Kato-Katz, 2 - HTX
    worms_lim_diag = 1 ###Number of couples
    #lambda_s::Float64 = 0.097
    #lambda_c::Float64 = 0.047
    prob_binomial::Float64 = 0.5
    max_worms::Int64 = 1
    tau::Float64 = 0.083
    latent_period::Float64 = tau*365
    mu_2line_s::Float64 = 21.6
    mu_line_s::Float64 = 10.8
    #z_star::Float64 = exp(-mu_line_s*tau)/(mu_2line_s/mu_line_s-(mu_2line_s/mu_line_s-1)*exp(-mu_line_s*tau))
    max_trials = 10
    max_age = 100
    mu_h::Float64 = 0.0137

    oviposition_rate::Float64 = 100*365 #annual 
    shedding_rate::Float64 = 58400
    miracidium_mortality::Float64 = (1/7)*24*365 ##all in years
    cercaria_mortality::Float64 = 365
    
    development_probability::Float64 = 0.83

    lambda_::Float64 =  0.120 #0.120 #0.221
    lambda_line::Float64 = 0.047 #0.047 #0.042

    infection_human::Float64 = 0.1 #estimar?? talvez de para tirar do artigo
    infection_snail::Float64 = 0.4 #estimar
    mu_w::Float64 = 0.09 #0.090 #0.079 ##Estimar
    immunity_parameter::Float64 = 8.07 #8.07 #8.96 ##Estimar
    beta::Float64 = 1.0 #lambda_line/lambda_ #0.3917 #0.190045 #estimar (do artigo?)
   
    immunity_period::Float64 = immunity_parameter*365
    mutation_prob::Float64 = 0.04

    a::Float64 = 3.4e-7
    b::Float64 = 1.2239e-2
    #T1::Float64 = 0.0719 #0.0719  #0.0537
    #T2::Float64 = 69.77 #114.12  #235.09       #69.77    #144.53
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
