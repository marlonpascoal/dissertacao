function poles = generate_Ksf_wrapper(popsize, system)
    ratio = .5;
    
    poles1 = generate_stable_poles_h2(floor(popsize * ratio), system);
    poles2 = generate_random_Ksf(floor(popsize * ratio), system);
    
    poles = [poles1 poles2];
    
    if (length(poles) < popsize)
        poles_fill = generate_random_Ksf(floor(popsize - length(poles)), system); 
        
        poles = [poles poles_fill];
    end

end

