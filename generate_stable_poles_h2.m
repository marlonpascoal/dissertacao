function poles = generate_stable_poles_h2(popsize, system)
    poles = {};
%     ju = 1;
%     sizeKSF = size(Z * F);
    vertices = length(system.Ai);
    param = struct();
    param.degW = 0;
    
    Info = parser_cont_realim_estado_h2(system.Ai,system.B1i,system.B2i,system.C1i,system.D2i, param);
    h2_ = Info.h2;
    
    Ksf_h2 = Info.Kst;
%     Ksf_h2 = cell(vertices,1);
%     for (k_i = 1:vertices)
%         Ksf_h2{k_i} = Info.Zi{k_i} * inv(Info.Gi{k_i});
%     end
    
%     poles{1} = Ksf_h2;
%     poles{1} = stats_hifoo.system_2.K.d * system.C2i{vertice};

%     for(i=2:popsize)
    count = 0;
    while (size(poles,2) < popsize)
        h2_temp = h2_ + 50 * rand;
        
        param = struct();
        param.h2 = h2_temp;
        param.degW = 0;
        
        Info = parser_cont_realim_estado_h2(system.Ai,system.B1i,system.B2i,system.C1i,system.D2i, param);
%         Ksf = Info.Zi{1} * inv(Info.Gi{1});

%         Ksf = cell(vertices,1);
        tempKsf = cell(vertices,1);
        Ksf_ = [];
        
%         for (k_i = 1:vertices)
%             Ksf{k_i} = Info.Zi{k_i} * inv(Info.Gi{1});
%         end
        Ksf = Info.Kst;
        
        size_ = size(Ksf{1});
        Ksf_ = do_reshape(Ksf, size_(1), size_(2), vertices);
        
        
        
        
%         tempKsf = reshape(Ksf, 1, size_(1) * size_(2));
        res_lmi = second_stage_agulhari_optimized(system, Ksf_, 0);
        
        if (res_lmi.gamas ~= Inf)
%         if (real(eig(system.Ai + system.B2i * Ksf)) < 0)
            poles{size(poles,2) + 1} =  Ksf;
            
            disp(sprintf('Generating Ksfs...%d%%',int8(100*(size(poles,2)/popsize))));
        else
            disp('Generated Ksf not feasible on second stage...');
        end
        
        count = count + 1;
        if (size(poles,2) == 0 && count > 50)
            break;
        end
        
    end
    
    