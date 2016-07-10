function poles = generate_random_Ksf(popsize, system)
%(sizeKSF, popsize, Atil, B2til)
%     [nx,nw]=size(system.B1i{1}); [nz,nx]=size(system.C1i{1}); [nz,nu]=size(D2i{1}); 
%     ny = size(C2i{1},1);
    
%     F = sdpvar(nx + nc,nx + nc,'full');
%         Z = sdpvar(nu + nc,nx + nc,'full');

%     sizeKSF = [system.nu, system.nx];
    
    poles = cell(0,0);
%     ju = 1;
%     sizeKSF = size(Z * F);
    vertices = length(system.Ai);

    while (size(poles, 2) < popsize)
%         qtdrows = 100;       
        
        rrr = [];
        for(jx = 1:system.Ksf_cols)
%             rrr = [rrr complex(-abs(50*randn(qtdrows,1)),(50*randn(qtdrows,1)))];
            rrr = [rrr -abs(50*randn)];%complex(-abs(50*randn(qtdrows,1)),(50*randn(qtdrows,1)))];
        end

    
%         for (jx = 1:size(rrr,1))
            Ksf_ = -place(system.Ai{1},system.B2i{1},rrr);
            Ksf = cell(vertices,1);
            if (real(eig(system.Ai{1} + system.B2i{1} * Ksf_)) < 0)
%                 kkk(ju,:) = (acc);
                for (i = 1:vertices)
                    Ksf{i} = Ksf_;
                end
                
                poles{size(poles,2) + 1} = Ksf;
%                 ju = ju + 1;
            end
%         end
    end
    
    
