function out =pso_uncertain_2(system, nc)
    tic;
    Ai = system.Ai;
    B1i = system.B1i;
    B2i = system.B2i;
    C1i = system.C1i;
    C2i = system.C2i;
    D1i = system.D1i;
    D2i = system.D2i;
    Dyi = system.Dyi;
    
    [nx,nw]=size(B1i{1}); [nz,nx]=size(C1i{1}); [nz,nu]=size(D2i{1}); 
    ny = size(C2i{1},1);
    vertices = length(Ai);
    popsize = 20;
    
    population = [];
    
    augSystem = {};
    augSystem.Ai = cell(vertices,1);
    augSystem.B1i = cell(vertices,1);
    augSystem.B2i = cell(vertices,1);
    augSystem.C1i = cell(vertices,1);
    augSystem.C2i = cell(vertices,1);
    augSystem.D1i = cell(vertices,1);
    augSystem.D2i = cell(vertices,1);
    augSystem.Dyi = cell(vertices,1);
    
    for(vertice = 1:vertices)
        augSystem.Ai{vertice} = [Ai{vertice} zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)];
        augSystem.B1i{vertice} = [B1i{vertice}; zeros(nc,nw)];
        augSystem.B2i{vertice} = [zeros(nx,nc) B2i{vertice}; eye(nc) zeros(nc,nu)];
        augSystem.C1i{vertice} = [C1i{vertice} zeros(nz,nc)];
        augSystem.C2i{vertice} = [zeros(nc,nx) eye(nc); C2i{vertice} zeros(ny,nc)];
        augSystem.D1i{vertice} = D1i{vertice};
        augSystem.D2i{vertice} = [zeros(nz,nc) D2i{vertice}];
        augSystem.Dyi{vertice} = [zeros(nc,nw); Dyi{vertice}];
    end
    
    [nx,nw]=size(augSystem.B1i{1}); [nz,nx]=size(augSystem.C1i{1}); [nz,nu]=size(augSystem.D2i{1}); 
    ny = size(augSystem.C2i{1},1);
    
    augSystem.nx = nx;
    augSystem.nw = nw;
    augSystem.nz = nz;
    augSystem.nu = nu;
    augSystem.ny = ny;
    augSystem.Ksf_rows = nu;
    augSystem.Ksf_cols = nx;
    
    
    system = augSystem;
    
%     for (vertice = 1:vertices)
%         clear *til P F Z temp_pop;
%         
%         ctr = (ctrb(Ai{vertice},B2i{vertice}));
%     
%         if rank(ctr) ~=  min(size(ctr))
%             out.status = [out.status '' disp('Vertice %d is Uncontrollable.', d)];
%             return;
%         end
%         
%         Atil = [Ai{vertice} zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)];
%         B1til = [B1i{vertice}; zeros(nc,nw)];
%         B2til = [zeros(nx,nc) B2i{vertice}; eye(nc) zeros(nc,nu)];
%         C1til = [C1i{vertice} zeros(nz,nc)];
%         C2til = [zeros(nc,nx) eye(nc); C2i{vertice} zeros(ny,nc)];
%         D1til = D1i{vertice};
%         D2til = [zeros(nz,nc) D2i{vertice}];
%         Dytil = [zeros(nc,nw); Dyi{vertice}];
%         
%         P = sdpvar(nx + nc,nx + nc,'sym');
%         F = sdpvar(nx + nc,nx + nc,'full');
%         Z = sdpvar(nu + nc,nx + nc,'full');
%         
%         popsize = 20;
%         sizeKSF = size(Z * F);
% 
% %         system_vertice = struct();
% %         system_vertice.Ai = Atil;
% %         system_vertice.B1i = B1til;
% %         system_vertice.B2i = B2til;
% %         system_vertice.C1i = C1til;
% %         system_vertice.C2i = C2til;
% %         system_vertice.D1i = D1til;
% %         system_vertice.D2i = D2til;
% %         system_vertice.Dyi = Dytil;
%         
% %         stable_poles = generate_stable_poles_h2(popsize, system, vertice);
%         stable_poles = generate_stable_poles(sizeKSF, popsize, Atil, B2til);
%         temp_pop = [];
%         for(kidx = 1: popsize)
%             
%             row_ = stable_poles{kidx};
%             size_ = size(row_);
%             temp_pop(kidx,:) = reshape(row_, 1, size_(1) * size_(2));
%         end
%         
%         population = [population temp_pop];
%     end


    stable_poles = generate_Ksf_wrapper(popsize, system);
%     stable_poles = generate_random_Ksf(popsize, system);
%     stable_poles = generate_stable_poles_h2(popsize, system);
%     size_ = size(stable_poles{1});
    
    for(kidx = 1: popsize)
        population(kidx,:) = do_reshape(stable_poles{kidx}, system.Ksf_rows, system.Ksf_cols, vertices);
    end
    
    history = zeros(1000);
    
    func = 'second_stage_agulhari_optimized';
    numInd = popsize;
    range = [1 3];

    tolerance = 1e-9;
    numIter = 1000; 
    pesoStoc = 0.5;
    c1 = .7;
    cmax = 1.43;

    range_min=range(1); % Range for initial swarm's elements
    range_max=range(2);
    fOb=func; % Objective function
    iter=0; % Number of iteration

%     disp('populacao inicial')

    ind = population;

    n_var = length(population(1,:));
    numVar=n_var; % Number of variables

    k = pesoStoc; % weight of stocastic element

    v = zeros(numInd,numVar); % Vector of swarm's velocity
    pd = ind;
    pd_costs = ones(numInd,1) * Inf;
    gd = zeros(1,numVar);
    gd_cost = Inf;
    ret_best = struct();

    radius=1000; % Initial radius for stop criterion
    costs = [];
    use_local_optimization = true;
    use_refresh_operator = true;

    while (true) %iter<numIter && radius>tolerance
        retLMI = {};
        for l=1:numInd
%             param = struct();
%             param.nc = nc;
%             param.K = ind(l,:);
            
%             ret = feval(fOb, Ai,B1i,B2i,C1i,C2i,D1i,D2i,Dyi, param)
            ret = feval(fOb, system, ind(l,:), nc);
            retLMI{l} = ret;
            
            if (use_local_optimization)
                %LOCAL OPTIMIZATION
                lastGama = ret.gamas;
                klocal = 2;
                klocal_max = 4;
                newGama = Inf;
                newGamas = (Inf * ones(klocal_max,1));
                newKdofs = cell(klocal_max,1);
                newKsfs = cell(klocal_max,1);

                newKsfs{1} = undo_reshape(ind(l,:), system.Ksf_rows, system.Ksf_cols, vertices);
%                 newKsfs{1} = ind(l,:);
                newGamas(1) = lastGama;
                newKdofs{1} = retLMI{l}.Kdofs{1};

    %             Kdof = ;
                while (klocal <= klocal_max && (isnan(abs(newGamas(klocal) - newGamas(klocal - 1))/newGamas(klocal)) || abs(newGamas(klocal) - newGamas(klocal - 1))/newGamas(klocal) > 1e-3))

                    if (isempty(newKdofs{klocal - 1}))
                        break
                    end
                    
                    tempKsf = cell(vertices,1);
                    for(i = 1:vertices)
                        tempKdof = newKdofs{klocal - 1};
                        tempKsf{i} = tempKdof * system.C2i{i};
                    end
                    
                    newKsfs{klocal} = tempKsf;

                    newInd = do_reshape(newKsfs{klocal}, system.Ksf_rows, system.Ksf_cols, vertices);
                    temp_res = feval(fOb, system, newInd, nc);

    %                 lastGama = newGama;

                    newGamas(klocal) = temp_res.gamas;
                    newKdofs{klocal} = temp_res.Kdofs{1};

                    klocal = klocal + 1;



                end

                [gama_ord,index] = sort(newGamas);
                ind(l,:) = do_reshape(newKsfs{index(1)}, system.Ksf_rows, system.Ksf_cols, vertices);

                if (index(1) ~= 1)
                    disp('local optimization performed...');
                end
                %FIM LOCAL OPTIMIZATION
                valF(l,1) = gama_ord(1);
            else
                valF(l,1) = ret.gamas;
            end
            
            
            

            if (valF(l,1) < pd_costs(l))
                pd(l,:) = ind(l,:);
                pd_costs(l) = valF(l,1);
                ret_best = ret;
            end
            
            disp(sprintf('performing fitness function... %d%%',100 * (l/popsize)));
        end
        
        [valF_ord,index]=sort(valF); % Sort the objective function's values for the swarm and identify the leader
        leader=ind(index(1),:);

        if (valF_ord(1) < gd_cost)
            gd_cost = valF_ord(1);
            gd = leader;
        end
        
        alpha_ = 1e-2;
        wi = 0.5 * (1 + tanh(alpha_ * pd_costs));
        wi = (kron(wi,ones(1,numVar)));
        v = wi .* v + kron(cmax*rand(numInd,1),ones(1,numVar)) .* (pd-ind) + kron(cmax*rand(numInd,1),ones(1,numVar)) .* (kron(gd, ones(numInd,1))-ind);
        ind = ind + v;
        
        radius = var(pd_costs);
        
        % REFRESH OPERATOR
        if (use_refresh_operator)
            should_reset = sum(abs(v) < 1e-1,2) == n_var;
            new_inds = zeros(popsize,n_var);
            for(reset_idx = 1:length(should_reset))
                if (should_reset(reset_idx) == 1)
                    disp('performing refresh operator...');
    %                 newInd = generate_stable_poles(sizeKSF, 1, system, vertice);

                    stable_poles = generate_random_Ksf(1, system);
%                     size_ = size(stable_poles{1});

    %                 for(kidx = 1: popsize)
                        new_inds(reset_idx,:) = do_reshape(stable_poles{1}, system.Ksf_rows, system.Ksf_cols, vertices);
    %                 end

    %                 new_inds(reset_idx,:) = do_reshape(cell2mat(newInd), size_(1), size_(2), vertices);
                end
            end

            ind(should_reset,:) = new_inds(should_reset,:);
            pd_costs(should_reset,:) = pd_costs(should_reset,:) * Inf ; % reset o p_d cost dos indiv�duos
            valF(should_reset,1) = valF(should_reset,1) * Inf;
        end
        % FIM REFRESH OPERATOR

        
        
        
%         radius=norm(leader-ind(index(end),:)); % Calculates the new radius
        
        iter=iter+1; % Increases the number of iteration

        hx = {};
        hv = {};
        for(m=1:length(ind(1,:)))
            hx{m} = ['x',num2str(m) ];
            hv{m} = ['v',num2str(m) ];
        end
        
        costs(iter) = gd_cost;
        start_variance = max([1, length(costs) - 30]);
        
        cost_variance = var(costs(start_variance:end));

    %     header = {hx, 'cost', 'best particle cost', hv};
        header = [hx 'cost', 'best particle cost' hv];
        disp(dataset({[ind valF pd_costs v], header{:}}))

    %     [pd pd_costs]
    %     [ind valF]
        tout = toc;
        header2 = {'iterations','iteration best cost', 'global best cost', 'stopping Criteria Value', 'Stopping Criteria', 'RunningTime(min)'};
        disp(dataset({[iter valF_ord(1) gd_cost, cost_variance, tolerance (tout/60)], header2{:}}))
        
        
        
        if (iter > 30 && (sum(costs(start_variance:end) == Inf) >= 30 || cost_variance < tolerance || (tout/60) > 3000) )
%         if iga>maxit | cost(1)<mincost
            break
        end
        
%         history(iter) = valF_ord(1);
        
%         figure(1)
%         plot(1:iter, history(1:iter));
    end

    if (sum(costs(start_variance:end) == Inf) >= 30 && ~isfield(ret_best, 'Kdofs'))
        out.status = 'Not feasible';
        out.iterations = iter;
        out.total_time = (toc/60);
    else
        out = struct('Ksf', gd, 'gama', gd_cost, 'Kdof', ret_best.Kdofs);
        out.status = 'OK';
        out.iterations = iter;
        out.total_time = (toc/60);
    end
