function out = ga_uncertain_2(system, nc)
    tic;
    out = struct();

    Ai = system.Ai;
    B1i = system.B1i;
    B2i = system.B2i;
    C1i = system.C1i;
    C2i = system.C2i;
    D1i = system.D1i;
    D2i = system.D2i;
    Dyi = system.Dyi;
    
    %_______________________________________________________
    % II Stopping criteria
    maxit=100; % max number of iterations
    mincost=-9999999; % minimum cost
    %_______________________________________________________
    % III GA parameters
    popsize=24; % set population size
    mutrate=.4; % set mutation rate
    selection=0.5; % fraction of population kept
	tolerance = 1e-9;
    use_local_optimization = true;
    
    population = [];
    
    [nx,nw]=size(B1i{1}); [nz,nx]=size(C1i{1}); [nz,nu]=size(D2i{1}); 
    ny = size(C2i{1},1);
    
    vertices = length(Ai);
    
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
	
	stable_poles = generate_Ksf_wrapper(popsize, system);
	
	for(kidx = 1: popsize)
        population(kidx,:) = do_reshape(stable_poles{kidx}, system.Ksf_rows, system.Ksf_cols, vertices);
    end
    
    
    ff='second_stage_agulhari_optimized'; % objective function
    npar= length(population(1,:)); % number of optimization variables
    Nt=npar; % continuous parameter GA Nt=#variables
    
    keep=floor(selection*popsize); % #population
    % members that survive
    nmut=ceil((popsize-1)*Nt*mutrate); % total number of
    % mutations
    M=ceil((popsize-keep)/2); % number of matings

    iga=0; % generation counter initialized

    % par=(varhi-varlo)*rand(popsize,npar)+varlo; % random
    ret = feval(ff,system, population, nc);
    cost= ret.gamas; % calculates population cost

    hx = {};
%     hv = {};
    for(m=1:length(population(1,:)))
        hx{m} = ['x',num2str(m) ];
%         hv{m} = ['v',num2str(m) ];
    end
    
%      [population cost]
     header = [hx 'cost'];
     disp(dataset({[population cost], header{:}}))
     
     rangeref = 1;
     
     minrange = -rangeref;
     maxrange = rangeref;
     
    % using ff
    [cost,ind]=sort(cost); % min cost in element 1
    population=population(ind,:); % sort continuous
    minc(1)=min(cost); % minc contains min of
    meanc(1)=mean(cost); % meanc contains mean of population

    M=ceil((popsize-keep)/2); % number of matings
    
    
    means___ = [];
    costs = [];
    
   
    
    %_______________________________________________________
    % Iterate through generations
    while (true) %iga<maxit
        
        
        iga=iga+1; % increments generation counter
        %_______________________________________________________
        
        means___(iga,:) = median(population,1);
        costs(iga) = cost(1);
        
       
%         M = ceil(1 / (1 + exp_help)) * keep
        
        % Pair and mate

        prob=flipud([1:keep]'/sum([1:keep])); % weights
        % chromosomes
        odds=[0 cumsum(prob(1:keep))']; % probability
        % distribution
        % function
        pick1=rand(1,M); % mate #1
        pick2=rand(1,M); % mate #2

        % ma and pa contain the indicies of the chromosomes
        % that will mate

        ic=1;
        while ic<=M
            for id=2:keep+1
                if pick1(ic)<=odds(id) & pick1(ic)>odds(id-1)
                    ma(ic)=id-1;
                end
                if pick2(ic)<=odds(id) & pick2(ic)>odds(id-1)
                    pa(ic)=id-1;
                end
            end
            ic=ic+1;
        end

        %_______________________________________________________
        % Performs mating using single point crossover
        ix=1:2:keep; % index of mate #1


        for ic=1:M
            beta_ = rand;
            population(keep+ix(ic),:)= beta_ * population(ma(ic),:) + (1 - beta_) * population(pa(ic),:); % 1st offspring
            population(keep+ix(ic)+1,:)= beta_ * population(pa(ic),:) + (1 - beta_) * population(ma(ic),:); % 1st offspring
        end

        %_______________________________________________________
        % Mutate the population
        mrow=sort(ceil(rand(1,nmut)*(popsize-1))+1);
        mcol=ceil(rand(1,nmut)*Nt);

        for ii=1:nmut
            increment = 1 + ( (maxrange-minrange)*rand+minrange);
            population(mrow(ii),mcol(ii))= population(mrow(ii),mcol(ii)) * increment;
            % mutation
        end % ii
        %_______________________________________________________
        % The new offspring and mutated chromosomes are
        % evaluated

        ret_full = feval(ff,system, population, nc);
        cost= ret_full.gamas; % calculates population cost
		
		for(l = 1:popsize)
			if (use_local_optimization)
				ret = struct('gamas', ret_full.gamas(l,:), 'Kdofs', ret_full.Kdofs{l});
			
                %LOCAL OPTIMIZATION
                lastGama = ret.gamas;
                klocal = 2;
                klocal_max = 4;
                newGama = Inf;
                newGamas = (Inf * ones(klocal_max,1));
                newKdofs = cell(klocal_max,1);
                newKsfs = cell(klocal_max,1);

                newKsfs{1} = undo_reshape(population(l,:), system.Ksf_rows, system.Ksf_cols, vertices);
%                 newKsfs{1} = ind(l,:);
                newGamas(1) = lastGama;
                newKdofs{1} = ret.Kdofs;

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
                    temp_res = feval(ff, system, newInd, nc);

    %                 lastGama = newGama;

                    newGamas(klocal) = temp_res.gamas;
                    newKdofs{klocal} = temp_res.Kdofs{1};

                    klocal = klocal + 1;



                end

                [gama_ord,index] = sort(newGamas);
                population(l,:) = do_reshape(newKsfs{index(1)}, system.Ksf_rows, system.Ksf_cols, vertices);

                if (index(1) ~= 1)
                    disp('local optimization performed...');
                end
                %FIM LOCAL OPTIMIZATION
                cost(l) = gama_ord(1);
                ret_full.Kdofs{l} = newKdofs{index(1)};
                
%            else
 %               valF(l,1) = ret.gamas;
            end
		end

        hx = {};
        for(m=1:length(population(1,:)))
            hx{m} = ['x',num2str(m) ];
        end
        
        header = [hx 'cost'];
        disp(dataset({[population cost], header{:}}))

        %_______________________________________________________
        % Sort the costs and associated parameters
        [cost,ind]=sort(cost);
        population=population(ind,:);
        %_______________________________________________________
        % Do statistics for a single nonaveraging run
        minc(iga+1)=min(cost);
        meanc(iga+1)=mean(cost);

%         header2 = {'iterations','iteration best cost'};
%         disp(dataset({[iga cost(1)], header2{:}}))
        
        
         % adaptative mutation control
        DL = cost(1) / mean(costs);
        exp_help = exp(1-DL);
        mutrate= exp_help / (1 + exp_help);
        
		costs(iga) = cost(1);
         
        start_variance = max([1, length(costs) - 30]);
        
        cost_variance = var(costs(start_variance:end));
        
        tout = toc;
        header2 = {'iterations','iteration best cost', 'stopping Criteria Value', 'Stopping Criteria', 'RunningTime(min)'};
        disp(dataset({[iga cost(1) cost_variance, tolerance (tout/60)], header2{:}}))
        
        
        %_______________________________________________________
        % Stopping criteria
        
		if (iga > 30 && (sum(costs(start_variance:end) == Inf) >= 30 || cost_variance < tolerance || (tout/60) > 3000))
%        if (iga > 30 && (cost_variance == NaN || cost_variance < 1e-9))
%         if iga>maxit | cost(1)<mincost
            break
        end
    end %iga
    
%     Ksf = par(1,:);
%     gama = cost(1);
%     Kdof = Kdofs(1);
    
    Kdofs_ = (ret_full.Kdofs);
%     Kdofs_ = (ret.Kdofs);
    out = struct('Ksf', population(1,:), 'gama', cost(1), 'Kdof', Kdofs_{1});
    out.status = 'OK';
    out.iterations = iga;
