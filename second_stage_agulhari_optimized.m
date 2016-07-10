function ret = second_stage_agulhari_optimized(system, Ksfs, nc)
%     system, Ksfs, nc
    
    
    warning('off','YALMIP:strict');
    gamas = [];
    Kdofs = {};
    
    [nx,nw]=size(system.B1i{1}); 
    [nz,nx]=size(system.C1i{1}); 
    [nz,nu]=size(system.D2i{1}); 
    ny = size(system.C2i{1},1);
    
    total = size(Ksfs, 1);
    for (qtd = 1:total) 
        param = struct();
        
        Ksf_ = Ksfs(qtd,:);
        vertices = length(system.B1i);
        l = length(Ksf_) / vertices;
        
        param.K = {};

        hash_system = genvarname(DataHash(system));
        hash_args = genvarname(DataHash({Ksf_, nc}));
        
        cache_fname = ['cache/' hash_system '.mat'];
        
        if (exist('secondstage_cache') == 0)
            if (exist(cache_fname, 'file') == 2)
                load(cache_fname);
            else
                secondstage_cache = struct();
            end
            
        end
        
        if (isfield(secondstage_cache, hash_args))
            gamas(qtd,:) = secondstage_cache.(hash_args).gama;
            Kdofs{qtd} = secondstage_cache.(hash_args).Kdof;
            
            disp('fitness function results loaded from cache...');
        else
%             for (i = 0:(vertices - 1)) 
%                 start = ((l * i) + 1);
    %             param.K{i+1} = ;
%                 param.K{i+1} = undo_reshape(Ksf_(start:(start+(l - 1))), nu + nc, nx + nc, vertices);
%             end

            param.K = undo_reshape(Ksf_, system.Ksf_rows, system.Ksf_cols, vertices);

    %         Ksf = reshape(Ksfs(qtd,:), nu + nc, nx + nc);


            param.nc = 0; %nc is set to 0 because the system is already augmented.

            Info = parser_cont_dynamic_hinf_K(system.Ai,system.B1i,system.B2i,system.C1i,system.C2i,system.D1i,system.D2i,system.Dyi,param);

            if (Info.feas == 1)
                gamas(qtd,:) = Info.hinf;
                Kdofs{qtd} = Info.Kout;
            else
                gamas(qtd,:) = Inf;
                Kdofs{qtd} = [];
            end

            secondstage_cache = setfield(secondstage_cache, (hash_args), struct('gama', gamas(qtd,:), 'Kdof', Kdofs{qtd}));
            save(cache_fname,'secondstage_cache');
        end
        
        if (total > 1)
            disp(sprintf('performing fitness function... %d%%',int8(100 * (qtd/total))));
        end
    end
    
    ret = struct('gamas',gamas, 'Kdofs', {Kdofs});
    