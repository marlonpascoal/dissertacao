% normas_ia_uncertain = struct();
load normas_ia_uncertain_ga_ncnx_part4;
load selected_uncertain;
load politopos_estabilizaveis_sof_comganho_vertices_Bprecknown_teste;

for (i = 16:20)%[5, 6, 7, 10, 13, 14, 15, 16, 21, 24, 26, 30, 32, 33, 34, 35, 36, 39, 40, 41, 43, 44, 45, 49])
    selecionado = selected_uncertain{i};
    
    system = sys{selecionado(1),selecionado(2),selecionado(3),selecionado(4),selecionado(5)};
    
    nx = selecionado(1);
    nc = nx;
    vertices = selecionado(4);
    system_temp = system;
    
    if (~isempty(system))
        system.Ai = cell(vertices,1);
        system.B1i = cell(vertices,1);
        system.B2i = cell(vertices,1);
        system.C1i = cell(vertices,1);
        system.C2i = cell(vertices,1);
        system.D1i = cell(vertices,1);
        system.D2i = cell(vertices,1);
        system.Dyi = cell(vertices,1);
        
%         for (i = 0:(vertices - 1)) 
%             start = ((l * i) + 1);
%     %             param.K{i+1} = ;
%             out{i+1} = reshape(row(start:(start+(l - 1))), nrows, ncols);
%         end
        
%                     B1i = rand(nx, nw);
%                     C1i = rand(nz, nx);
%                     D1i = rand(nz, nw);
%                     D2i = rand(nz, nu);
%                     Dyi

        tamanhos = struct();
        tamanhos.Ai = (size(system_temp.Ai,2) / vertices);
%         tamanhos.B1i = (size(system_temp.B1i,2) / vertices);
        tamanhos.B2i = (size(system_temp.B2i,2) / vertices);
%         tamanhos.C1i = (size(system_temp.C1i,2) / vertices);
        tamanhos.C2i = (size(system_temp.C2i,2) / vertices);
%         tamanhos.D1i = (size(system_temp.D1i,2) / vertices);
%         tamanhos.D2i = (size(system_temp.D2i,2) / vertices);
%         tamanhos.Dyi = (size(system_temp.Dyi,2) / vertices);
        for (o = 0:(vertices - 1))
            idx = o + 1;
            system.Ai{idx} = system_temp.Ai(:,(tamanhos.Ai * o + 1):((tamanhos.Ai * o + 1)+(tamanhos.Ai - 1)));
            system.B1i{idx} = system_temp.B1i;
%             system.B1i{idx} = system_temp.B1i(:,(tamanhos.B1i * o + 1):((tamanhos.B1i * o + 1)+(tamanhos.B1i - 1)));
            system.B2i{idx} = system_temp.B2i(:,(tamanhos.B2i * o + 1):((tamanhos.B2i * o + 1)+(tamanhos.B2i - 1)));
            system.C1i{idx} = system_temp.C1i;
%             system.C1i{idx} = system_temp.C1i(:,(tamanhos.C1i * o + 1):((tamanhos.C1i * o + 1)+(tamanhos.C1i - 1)));
            system.C2i{idx} = system_temp.C2i(:,(tamanhos.C2i * o + 1):((tamanhos.C2i * o + 1)+(tamanhos.C2i - 1)));
            system.D1i{idx} = system_temp.D1i;
%             system.D1i{idx} = system_temp.D1i(:,(tamanhos.D1i * o + 1):((tamanhos.D1i * o + 1)+(tamanhos.D1i - 1)));
            system.D2i{idx} = system_temp.D2i;
%             system.D2i{idx} = system_temp.D2i(:,(tamanhos.D2i * o + 1):((tamanhos.D2i * o + 1)+(tamanhos.D2i - 1)));
            system.Dyi{idx} = system_temp.Dyi;
%             system.Dyi{idx} = system_temp.Dyi(:,(tamanhos.Dyi * o + 1):((tamanhos.Dyi * o + 1)+(tamanhos.Dyi - 1)));
        end
        
        
%         res_ga = ga_uncertain(system,nc);
        res_ga = ga_uncertain_2(system,nc);
%         stats_hifoo.(['system_' num2str(i)]);

%         if (res_ga.gama < res_pso.gama)
%             Kdof = res_ga.Kdof;
%             ia = ('GA');
%         else
            
            ia = ('GA');
%         end
        
        normas = struct();
        
        if (strcmp(res_ga.status,'OK') && res_ga.gama ~= Inf)
            Kdof = res_ga.Kdof;
            normas = calc_normas(system, Kdof, nc);
        end
        normas = setfield(normas, 'ia', ia);
        normas.res_ga = res_ga;
        normas_ia_uncertain_ga_ncnx = setfield(normas_ia_uncertain_ga_ncnx, ['system_' num2str(i)] ,normas);
        
        
        save normas_ia_uncertain_ga_ncnx_part4 normas_ia_uncertain_ga_ncnx;
    end
end