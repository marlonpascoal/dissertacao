load ga_precise_known_ncnx_part3;
load politopos_estabilizaveis_sof_comganho_vertices_Bprecknown;
load sistemas_selecionados_estatistica;

for (i = 41:50)%[5, 6, 7, 10, 13, 14, 15, 16, 21, 24, 26, 30, 32, 33, 34, 35, 36, 39, 40, 41, 43, 44, 45, 49])
    selecionado = selected{i};
    
    system = sys{selecionado(1),selecionado(2),selecionado(3),selecionado(4),selecionado(5)};
    
    nx = selecionado(1);
    nc = nx;
    
    
    if (~isempty(system))
        system.Ai = {system.Ai};
        system.B1i = {system.B1i};
        system.B2i = {system.B2i};
        system.C1i = {system.C1i};
        system.C2i = {system.C2i};
        system.D1i = {system.D1i};
        system.D2i = {system.D2i};
        system.Dyi = {system.Dyi};
        
%         res_ga = ga_uncertain(system,nc);
        res_ga = ga_uncertain_2(system,nc);
%         stats_hifoo.(['system_' num2str(i)]);

%         if (res_ga.gama < res_pso.gama)
%             Kdof = res_ga.Kdof;
%             ia = ('GA');
%         else
            ia = ('PSO');
%         end

		normas = struct();

		if (strcmp(res_ga.status,'OK') && res_ga.gama ~= Inf)
            Kdof = res_ga.Kdof;
            normas = calc_normas(system, Kdof, nc);
        end
    
        normas = setfield(normas, 'ia', ia);
        normas.pso_res = res_ga;
        normas_ia_v3_ncnx = setfield(normas_ia_v3_ncnx, ['system_' num2str(i)] ,normas);
        
        
        save ga_precise_known_ncnx_part3 normas_ia_v3_ncnx;
    end
end