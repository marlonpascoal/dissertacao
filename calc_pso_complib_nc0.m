load pso_complib_nc0;
addpath('COMPlib_r1_0');

sistemas = obtem_labels_compleib();

for (idx = 11:12)
    name = cell2mat(sistemas(idx));
    
    [A,B1,B,C1,C,D11,D12,D21,nx,nw,nu,nz,ny] = COMPleib(name);
    
%     system = 'AC1';
    nc = 0;
%     tic;

    system = struct();
    system.Ai = {A};
        system.B1i = {B1};
        system.B2i = {B};
        system.C1i = {C1};
        system.C2i = {C};
        system.D1i = {D11};
        system.D2i = {D12};
        system.Dyi = {D21};

    res_pso = pso_uncertain_2(system,nc);
    
    ia = ('PSO');
%         end

    normas = struct();

    if (strcmp(res_pso.status,'OK'))
        Kdof = res_pso.Kdof;
        normas = calc_normas(system, Kdof, nc);
    end

    normas = setfield(normas, 'ia', ia);
    normas.pso_res = res_pso;
    pso_complib_nc0 = setfield(pso_complib_nc0, [name] ,normas);


    save pso_complib_nc0 pso_complib_nc0;
    
end