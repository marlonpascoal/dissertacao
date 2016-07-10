function Info = parser_cont_dynamic_hinf_K(Ai,B1i,B2i,C1i,C2i,D11i,D12i,D21i,param)
% Given a stabilizing parameter-dependent state-feedback controller, the
% conditions for obtaining the robust output-feedback controller are stated,
% considering the H-infty criterium
%input:  Mi         -> vertices of the polytope (Mi = [M1 ... MN]) or
%                     Mi{1} = M1 ... Mi{N} = MN
%        param.precision -> precision of the solution
%        param.K    -> State-feedback gain
%        param.hinf -> Target H-infty norm, 0 if one wants the minimization
%        param.nc   -> Order of the controller
%        param.degP -> Degree of polynomial P
%        param.polyaP -> Degree of Polya relaxations on P
%        param.degF -> Degree of polynomial F
%        param.polyaF -> Degree of Polya relaxations on F
%        param.degV -> Degree of polynomial V
%        param.polyaV -> Degree of Polya relaxations on V
%        param.degH -> Degree of polynomial H
%        param.polyaH -> Degree of Polya relaxations on H
%output: Info.feas     -> feasible (1) or not (0)
%        Info.cpusec   -> cpu time to solve the LMIs (seconds)
%        Info.cpusec_m -> cpu time to mount the LMIs (seconds)
%        Info.Pi       -> Solution variables P
%        Info.Vi       -> Solution variables V
%        Info.Fi       -> Solution variables F
%        Info.Hi       -> Solution variables H
%        Info.Rm       -> Solution variable R
%        Info.Lm       -> Solution variable L
%        Info.Kout    -> Output feedback gain
%        Info.K        -> number of scalar variables
%        Info.L        -> number of LMI rows
%
%
% Date: 11/01/2010
% Author: agulhari@dt.fee.unicamp.br

if nargin == 9
    if isfield(param,'precision')
        precision = param.precision;
    else
        precision = -7;
    end
    
    if isfield(param,'nc')
        nc = param.nc;
    else
        nc = 0;
    end
    
    if isfield(param,'degP')
        degP = param.degP;
    else
        degP = 1;
    end

    if isfield(param,'polyaP')
        polyaP = param.polyaP;
    else
        polyaP = 0;
    end

    if isfield(param,'degF')
        degF = param.degF;
    else
        degF = 1;
    end

    if isfield(param,'polyaF')
        polyaF = param.polyaF;
    else
        polyaF = 0;
    end

    if isfield(param,'degV')
        degV = param.degV;
    else
        degV = 1;
    end

    if isfield(param,'polyaV')
        polyaV = param.polyaV;
    else
        polyaV = 0;
    end
    
    if isfield(param,'degH')
        degH = param.degH;
    else
        degH = 1;
    end

    if isfield(param,'polyaH')
        polyaH = param.polyaH;
    else
        polyaH = 0;
    end

    if isfield(param,'K')
        Kst = param.K;
    else
        Kst = [];
    end

    if isfield(param,'hinf')
        hinf = param.hinf;
    else
        hinf = 0;
    end
    
    if isfield(param,'xi')
        xi = param.xi;
    else
        xi = 1;
    end
    
    
    if isfield(param,'statedegP')
        paramstate.degW = param.statedegP;
    else
        paramstate.degW = 1;
    end
    
    if isfield(param,'statepolyaP')
        paramstate.polyaW = param.statepolyaP;
    else
        paramstate.polyaW = 0;
    end
    
    if isfield(param,'statedegF')
        paramstate.degF = param.statedegF;
    else
        paramstate.degF = 1;
    end
    
    if isfield(param,'statepolyaF')
        paramstate.polyaF = param.statepolyaF;
    else
        paramstate.polyaF = 0;
    end
    
    if isfield(param,'statedegZ')
        statedegZ = param.statedegZ;
    else
        statedegZ = 1;
    end
    
    if isfield(param,'statepolyaZ')
        statepolyaZ = param.statepolyaZ;
    else
        statepolyaZ = 0;
    end
    
    if isfield(param,'statedegG')
        statedegG = param.statedegG;
    else
        statedegG = 1;
    end
    
    if isfield(param,'statepolyaG')
        statepolyaG = param.statepolyaG;
    else
        statepolyaG = 0;
    end
else
    precision = -7;
    nc = 0;
    Kst = [];
    hinf = 0;
    degP = 1;
    polyaP = 0;
    degF = 1;
    polyaF = 0;
    degV = 1;
    polyaV = 0;
    degH = 1;
    polyaH = 0;
    xi = 1;
    
    paramstate.degW = 1;
    paramstate.polyaW = 0;
    paramstate.degF = 1;
    paramstate.polyaF = 0;
    statedegZ = 1;
    statepolyaZ = 0;
    statedegG = 1;
    statepolyaG = 0;
end

%Determina��o da ordem e do n�mero de v�rtices
if (iscell(Ai))
    order = length(Ai{1});
    vertices = length(Ai);
    inputs_w = size(B1i{1},2);
    inputs_u = size(B2i{1},2);
    outputs_w = size(C1i{1},1);
    outputs_u = size(C2i{1},1);
else
    [order, vertices] = size(Ai);
    vertices = vertices/order;
    inputs_w = size(B1i,2)/vertices;
    inputs_u = size(B2i,2)/vertices;
    outputs_w = size(C1i,1);
    outputs_u = size(C2i,1);
    
    
    for cont=1:vertices
        Aaux{cont} = Ai(:,(cont-1)*order+1:cont*order);
        B1aux{cont} = B1i(:,(cont-1)*inputs_w+1:cont*inputs_w);
        B2aux{cont} = B2i(:,(cont-1)*inputs_u+1:cont*inputs_u);
        C1aux{cont} = C1i(:,(cont-1)*order+1:cont*order);
        C2aux{cont} = C2i(:,(cont-1)*order+1:cont*order);
        D11aux{cont} = D11i(:,(cont-1)*inputs_w+1:cont*inputs_w);
        D12aux{cont} = D12i(:,(cont-1)*inputs_u+1:cont*inputs_u);
        D21aux{cont} = D21i(:,(cont-1)*inputs_w+1:cont*inputs_w);
    end
    clear Ai B1i B2i C1i C2i D11i D12i D21i;
    Ai = Aaux;
    B1i = B1aux;
    B2i = B2aux;
    C1i = C1aux;
    C2i = C2aux;
    D11i = D11aux;
    D12i = D12aux;
    D21i = D21aux;
    clear Aaux B1aux B2aux C1aux C2aux D11aux D12aux D21aux;
end

%Define the augmented system
%C3 = randn(outputs_u,nc);%%%%%%%%%%
%Info.C3 = C3;%%%%%%%%%
for cont=1:vertices
    Atil{cont} = [Ai{cont} zeros(order,nc); zeros(nc,order) zeros(nc,nc)];
    B1til{cont} = [B1i{cont}; zeros(nc,inputs_w)];
    B2til{cont} = [zeros(order,nc) B2i{cont}; eye(nc) zeros(nc,inputs_u)];
    C1til{cont} = [C1i{cont} zeros(outputs_w,nc)];
    %C2til{cont} = [zeros(nc,order) eye(nc); C2i{cont} -C3];
    C2til{cont} = [zeros(nc,order) eye(nc); C2i{cont} zeros(outputs_u,nc)];
    D11til{cont} = [D11i{cont}];
    D12til{cont} = [zeros(outputs_w,nc) D12i{cont}];
    D21til{cont} = [zeros(nc,inputs_w); D21i{cont}];
end

Info.cpusec_m = clock;

%sistema de LMIs
LMIs = set([]);
Info.L = 0;
Info.feas = 0;


% if (isempty(Kst))
%     cd ../Dynamic_Continuo_Estab/
%     %K(alpha) = Z*inv(G(alpha))
%     paramstate.degZ = 0;
%     paramstate.degG = statedegG;
%     paramstate.nc = nc;
%     contxi = 1;
%     InfoState.feas = 0;
%     degK = statedegG;
%     while ((contxi <= length(xi)) && (InfoState.feas == 0))
%         paramstate.xi = xi(contxi);
%         InfoState = parser_cont_dynamic_estado_stab_finsler(Ai,B2i,paramstate);
%         contxi = contxi + 1;
%     end
%     if (InfoState.feas == 0)
%         %K(alpha) = Z(alpha)*inv(G)
%         paramstate.degZ = statedegZ;
%         paramstate.degG = 0;
%         degK = statedegZ;
%         contxi = 1;
%         while ((contxi <= length(xi)) && (InfoState.feas == 0))
%             param.xi = xi(contxi);
%             InfoState = parser_cont_dynamic_estado_stab_finsler(Ai,B2i,paramstate);
%             contxi = contxi + 1;
%         end
%         InfoState = parser_cont_dynamic_estado_stab_finsler(Ai,B2i,paramstate);
%     end
%     cd ../Dynamic_Continuo_Hinf/
%     
%     if (InfoState.feas == 1)
%         polyKst = poly_struct(InfoState.Kst,'K',vertices,max([paramstate.degZ paramstate.degG]));
%     end
% else
    InfoState.feas = 1;
    if (iscell(Kst))
        numcoefsK = length(Kst);
    else
        numcoefsK = size(Kst,2)/(order+nc);
    end
    if (numcoefsK == vertices)
        degK = 1;
    elseif(numcoefsK == 1)
        degK = 0;
    else
        auxcoefsK = vertices;
        degK = 1;
        while (auxcoefsK < numcoefsK)
            degK = degK + 1;
            auxcoefsK = (factorial(vertices + degK - 1))/(factorial(degK)*factorial(vertices-1));
        end
    end
    
    polyKst = poly_struct(Kst,'K',vertices,degK);
% end

% fid = lmifiles('o','parser_cont_dynamic_hinf_K_deg2all');%%%%%%%%%

if (InfoState.feas==1)
    %Defining the decision variables
    numcoefsP = (factorial(vertices + degP - 1))/(factorial(degP)*factorial(vertices-1));
    numcoefsF = (factorial(vertices + degF - 1))/(factorial(degF)*factorial(vertices-1));
    numcoefsV = (factorial(vertices + degV - 1))/(factorial(degV)*factorial(vertices-1));
    numcoefsH = (factorial(vertices + degH - 1))/(factorial(degH)*factorial(vertices-1));
%     for cont=1:numcoefsP
%         P{cont} = sdpvar(order+nc,order+nc,'symmetric');
%     end
%     for cont=1:numcoefsF
%         F{cont} = sdpvar(order+nc,order+nc,'full');
%     end
%     for cont=1:numcoefsV
%         V{cont} = sdpvar(order+nc,order+nc,'full');
%     end
%     for cont=1:numcoefsH
%         H{cont} = sdpvar(outputs_w,outputs_w,'full');
%     end
    
    %R = [sdpvar(nc,nc,'full') zeros(nc,inputs_u); sdpvar(inputs_u,nc,'full') sdpvar(inputs_u,inputs_u,'full')];
    R = sdpvar(inputs_u+nc,inputs_u+nc,'full');
    %L = [zeros(nc,nc+outputs_u); zeros(inputs_u,outputs_u) sdpvar(inputs_u,outputs_u,'full')];
    L = sdpvar(inputs_u+nc,outputs_u+nc,'full');

    if (hinf==0)
        gama_quad = sdpvar(1,1);
        optimize = true;
    else
        gama_quad = hinf^2;
        optimize = false;
    end
    
    %Constructing the polynomials
    [vecpoly{1}] = poly_struct(Atil,'A',vertices,1);
    [vecpoly{2}] = poly_struct(B1til,'B1',vertices,1);
    [vecpoly{3}] = poly_struct(B2til,'B2',vertices,1);
    [vecpoly{4}] = poly_struct(C1til,'C1',vertices,1);
    [vecpoly{5}] = poly_struct(C2til,'C2',vertices,1);
    [vecpoly{6}] = poly_struct(D11til,'D11',vertices,1);
    [vecpoly{7}] = poly_struct(D12til,'D12',vertices,1);
    [vecpoly{19}] = poly_struct(D21til,'D21',vertices,1);
    [vecpoly{8}] = poly_struct(order+nc,order+nc,'P','symmetric',vertices,degP);
    [vecpoly{9}] = poly_struct(order+nc,order+nc,'F','full',vertices,degF);
    [vecpoly{10}] = poly_struct(order+nc,order+nc,'V','full',vertices,degV);
    [vecpoly{11}] = poly_struct(outputs_w,outputs_w,'H','full',vertices,degH);
    [vecpoly{12}] = poly_struct(R,'R');
    [vecpoly{13}] = poly_struct(L,'L');
    

    vecpoly{14} = polyKst;
    
    [vecpoly{15}] = poly_struct(zeros(order+nc,outputs_w),'0np');
    [vecpoly{16}] = poly_struct(eye(outputs_w),'Ip');
    [vecpoly{17}] = poly_struct(eye(inputs_w),'Ir');
    [vecpoly{18}] = poly_struct(gama_quad,'gama','scalar');

    %Polya relaxations
    poly1 = poly_struct(1,'1',vertices,0);
    for cont=1:polyaP
        vecpoly{8} = operation_poly(vecpoly{8},poly1,'*');
    end
    for cont=1:polyaF
        vecpoly{9} = operation_poly(vecpoly{9},poly1,'*');
    end
    for cont=1:polyaV
        vecpoly{10} = operation_poly(vecpoly{10},poly1,'*');
    end
    for cont=1:polyaH
        vecpoly{11} = operation_poly(vecpoly{11},poly1,'*');
    end
    
    
    
        
    tamlmis = 0;
    %Constructing the LMIs
    TermP{1,1} = parser_poly('P',vecpoly,vertices);
    
%     lmifiles('i',fid,TermP,'>');%%%%%%%%%
    
    LMIs = construct_lmi_terms(TermP,'>');
    Info.L = Info.L + (length(LMIs)-tamlmis)*(size(TermP{1,1}.data(1).value,1));
    tamlmis = length(LMIs);

    
%     T = parser_poly('F*A + F*B2*K',vecpoly,vertices);
%     T.label = 'T';
%     vecpoly{20} = T;
    Term{1,1} = parser_poly('F*A + F*B2*K + (F*A + F*B2*K)''',vecpoly,vertices);
    Term{1,2} = parser_poly('P - F + A''*V'' + K''*B2''*V''',vecpoly,vertices);
    Term{1,3} = parser_poly('F*B1',vecpoly,vertices);
    Term{1,4} = parser_poly('C1''*H + K''*D12''*H',vecpoly,vertices);
    Term{1,5} = parser_poly('F*B2 + C2''*L'' - K''*R''',vecpoly,vertices);
    Term{2,2} = parser_poly('-V-V''',vecpoly,vertices);
    Term{2,3} = parser_poly('V*B1',vecpoly,vertices);
    Term{2,4} = parser_poly('0np',vecpoly,vertices);
    Term{2,5} = parser_poly('V*B2',vecpoly,vertices);
    Term{3,3} = parser_poly('-gama*Ir',vecpoly,vertices);
    Term{3,4} = parser_poly('D11''*H',vecpoly,vertices);
    Term{3,5} = parser_poly('D21''*L''',vecpoly,vertices);
    Term{4,4} = parser_poly('Ip - H - H''',vecpoly,vertices);
    Term{4,5} = parser_poly('H''*D12',vecpoly,vertices);
    Term{5,5} = parser_poly('-R-R''',vecpoly,vertices);

%     lmifiles('i',fid,Term,'<');%%%%%%%%%
    
    LMIs = LMIs + construct_lmi_terms(Term,'<');
    soma = 0;
    for cont=1:5
        soma = soma + size(Term{cont,cont}.data(1).value,1);
    end
    Info.L = Info.L + (length(LMIs)-tamlmis)*soma;
    tamlmis = length(LMIs);
    
    
    
    Info.cpusec_m = etime(clock,Info.cpusec_m);

    Info.K = size(getvariables(LMIs),2);

    ops = sdpsettings;
    ops.verbose = 0;
    ops.solver = 'sedumi';
    
%     lmifiles('i',fid,'gama');%%%%%%%%%%%%%%%%
%     lmifiles('i',fid,'sdpsettings','verbose',0,'solver','sedumi');%%%%%%%

    
    if (optimize)
        sol = solvesdp(LMIs,gama_quad,ops);
    else
        sol = solvesdp(LMIs,[],ops);
    end
    Info.cpusec = sol.solvertime;
    [p,d]=checkset(LMIs);
%     p = -10;

%     lmifiles('c',fid);%%%%%%%
%     disp('Terminou');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Info.feas = 0;
    %Isola as solu��es, se existirem
    if ((sum(roundn(p,precision) < 0) == 0) && ~isnan(sum(p)))
        Info.Atil = Atil;
        Info.B1til = B1til;
        Info.B2til = B2til;
        Info.C1til = C1til;
        Info.C2til = C2til;
        Info.D11til = D11til;
        Info.D12til = D12til;
        Info.D21til = D21til;
        P = rolmip('getvar','P');
        F = rolmip('getvar','F');
        V = rolmip('getvar','V');
        H = rolmip('getvar','H');
        for cont=1:numcoefsP
            Info.Pi{cont} = double(P.data(cont).value);
        end
        for cont=1:numcoefsF
            Info.Fi{cont} = double(F.data(cont).value);
        end
        for cont=1:numcoefsV
            Info.Vi{cont} = double(V.data(cont).value);
        end
        for cont=1:numcoefsH
            Info.Hi{cont} = double(H.data(cont).value);
        end
        Info.Lm = double(L);
        Info.Rm = double(R);
        Info.Kout = inv(double(R))*double(L);
        Info.hinf = sqrt(double(gama_quad));
        Info.feas = 1;

    end
else
    Info.cpusec_m = etime(clock,Info.cpusec_m);
end
    
return