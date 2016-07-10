function Info = parser_cont_realim_estado_h2(Ai,B1i,B2i,Ci,Di,param)
% Calculates a stabilizing parameter-dependent state-feedback gain.
%input:  Ai       -> vertices of the polytope (Ai = [A1 A2 ... AN])
%        Bi       -> vertices of the polytope (Bi = [B1 B2 ... BN])
%        Ci       -> vertices of the polytope (Ci = [C1 C2 ... CN])
%        param.degW -> Degree of polynomial W
%        param.polyaW -> Degree of Polya relaxations on W
%        param.degZ -> Degree of polynomial Z
%        param.polyaZ -> Degree of Polya relaxations on Z
%        param.precision -> Precision for the tolerance on the primal
%                           value.
%output: Info.feas     -> feasible (1) or not (0)
%        Info.cpusec   -> cpu time to solve the LMIs (seconds)
%        Info.cpusec_m -> cpu time to mount the LMIs (seconds)
%        Info.Wi       -> Solution variables W
%        Info.Zi       -> Solution variables Z
%        Info.Gi       -> Equal to Info.Wi for compatibility reasons
%        Info.K        -> number of scalar variables
%        Info.L        -> number of LMI rows
%
%example: 2 states and 2 vertices
% A1 = [1.0  -0.5;
%       -0.3  -0.1 ];
% A2 = [-0.6  0.6;
%       0.7    -0.3];
% B1 = [1.3; -2.6];
% B2 = [1.3; -2.6];
% C1 = eye(2);
% C2 = eye(2);
% Ai = [A1 A2];
% Bi = [B1 B2];
% Ci = [C1 C2];
% Info = parser_cont_realim_estado_h2(Ai,Bi,Ci,param);
%
% Date: 18/01/2010
% Author: agulhari@dt.fee.unicamp.br

if nargin == 6    
    if isfield(param,'degW')
        degW = param.degW;
    else
        degW = 1;
    end
    
    if isfield(param,'polyaW')
        polyaW = param.polyaW;
    else
        polyaW = 0;
    end

    if isfield(param,'degZ')
        degZ = param.degZ;
    else
        degZ = 1;
    end
    
    if isfield(param,'polyaZ')
        polyaZ = param.polyaZ;
    else
        polyaZ = 0;
    end

    if isfield(param,'h2')
        h2 = param.h2;
    else
        h2 = 0;
    end
    
    if isfield(param,'precision')
        precision = param.precision;
    else
        precision = -7;
    end
else
    degW = 1;
    polyaW = 0;
    degZ = 1;
    polyaZ = 0;
    h2 = 0;
    precision = -7;
end

%Determina��o da ordem e do n�mero de v�rtices
if (iscell(Ai))
    order = length(Ai{1});
    vertices = length(Ai);
    inputs_w = size(B1i{1},2);
    inputs_u = size(B2i{1},2);
    outputs = size(Ci{1},1);
else
    [order, vertices] = size(Ai);
    vertices = vertices/order;
    inputs_w = size(B1i,2)/vertices;
    inputs_u = size(B2i,2)/vertices;
    outputs = size(Ci,1);
end

Info.cpusec_m = clock;

%sistema de LMIs
LMIs = set([]);
Info.L = 0;
Info.feas = 0;


%Defining the decision variables
numcoefsW = (factorial(vertices + degW - 1))/(factorial(degW)*factorial(vertices-1));
numcoefsZ = (factorial(vertices + degZ - 1))/(factorial(degZ)*factorial(vertices-1));
for cont=1:numcoefsW
    W{cont} = sdpvar(order,order,'symmetric');
end
for cont=1:numcoefsZ
    Z{cont} = sdpvar(inputs_u,order,'full');
end
if (h2==0)
    h2_quad = sdpvar(1,1);
    optimize = true;
else
    h2_quad = h2*h2;
    optimize = false;
end
X = sdpvar(outputs,outputs,'symmetric');
trX = trace(X);

%Constructing the polynomials
[vecpoly{1}] = poly_struct(Ai,'A',vertices,1);
[vecpoly{2}] = poly_struct(B1i,'B1',vertices,1);
[vecpoly{3}] = poly_struct(B2i,'B2',vertices,1);
[vecpoly{4}] = poly_struct(Ci,'C',vertices,1);
[vecpoly{5}] = poly_struct(Di,'D',vertices,1);
[vecpoly{6}] = poly_struct(W,'W',vertices,degW);
[vecpoly{7}] = poly_struct(Z,'Z',vertices,degZ);
[vecpoly{8}] = poly_struct(X,'X',vertices,0);
[vecpoly{9}] = poly_struct(trX,'trX',vertices,0);
[vecpoly{10}] = poly_struct(h2_quad,'h2',vertices,0);
[vecpoly{11}] = poly_struct(eye(inputs_w),'Im',vertices,0);

%Polya relaxations
poly1 = poly_struct(1,'1',vertices,0);
for cont=1:polyaW
    vecpoly{6} = operation_poly(vecpoly{6},poly1,'*');
end
for cont=1:polyaZ
    vecpoly{7} = operation_poly(vecpoly{7},poly1,'*');
end



tamlmis = 0;
%Constructing the LMIs
Term{1,1} = parser_poly('W',vecpoly,vertices);

LMIs = construct_lmi(Term,vertices,vecpoly,'>');
Info.L = Info.L + (length(LMIs)-tamlmis)*(size(Term{1,1}.data(1).value,1));
tamlmis = length(LMIs);

clear Term;
Term{1,1} = parser_poly('trX - h2',vecpoly,vertices);
LMIs = LMIs + construct_lmi(Term,vertices,vecpoly,'<');
Info.L = Info.L + 1;
tamlmis = length(LMIs);

clear Term;
Term{1,1} = parser_poly('X',vecpoly,vertices);
Term{1,2} = parser_poly('C*W + D*Z',vecpoly,vertices);
Term{2,2} = parser_poly('W',vecpoly,vertices);
LMIs = LMIs + construct_lmi(Term,vertices,vecpoly,'>=');
Info.L = Info.L + (length(LMIs)-tamlmis)*(size(Term{1,1}.data(1).value,1) + size(Term{2,2}.data(1).value,1));

clear Term;
Term{1,1} = parser_poly('A*W + B2*Z + W*A'' + Z''*B2''',vecpoly,vertices);
Term{1,2} = parser_poly('B1',vecpoly,vertices);
Term{2,2} = parser_poly('-Im',vecpoly,vertices);

LMIs = LMIs + construct_lmi(Term,vertices,vecpoly,'<=');
soma = 0;
for cont=1:size(Term,2)
    soma = soma + size(Term{cont,cont}.data(1).value,1);
end
Info.L = Info.L + (length(LMIs)-tamlmis)*soma;
tamlmis = length(LMIs);



Info.cpusec_m = etime(clock,Info.cpusec_m);

Info.K = size(getvariables(LMIs),2);

ops = sdpsettings;
ops.verbose = 0;
ops.solver = 'sedumi';
ops.sedumi.eps = 1e-12; 
if (optimize)
    sol = solvesdp(LMIs,h2_quad,ops);
else
    sol = solvesdp(LMIs,[],ops);
end
Info.cpusec = sol.solvertime;
[p,d]=checkset(LMIs);

Info.feas = 0;
%Isola as solu��es, se existirem
if (((~optimize && (sum(p <= 0) == 0)) || (optimize && sum(roundn(p,precision) < 0) == 0)) && ~isnan(sum(p)))
    for cont=1:numcoefsW
        Info.Wi{cont} = double(W{cont});
        Info.Gi{cont} = double(W{cont});
    end
    for cont=1:numcoefsZ
        Info.Zi{cont} = double(Z{cont});
    end
    if (degW == 0)
        for cont=1:numcoefsZ
            Info.Kst{cont} = Info.Zi{cont}*inv(Info.Gi{1});
        end
    end
    Info.X = double(X);
    Info.h2 = sqrt(double(h2_quad));
    Info.feas = 1;
end

    
return