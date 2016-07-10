system = struct();

% k1 = 2 % 1 <= k1 <= 4
k2 = .5;
m1 = 2;
m2 = 1;
% c0 = 2 % 1 <= c0 <= 4

% A = [0 0 1 0; 0 0 0 1; (-k1 + k2) / m1 k2/m1 -c0/m1 0; k2/m2 -k2/m2 0 -c0/m2]

Ai1 = [0 0 1 0; 0 0 0 1; -(1 + k2) / m1 k2/m1 -1/m1 0; k2/m2 -k2/m2 0 -1/m2];
Ai2 = [0 0 1 0; 0 0 0 1; -(4 + k2) / m1 k2/m1 -1/m1 0; k2/m2 -k2/m2 0 -1/m2];
Ai3 = [0 0 1 0; 0 0 0 1; -(1 + k2) / m1 k2/m1 -4/m1 0; k2/m2 -k2/m2 0 -4/m2];
Ai4 = [0 0 1 0; 0 0 0 1; -(4 + k2) / m1 k2/m1 -4/m1 0; k2/m2 -k2/m2 0 -4/m2];

system.Ai = {Ai1 Ai2 Ai3 Ai4}

B1 = [0 1 1 0]';
B2 = [0 0 m1 0]';
C1 = [0 1 0 0];
C2 = [0 0 1 0; 0 0 0 1];
D1 = 0;
D2 = 0;
Dy = [0 0]';
Ze = [0 0]';

system.B1i = {B1 B1 B1 B1};
system.B2i = {B2 B2 B2 B2};
system.C1i = {C1 C1 C1 C1};
system.C2i = {C2 C2 C2 C2};
system.D1i = {D1 D1 D1 D1};
system.D2i = {D2 D2 D2 D2};
system.Dyi = {Dy Dy Dy Dy};

[nx,nw]=size(system.B1i{1});

nc = nx

% pso_massa_mola_nc0 = struct();
load ga_massa_mola_ncnx;
for (iter = 1:10)
    % pso_uncertain(system, nc)
    res_ga = ga_uncertain_2(system, nc)
    % ga_uncertain(system, nc)
    
    ga_massa_mola_ncnx = setfield(ga_massa_mola_ncnx, ['iter_' num2str(iter)],res_ga);
    save ga_massa_mola_ncnx ga_massa_mola_ncnx;
end
