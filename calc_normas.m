function out = calc_normas(system, Kdof, nc)

A = system.Ai{1};
B1 = system.B1i{1};
B2 = system.B2i{1};
C1 = system.C1i{1};
C2 = system.C2i{1};
D1 = system.D1i{1};
D2 = system.D2i{1};
Dy = system.Dyi{1};

% Kdof = 
%nc = 0;

[nx,nw]=size(B1); [nz,nx]=size(C1); [nz,nu]=size(D2); 
ny = size(C2,1);

Atil = [A zeros(nx,nc); zeros(nc,nx) zeros(nc,nc)];
B1til = [B1; zeros(nc,nw)];
B2til = [zeros(nx,nc) B2; eye(nc) zeros(nc,nu)];
C1til = [C1 zeros(nz,nc)];
C2til = [zeros(nc,nx) eye(nc); C2 zeros(ny,nc)];
D1til = D1;
D2til = [zeros(nz,nc) D2];
Dytil = [zeros(nc,nw); Dy];

Acl = Atil + B2til * Kdof * C2til;
Bcl = B1til + B2til * Kdof * Dytil;
Ccl = C1til + D2til * Kdof * C2til;
Dcl = D1til + D2til * Kdof * Dytil;
sss = size(Dcl);
Dcl_h2 = zeros(sss(1),sss(2));

scl = ss(Acl, Bcl, Ccl, Dcl);
scl_h2 = ss(Acl, Bcl, Ccl, Dcl_h2);
Ksf = Kdof * C2til;

[sv, w] = sigma(ss(Acl,Bcl,Ccl,Dcl));
Hm = min(min(sv));

out = struct();
out = setfield(out, 'Kdof', Kdof);
out = setfield(out, 'hinf', norm(scl, 'inf'));
out = setfield(out, 'h2', norm(scl_h2));
out = setfield(out, 'Ksf', Ksf);
out = setfield(out, 'norm_Ksf', norm(Ksf));
out = setfield(out, 'eig_min', min(min(eig(Acl))));
out = setfield(out, 'eig_max', max(max(eig(Acl))));
out = setfield(out, 'h_', Hm);