%%%% MMA%%%%
function topMMA(nelx,nely,volfrac,penal,rmin)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
 
%F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
%fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
 
%F = sparse(2*(nelx*(nely+1)+fix(nely/2)+1),1,-1,2*(nely+1)*(nelx+1),1);
%fixeddofs = 1:2*(nely+1);
 
F = sparse(2*(fix(nely/2)+1)-1,1,1,2*(nely+1)*(nelx+1),1);
fixeddofs = [1 2 2*nely+1 2*nely+2];
 
U = zeros(2*(nely+1)*(nelx+1),1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
x = repmat(volfrac,nely,nelx);
n = nelx*nely;
xold1 = repmat(0,n,1);
xold2 = xold1; 
low = repmat(0,n,1);
upp = repmat(1,n,1);
loop = 0;
change = 1;
%% START ITERATION
while change > 0.01
  loop = loop + 1;
  %% FE-ANALYSIS
  sK = reshape(KE(:)*(Emin+x(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
  f = sum(sum((Emin+x.^penal*(E0-Emin)).*ce));
  dc = -penal*(E0-Emin)*x.^(penal-1).*ce;
  %% FILTERING/MODIFICATION OF SENSITIVITIES
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
  xval = reshape(x,n,1);
  xmin = zeros(n,1);
  xmax = ones(n,1);
  f0val = f;
  df0dx = reshape(dc,n,1);
  fval = sum(xval)-volfrac*n;
  dfdx = ones(1,n);
  a0 =1;
  a = 0;
  c = 100000;
  d = 0;
  
  
  [xmma,~,~,~,~,~,~,~,~,low1,upp1] = mmasub(1,n,loop,xval,xmin,xmax,xold1,xold2,f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d);
 
  xold2 = xold1; 
  xold1 = xmma;
  low = low1;
  upp = upp1;
  
  
  change = max(abs(xmma(:)-x(:)));
  x = reshape(xmma,nely,nelx);
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,f, ...
    mean(x(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-x); caxis([0 1]); axis equal; axis off; drawnow;
end


