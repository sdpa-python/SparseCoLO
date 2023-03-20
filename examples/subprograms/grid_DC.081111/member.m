function [val1,val2,val3,val4,val5,val6,val7] = member(nn)
% 
%%%%% val1 = dll
%%%%% val2 = vec_g
%%%%% val3 = nd
%%%%% val4 = nm
% 
%%% coord_x (nk, 2)
nk = (nn + 1)^2;
nm = (2 * nn * (nn+1)) + (2 * nn * nn);
% 
x0 = [0:nn];
y0 = [0:nn];
y0 = repmat(y0,(nn+1),1);
yy = size(y0,1) * size(y0,2);
coord_x = [repmat(x0',(nn+1),1), reshape(y0,yy,1), zeros(nk,1)];
coord_x = 2 * coord_x;
% 
%%% irr (nm, 2)
irr = [];
x0 = [1:nn]';
xn = -1;
for j=1:(nn+1)
    irr = [irr;...
            x0 + ((xn+1) * ones(nn,1)), x0 + ((xn+2) * ones(nn,1))];
    xn = xn + (nn+1);
end
y0 = [1:(nn+1)]';
yn = 0;
for j=1:nn
    irr = [irr;...
            y0 + yn * ones(nn+1,1), y0 + (yn+nn+1) * ones(nn+1,1)];
    yn = yn + (nn+1);
end
clear x0 y0 xn yn
% 
x0 = [1:nn]';
y0 = [(nn+3):((2*nn)+2)]';
yn = 0;
for j=1:nn
    irr = [irr;...
            x0 + yn * ones(nn,1), y0 + yn * ones(nn,1)];
    yn = yn + (nn+1);
end
clear x0 y0 yn
% 
x0 = [2:(nn+1)]';
y0 = [(nn+2):((2*nn)+1)]';
yn = 0;
for j=1:nn
    irr = [irr;...
            x0 + yn * ones(nn,1), y0 + yn * ones(nn,1)];
    yn = yn + (nn+1);
end
clear x0 y0 yn
%
%%% ird (nk, 3)
ird = ones(nk,3);
ird(1,:) = zeros(1,3);
ird(nn+1,:) = zeros(1,3);
%
nd = sum(ird); nd = sum(nd);
%
ii = 0;
for j=1:nk
    for k=1:3
      if ird(j,k) == 1
      ii = ii + 1;
      ird(j,k) = ii;
      else ird(j,k) = nd +1;
      end
    end
end
% 
for i=1:nm
    for j=1:3
        ir(i,j) = ird(irr(i,1), j);
        ir(i,j+3) = ird(irr(i,2), j);
    end
end
% 
for i=1:nm
    j1 = irr(i,1);
    j2 = irr(i,2);
    dx = coord_x(j2,1) - coord_x(j1,1);
    dy = coord_x(j2,2) - coord_x(j1,2);
    dll(i,1) = norm([dx,dy]);
    T0 = [dx/dll(i), dy/dll(i), 0;...
            -dy/dll(i), dx/dll(i), 0;...
            0, 0, 1];
    matT(:,:,i) = [T0, zeros(3,3);...
            zeros(3,3), T0];
end


for i=1:nm
    T0 = zeros(6,nd+1);
    for k=1:6
        for j=1:6
            T0(k,ir(i,j)) = matT(k,j,i);
        end
    end
    TT(:,:,i) = T0(:,1:nd);
end
% 
val1 = dll;
val2 = TT;
val3 = nd;
val4 = nm;
val5 = coord_x;
val6 = ir;
val7 = irr;