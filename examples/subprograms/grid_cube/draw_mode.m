function [dummy] = draw_mode(vec_u,coord_x,irr,ir,dll,nm,matT,ccs)
% 
dummy = 1;
% 
hold off;
for i=1:nm
    j1 = irr(i,1);
    j2 = irr(i,2);
    mem0(1,:,i) = coord_x(j1,:);
    mem0(2,:,i) = coord_x(j2,:);
end
% axes('DataAspectRatioMode','manual')
hold on;
for i=1:nm
    axis equal;
    plot(mem0(:,1,i), mem0(:,2,i), ':', 'LineWidth',2);
end
% 
uuu = [vec_u; 0];
for i=1:nm
    for j=1:6
        dis_m(j,i) = uuu(ir(i,j));
    end
end
clear uuu
for i=1:nm
    def(:,:,i) = mem0(:,:,i) + reshape(dis_m(:,i),3,2)';
end
% for i=1:nm
%     plot(def(:,1,i), def(:,2,i), ':', 'LineWidth',2);
% end
for i=1:nm
    dis_l(:,i) = matT(:,:,i) * vec_u;
end
for i=1:nm
    ll(i) = dis_l(4,i) - dis_l(1,i) + dll(i);
end
for i=1:nm
    bb(1,i) = dis_l(2,i);
    bb(2,i) = dis_l(3,i);
    bb(3,i) = -((2*dis_l(3,i)+dis_l(6,i))/dll(i))...
        + (3*(dis_l(5,i)-dis_l(2,i))/(dll(i)^(2)));
    bb(4,i) = ((dis_l(3,i)+dis_l(6,i))/(dll(i)^(2)))...
        - 2*(dis_l(5,i)-dis_l(2,i))/(dll(i)^(3));
end
for i=1:nm
    j1 = irr(i,1);
    j2 = irr(i,2);
    dx = coord_x(j2,1) - coord_x(j1,1);
    dy = coord_x(j2,2) - coord_x(j1,2);
    dir_cos(:,i) = [dx/dll(i); dy/dll(i); -dy/dll(i); dx/dll(i)];
end
hold on
nnn = 10;
for i=1:nm
    if ccs(i) > 10^(-7)
    xx(1,:) = 0:nnn;
    xx(1,:) = xx(1,:) * dll(i) / nnn;
    %%%% >>> for jj=1
    xx(2,1) = bb(1,i) + (bb(2,i)*xx(1,1))...
        + (bb(3,i)*(xx(1,1)^(2))) + (bb(4,i)*(xx(1,1)^(3)));
    xx(1,1) = xx(1,1) * (ll(i) / dll(i));
    dir_theta = mat(dir_cos(:,i)');
    xx(:,1) = dir_theta * xx(:,1);
    dis_tmp = (def(1,1:2,i)') - xx(:,1);
    xx(:,1) = xx(:,1) + dis_tmp;
    %%%% <<< for jj=1
    for jj=2:(nnn+1)
        xx(2,jj) = bb(1,i) + (bb(2,i)*xx(1,jj))...
            + (bb(3,i)*(xx(1,jj)^(2))) + (bb(4,i)*(xx(1,jj)^(3)));
        xx(1,jj) = xx(1,jj) * (ll(i) / dll(i));
        xx(:,jj) = (dir_theta * xx(:,jj)) + dis_tmp;
    end
    plot(xx(1,:), xx(2,:), 'r-', 'LineWidth',2);
end
end
clear xx


