clear all
clc

hor = 2.01;
L = linspace(1,5,100);
Go = .4;
K = 5.81E9;

pts1 = [];
p_pts1 = [];
for m = 1:length(L)
    hor = 20*L(m)+.01;
    c = 18*K/(pi*hor^4);
    lim1 = floor(hor/L(m));
    
    for z = 1:lim1
        ref = [0,0,-z*L(m)];
        for i = -lim1:1:lim1
            for j = -lim1:1:lim1
                for k = -lim1:1:lim1
                    pts1(end+1,1) = i*L(m);
                    pts1(end+1,2) = j*L(m);
                    pts1(end+1,3) = k*L(m);
                    pts1(end+1,4) = sqrt((i*L(m)-ref(1))^2+(j*L(m)-ref(2))^2+(k*L(m)-ref(3))^2);
                end
            end
        end

        for i = 1:length(pts1)
            if pts1(i,4) < hor && pts1(i,3) > -1
                p_pts1(end+1,:) = pts1(i,:);
            end
        end
    end
    
    s0_1(m) = sqrt((2*Go)/(c*sum(p_pts1(:,4))));
end

% for m = 1:length(L)
%     
% end
plot(L,s0_1)
axis([0,max(s0_1)],[0,max(L)])
%p_pts(end+1,:) = [0,0,-1,0];

% figure
% scatter3(p_pts(:,1),p_pts(:,2),p_pts(:,3))
% hold on
% for i = 1:length(p_pts)-1
%     plot3([p_pts(i,1),ref(1)],[p_pts(i,2),ref(2)],[p_pts(i,3),ref(3)],'k')
% end