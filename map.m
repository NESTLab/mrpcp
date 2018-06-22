%% Environment map

% Targets
xt = randi(25,targets,1);
yt = randi(25,targets,1); 
% xt = [5;7;8;2;7];
% yt = [3;9;1;4;8];

% Depots
xd = randi(25,depots,1);
yd = randi(25,depots,1);
% xd = [6;9;5];
% yd = [5;3;7];

x_pos = [xt',xd'];
y_pos = [yt',yd'];
nodes = [x_pos',y_pos'];

A = zeros(total_nodes);
G = digraph(A);

%% Visualization

figure;
h = plot(G,'XData', x_pos, 'YData', y_pos);
axis square
hold on
title("Map");
highlight(h,[1:1:targets])

%% Parameters

% Distance between two nodes
e_dist = pdist2(nodes,nodes);
cij_per_robot = reshape(e_dist',total_nodes^2,1);
cij=zeros(total_nodes^2 *K,1);
for k=1:K
    cij(1+(k-1)*total_nodes^2:(k*total_nodes^2),1)=cij_per_robot;
end
% e_dist = distances(G);
% cij_per_robot = reshape(e_dist',total_nodes^2,1);
% cij=zeros(total_nodes^2 *K,1);
% for k=1:K
%     cij(1+(k-1)*total_nodes^2:(k*total_nodes^2),1)=cij_per_robot;
% end

% Fuel cost between two nodes
fij = ones(total_nodes);






