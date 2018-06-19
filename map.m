%% Environment map

% Targets
xt = randi(10,targets,1);
yt = randi(10,targets,1);

% Depots
xd = randi(10,depots,1);
yd = randi(10,depots,1);

nodes = [xt,yt;xd,yd];

A = ones(total_nodes);
t = diag(ones(1,total_nodes));
A = A -t;
G = graph(A);

%% Visualization
% figure;
% axis equal;
% scatter(xt,yt);
% xlim([0 10]);
% ylim([0 10]);
% hold on
% scatter(xd,yd);
% title("Map");
% legend('Targets', 'Depots');

%plot(G);

%% Parameters

% Distance between two nodes
e_dist = pdist2(nodes,nodes);
cij_per_robot = reshape(e_dist',total_nodes^2,1);
cij=zeros(total_nodes^2 *K,1);
for k=1:K
    cij(1+(k-1)*total_nodes^2:(k*total_nodes^2),1)=cij_per_robot;
end

% Fuel cost between two nodes
fij = ones(total_nodes);






