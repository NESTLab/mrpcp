
%% MILP Formulation of Multi-Robot Long-Term Persistent Coverage Problem
%
% Reference: D. Mitchell, M. Corah, N. Chakraborty, K. Sycara and
% N. Michael, "Multi-robot long-term persistent coverage with fuel
% constrained robots," 2015 IEEE (ICRA)
%
% $Author Dharini Dutia     `        $Created June 2018
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cleanup
clc; clear all; close all;
% Mersenne-Twister with seed 0
% This makes the current execution repeatable
rng('default');


% Defining the environment
% Nodes
T = 5; targets = 1:1:T;
D = 3;  depots = T+1:1:T+D;
N = T+D; total_nodes = [targets,depots];

% Number of robots
K = 2;

% Fuel capacity
L = 2*T; %todo :revise

% Ratio of time needed to refuel and time spent traversing the tour
qk = 0.1*ones(K,1);%rand(K,1)/2;

% Define the starting point of the robots
Bk = randi([depots(1),depots(D)],K,1); %D(1:K); %starting from depots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Node Generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Targets
xt = randi(25,T,1);
yt = randi(25,T,1); 
% xt = [11;3;7;11;15;7;16;18;6;3];  % For testing
% yt = [8;8;11;13;3;7;21;1;24;19];

% Depots
xd = randi(25,D,1);
yd = randi(25,D,1);
% xd = [13;15;6];
% yd = [12;25;14];

% Combining the co-ordinates
x_pos = [xt',xd'];
y_pos = [yt',yd'];
nodes = [x_pos',y_pos'];

% Shortest distance between the nodes
e_dist = pdist2(nodes,nodes);
% cij = time required to traverse each edges
cij_per_robot = reshape(e_dist',N^2,1);
% if robots' velocities are not equivalent:
% use next block to obtain different cij for each robot
% cij=zeros(total_nodes^2 *K,1);
% for k=1:K
%     cij(1+(k-1)*total_nodes^2:(k*total_nodes^2),1)=cij_per_robot;
% end

% Fuel cost between two nodes
fij = [e_dist(1:T,1:T),0.5*L*ones(T,D); ...
       0.5*L*ones(D,T),0.8*L*ones(D,D)];
%%%%%%%%%%%%%%%%%%%%%%%%%% Decision variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SYMBOL : DESCRIPTION [MULTIPLICITY]
%
% Pmax   : max path length among robots     [1]
% x_kij  : {0,1}       for k in K, i,j in N [K * N * N]
% p_kij  : {0,1,...,T} for k in K, i,j in N [K * N * N]
% ri(k)  : {0,1,...,L} for i in T [T * K]

% Total number of variables
total_vars = 1 + (N^2 *K)*2 + T*K;

% Obtaining ctype: for CPLEX MILP solver
ctype(1) = char('C');
for i =2:total_vars
    ctype(i) = char('I');
end
for k=1:K
    for i=1:N
        for j=1:T
           ctype(j+(i-1)*N+(k-1)*N^2 +1) = char('B');
        end
    end
end

% Objective function to be minimized
f = [1;zeros((N^2 *K)*2 + T*K,1)];


%% Integer (bound) constraints

% Equation 4
lb1 = zeros(N^2 *K,1); %size of X
% upper bound 1 for targets and |T| for depots
ub1 = T*ones(N^2 *K,1);
for k=1:K
    for i=1:N
        for j=1:T
            % Equation 5
            ub1(j+(i-1)*N+(k-1)*N^2,1) = 1;
        end
    end
end

% Equation 14 - Part 1/2
lb2 = zeros(N^2 *K,1);
ub2 = T*ones(N^2 *K,1); 

% Equation 20
lb3 = zeros(T*K,1);
ub3 = L*ones(T*K,1);

lb = [0;lb1;lb2;lb3];
ub = [Inf;ub1;ub2;ub3];


%% Including the constraint of Pmax
% minmax -> min constraint
% Equation 3
bineq_pmax = zeros(K,1);

temp = zeros(K,N^2 *K);
for k=1:K
    for i=1:N
        for j=1:N
            temp(k,j+(i-1)*N+(k-1)*N^2) = (1+qk(k))*cij_per_robot(j+(i-1)*N);
        end
    end
end

Aineq_pmax = [-1*ones(K,1),temp,zeros(K,(N^2 *K) + T*K)];


%% Degree Constraints

% To ensure only one robot arrives and depart from each target
Aeq6_7 = zeros(T*2,N^2 *K);
for k=1:K
    for i=1:T
        for j=1:N
            % Equation 6
            Aeq6_7(i,j+(i-1)*N+(k-1)*N^2) = 1;
            
            % Equation 7
            Aeq6_7(i+T, i+(j-1)*N+(k-1)*N^2) = 1;
        end        
    end
end
Aeq6_7 = [zeros(T*2,1),Aeq6_7,zeros(T*2,N^2 *K+T*K)];
beq6_7 = ones(T*2,1);

% Robot begins and end at starting position
Aineq8_9 = zeros(K*2,N^2 *K);
for k=1:K
    for i=1:N
        % Equation 8
        Aineq8_9(k,i+(Bk(k)-1)*N+(k-1)*N^2) = 1;
        
        % Equation 9
        Aineq8_9(k+K,Bk(k)+(i-1)*N+(k-1)*N^2) = 1;        
    end
end
Aineq8_9 = [zeros(K*2,1),Aineq8_9, zeros(K*2,N^2 *K +T*K)];
bineq8_9 = ones(K*2,1);

% Every robot visits a target, leaves it
% Equation 10
Aeq10 = zeros(N*K,N^2 *K);
for k=1:K
    for j=1:N
        for i=1:N
            Aeq10(j+(k-1)*N,j +(i-1)*N +(k-1)*N^2) = 1;
            Aeq10(j+(k-1)*N,i +(j-1)*N +(k-1)*N^2) = -1;
        end
    end
end
Aeq10 = [zeros(N*K,1),Aeq10,zeros(N*K,N^2 *K+T*K)];
beq10 = zeros(N*K,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Block test 1/2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bineq = [bineq_pmax;bineq8_9];
% Aineq = [Aineq_pmax;Aineq8_9];
% 
% beq = [beq6_7;beq10];
% Aeq = [Aeq6_7;Aeq10];
% 
% X = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype)
% 
% % Visualization
% if ~isempty(X)
%     adaj = zeros(N,N);
%     for k=1:K
%         start = 2+(k-1)*N^2; 
%         A(:,:,k) = transpose(reshape(X(start:start+N^2 -1),[N,N]));
%         
%         map(A(:,:,k), k, T, N, x_pos, y_pos);
%     end
% end


%% Capacity & Flow Constraints

% Equation 11: flow through the starting node
x_Aeq11 = zeros(K,(N^2 *K));
p_Aeq11 = zeros(K,(N^2 *K));
for k=1:K
    for i=1:N
        for j=1:N
            if i<=T
                x_Aeq11(k,j+(i-1)*N+(k-1)*N^2) = -1;
            end 
        end
        if i~= Bk(k)
            p_Aeq11(k,i+(Bk(k)-1)*N+(k-1)*N^2) = 1;
            p_Aeq11(k,Bk(k)+(i-1)*N+(k-1)*N^2) = -1;
        end
    end
end
Aeq11 = [zeros(K,1),x_Aeq11,p_Aeq11,zeros(K,T*K)];
beq11(K,1) = 0;

% Equation 12 : Capacity updated after visiting each node
x_Aeq12 = zeros(T*K,(N^2 *K));
p_Aeq12 = zeros(T*K,(N^2 *K));
for k=1:K
    for i=1:T
        for j=1:N
            x_Aeq12(i+(k-1)*T,j+(i-1)*N+(k-1)*N^2) = -1;
            if j~=i
                p_Aeq12(i+(k-1)*T,j +(i-1)*N +(k-1)*N^2) = -1;
                p_Aeq12(i+(k-1)*T,i +(j-1)*N +(k-1)*N^2) = 1;
            end                      
        end
    end
end
Aeq12 = [zeros(T*K,1),x_Aeq12,p_Aeq12,zeros(T*K,T*K)];
beq12(T*K,1) = 0; 

% Equation 13 : Capacity remains the same after passing a depot
x_Aeq13 = zeros(D*K,(N^2 *K));
p_Aeq13 = zeros(D*K,(N^2 *K));
count = 1;
for k=1:K
    for i=T+1:N
        for j=1:N
           if j~=i && i~=Bk(k)
            p_Aeq13(count,j +(i-1)*N +(k-1)*N^2) = -1;
            p_Aeq13(count,i +(j-1)*N +(k-1)*N^2) = 1;
           end
        end
        count = count +1;
    end
end
Aeq13 = [zeros(D*K,1),x_Aeq13, p_Aeq13,zeros(D*K,T*K)];
beq13 = zeros(D*K,1);

%Equation 14 - Part 2/2
% target capacity should not exceed T
count = 1;
x_Aineq14 = zeros(N^2 *K,(N^2 *K));
p_Aineq14 = zeros(N^2 *K,(N^2 *K));
for k=1:K
    for i=1:N
        for j=1:N
            p_Aineq14(count,j+(i-1)*N+(k-1)*N^2) = 1;
            x_Aineq14(count,j+(i-1)*N+(k-1)*N^2) = -T;
            count = count+1;
        end
    end
end
Aineq14 = [zeros(N^2*K,1),x_Aineq14,p_Aineq14,zeros(N^2 *K,T*K)];
bineq14 = zeros(N^2 *K,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Block test 2/2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bineq = [bineq_pmax;bineq8_9;bineq14];
% Aineq = [Aineq_pmax;Aineq8_9;Aineq14];
% 
% beq = [beq6_7;beq10;beq11;beq12;beq13];
% Aeq = [Aeq6_7;Aeq10;Aeq11;Aeq12;Aeq13];
% 
% X = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype)
% 
% % Visualization
% if ~isempty(X)
%     adaj = zeros(N,N);
%     for k=1:K
%         start = 2+(k-1)*N^2; 
%         A(:,:,k) = transpose(reshape(X(start:start+N^2 -1),[N,N]));
%         
%         map(A(:,:,k), k, T, N, x_pos, y_pos);
%     end
% end



%% Fuel Constraints

% Large constant
M = (L + max(fij(:)));

% Equation 15 & 16
% fuel lost between two nodes = fuel cost of travelling between them
r_Aineq15_16 = zeros(T^2 *K *2,T*K);
x_Aineq15_16 = zeros(T^2 *K *2,(N^2 *K));
count = 1;
for k=1:K
    for i=1:T
        for j=1:T
            % Equation 15
            %if j~= i   %not sure
            r_Aineq15_16(count,j+(k-1)*T) = 1;
            r_Aineq15_16(count,i+(k-1)*T) = -1;
            %end
            x_Aineq15_16(count,j+(i-1)*N+(k-1)*N^2) = M;
            
            % Equation 16
            %if j~= i   %not sure
            r_Aineq15_16(T^2 *K+count,j+(k-1)*T) = -1;
            r_Aineq15_16(T^2 *K+count,i+(k-1)*T) = 1;
            %end
            x_Aineq15_16(T^2 *K+count,j+(i-1)*N+(k-1)*N^2) = M;
            
            count = count+1;
        end
    end
end
Aineq15_16 = [zeros(T^2 *K *2,1),x_Aineq15_16,zeros(T^2 *K*2,(N^2 *K)),r_Aineq15_16];
bineq15_16 = [(M-fij(i,j))*ones(T^2 *K,1); (M+fij(i,j))*ones(T^2 *K,1)];
                      
% Equation 17 & 18
% fuel level at target visited after leaving a depot = fuel capacity - fuel
% cost of traversal
r_Aineq17_18 = zeros(D*T*K *2,T*K);
x_Aineq17_18 = zeros(D*T*K *2,(N^2 *K));
count = 1;
for k=1:K
    for i=T+1:N
        for j=1:T
            % Equation 17
            r_Aineq17_18(count,j+(k-1)*T) = -1;
            x_Aineq17_18(count,j+(i-1)*N+(k-1)*N^2) = M;
            
            % Equation 18
            r_Aineq17_18(D*T*K+count,j+(k-1)*T) = 1;
            x_Aineq17_18(D*T*K+count,j+(i-1)*N+(k-1)*N^2) = M;
            
            count = count+1;
        end
    end
end
Aineq17_18 = [zeros(D*T*K *2,1),x_Aineq17_18,zeros(D*T*K*2,(N^2 *K)),r_Aineq17_18];
bineq17_18 = [(M-L+fij(i,j))*ones(D*T*K,1); (M+L-fij(i,j))*ones(D*T*K,1)];
            
% Equation 19
% restricts fuel lost in approaching a depot to being most the cost to 
% travel from the preceding target
r_Aineq19 = zeros(T*D*K,T*K);
x_Aineq19 = zeros(T*D*K,(N^2 *K));
count = 1;
for k=1:K
    for i=1:T
        for j=T+1:N
            r_Aineq19(count,i+(k-1)*T) = -1;
            x_Aineq19(count,j+(i-1)*N+(k-1)*N^2) = M;
            count = count+1;
        end
    end
end
Aineq19 = [zeros(T*D*K,1),x_Aineq19,zeros(T*D*K,(N^2 *K)),r_Aineq19];
bineq19 = (M-fij(i,j))*ones(T*D*K,1);

%% CPLEX optimization

% Combining the matrices
bineq = [bineq_pmax;bineq8_9;bineq14;bineq15_16;bineq17_18;bineq19];
Aineq = [Aineq_pmax;Aineq8_9;Aineq14;Aineq15_16;Aineq17_18;Aineq19];

beq = [beq6_7;beq10;beq11;beq12;beq13];
Aeq = [Aeq6_7;Aeq10;Aeq11;Aeq12;Aeq13];

% To get the total number of in/equality equations
[m,~] = size(beq);
[n,~] = size(bineq);
eq_count = m;
ineq_count = n;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Final formulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
X = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype)
toc

X = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype)

% Visualization
if ~isempty(X)
    adaj = zeros(N,N);
    for k=1:K
        start = 2+(k-1)*N^2; 
        A(:,:,k) = transpose(reshape(X(start:start+N^2 -1),[N,N]));
        
        map(A(:,:,k), k, T, N, x_pos, y_pos);
    end
end


