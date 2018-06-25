
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
T = 10; targets = 1:1:T;
D = 3;  depots = T+1:1:T+D;
N = T+D; total_nodes = [targets,depots];

% Number of robots
K = 2;

% Fuel capacity
L = 1000; %todo :revise

% Ratio of time needed to refuel and time spent traversing the tour
qk = randi([0,1],K,1); %todo: revise

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
fij = ones(N);    %todo: revise

%%%%%%%%%%%%%%%%%%%%%%%%%% Decision variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SYMBOL : DESCRIPTION [MULTIPLICITY]
%
% Pmax   : max path length among robots     [1]
% x_kij  : {0,1}       for k in K, i,j in N [K * N * N]
% p_kij  : {0,1,...,T} for k in K, i,j in N [K * N * N]
% ri     : {0,1,...,L} for i in N [T]

% Total number of variables
total_vars = 1 + (N^2 *K)*2 + T;

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
f = [1;zeros((N^2 *K)*2 + T,1)];


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
lb3 = zeros(T,1);
ub3 = [L*ones(T,1)];

lb = [0;lb1;lb2;lb3];
ub = [Inf;ub1;ub2;ub3];


%% Including the constraint of Pmax
% minmax -> min constraint

bineq_pmax = zeros(K,1);

temp = zeros(K,N^2 *K);
for k=1:K
    for i=1:N
        for j=1:N
            temp(k,j+(i-1)*N+(k-1)*N^2) = (1+qk(k))*cij_per_robot(j+(i-1)*N);
        end
    end
end

Aineq_pmax = [-1*ones(K,1),temp,zeros(K,(N^2 *K) + T)];


%% Degree Constraints

% To ensure only one robot arrives and depart from each target
Aeq6_7 = zeros(T*2,N^2 *K);
beq6_7 = ones(T*2,1);
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

% Robot begins and end at starting position
Aineq8_9 = zeros(K*2,N^2 *K);
bineq8_9 = ones(K*2,1);
for k=1:K
    for i=1:N
        % Equation 8
        Aineq8_9(k,i+(Bk(k)-1)*N) = 1;
        
        % Equation 9
        Aineq8_9(k+K,Bk(k)+(i-1)*N) = 1;        
    end
end

% Every robot visits a target, leaves it
% Equation 10
Aeq10 = zeros(N*K,N^2 *K);
beq10 = zeros(N*K,1);
for k=1:K
    for j=1:N
        for i=1:N
            if j==i
                flag=0;
            else
                flag=1;
            end
            flag = 1; % todo
            Aeq10(j+(k-1)*N,j +(i-1)*N +(k-1)*N^2) = 1*flag;
            Aeq10(j+(k-1)*N,i +(j-1)*N +(k-1)*N^2) = -1*flag;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Block test 1/2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,~] = size(bineq8_9);
bineq = [bineq_pmax;bineq8_9];
Aineq = [Aineq_pmax; zeros(n,1),Aineq8_9, zeros(n,N^2 *K +T)];

beq = [beq6_7;beq10];
[m,~] = size(beq);
Aeq = [zeros(m,1),[Aeq6_7;Aeq10], zeros(m,N^2 *K +T)];

X = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype)

% Visualization
if ~isempty(X)
    adaj = zeros(N,N);
    for k=1:K
        start = 2+(k-1)*N^2; 
        A(:,:,k) = reshape(X(start:start+N^2 -1),[N,N]);
        
        map(A(:,:,k), k, T, N, x_pos, y_pos);
    end
end


%% Capacity & Flow Constraints

x_Aeq2 = zeros(1,(N^2 *K));
p_Aeq2 = zeros(1,(N^2 *K));
temp = zeros(1,(N^2 *K));

x_Aineq2 = zeros(1,(N^2 *K));
p_Aineq2 = zeros(1,(N^2 *K));
bineq2 = zeros(1,1);

% flow through the starting node
for k=1:K
    for i=1:N
        for j=1:N
            if i<=T
                % Equation 11 - Part 1/2
                x_Aeq2(k,j+(i-1)*N+(k-1)*N^2) = -1;
                beq2(k,1) = 0;
                
                % Equation 12 - Part 1/2
                x_Aeq2(K+i+(k-1)*T,j+(i-1)*N+(k-1)*N^2) = -1;
                beq2(K+i+(k-1)*T,1) = 0;
            end
            
        end
        % Equation 11 - Part 2/2
        p_Aeq2(k,i+(Bk(k)-1)*N+(k-1)*N^2) = 1;
        temp(1,Bk(k)+(i-1)*N+(k-1)*N^2) = -1;
        p_Aeq2(k,:) = p_Aeq2(k,:) + temp(1,:);
        temp = zeros(1,(N^2 *K));
    end
end

temp = zeros(1,(N^2 *K));
[y,~] = size(p_Aeq2);
count =1;
for k=1:K
    for i=1:N
        for j=1:N
            % Equation 12 - Part 2/2
            if i<=T
                if j==i
                    flag=0;
                else
                    flag=1;
                end
                p_Aeq2(y+i+(k-1)*T,j +(i-1)*N +(k-1)*N^2) = -1*flag;
                temp(1,i +(j-1)*N +(k-1)*N^2) = 1*flag;
                p_Aeq2(y+i+(k-1)*T,:) = p_Aeq2(y+i+(k-1)*T,:) + temp(1,:);
                temp = zeros(1,(N^2 *K));
            end
        end
    end
end

[y,~] = size(p_Aeq2);
[m,~] = size(beq2);
count =1;
for k=1:K
    for i=T+1:N
        for j=1:N
            % Equation 13
            if i>T
                if j==i
                    flag=0;
                else
                    flag=1;
                end
                p_Aeq2(y+count,j +(i-1)*N +(k-1)*N^2) = -1*flag;
                temp(1,i +(j-1)*N +(k-1)*N^2) = 1*flag;
                p_Aeq2(y+count,:) = p_Aeq2(y+count,:) + temp(1,:);
                temp = zeros(1,(N^2 *K));
               
                beq2(m+count,1) = 0;
            end
        end
        count = count +1;
    end
end

% %Equation 14 - Part 2/2
% count = 1;
% for k=1:K
%     for i=1:N
%         for j=1:N
%             p_Aineq2(count,j+(i-1)*N+(k-1)*N^2) = 1;
%             x_Aineq2(count,j+(i-1)*N+(k-1)*N^2) = -T;
%             bineq2(count,1) = 0;
%             count = count+1;
%         end
%     end
% end

[x,~] = size(x_Aeq2);
[y,~] = size(p_Aeq2);
[m,~] = size(beq2);
Aeq2 = [[x_Aeq2;zeros(y-x,(N^2 *K))],p_Aeq2,zeros(m,T)];

[x,~] = size(x_Aineq2);
[n,~] = size(bineq2);
Aineq2 = [x_Aineq2,p_Aineq2,zeros(x,T)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Block test 2/2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Aineq = [Aineq; zeros(n,1), Aineq2];
bineq = [bineq;bineq2];
Aeq = [Aeq;zeros(m,1), Aeq2];
beq = [beq;beq2];

X = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype)

% Visualization
if ~isempty(X)
    adaj = zeros(N,N);
    for k=1:K
        start = 2+(k-1)*N^2; 
        A(:,:,k) = reshape(X(start:start+N^2 -1),[N,N]);
        
        map(A(:,:,k), k, T, N, x_pos, y_pos);
    end
end



%% Fuel Constraints

% Large constant
M = (L + max(fij(:)));

r_Aineq2 = zeros(1,N);
x_Aineq2f = zeros(1,(N^2 *K));

count = 1;
[m,~] = size(bineq2);
[r,~] = size(r_Aineq2);
for k=1:K
    for i=1:T
        for j=1:T
            % Equation 15
            r_Aineq2(count,j) = 1;
            r_Aineq2(count,i) = -1;
            x_Aineq2f(count,j+(i-1)*N+(k-1)*N^2) = M;
            
            bineq2(m+count,1) = M-fij(i,j);
            
            % Equation 16
            r_Aineq2(T^2 *K+count,j) = -1;
            r_Aineq2(T^2 *K+count,i) = 1;
            x_Aineq2f(T^2 *K+count,j+(i-1)*N+(k-1)*N^2) = M;
            
            bineq2(T^2 *K+m+count,1) = M+fij(i,j);
            
            count = count+1;
        end
    end
end

count = 1;
[m,~] = size(bineq2);
[r,~] = size(r_Aineq2);
for k=1:K
    for i=T+1:N
        for j=1:T
            % Equation 17
            r_Aineq2(r+count,j) = -1;
            x_Aineq2f(r+count,j+(i-1)*N+(k-1)*N^2) = M;
            bineq2(m+count,1) = M-L+fij(i,j);
            
            % Equation 18
            r_Aineq2(D*T*K+r+count,j) = 1;
            x_Aineq2f(D*T*K+r+count,j+(i-1)*N+(k-1)*N^2) = M;
            bineq2(D*T*K+m+count,1) = M+L-fij(i,j);
            
            count = count+1;
        end
    end
end

count = 1;
[m,~] = size(bineq2);
[r,~] = size(r_Aineq2);
for k=1:K
    for i=1:T
        for j=T+1:N
            % Equation 19
            r_Aineq2(r+count,j) = -1;
            x_Aineq2f(r+count,j+(i-1)*N+(k-1)*N^2) = M;
            bineq2(m+count,1) = M-fij(i,j);
            count = count+1;
        end
    end
end

[x,~] = size(x_Aineq2f);
Aineq2 = [Aineq2; x_Aineq2f,zeros(x,N^2 *K),r_Aineq2];

% To get the total number of in/equality equations
[m,~] = size(beq2);
[n,~] = size(bineq2);
eq_count2 = m;
ineq_count2 = n;


%% CPLEX optimization

% total number of in/equality equations
eq_count = eq_count1 + eq_count2;
ineq_count = ineq_count1 + ineq_count2 + K;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Final formulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combining the matrices
Aineq = [Aineq_pmax;zeros(ineq_count1,1),Aineq1, zeros(ineq_count1,N^2 *K +N);zeros(ineq_count2,1),Aineq2];
bineq = [bineq_pmax;bineq1;bineq2];
Aeq = [zeros(eq_count1,1),Aeq1,zeros(eq_count1,N^2 *K+N);zeros(eq_count2,1),Aeq2];
beq = [beq1;beq2];

X = intlinprog(f,(N^2 *K)*2 + N +1,Aineq,bineq,Aeq,beq,lb,ub)

tic
Y = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype)
toc

if ~isempty(X)
    adaj = zeros(N,N);
    for k=1:K
        start = 2+(k-1)*N^2; 
        A(:,:,k) = reshape(X(start:start+N^2 -1),[N,N]);
        adaj = adaj + A(:,:,k);
    end
    plot(graph(adaj),'XData', x_pos, 'YData', y_pos);
end


