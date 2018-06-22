
%% MILP Formulation of Multi-Robot Long-Term Persistent Coverage Problem
%
% Reference: D. Mitchell, M. Corah, N. Chakraborty, K. Sycara and
% N. Michael, "Multi-robot long-term persistent coverage with fuel
% constrained robots," 2015 IEEE (ICRA)
%
% $Author Dharini Dutia     `        $Created June 2018
%
% The indices are modular in terms of the changes in environmental
% parameters. But not with placing, it works sequentially.


clc; clear all; close all;

%% Defining the environment

%Fuel capacity
L = 1000; %todo

%Nodes
targets = 10; T = 1:1:targets;
depots = 3;  D = targets+1:1:targets+depots;
total_nodes = targets+depots; N = [T,D];

%Number of robots
K = 2;

%Ratio of time needed to refuel and time spent traversing the tour
qk = 0.1*ones(K,1); %randi([0,1],K,1); %todo: revise

%Define the starting point of the robots
Bk = randi([D(1),D(depots)],K,1); %D(1:K); %starting from depots

%Shortest distance between the nodes
map

total_vars = (total_nodes^2 *K)*2 + total_nodes + 1; 


%% Integer (bound) constraints
% Equation 4
lb1 = zeros(total_nodes^2 *K,1); %size of X
% upper bound 1 for targets and |T| for depots
ub1 = targets*ones(total_nodes^2 *K,1);
for k=1:K
    for i=1:total_nodes
        for j=1:targets
            % Equation 5
            ub1(j+(i-1)*total_nodes+(k-1)*total_nodes^2,1) = 1;
        end
    end
end

% Equation 14 - Part 1/2
lb2 = zeros(total_nodes^2 *K,1);
ub2 = targets*ones(total_nodes^2 *K,1); %Inf*ones(total_nodes^2 *K,1); %

% Equation 20
lb3 = zeros(total_nodes,1);
ub3 = [L*ones(targets,1); Inf*ones(depots,1)];

lb = [0;lb1;lb2;lb3];
ub = [Inf;ub1;ub2;ub3];

% Obtaining ctype
for i =1:total_vars
    ctype(i) = char('I');
end
for k=1:K
    for i=1:total_nodes
        for j=1:targets
           ctype(j+(i-1)*total_nodes+(k-1)*total_nodes^2 +1) = char('B');
        end
    end
end

%% Including the constraint of Pmax

Pmax_ineq = zeros(1,(total_nodes^2 *K)*2 + total_nodes +1);
bineq_pmax = zeros(K,1);

temp = zeros(K,total_nodes^2 *K);
for k=1:K
    for i=1:total_nodes
        for j=1:total_nodes
            temp(k,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = (1+qk(k))*cij(j+(i-1)*total_nodes+(k-1)*total_nodes^2);
        end
    end
end

Pmax_ineq = [-1*ones(K,1),temp,zeros(K,(total_nodes^2 *K) + total_nodes)];


 %%%%%%%%%%%%%%%%%%%% Block test 1/3  %%%%%%%%%%%%%%%%%%%%
f_pmax = [1;zeros((total_nodes^2 *K)*2 + total_nodes,1)];
% X = intlinprog(f_pmax,total_vars,[],[],[],[],lb,ub)
% Y = cplexmilp(f_pmax,Pmax_ineq,bineq_pmax,[],[], [],[],[],lb,ub,ctype)
Pmax_ineq2 = [-1*ones(K,1),temp];


%% Degree Constraints

Aeq1 = zeros(1,total_nodes^2 *K);
beq1 = zeros(1,1);
Aineq1 = zeros(1,total_nodes^2 *K);
bineq1 = zeros(1,1);

for k=1:K
    %set without starting node
    %woB = N(N~=Bk(k));
    for i=1:total_nodes
        for j=1:total_nodes
            if i<=targets
                % Equation 6
                Aeq1(i,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = 1;
                beq1(i,1) = 1;
                
                % Equation 7
                Aeq1(i+targets, i+(j-1)*total_nodes+(k-1)*total_nodes^2) = 1;
                beq1(i+targets,1) = 1;
                
            end
        end
        % Equation 8
        Aineq1(k,i+(Bk(k)-1)*total_nodes) = 1;
        bineq1(k,1) = 1;
        % Equation 9
        Aineq1(k+K,Bk(k)+(i-1)*total_nodes) = 1;
        bineq1(k+K,1) = 1;
    end
end

% Equation 10
for k=1:K
    for j=1:total_nodes
        [x,~] = size(Aeq1);
        [m,~] = size(beq1);
        for i=1:total_nodes
            if j==i
                flag=0;
            else
                flag=1;
            end
            flag = 1;
            Aeq1(x+1,j +(i-1)*total_nodes +(k-1)*total_nodes^2) = 1*flag;
            Aeq1(x+1,i +(j-1)*total_nodes +(k-1)*total_nodes^2) = -1*flag;
            beq1(m+1) = 0;
        end
    end
end


% number of in/equality equations
[m,~] = size(beq1);
[n,~] = size(bineq1);
eq_count1 = m;
ineq_count1 = n;

 %%%%%%%%%%%%%%%%%%%% Block test 2/3  %%%%%%%%%%%%%%%%%%%%
% f_deg = [1;zeros(total_nodes^2 *K,1)];
% Aineq_deg = [Pmax_ineq2; zeros(n,1),Aineq1];
% bineq_deg = [bineq_pmax;bineq1];
% Aeq_deg = [zeros(m,1),Aeq1]; 
% X = cplexmilp(f_deg,Aineq_deg,bineq_deg,Aeq_deg,beq1, [],[],[],[0;lb1],[Inf;ub1],ctype(1:total_nodes^2 *K +1))

Aineq_deg = [Pmax_ineq; zeros(n,1),Aineq1, zeros(n,total_nodes^2 *K +total_nodes)];
bineq_deg = [bineq_pmax;bineq1];
Aeq_deg = [zeros(m,1),Aeq1, zeros(m,total_nodes^2 *K +total_nodes)];
X = cplexmilp(f_pmax,Aineq_deg,bineq_deg,Aeq_deg,beq1, [],[],[],lb,ub,ctype)

if ~isempty(X)
    adaj = zeros(total_nodes,total_nodes);
    for k=1:K
        start = 2+(k-1)*total_nodes^2; 
        A(:,:,k) = reshape(X(start:start+total_nodes^2 -1),[total_nodes,total_nodes]);
        adaj = adaj + A(:,:,k);
    end
    plot(digraph(adaj),'XData', x_pos, 'YData', y_pos);
end



%% Capacity & Flow Constraints

beq2 = zeros(1,1);
x_Aeq2 = zeros(1,(total_nodes^2 *K));
p_Aeq2 = zeros(1,(total_nodes^2 *K));
temp = zeros(1,(total_nodes^2 *K));

Aineq2 = zeros(1,(total_nodes^2 *K) *2 +total_nodes);
x_Aineq2 = zeros(1,(total_nodes^2 *K));
p_Aineq2 = zeros(1,(total_nodes^2 *K));
bineq2 = zeros(1,1);

for k=1:K
    for i=1:total_nodes
        for j=1:total_nodes
            if i<=targets
                % Equation 11 - Part 1/2
                x_Aeq2(k,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = -1;
                beq2(k,1) = 0;
                
                % Equation 12 - Part 1/2
                x_Aeq2(K+i+(k-1)*targets,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = -1;
                beq2(K+i+(k-1)*targets,1) = 0;
            end
            
        end
        % Equation 11 - Part 2/2
        p_Aeq2(k,i+(Bk(k)-1)*total_nodes+(k-1)*total_nodes^2) = 1;
        temp(1,Bk(k)+(i-1)*total_nodes+(k-1)*total_nodes^2) = -1;
        p_Aeq2(k,:) = p_Aeq2(k,:) + temp(1,:);
        temp = zeros(1,(total_nodes^2 *K));
    end
end

temp = zeros(1,(total_nodes^2 *K));
[y,~] = size(p_Aeq2);
count =1;
for k=1:K
    for i=1:total_nodes
        for j=1:total_nodes
            % Equation 12 - Part 2/2
            if i<=targets
                if j==i
                    flag=0;
                else
                    flag=1;
                end
                p_Aeq2(y+i+(k-1)*targets,j +(i-1)*total_nodes +(k-1)*total_nodes^2) = -1*flag;
                temp(1,i +(j-1)*total_nodes +(k-1)*total_nodes^2) = 1*flag;
                p_Aeq2(y+i+(k-1)*targets,:) = p_Aeq2(y+i+(k-1)*targets,:) + temp(1,:);
                temp = zeros(1,(total_nodes^2 *K));
            end
        end
    end
end

[y,~] = size(p_Aeq2);
[m,~] = size(beq2);
count =1;
for k=1:K
    for i=targets+1:total_nodes
        for j=1:total_nodes
            % Equation 13
            if i>targets
                if j==i
                    flag=0;
                else
                    flag=1;
                end
                p_Aeq2(y+count,j +(i-1)*total_nodes +(k-1)*total_nodes^2) = -1*flag;
                temp(1,i +(j-1)*total_nodes +(k-1)*total_nodes^2) = 1*flag;
                p_Aeq2(y+count,:) = p_Aeq2(y+count,:) + temp(1,:);
                temp = zeros(1,(total_nodes^2 *K));
               
                beq2(m+count,1) = 0;
            end
        end
        count = count +1;
    end
end

%Equation 14 - Part 2/2
count = 1;
for k=1:K
    for i=1:total_nodes
        for j=1:total_nodes
            p_Aineq2(count,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = 1;
            x_Aineq2(count,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = -targets;
            bineq2(count,1) = 0;
            count = count+1;
        end
    end
end

[x,~] = size(x_Aeq2);
[y,~] = size(p_Aeq2);
[m,~] = size(beq2);
Aeq2 = [[x_Aeq2;zeros(y-x,(total_nodes^2 *K))],p_Aeq2,zeros(m,total_nodes)];

[x,~] = size(x_Aineq2);
[n,~] = size(bineq2);
Aineq2 = [x_Aineq2,p_Aineq2,zeros(x,total_nodes)];

 %%%%%%%%%%%%%%%%%%%% Block test 3/3  %%%%%%%%%%%%%%%%%%%%
Aineq_cap = [Pmax_ineq; zeros(ineq_count1,1),Aineq1, zeros(ineq_count1,total_nodes^2 *K +total_nodes); zeros(n,1), Aineq2];
bineq_cap = [bineq_pmax;bineq1;bineq2];
Aeq_cap = [zeros(eq_count1,1),Aeq1, zeros(eq_count1,total_nodes^2 *K +total_nodes); zeros(m,1), Aeq2];
beq_cap = [beq1;beq2];
X = cplexmilp(f_pmax,Aineq_cap,bineq_cap,Aeq_cap,beq_cap, [],[],[],lb,ub,ctype)

if ~isempty(X)
    adaj = zeros(total_nodes,total_nodes);
    for k=1:K
        start = 2+(k-1)*total_nodes^2; 
        A(:,:,k) = reshape(X(start:start+total_nodes^2 -1),[total_nodes,total_nodes]);
        adaj = adaj + A(:,:,k);
    end
    plot(graph(adaj),'XData', x_pos, 'YData', y_pos);
end



%% Fuel Constraints

% Large constant
M = (L + max(fij(:)));

r_Aineq2 = zeros(1,total_nodes);
x_Aineq2f = zeros(1,(total_nodes^2 *K));

count = 1;
[m,~] = size(bineq2);
[r,~] = size(r_Aineq2);
for k=1:K
    for i=1:targets
        for j=1:targets
            % Equation 15
            r_Aineq2(count,j) = 1;
            r_Aineq2(count,i) = -1;
            x_Aineq2f(count,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = M;
            
            bineq2(m+count,1) = M-fij(i,j);
            
            % Equation 16
            r_Aineq2(targets^2 *K+count,j) = -1;
            r_Aineq2(targets^2 *K+count,i) = 1;
            x_Aineq2f(targets^2 *K+count,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = M;
            
            bineq2(targets^2 *K+m+count,1) = M+fij(i,j);
            
            count = count+1;
        end
    end
end

count = 1;
[m,~] = size(bineq2);
[r,~] = size(r_Aineq2);
for k=1:K
    for i=targets+1:total_nodes
        for j=1:targets
            % Equation 17
            r_Aineq2(r+count,j) = -1;
            x_Aineq2f(r+count,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = M;
            bineq2(m+count,1) = M-L+fij(i,j);
            
            % Equation 18
            r_Aineq2(depots*targets*K+r+count,j) = 1;
            x_Aineq2f(depots*targets*K+r+count,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = M;
            bineq2(depots*targets*K+m+count,1) = M+L-fij(i,j);
            
            count = count+1;
        end
    end
end

count = 1;
[m,~] = size(bineq2);
[r,~] = size(r_Aineq2);
for k=1:K
    for i=1:targets
        for j=targets+1:total_nodes
            % Equation 19
            r_Aineq2(r+count,j) = -1;
            x_Aineq2f(r+count,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = M;
            bineq2(m+count,1) = M-fij(i,j);
            count = count+1;
        end
    end
end

[x,~] = size(x_Aineq2f);
Aineq2 = [Aineq2; x_Aineq2f,zeros(x,total_nodes^2 *K),r_Aineq2];

% To get the total number of in/equality equations
[m,~] = size(beq2);
[n,~] = size(bineq2);
eq_count2 = m;
ineq_count2 = n;


%% Objective Function & Co-efficient matrices

% total number of in/equality equations
eq_count = eq_count1 + eq_count2;
ineq_count = ineq_count1 + ineq_count2 + K;

% Objective function to be minimized
f = [1;zeros((total_nodes^2 *K)*2 + total_nodes,1)];


% Combining the matrices
Aineq = [Pmax_ineq;zeros(ineq_count1,1),Aineq1, zeros(ineq_count1,total_nodes^2 *K +total_nodes);zeros(ineq_count2,1),Aineq2];
bineq = [bineq_pmax;bineq1;bineq2];
Aeq = [zeros(eq_count1,1),Aeq1,zeros(eq_count1,total_nodes^2 *K+total_nodes);zeros(eq_count2,1),Aeq2];
beq = [beq1;beq2];


%% CPLEX optimization

% State variables = Pmax; xij of robot1; ... xij robot k; pij of robot1; ... robotk; ri..N
% total variables = (total_nodes^2 *K)*2 + total_nodes +1

X = intlinprog(f,(total_nodes^2 *K)*2 + total_nodes +1,Aineq,bineq,Aeq,beq,lb,ub)

tic
Y = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype)
toc

if ~isempty(X)
    adaj = zeros(total_nodes,total_nodes);
    for k=1:K
        start = 2+(k-1)*total_nodes^2; 
        A(:,:,k) = reshape(X(start:start+total_nodes^2 -1),[total_nodes,total_nodes]);
        adaj = adaj + A(:,:,k);
    end
    plot(graph(adaj),'XData', x_pos, 'YData', y_pos);
end


