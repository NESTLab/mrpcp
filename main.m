
%% MILP Formulation of Multi-Robot Long-Term Persistent Coverage Problem

clc; clear all; close all;

%% Defining the environment

%Fuel capacity
L = 10; %todo

%Nodes
targets = 5; T = 1:1:targets;
depots = 3;  D = targets+1:1:targets+depots;
total_nodes = targets+depots; N = [T,D];

%Number of robots
K = 2;

%Ratio of time needed to refuel and time spent traversing the tour
qk = 0.35; %todo

%Define the starting point of the robots
Bk = [targets+1, targets+2]; %starting depot

%Shortest distance between the nodes
map

%% Constraints

% Integer (bound) constraints
lb1 = zeros(total_nodes^2 *K,1); %size of X
ub1 = targets*ones(total_nodes^2 *K,1);
for k=1:K
    for i=1:total_nodes
        for j=1:targets
            ub1(j+(i-1)*total_nodes+(k-1)*total_nodes^2,1) = 1;
        end
    end
end

% Degree constraints
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
                 Aeq1(i,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = 1;
                 beq1(i,1) = 1;
             
                 Aeq1(i+targets, i+(j-1)*total_nodes) = 1;
                 beq1(i+targets,1) = 1;
             
             end
         end
         
         Aineq1(k,i+(Bk(k)-1)*total_nodes) = 1;
         bineq1(k,1) = 1;
         Aineq1(k+K,Bk(k)+(i-1)*total_nodes) = 1;
         bineq1(k+K,1) = 1;
     end
end


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


%% Capacity & Flow Constraints

lb2 = zeros(total_nodes^2 *K,1);
ub2 = targets*ones(total_nodes^2 *K,1);


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
            if j<=targets
                x_Aeq2(k,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = -1;
                beq2(k,1) = 0;
                
             
                x_Aeq2(K+i,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = -1;
                beq2(K+i,1) = 0;
            end
            
        end
        p_Aeq2(k,i+(Bk(k)-1)*total_nodes+(k-1)*total_nodes^2) = 1;
        temp(1,Bk(k)+(i-1)*total_nodes+(k-1)*total_nodes^2) = -1; 
        p_Aeq2(k,:) = p_Aeq2(k,:) + temp(1,:);
        temp = zeros(1,(total_nodes^2 *K));
    end
end

temp = zeros(1,(total_nodes^2 *K));
for k=1:K
    for i=1:total_nodes
        for j=1:total_nodes
            if j<=targets
                if j==i
                    flag=0;
                else
                    flag=1;
                end
                p_Aeq2(K+i,j +(i-1)*total_nodes +(k-1)*total_nodes^2) = 1*flag;
                temp(1,i +(j-1)*total_nodes +(k-1)*total_nodes^2) = -1*flag;
                p_Aeq2(K+i,:) = p_Aeq2(K+i,:) + temp(1,:);
                temp = zeros(1,(total_nodes^2 *K));
            end
            
            if j>targets 
                if j==i
                    flag=0;
                else
                    flag=1;
                end
                p_Aeq2(K+total_nodes+i,j +(i-1)*total_nodes +(k-1)*total_nodes^2) = 1*flag;
                temp(1,i +(j-1)*total_nodes +(k-1)*total_nodes^2) = -1*flag;
                p_Aeq2(K+total_nodes+i,:) = p_Aeq2(K+total_nodes+i,:) + temp(1,:);
                temp = zeros(1,(total_nodes^2 *K));
                
                beq2(K+total_nodes+i,1) = 0;
            end
            
        end 
    end
end

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
Aineq2 = [x_Aineq2,p_Aineq2,zeros(x,total_nodes)]; 

%% Fuel Constraints

lb3 = zeros(total_nodes,1);
ub3 = L*ones(total_nodes,1);

% Large constant
M = L + max(fij(:));

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
            count = count+1;
            
            % Equation 16
            r_Aineq2(targets^2*k+count,i) = 1;
            r_Aineq2(targets^2*k+count,j) = -1;
            x_Aineq2f(targets^2*k+count,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = M;
            
            bineq2(targets^2*k+m+count,1) = M+fij(i,j);
            
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
            r_Aineq2(depots*targets*2+r+count,j) = 1;
            x_Aineq2f(depots*targets*2+r+count,j+(i-1)*total_nodes+(k-1)*total_nodes^2) = M;
            bineq2(depots*targets*2+m+count,1) = M+L-fij(i,j);
            
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
            bineq2(m+count,1) = M+fij(i,j);
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



%% Objective Function

f = [cij;zeros((total_nodes^2 *K) + total_nodes,1)];

 
%Pmax = (1+qk)* cij* xijk;
%f = (1+qk)*cij;

Aineq = [Aineq1, zeros(ineq_count1,total_nodes^2 *K +total_nodes);Aineq2];
bineq = [bineq1;bineq2];
Aeq = [Aeq1,zeros(eq_count1,total_nodes^2 *K+total_nodes);Aeq2];
beq = [beq1;beq2];
lb = [lb1;lb2;lb3];
ub = [ub1;ub2;ub3];

%% CPLEX optimization

% X = xij of robot1; ... xij robot k; pij of robot1; ... robotk; ri..N

%X = intlinprog(f,total_nodes^2 *K,Aineq,bineq,Aeq,beq,lb,ub)

%doesn't work with cplexmilp  %todo
tic
Y = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub)
toc
