
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
                 Aeq1(i,j+(i-1)*total_nodes) = 1;
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
           Aeq1(x+1,j +(i-1)*total_nodes +(k-1)*total_nodes^2) = 1;
           Aeq1(x+1,i +(j-1)*total_nodes +(k-1)*total_nodes^2) = -1;
           beq1(m+1) = 0;
        end 
    end
end


% number of in/equality equations
[m,~] = size(beq1);
[n,~] = size(bineq1);
eq_count1 = m;
ineq_count1 = n;

Aeq1
beq1 
Aineq1
bineq1
eq_count1
ineq_count1

% Capacity & Flow Constraints

lb2 = zeros(total_nodes^2 *K,1);
ub2 = targets*ones(total_nodes^2 *K,1);


beq2 = zeros(1,1);
x_Aeq2 = zeros(1,(total_nodes^2 *K));
p_Aeq2 = zeros(1,(total_nodes^2 *K));
temp = zeros(1,(total_nodes^2 *K));

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

[x,~] = size(x_Aeq2);
[y,~] = size(p_Aeq2);
[m,~] = size(beq2);
Aeq2 = [[x_Aeq2;zeros(y-x,(total_nodes^2 *K))],p_Aeq2,zeros(m,total_nodes)];


% Fuel Constraints

Aineq2 = zeros(1,(total_nodes^2 *K) *2 +total_nodes);
bineq2 = zeros(1,1);



% To get the total number of in/equality equations
[m,~] = size(beq2);
[n,~] = size(bineq2);
eq_count2 = m;
ineq_count2 = n;

Aeq2
beq2 
Aineq2
bineq2
eq_count2
ineq_count2

%% Objective Function

%Shortest distance between the nodes
map

%f = (ones(total_nodes^2 *K,1)+cij);
f = cij;
 
%Pmax = (1+qk)* cij* xijk;
%f = (1+qk)*cij;

Aineq = [Aineq1, zeros(ineq_count1,total_nodes^2 *K +total_nodes);Aineq2];
bineq = [bineq1;bineq2];
Aeq = [Aeq1,zeros(eq_count1,total_nodes^2 *K+total_nodes);Aeq2];
beq = [beq1;beq2];
lb = [lb1;lb2];
ub = [ub1;ub2];

%% CPLEX optimization

% X = xij of robot1; ... xij robot k; pij of robot1; ... robotk; ri..N

%X = intlinprog(f,total_nodes^2 *K,Aineq,bineq,Aeq,beq,lb,ub)

%doesn't work with cplexmilp  %todo
Y = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub)


