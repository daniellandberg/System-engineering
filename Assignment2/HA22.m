% HA2 - Systems Engineering
% Birgir Steinn Hermannsson and Daniel Landberg
close all; clear all;

%% Parameters

n = 20; %no. of expected efficient solutions before Cmax breaks

Cmax = 500;                                     %Budget
lambda = (1/1000)*[50 40 45 51 25 48 60 35 15]; %Malfunctioning intensities for LRUs
T = [4 7 14 5 10 18 24 8 12];                   %Repair times for LRUs
c = [14 19 25 15 10 45 80 33 30];               %Purchase cost for LRUs

%% Generate table of quotients

R = [];
Rstore = [];
z = [];
for sj = 0:n
    for j = 1:9
        p = [0 0 0 0 0 0 0 0 0];
        for k =(sj+1):10
            p(j) = p(j)+((lambda(j)*T(j))^k/(factorial(k)))*exp(-lambda(j)*T(j));
        end
        R(j) = p(j);
        Rstore(sj+1,j) = R(j);
        z(sj+1,j) = R(j)/(c(j));
    end  
end
zbackup = z;

%% Marginal allocation algorithm

%Initiation of vectors and matrices 
s = zeros(n,9);                       %Initiate spare allocation
C = zeros(n,1);                       %Initiate total cost  
EBO = zeros(n,1);                     %Initiate expected no. of backorders
C(1) = 0;                             %Set C^(0)=0, because no spare parts  
EBO(1) = 0;
for j = 1:9
    EBO(1) = EBO(1)+lambda(j)*T(j);   %Set EBO^(0) accordingly
end 

%Selects the largest uncanceled quotient in the z table.
%If several equally large, one is chosen arbitrarily.
%The largest quotient is set to 0 and its column number saved.
k = 1;
for j = 1:n
    zmax = 0;                   %Initate variable zmax
    l = 0;                      %Initate index l
    for i = 1:9                 %Iterate for every column i
        t = 1;                  %S        tart at row t=1
        q = z(t,i);             %Find value at z(1,i)
        if q == 0               %Check if quotient has been cancelled
            while q == 0        %While value is 0...
                t = t+1;        %Increase row index
                q = z(t,i);     %Update value from z
            end
        end
        if q>zmax 
            zmax = q;           %Update max value if larger than before
            tmax = t;           %Update row index for max value
            lmax = i;           %Update column index for max value
        end
    end
    z(tmax,lmax) = 0;           %After all columns have been checked, max value set to 0
    k = k+1;                    %Update counter
    s(k,lmax) = s(k-1,lmax)+1;  %Increase s_j^(k) by 1 (i.e. part s_j that 
    for i = 1:9
        if i ~= lmax
            s(k,i) = s(k-1,i);  %Set all other vales in s to previous values
        end
    end  
    C(k) = C(k-1)+c(lmax);                  %Update total cost
    EBO(k) = EBO(k-1)-Rstore(tmax,lmax);    %Update EBO
    if C(k) >= Cmax
    break
    end
end

%% Plot the results

EBOPlot = nonzeros(EBO);
range = length(EBOPlot);
plot(C(1:range),EBOPlot,'-o')
xlabel('C(s)')

ylabel('EBO(s)')

%% Dynamic programming part

% States: if we have budget for the part or not
% Stages: 500?
% decisions: if we buy the part or keeps it as a EBO
% Sn+1 = min(ebo(s)) while c(s) <= current budget) 
% Value-functions: Minimal EBO
% The realation between recusion and Value-function:
% The value function is finding the optimal for every subproblem. The
% recurssion then updates the matrix accordingly to the optimal for every
% subproblem and the calls the value function again to find the next
% optimal.


% Setup

% create EBO matrix
EBO = zeros(51,length(T));
for i = 1:length(T)
    EBO(1,i) = T(i)*lambda(i);
    for s = 1:50
        prob = 1 - poisscdf(s-1, EBO(1,i)); 
        EBO(s+1,i) = EBO(s,i) - prob;
    end
end

% crete a matrix with all possible purchases for every budget 
% ex: column 5 (with part cost 10) will have the number 50 at the last row
LRU = zeros(Cmax+1, length(T));
for i = 1 : Cmax+1
    row = ones(1,length(T)) * i-1;
    LRU(i,:) = fix(row./c); % rounding to lower integer 
end

% Setting up recursion, Starting from LRU9
Matrix = zeros(Cmax+1, length(T)+1);
Matrix(:,length(T)) = LRU(:,length(T));
for i = 1 : Cmax+1
    idx = LRU(i,length(T))+1;
    Matrix(i,length(T)+1) = EBO(idx,length(T));
end

% end of setup

% recursion for every LRU
for i = length(T)-1 :-1:1
    Matrixtemp = Matrix; % Need a copy of Matrix in recursion
    Matrix = allVals(Matrix,Matrixtemp, i, Cmax, LRU, T, c, EBO); % Calling allVals which loops over every budget
    
end

%% Plot Dynamic programming problem 2
plot(0:Cmax,Matrix(1:end,length(T)+1))
hold on
xlabel('C(s)')
ylabel('EBO(s)')
title('Dynamic programming')
%% Plotting problem 1 and 2 together

plot(C(1:range),EBOPlot,'-o')
hold on
plot(0:Cmax,Matrix(1:end,length(T)+1))
xlabel('C(s)')
ylabel('EBO(s)')
title('Dynamic Programming vs Marginal Allocation')
legend('Marginal allocation','Dynamic programming')
%%
% Looping over evry single budget
function Matrix = allVals(Matrix,Matrixtemp, i, Cmax,LRU, T, c, EBO)
    for bud = 1 : Cmax+1
        NrParts = LRU(bud,i); % Nr of parts with current budget of LRUi
        vec = partloop(NrParts,i,bud, T, c, EBO, Matrixtemp); % Partloop: bf over all setups
        Matrix(bud,:) = vec;
    end
    
end

%This function is finding the optimal for every subproblem (value function)
function ReturnArray = partloop(NrParts, i, bud, T, c, EBO, Matrixtemp)
    MiniMatrix = zeros(NrParts+1,length(T)+1);
    for part = 0 : NrParts
        MoneyLeft = bud - part*c(i); % this is current budget - how many parts i we have purchased
        LRUtemp = Matrixtemp(MoneyLeft,:); % Picking out the row matching our budget (extra cash)
        vec = LRUtemp(1,:); 
        vec(1,i) = part;
        MiniMatrix(part+1,:) =  vec; %testing different amount of parts
        MiniMatrix(part+1, length(T)+1) = EBO(part+1, i) + LRUtemp(length(T)+1); % EBO for this setup
    end
    [~, idx] = min(MiniMatrix(:,length(T)+1)); % taking the one with the lowest EBO's
    ReturnArray = MiniMatrix(idx,:); %returning to uppdate the solution matrix with the best row for this budget
end

