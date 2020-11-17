
d12 = 130;
d1 = 40;
d2 = 90;
d0 = 20;

lam1 = 0.3;
lam2 = 0.6;

mx = 0.7;
my = 0.7;
mz = 0.5;
mxy = 1.8;
mxz = 1.4;
myz = 1.4;


%%
% Question 1-7
%Setting up matrix Q
Q1 = [-(lam1+lam2) lam1 lam2 0; mxy -(mxy + lam2) 0 lam2; mxy 0 -(mxy+lam1) lam1; 0 mz mxy -(mxy+mz)];
vec = [0 0 0 0 1]';
Q = Q1';
Q1 = zeros(4,5)';
for i = 1 : 4
    Q1(i,:)=Q(i,:);
end
Q1(end,:) = ones(4,1);
Q2 = Q1;
Q2(:,4)=[0 myz mx -(myz+mx),1];
Q3 = Q1;
Q3(:,4)=[0 mxy mz -(mxy+mz),1];
Q4 = Q1;
Q4(:,4)=[0 mx myz -(mx+myz),1];

%strategy 1
p1=Q1\vec;
%strategy 2
p2=Q2\vec;
%strategy 3
p3=Q3\vec;
%strategy 4
p4=Q4\vec;

prod = [d12,d2,d1,d0];
% Calculating the production
pr1 = prod*p1;
pr2 = prod*p2;
pr3 = prod*p3;
pr4 = prod*p4;

%%
% Question 8-10
Q = [-(lam1+lam2) lam1 lam2 0; mxy -(mxy + lam2) 0 lam2; mxy 0 -(mxy+lam1) lam1; 0 0 0 0];
state = [1 2 3 4];
prod = [d12 d2 d1 d0];
Results = zeros(4,4);
LastRow =[0 mz mxy -(mxy+mz); 0 myz mx -(myz+mx); 0 mxy mz -(mxy+mz); 0 mx myz -(mx+myz)];
for j = 1 : 4
    Q(4,:) = LastRow(j,:);
    maxT = 100000; 
    currentT = 0;
    St = 1;
    flag = 1; 
    T = zeros(4,1);
    pr = 0;
    Ttotal = zeros(4,1);

    while flag == 1 % running until we exceeds time
        % Taking the random time-steps with mean set to Q(n,m)
        stateVec = Q(St,:);
        T(1,1) = exprnd(stateVec(1,1)^(-1));
        T(2,1) = exprnd(stateVec(1,2)^(-1));
        T(3,1) = exprnd(stateVec(1,3)^(-1));
        T(4,1) = exprnd(stateVec(1,4)^(-1));
        T(St,1) = 0;

        Ti = min(T(T>0)); % taking the min time

        for i = 1:4
            if T(i) == Ti
                newstate = i; %uppdating state
            end
        end

        Ttotal(St) = Ttotal(St) + Ti; 
        pr = Ti*prod(St)+pr;

        jumpTo = newstate;
        St = jumpTo;
        currentT = currentT + Ti;

        if currentT >= maxT
            flag = 0;
        end
    end
    Results(:,j) = Ttotal;
end
% Calculating the production
prc1 = prod*Results(:,1)/maxT
prc2 = prod*Results(:,2)/maxT
prc3 = prod*Results(:,3)/maxT
prc4 = prod*Results(:,4)/maxT

error1 = prc1-pr1;
error2 = prc2-pr2;
error3 = prc3-pr3;
error4 = prc4-pr4;

%%
% Question 11-13
% creating our P matrix
h = 0.001;
Q = [-(lam1+lam2) lam1 lam2 0; mxy -(mxy + lam2) 0 lam2; mxy 0 -(mxy+lam1) lam1; 0 0 0 0];
Q = Q*h;
state = [1 2 3 4];
prod = [d12 d2 d1 d0];
LastRow =[0 mz mxy -(mxy+mz); 0 myz mx -(myz+mx); 0 mxy mz -(mxy+mz); 0 mx myz -(mx+myz)];
LastRow = h*LastRow;
Q(1,1) = 1+Q(1,1); % Adding 1 to diagonal
Q(2,2) = 1+Q(2,2);
Q(3,3) = 1+Q(3,3);
LastRow(:,4) = 1 + LastRow(:,4);
Pim = zeros(4,4);

n= 10000000;
for i = 1 : 4
    Q(4,:) = LastRow(i,:);
    P = Q;
    state = 1;
    for k = 1: n
        x1 = P(state,1);
        x2 = P(state,2);
        x3 = P(state,3);
        jump = rand(); % if this value is close to 1, we jump
        if jump > x1+x2+x3              % State 4
            state = 4;
            Pim(i,4) = Pim(i,4) + 1; 
        elseif jump > x1+x2             % State 3
            state = 3;
            Pim(i,3) = Pim(i,3) + 1; 
        elseif jump > x1                % State 2
            state = 2;
            Pim(i,2) = Pim(i,2) + 1;
        else                            % State 1
            state = 1; 
            Pim(i,1) = Pim(i,1) + 1;
        end
    end
end
% Calcuating the production for ea strategy
Pim=Pim/n;
prd1 = Pim(1,:)*prod';
prd2 = Pim(2,:)*prod';
prd3 = Pim(3,:)*prod';
prd4 = Pim(4,:)*prod';

%%
% Question 14
% We only need to count the times we are jumping away from 
% state 1
broken = 0;
broken2 = 0;
b=0.4;
for i = 1 : 10000
    xtime = exprnd(1/lam1);
    ytime = exprnd(1/lam2);
    if xtime < b || ytime < b % b = 0.4
        broken = broken +1;
    end
    if xtime < 2*b/3 || ytime < 2*b/3 % b = 2*0.4/3
        broken2 = broken2 +1;
    end
end
Prob_turbine_works= 1-(broken/10000) %Probability that they work b = 0.4
Prob_turbine_works_anotherB= 1-(broken2/10000) %Probability that they work b = 2*0.4/3


%%
% Question 14
% Only counting the times we are jumping away from state 1
h = 0.001;
jumping = 0;
jumping2 = 0;
for j = 1 :10000
    for i = 1 : b/h % b = 0.4
        jump = rand();
        if jump > P(1,1)
            jumping = jumping + 1;
            break
        end    
    end
    for i = 1 : (2*b/3)/h % b = 2*0.4/3
        jump = rand();
        if jump > P(1,1)
            jumping2 = jumping2 + 1;
            break
        end    
    end
end
Prob_turbine_worksDisc= 1-(jumping/10000) %Probability that they work b = 0.4
Prob_turbine_works_anotherBDisc= 1-(jumping2/10000) %Probability that they work b = 2*0.4/3

