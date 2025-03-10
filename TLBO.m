function [X, FVAL, BestFVALIter, pop] = TLBO(x0,FITNESSFCN,lb,ub,T,NPop)
% Teaching Learning Based optimization (TLBO)
% SanitizedTLBO attempts to solve problems of the following forms:
%         min F(X)  subject to: lb <= X <= ub
%          X        
%                           
%  [X,FVAL,BestFVALIter, pop] = SanitizedTLBO(FITNESSFCN,lb,ub,T,NPop)
%  FITNESSFCN   - function handle of the fitness function
%  lb           - lower bounds on X
%  ub           - upper bounds on X
%  T            - number of iterations
%  NPop         - size of the population (class size)  
%  X            - minimum of the fitness function determined by SanitizedTLBO
%  FVAL         - value of the fitness function at the minima (X)
%  BestFVALIter - the best fintess function value in each iteration
%  pop          - the population at the end of the specified number of iterations 
% preallocation to store the best objective function of every iteration
% and the objective function value of every student
count=0;
BestFVALIter = NaN(T,1);
obj = NaN(NPop,1);
% Determining the size of the problem
D = length(lb);
% Generation of initial population
pop = x0;
%  Evaluation of objective function
%  Can be vectorized 
for p = 1:NPop
    obj(p) = FITNESSFCN(pop(p,:));
    count=count+1;
end
for gen = 1: T
    
    % Partner selection for all students
    % Note that randperm has been used to speedup the partner selection.
    Partner = randperm(NPop);
    % There is a remote possibility that the ith student will have itself as its partner
    % No experiment is available in literature on the disadvantages of
    % a solution having itself as partner solution.
    
    for i = 1:NPop
        
        % ----------------Begining of the Teacher Phase for ith student-------------- %
        mean_stud = mean(pop);
        
        % Determination of teacher
        [~,ind] = min(obj);
        best_stud = pop(ind,:);
        
        % Determination of the teaching factor
        TF = randi([1 2],1,1);
        
        % Generation of a new solution
        NewSol = pop(i,:) + rand(1,D).*(best_stud - TF*mean_stud);
        
        % Bounding of the solution
        for z=1:length(ub)
        NewSol(z) = max(min(ub(z), NewSol(z)),lb(z));
        end
        % Evaluation of objective function
        NewSolObj = FITNESSFCN(NewSol);
        count=count+1;
        
        % Greedy selection
        if (NewSolObj < obj(i))
            pop(i,:) = NewSol;
            obj(i) = NewSolObj;
        end
        % ----------------Ending of the Teacher Phase for ith student-------------- %
        
        
        % ----------------Begining of the Learner Phase for ith student-------------- %
        % Generation of a new solution
        if (obj(i)< obj(Partner(i)))
            NewSol = pop(i,:) + rand(1, D).*(pop(i,:)- pop(Partner(i),:));
        else
            NewSol = pop(i,:) + rand(1, D).*(pop(Partner(i),:)- pop(i,:));
        end
        
        % Bounding of the solution
        for z=1:length(ub)
        NewSol(z) = max(min(ub(z), NewSol(z)),lb(z));
        end
        
        % Evaluation of objective function
        NewSolObj =  FITNESSFCN(NewSol);
        count=count+1;
        
        % Greedy selection
        if(NewSolObj< obj(i))
            pop(i,:) = NewSol;
            obj(i) = NewSolObj;
        end
        % ----------------Ending of the Learner Phase for ith student-------------- %
        
    end
    
    % This is not part of the algorithm but is used to keep track of the
    % best solution determined till the current iteration
    [BestFVALIter(gen),ind] = min(obj);
end
% Extracting the best solution
X = pop(ind,:);
FVAL = BestFVALIter(gen);
count