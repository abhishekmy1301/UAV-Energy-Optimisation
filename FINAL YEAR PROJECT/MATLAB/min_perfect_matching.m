function [ indices, cost ] = min_perfect_matching( G )
% Minimum Perfect Matching Tool
%
% Synopsis
%   [indices, cost] = min_perfect_matching( G ) 
%   
% Description
%   Function to solve the Minimum Perfect Matching on non-biparite graphs
%   problem using Integer linear programming.
%
%   Returns vector of matched indices and cost of the match. Requires
%   symmetric adjacent matrix of even rank.
% 
% Author: Vojtech Knyttl, knyttvoj@fel.cvut.cz

    %% dimensions
    if( ~isequal(G, G') )
      e = MException('PerfectMatching:MxNotSymetric', 'Matrix is not symmetric.' );
      throw(e);
    end;

    len = uint16(length( G ));
    if( mod(len,2) == 1 )
      e = MException('PerfectMatching:MxNotEvenRank', 'Matrix must be of even rank.' );
      throw(e);
    end;
    
    llen = (len*len)-len;
 
    %% function to minimize
    f = zeros(1,llen);
    for i=1:len-1
      f((i-1)*len+1:i*len)=G(i,:);
    end;
       
    %% respect matrix
    A = zeros(len,llen);
    for i=1:len
      for j=1:len*len-len
        idiv = idivide(j-1,len,'floor')+1;
        if idiv == i
          if and(mod(j-1,len)>=i,mod(j,len)~=i)
            A(i,j) = 1;
          end;
        elseif idiv < i
          if mod(j-1,len)+1==i
            A(i,j) = 1;  
          end;
        end;
      end;
    end;
    b = ones( len, 1 );
    
    %% removing zero columns
    remove_cols = find(all(A==0));
    f(:,remove_cols)=[];
    A(:,remove_cols)=[];
    
    %% glpk ilinprog settings
    sense=1;                            % minimization
    ctype = repmat( 'S', 1, len );      % equalities
    lb = zeros( llen/2, 1 );            % lower bound
    ub = ones( llen/2, 1 );             % upper bound
    i = 1:llen/2; 
    e=2^-24;
    % vartype = repmat( 'I', 1, llen/2 ); % integral vals
   % param.msglev = 1;  
   % param.itlim = 100;
    
    xmin = IP1(f,[],[],A,b,lb,ub,i,e);

    %xmin = glpk(f,A,b,lb,ub,ctype,vartype,sense,param);
    
    %% adding remove columns
    match = find(xmin==1)';
    
    for i=remove_cols
      match(match>=i)=match(match>=i)+1;
    end;
    
    %% forming the result
    indices = zeros(1,len); cost = 0;
    for i=match
      x = idivide(i-1,len,'floor')+1;
      y = mod(i-1,len)+1;
      indices(x)=y;
      indices(y)=x;
      cost = cost + G(x,y);
    end;
end

% By Sherif A. Tawfik, Faculty of Engineering, Cairo University
% [x,val,status]=IP1(f,A,b,Aeq,beq,lb,ub,M,e)
% this function solves the following mixed-integer linear programming problem
%   min f*x
%  subject to
%        A*x <=b
%        Aeq * x = beq   
%        lb <= x <= ub
%        M is a vector of indeces for the variables that are constrained to be integers
%        e is the integarilty tolerance
% the return variables are :
% x : the solution
% val: value of the objective function at the optimal solution
% status =1 if successful
%        =0 if maximum number of iterations reached in he linprog function
%        =-1 if there is no solution
% Example:
%        maximize 17 x1 + 12 x2 
%        subject to
% 	             10 x1 + 7 x2 <=40
%                   x1 +   x2 <= 5
%                   x1, x2 >=0 and are integers
% f=[-17, -12]; %take the negative for maximization problems
% A=[ 10  7; 1 1];
% B=[40; 5];
% lb=[0 0];
% ub=[inf inf];
% M=[1,2];
% e=2^-24;
% [x v s]= IP(f,A,B,[],[],lb,ub,M,e)
function [x,val,status]=IP1(f,A,b,Aeq,beq,lb,ub,M,e)
options = optimset('display','off');
bound=inf; % the initial bound is set to +ve infinity
[x0,val0]=linprog(f,A,b,Aeq,beq,lb,ub,[],options); 
[x,val,status,b]=rec(f,A,b,Aeq,beq,lb,ub,x0,val0,M,e,bound);
end % a recursive function that processes the BB tree 
function [xx,val,status,bb]=rec(f,A,b,Aeq,beq,lb,ub,x,v,M,e,bound) 
options = optimset('display','off');
% x is an initial solution and v is the corressponding objective function value
% solve the corresponding LP model with the integarily constraints removed
[x0,val0,status0]=linprog(f,A,b,Aeq,beq,lb,ub,[],options); 

% if the solution is not feasible or the value of the objective function is
% higher than the current bound return with the input intial solution
if status0<=0 | val0 > bound  
    xx=x; val=v; status=status0; bb=bound;
    return;
end
% if the integer-constraint variables turned to be integers within the
% input tolerance return
ind=find( abs(x0(M)-round(x0(M)))>e ); 
if isempty(ind)
    status=1;        
    if val0 < bound    % this solution is better than the current solution hence replace
        x0(M)=round(x0(M));
        xx=x0;        
        val=val0;
        bb=val0;
    else
        xx=x;  % return the input solution
        val=v;
        bb=bound;
    end
    return
end
% if we come here this means that the solution of the LP relaxation is
% feasible and gives a less value than the current bound but some of the
% integer-constraint variables are not integers. 
% Therefore we pick the first one that is not integer and form two LP problems
% and solve them recursively by calling the same function (branching)
% first LP problem with the added constraint that Xi <= floor(Xi) , i=ind(1)
br_var=M(ind(1));
br_value=x(br_var);
if isempty(A)
    [r c]=size(Aeq);
else
    [r c]=size(A);
end
A1=[A ; zeros(1,c)];
A1(end,br_var)=1;
b1=[b;floor(br_value)];
% second LP problem with the added constraint that Xi >= ceil(Xi) , i=ind(1)
A2=[A ;zeros(1,c)];
A2(end,br_var)=-1;
b2=[b; -ceil(br_value)];
% solve the first LP problem
[x1,val1,status1,bound1]=rec(f,A1,b1,Aeq,beq,lb,ub,x0,val0,M,e,bound);
status=status1;
if status1 >0 & bound1<bound % if the solution was successfull and gives a better bound
   xx=x1;
   val=val1;
   bound=bound1;
   bb=bound1;
else
    xx=x0;
    val=val0;
    bb=bound;
end
    
% solve the second LP problem
[x2,val2,status2,bound2]=rec(f,A2,b2,Aeq,beq,lb,ub,x0,val0,M,e,bound);
if status2 >0 & bound2<bound % if the solution was successfull and gives a better bound
    status=status2;
    xx=x2;
    val=val2;
    bb=bound2;
end
end
