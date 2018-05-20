function x = thomas(A,f) 
% solves the problem Ax = f where A is an n x n tridiagonal (and hence
% square) matrix; f is a column vector. Source: Datta, Linear Algebra and
% its Applications

cmat = diag(A,-1); %the lower diagonal
c = [0; cmat]; %in the book actually it says [c2,...,cm] 

a = diag(A,0); %the main diagonal

bmat = diag(A,1); %the upper diagonal
b = [bmat; 0]; %because in the book it says [b1,...,b_{n-1}

%if A is not square we want the code to crash 
[m,n] = size(A); 
if m~=n
    x = NaN;
    disp('A is not square')
    return 
end

u = zeros(1,n); 
l = zeros(1,n); 

u(1) = a(1); 
%Ax = LUx
for ii = 2:n 
    l(ii) = c(ii)/u(ii-1);
    u(ii) = a(ii) - l(ii) * b(ii-1); 
end

%now need to solve the problem Ax = LUx = f. So we first solve Ly = f. Note
%that this problem is easy since L = diag(l(2:n),-1) + diag(ones(1,n))

y = zeros(1,n); 
y(1) = f(1); 
for ii = 2:n
    y(ii) = f(ii) - l(ii)*y(ii-1); 
end

%now we need to solve the problem Ux = y. Note that this problem is easy
%since U = diag(u) + diag(b(1:n-1),-1)

x = zeros(n,1); %note that we want this to be a column vector 
x(n) = y(n)/u(n); 
for ii = n-1:-1:1
    x(ii) = (y(ii) - b(ii)*x(ii+1))/u(ii); 
end