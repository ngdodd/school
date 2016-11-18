function laplace(N)
%   laplace(N) solves Laplace's equation u_xx + u_yy = 0 on the unit square
%   with (N+1)x(N+1) gridpoints i = 1,...,N+1 and j = 1,...,N+1
%   In u(i,j), i labels rows y and j labels columns x
%   i = 1 labels y = 1, i = N+1 labels y = 0
%   j = 1 labels x = 0, j = N+1 labels x = 1
%   BCs are Dirichlet with u = 0 on boundary except u(x=1,y) = 1
%   Uses Preconditioned Conjugate Gradient Method with
%   M^-1 = Minv = (I - L*Dinv)*(I - Dinv*L') due to
%   Ament, Knittel, Weiskopf, & Strasser (2010)

tic
h = 1/(N+1); % dx = dy = h
i = (1:N+1);
xcoord = (i-1)/N;
ycoord = (i-1)/N;
ufull = zeros(N+1,N+1);
ufull(i,N+1) = 1; % nonzero BC u(x=1,y) = 1
u = zeros(N-1,N-1); % interior unknown values of u
u = ufull(2:N,2:N); 
u = reshape(u,(N-1)^2,1);

b = zeros((N-1)^2,1);
j = 1:N-1;
n = (N-2)*(N-1) + j;
b(n) = b(n)-1; % East BC 

% number grid points in (N-1)*(N-1) interior
G = numgrid('S',N+1);
% calculate Laplacian
A = -delsq(G)
% Matrix_A = full(A) % Set N = 4
% Vector_b = b

L = tril(A,-1) ;
D = diag(diag(A)) ;
Dinv = inv(D) ;
Minv = (speye((N-1)^2)-L*Dinv)*(speye((N-1)^2)-Dinv*L') ;

residual = b - A*u;
z = Minv*residual;
d = z;
zr_old = z'*residual;
k = 1; EPSILON = 10^-6;
while k<length(b) && norm(residual,1)>EPSILON
    k = k+1;
    Ad = A*d;
    alpha = zr_old/(d'*Ad);
    u = u + alpha*d;
    residual = residual - alpha*Ad;
    z = Minv*residual;
    zr_new = z'*residual;
    d = z + (zr_new/zr_old)*d;
    zr_old = zr_new;
end
iterations = k
toc

u = reshape(u,N-1,N-1);
ufull(2:N,2:N) = u;

figure
surf(xcoord,ycoord,ufull)
colormap('jet')
shading flat
xlabel('x','FontSize',16); 
ylabel('y','FontSize',16); 
zlabel('u','FontSize',16)

end
