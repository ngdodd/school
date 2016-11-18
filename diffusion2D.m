function diffusion2D(N,tf)
%   diffusion2D(N,tf) solves the 2D diffusion (heat) equation 
%     u_t = u_xx + u_yy
%   on the unit square with (N+1)*(N+1) gridpoints (1,...,N+1)^2 to time tf
%   with homogeneous Dirichlet BCs and a Gaussian IC.
%   Uses backward Euler or TRBDF2 method and fixed dt for movie
% 
%   Try: diffusion2D(100,0.01)

tic
GAMMA = 2-sqrt(2);
h = 1/N;
k = (1:N+1)';
x = (k-1)*h; y = (k-1)*h;  
% dt = h^2/2; % forward Euler stability limit
dt = 10*h^2/2; % to watch movie; dt could be set much larger
steps = tf/dt+1;

u = zeros(N-1,N-1); umid = zeros((N-1)^2,1); % interior
U = zeros(N+1,N+1); % full region
i = 1:N-1; j = 1:N-1;
u(i,j) = exp(-40*((y(i+1)-0.5).^2))*exp(-40*((x(j+1)'-0.5).^2));

% number grid points in (N-1)*(N-1) interior
G = numgrid('S',N+1);
% calculate Laplacian D2
D2 = -delsq(G);
%full(D2) % Set N = 4

CONST = GAMMA/(2*h^2);
CONST1 = (1-GAMMA)/((2-GAMMA)*h^2);
CONST2 = 1/(GAMMA*(2-GAMMA));
CONST3 = (1-GAMMA)^2/(GAMMA*(2-GAMMA));
for n = 1:steps % timestep loop
    u = reshape(u,(N-1)^2,1);
    % u = (speye((N-1)^2,(N-1)^2)-dt*D2/h^2)\u; % backward Euler
    umid = (speye((N-1)^2,(N-1)^2)-CONST*dt*D2)\...
        ((speye((N-1)^2,(N-1)^2)+CONST*dt*D2)*u); % TR
    u = (speye((N-1)^2,(N-1)^2)-CONST1*dt*D2)\(CONST2*umid-CONST3*u); % BDF2
    u = reshape(u,N-1,N-1);
    
    % For Movie:
    U(i+1,j+1) = u(i,j);
    surf(x,y,U,'FaceColor','interp','EdgeColor','none','FaceLighting',...
        'gouraud')
    view(-50,30);
    camlight left;
    xlabel('x','FontSize',14); ylabel('y','FontSize',14);
    title('u','FontSize',16);
    axis([0 1 0 1 0 1.1]);
    M(n) = getframe(gcf);
end
toc
end
