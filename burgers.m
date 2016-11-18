function burgers(N,steps,CFL,nu)
%   burgers(N,steps,CFL,nu) solves the inviscid Burgers' equation 
%       u_t + f_x = 0, f(u) = u^2/2
%   with N+1 gridpoints for 0 <= x <= 1, h = 1/N, dt = CFL*h, and
%   t_f = steps*dt 
%   Uses two-step LW with artificial viscosity, and through-flow BCs
%   NOTE: nu_A = nu h |u|, which -> 0 as h -> 0
%
%   Try: 
%        burgers(200,100,0.9,0.2) for N-wave
h = 1/N;
x = linspace(0,1,N+1); t = 0;
u = zeros(N+1,1); % N+1 dimensional column vector

% ICs: Riemann problem shock wave
%j = 1:N/2-1; u(j) = 2; 
%u(N/2) = 1.5; 
%j = N/2+1:N+1; u(j) = 1;

% ICs: Riemann problem rarefaction wave
% j = 1:N/2-1; u(j) = 1; 
% u(N/2) = 1.5; 
% j = N/2+1:N+1; u(j) = 2;

% ICs: N-wave
% j = N/4:3*N/4; u(j) = 4*(j-N/4)/N - 1;

% ICs: merging shocks
 j = 1:N/4-1; u(j) = 3; 
 u(N/4) = 2.5;
 j = N/4+1:N/2-1; u(j) = 2;
 u(N/2) = 1.5;
 j = N/2+1:N+1; u(j) = 1;

for n = 1:steps % timestep loop
    umax = max(abs(u));
    dt = CFL*h/umax; t = t + dt;
    % with thru-flow BCs
    % gridpts: 1 | 1 2 ... N N+1 | N+1
    uleft = [u(1);u(1:N)];
    uright = [u(2:N+1);u(N+1)];
    
    % LW: LF partial step
    umidright = (u+uright)/2 - dt*(uright.^2-u.^2)/(4*h);
    umidleft = (u+uleft)/2 - dt*(u.^2-uleft.^2)/(4*h);
    % LW: leapfrog partial step
    u = u - dt*(umidright.^2-umidleft.^2)/(2*h) + ...
        dt*nu*abs(u).*(uright-2*u+uleft)/h; % artificial viscosity
    
    plot(x,u,'r-')
    xlabel('x','FontSize',16)
    ylabel('u','FontSize',16)
    axis([0 1 0 3.2]);
    % axis([0 1 -1 1]); % for N-wave
    getframe;
end
tf = t

end
