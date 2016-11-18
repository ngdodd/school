function wave(N,steps,CFL)
%   wave3(N,steps,CFL) solves the wave equation w_tt - w_xx = 0 as a 
%   first-order system (u = w_x, v = w_t):
%      u_t = v_x 
%      v_t = u_x
%   with N+1 gridpoints for 0 <= x <= 1, h = 1/N, dt = CFL*h, and
%   t_f = steps*dt
%   For the ICs, N >= 20 should be an integer multiple of 10
%   Uses two-step LW method and periodic BCs
%
%   Try: wave(200,280,1)

h = 1/N;
dt = CFL*h;
u = zeros(N+1,1); v = zeros(N+1,1); % N+1 dimensional column vectors
j = N/2-N/20:N/2+N/20;
u(j+1) = -cos(10*j*pi/N); % ICs
x = linspace(0,1,N+1);
for n = 1:steps % timestep loop
    % with periodic BCs
    % gridpts: ... N+1 | 1 2 ... N N+1 | 1 ...
    uleft = [u(N+1);u(1:N)];
    uright = [u(2:N+1);u(1)];
    vleft = [v(N+1);v(1:N)];
    vright = [v(2:N+1);v(1)];
    % LF partial step
    umidright = (u+uright)/2 + dt*(vright-v)/(2*h);
    vmidright = (v+vright)/2 + dt*(uright-u)/(2*h);
    umidleft = (u+uleft)/2 + dt*(v-vleft)/(2*h);
    vmidleft = (v+vleft)/2 + dt*(u-uleft)/(2*h);
    % leapfrog partial step
    u = u + dt*(vmidright-vmidleft)/h;
    v = v + dt*(umidright-umidleft)/h;
    plot(x,u,'r-')
    xlabel('x','FontSize',16)
    ylabel('u','FontSize',16)
    axis([0 1 0 1]);
    getframe;
end

end
