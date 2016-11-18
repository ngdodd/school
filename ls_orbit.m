%% Planetary orbit
% The expression z = ax^2 + bxy + cy^2 + dx + ey + f is known
% as a quadratic form. The set of points (x,y) where z = 0, 
% is a conic section. It can be an ellipse, a parabola, or a 
% hyperbola, depending on the sign of the discriminant 
% b^2 - 4ac. Circles and lines are special cases. The equation 
% z = 0 can be normalized by dividing the quadratic form by any
% nonzero coefficient. For example, if f ~= 0, we can divide 
% all the other coefficients by f and obtain a quadratic form 
% with the constant term equal to one.
function ls_orbit()

% A planet follows an elliptical orbit. The following are ten
% observations of its position in the xy-plane.
x = [1.02 .95 .87 .77 .67 .56 .44 .30 .16 .01] ;
y = [0.39 .32 .27 .22 .18 .15 .13 .12 .13 .15] ;

%figure('units','normalized','outerposition',[0 0 1 1])
fprintf('Least squares fit of original orbit measurements: \n' ) ;
[ho1,hm1,coeff1] = planetary_orbit( x, y, 'b', 'k*' ) ;
hold on

% This least squares problem is nearly rank deficient. To see
% what effect this has on the solution, perturb the data 
% slightly by adding to each coordinate of each data point a 
% random number uniformly distributed in the interval 
% [-0.0005,0.0005]. Compute the new coefficients resulting from
% the perturbed data. Plot the new orbit on the same plot with
% the old orbit.
x_noise = -0.0005 + (0.0005-(-0.0005)).*rand(1,length(x)) ;
y_noise = -0.0005 + (0.0005-(-0.0005)).*rand(1,length(y)) ;

fprintf('Least squares fit of perturbed orbit measurements: \n' ) ;
[ho2,hm2,coeff2] = planetary_orbit( x + x_noise, y + y_noise, 'g', 'rx' ) ;
legend( [ho1,ho2,hm1,hm2], 'Orbit Fit', 'Perturbed Orbit Fit', 'Measurements', 'Perturbed Measurements' ) ;

% Relative error of fitted orbit
fprintf( 'Relative error of fitted coefficients using 2-norm: ' )
fprintf( '%f\n', norm( (coeff1-coeff2), 2)/norm(coeff1, 2) ) ;

fprintf( 'Relative error of fitted coefficients using infinity norm: ' )
fprintf( '%f\n', norm( (coeff1-coeff2), Inf)/norm(coeff1, Inf) ) ;
end

function [orbit,meas,coeff]=planetary_orbit( x, y, orbit_plt_options, data_plt_options )
% Create the linear system from the observation data
b = -1*ones( length(x), 1 ) ;
A = [] ;

% Use quadratic form and orbit data to create each row of A
for i = 1:length(x)
    % Form the ith row of A
    A = [ A ; [ x(i).^2, x(i)*y(i), y(i).^2, x(i), y(i) ] ] ;
end

% Determine the coefficients in the quadratic form that fits
% these data in the least squares sense by setting one of the 
% coefficients equal to one and solving a 10x5 overdetermined 
% system of linear equations for the other five coefficients. 
% Use the method of normal equations or the method based on 
% QR factorization. Plot the orbit with x on the x-axis and y 
% on the y-axis. Superimpose the ten data points on the plot.

% Compute the Cholesky factorization of A'A
G = Cholesky_Factorization( A'*A ) ;
% Solve the lower-triangular system Gw = A'b
w = Forward_Substitution( G, A'*b ) ;
% Solve G'x = w, where x is the vector of coefficients that 
% provides the fit
coeff = Back_Substitution( G', w )

% Define x-axis dicretization parameters
x_min = -0.75 ;
x_max = 1.5 ;
delta_x = 0.1 ;

% Define y-axis discretization parameters
y_min = 0 ;
y_max = 1.5 ;
delta_y = 0.1 ;

% Discretize 2-D grid
[X,Y] = meshgrid( x_min:delta_x:x_max, y_min:delta_y:y_max ) ;
Z = coeff(1)*X.^2 + coeff(2)*X.*Y + coeff(3)*Y.^2 + coeff(4)*X + coeff(5)*Y + 1 ;
fprintf( 'Elliptical orbit fitted via method of normal equations: \n' )
fprintf( '%fx^2 + %fxy + %fy^2 + %fx + %fy + 1 = 0\n\n', coeff(1), coeff(2), coeff(3), coeff(4), coeff(5) ) ;
[~,orbit] = contour( X, Y, Z, [0,0], orbit_plt_options ) ;
hold on
meas = scatter( x, y, data_plt_options ) ;
end


%% Cholesky Factorization for square matrix

% Given a symmetric positive definite matrix A, the following 
% algorithm computes a lower triangular G such that A = GG^T.
%
% Inputs
% ------------------------------------------------------------------------
% A: Square linear system to factorize.
%
% Outputs
% ------------------------------------------------------------------------
% At the end of the algorithm, the unique lower triangular 
% matrix G such that G*G^T = A will be returned.
function G=Cholesky_Factorization( A )
G = A ;
[~,n] = size( A ) ;

for j = 1:n
    if j > 1
        G(j:n,j) = G(j:n,j) - G(j:n,1:j-1)*G(j,1:j-1)' ;
    end
    
    G(j:n,j) = G(j:n,j)./sqrt(G(j,j)) ;
end

G = tril(G) ;
end

%% Forward substitution to solve the system Ly = b


% Computes the forward substitution step in solving the 
% linear system via LU decompostion.
%--------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------
% L: Unit lower triangular system of equations being solved.
% b: Vector of length n, where n is the number of rows and 
% columns of L.
%
% Outputs
% y: Vector of length n, such that y = inv(L)*b
% ------------------------------------------------------------------------
function y=Forward_Substitution(L, b)
%L is assumed to be nonsingular
[n,~] = size(b) ;
y = b ;
for j = 1:n-1
    y(j) = y(j)/L(j,j) ;
    y(j+1:n) = y(j+1:n) - y(j)*L(j + 1:n,j) ;
end
y(n) = y(n)/L(n,n) ;
end


%% Back Substitution to solve the system Ux = y
%------------------------------------------------------------------------

% Computes the back substitution step in solving the 
% linear system via LU decompostion.
%--------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------
% U: Upper triangular system of equations being solved.
% y: Vector of length n, solved for in the forward 
% substitution step.
%
% Outputs
% x: Solution to the original system Ax = b.
% ------------------------------------------------------------------------
function x=Back_Substitution(U,y)
%U is assumed to be nonsignular

[n,~] = size(y) ;
x = y ;
for j = n:-1:2
    x(j) = x(j)/U(j,j) ;
    x(1:j-1) = x(1:j-1) - x(j)*U(1:j-1,j) ;
end
x(1) = x(1)/U(1,1) ;
end