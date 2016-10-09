function c = ConvectionExactSolution(x, y, time)
%% parameters
u  = 1.0;
v  = 1.0;
Dx = 1e-2;
Dy = 1e-2;
x0 = .5;
y0 = .5;

cx = exp( -(x-x0-u*time).^2 /Dx/(4*time+1) );
cy = exp( -(y-y0-v*time).^2 /Dy/(4*time+1) );

c  = cx.*cy/(4*time+1);
