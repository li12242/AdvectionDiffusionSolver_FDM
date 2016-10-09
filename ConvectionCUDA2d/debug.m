Dx = 1e-2; 
Dy = 1e-2;
u  = 1.0;
v  = 1.0;

CFL = 0.1;
np = 201;
dx = 2/(np-1);

dt = CFL*min(dx/u, dx*dx/Dx);
finalTime = 0.5;
nt = floor(0.5/dt);

x = load('x.txt');
y = load('y.txt');
L2 = zeros(nt, 1);
for i=1:nt
    filename = ['result-', num2str(i), '.txt'];
    c = load(filename);
    time  = i*dt;
    c_ext = ConvectionExactSolution(x, y, time);
    L2(i) = sqrt(sum(sum((c - c_ext).^2))/np/np);
    fprintf('time step %d, L2 error %f\n', i, L2(i));
end