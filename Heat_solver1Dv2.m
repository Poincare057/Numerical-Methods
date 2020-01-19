#Heat Equation Solver 2 with Crank Nicolson
xmesh = 200
tmesh = 50
x = linspace(0, 1, xmesh);
t = linspace(0, 2, tmesh);
v = [;]
v(1, 1:xmesh) = sin(3.1415926*x);       #inital condition along t = 0
b = (1/(3.1415926))^2
h = 1/xmesh
dt = 2/tmesh
l = dt/(h^2)
f = b*l
disp(f)
v(2:tmesh, 1) = 0                  #boundary condition along x = 0
A = [;]
B = [;]
for i = 1:xmesh - 2
  A(i,i) = 1 + 2*f;
endfor  
for i = 2:xmesh - 3
  A(i, i-1) = -f;
endfor
for i = 1:xmesh - 3
  A(i, i+1) = -f;
endfor
for i = 1:xmesh - 2
  B(i,i) = 1 - 2*f;
endfor  
for i = 2:xmesh - 3
  B(i, i-1) = f;
endfor
for i = 1:xmesh - 3
  B(i, i+1) = f;
endfor

C = inv(A)*B;
for n = 1:tmesh - 1
  v(n+1, 2:xmesh-1) = C*transpose(v(n, 2:xmesh - 1));
  v(n+1, xmesh) = v(n+1, xmesh - 1);
  v(n+1, 1) = v(n+1, 2);
endfor

surf(x, t, v(1:tmesh, 1:xmesh))