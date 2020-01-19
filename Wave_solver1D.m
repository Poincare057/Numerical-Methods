#1st order Wave Equation Solvers
#1. Godunov Scheme
##a = 4
##v = []
##x = linspace(0, 1, 10)
##t = linspace(0, 2, 150)
##h = 1/10
##dt = 2/150
##l = dt/h
##v(1:10, 1) = sin(3.1416*x) + cos(3.1416*x);
##for n = 1:150
##  for i = 2:9
##    v(i, n+1) = v(i, n) - l*((1 + sign(a))/2)*a*(v(i,n) - v(i-1, n)) - l*((1 - sign(a))/2)*a*(v(i+1,n) - v(i, n)) ;
##endfor;
##endfor;
##surf(x, t, transpose(v(1:10, 1:150)))

#2. Lax-Friedrichs Scheme
##a = 4
##v = []
##x = linspace(0, 1, 10)
##t = linspace(0, 2, 150)
##h = 1/10
##dt = 2/150
##l = dt/h
##v(1:10, 1) = sin(3.1416*x) + cos(3.1416*x);
##v(1, 2:150) = 0#1 - sin(t(2:150))/100
##v(10, 2:150) = 0#1 - sin(t(2:150))/100
##for n = 1:150
##  for i = 2:9
##    v(i, n+1) = 0.5*(1 - a*l)*v(i+1, n) + 0.5*(1+a*l)*v(i-1, n);
##endfor;
##endfor;
##surf(x, t, transpose(v(1:10, 1:150)))

###3. Lax-Wendroff Scheme
##a = 4
##v = []
##x = linspace(0, 1, 10);
##t = linspace(0, 2, 150);
##h = 1/10
##dt = 2/150
##l = dt/h
##v(1:10, 1) = cos(3*3.1416*x);
##v(1, 2:150) = exp(-t(2:150));
##v(10, 2:150) = exp(-t(2:150));
##for n = 1:149
##  for i = 2:9
##    v(i, n+1) = 0.5*a*l*(v(i+1, n) - v(i-1, n))+ 0.5*(a^2)*(l^2)*(v(i+1, n) - 2*v(i, n) + v(i-1, n));
##endfor;
##endfor;
##surf(x, t, transpose(v(1:10, 1:150)))

#4. Crank-Nicolson Scheme
##a = 4
##v = [;]
##x = linspace(-10, 10, 1000)
##t = linspace(0, 2, 2000)
##h = 20/1000
##dt = 2/2000
##l = dt/(h)
##v(1:500 , 1) = 1 ;
##v(501:550, 1) = 2*(x(501:550).^3) - 3*(x(501:550).^2) + 1;
##v(551:1000, 1) = 0;
##plot(x, v)
##
##A = [;]
##B = [;]
##for i = 1:1000
##  A(i,i) = 1;
##endfor  
##for i = 2:1000
##  A(i, i-1) = -a*l*0.25;
##endfor
##for i = 1:999
##  A(i, i+1) = a*l*0.25;
##endfor
##for i = 1:1000
##  B(i,i) = 1;
##endfor  
##for i = 2:1000
##  B(i, i-1) = a*l*0.25;
##endfor
##for i = 1:999
##  B(i, i+1) = -a*l*0.25;
##endfor
##
##C = inv(A)*B;
##for n = 1:1999
##  v(1:1000, n+1) = C*(v(1:1000, n));
##endfor
##
##for p = 1:2000
##  plot(x, v(1:1000, p))
##  hold on
##endfor

###5. Backward-Time Forward-Space Scheme   (information moves from right to left)              #l2 stable for |al| >= 1
##xmesh = 200
##tmesh = 100
##x = linspace(-10, 10, xmesh)
##t = linspace(0, 10, tmesh)
##a = 1
##dt = 10/tmesh
##h = 20/xmesh
##l = a*dt/h
##v = [;]
##v(1:xmesh/2, 1) = 1;
##v((xmesh/2)+1: int32((11/20)*xmesh/1 + 1), 1) = 2*(x((xmesh/2)+1: int32((11/20)*xmesh/1 + 1)).^3) - 3*(x((xmesh/2)+1: int32((11/20)*xmesh/1 + 1)).^2) + 1;
##v(int32((11/20)*xmesh/1 + 1):xmesh, 1) = 0;
##v(1, 1:tmesh) = 1;
##A = [;]
##for i = 1:tmesh
##  A(i,i) = a*l;
##endfor
##
##for i = 1:tmesh-1
##  A(i,i+1) = 1- a*l;
##endfor
##B = A^(-1);
##for ti = 2:tmesh
##  C = B*(v(1:xmesh, ti-1))
##  v(1:xmesh, ti) = C(1:xmesh);
##endfor
##
##
##for ti = 1:10:50
##  plot(x, v(1:xmesh, ti))
##  hold on
##endfor

#surf(x,t,transpose(v(1:tmesh, 1:xmesh)))


#6. Backward-Space Forward-Time Scheme          (information moves from left to right)
xmesh = 150
tmesh = 100
x = linspace(-10, 10, xmesh);
t = linspace(0, 10, tmesh);
a = 1
dt = 10/tmesh
h = 20/xmesh
l = dt/h
v = [;]
v(1:xmesh/2, 1) = 1;
v((xmesh/2)+1: int32((11/20)*xmesh/1 + 1), 1) = 2*(x((xmesh/2)+1: int32((11/20)*xmesh/1 + 1)).^3) - 3*(x((xmesh/2)+1: int32((11/20)*xmesh/1 + 1)).^2) + 1;
v(int32((11/20)*xmesh/1 + 1):xmesh, 1) = 0;
v(1, 1:tmesh) = 1;
A = [;]
for i = 1:xmesh
  A(i,i) = a*l;
endfor
for i = 1:xmesh-1
  A(i,i+1) = 1 - a*l;
endfor
for ti = 1:tmesh-1
  B = A*(v(1:xmesh, ti));
  v(2:xmesh, ti+1) = B(1:xmesh-1);
  v(1, ti+1) = 1;
endfor
for ti = 1:10:50
  plot(x, v(1:xmesh, ti))
  hold on
endfor

# Analytic Solution
a = 1

for i = 1:xmesh
  for j = 1:tmesh
    if x(i) <= t(j)
      k(i, j) = 1;
    elseif x(i) - t(j) < 1
      k(i,j) = 2*(x(i) - t(j))^3 - 3*(x(i) - t(j))^2 + 1;
    elseif x(i) >= t(j) + 1
      k(i,j) = 0;
    endif;
  endfor;
endfor;
for ti = 1:10:50
  plot(x, k(1:xmesh, ti))
  hold on
endfor
