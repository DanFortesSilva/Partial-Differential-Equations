# This script solves the 2D Poisson Equation using Finite Differences Method

using PyPlot # Calling python functions and plotting function matplotlib

N     = 50;        # Number of Steps in space(x) and space(y)
Niter = 2^13       # Number of Iteration
dx  = 2/(N-1)      # Width of space step (x)
dy  = 2/(N-1)      # Width of space step (y)
dx2 = dx*dx
dy2 = dy*dy

u = zeros(N,N) # initializing solution u
u[1,:] .= 0 # boundary condition u=0 at y=0
u[N,:] .= 0 # boundary condition u=1 at y=1
u[:,1] .= 0 # boundary condition u=0 at x=0
u[:,N] .= 0 # boundary condition u=1 at x=1

# Source b
b = zeros(N,N) # initializing source b
b[Int(floor(N/4)),Int(floor(N/4))]  = 10; # changing location and value of source
b[Int(floor(3*N/4)),Int(floor(3*N/4))]  = -10; # changing location and value of sink

# Solving Poisson Equation with central differences discritization

function Poisson2D()
    for i = 1:Niter
        u[2:N-1, 2:N-1] = ((u[1:N-2, 2:N-1] + u[3:N, 2:N-1])*dy2 + (u[2:N-1,1:N-2] + u[2:N-1, 3:N])*dx2 + b[2:N-1,2:N-1]*dx2*dy2) * (1/(2*(dx2+dy2)));
    end
end

Poisson2D()

# Creating a meshgrid to plot solution u
x = range(0, length = N, stop = 1)
y = range(0, length = N, stop = 1)

# Plotting solution using PyPlot.jl functions
figure(1)
su = surf(x,y,u, cmap=ColorMap("jet"))
xlabel("X")
ylabel("Y")
zlabel("T")
colorbar(su)

figure(2)
cu = contourf(x,y,u, cmap=ColorMap("jet"))
xlabel("X")
ylabel("Y")
colorbar(cu)
