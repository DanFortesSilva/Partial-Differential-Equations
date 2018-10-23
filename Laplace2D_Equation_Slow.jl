# This script solves the 2D Laplace Equation using Finite Differences Method

using PyPlot, PyCall # Calling python functions and plotting function matplotlib
@pyimport numpy as num

N     = 50;        # Number of Steps in space(x) and space(y)
Niter = 2^15       # Number of Iteration
dx  = 2/(N-1)      # Width of space step (x)
dy  = 2/(N-1)      # Width of space step (y)
dx2 = dx*dx
dy2 = dy*dy

u = zeros(Float64,N,N) # initializing solution u
u[1,:] .= 0 # boundary condition u=0 at y=0
u[N,:] .= 1 # boundary condition u=1 at y=1
u[:,1] .= 0 # boundary condition u=0 at x=0
u[:,N] .= 0 # boundary condition u=1 at x=1


# Solving Laplace Equation
function Laplace2D_Slow(u::Array{Float64,2},N::Integer,Niter::Integer,dx::Float64,dx2::Float64,dy::Float64,dy2::Float64)
    for i = 1:Niter
        u[2:N-1, 2:N-1] = ((u[1:N-2, 2:N-1] + u[3:N, 2:N-1])*dy2 + (u[2:N-1,1:N-2] + u[2:N-1, 3:N])*dx2) * (1/(2*(dx2+dy2)));
    end
end

@time Laplace2D_Slow(u,N,Niter,dx,dx2,dy,dy2)

# Creating a meshgrid to plot solution u
x = range(0, length = N, stop = 1)
y = range(0, length = N, stop = 1)
X,Y = num.meshgrid(x,y)

# Plotting solution using PyPlot.jl functions
figure(1)
su = surf(X,Y,u, cmap=ColorMap("jet"))
xlabel("X")
ylabel("Y")
zlabel("T")
colorbar(su)

figure(2)
cu = contourf(X,Y,u, cmap=ColorMap("jet"))
xlabel("X")
ylabel("Y")
colorbar(cu)
