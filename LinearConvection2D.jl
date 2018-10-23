# 2D Linear Convection (Advection) (Coupled) Equation
# Timestep will be discretized as a forward difference and both spatial steps will be discretized as backward differences
using PyPlot

# PHYSICAL PARAMTERS
Lx = 2.0
Ly = 2.0
# NUMERICAL PARAMETERS
NT = 100
NX = 81
NY = 81
c = 1.0 # constant velocity
dx = Lx/(NX-1)
dy = Ly/(NY-1)
sigma = 0.2
dt = sigma*dx

x = range(0, length = NX, stop = Lx)
y = range(0, length = NY, stop = Ly)

u  = ones(Float64,NY,NX) # Create a 1xn vector of 1's
un = ones(Float64,NY,NX) #
v  = ones(Float64,NY,NX) # Create a 1xn vector of 1's
vn = ones(Float64,NY,NX) #

###Assign initial conditions
u[Int(0.5/dy):Int(1/dy),Int(0.5/dx):Int(1/dx)] .= 2.0 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
v[Int(0.5/dy):Int(1/dy),Int(0.5/dx):Int(1/dx)] .= 2.0 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

figure(1)
s1 = surf(x,y,u,cmap=ColorMap("coolwarm"))
xlabel("X [m]")
ylabel("Y [m]")
zlabel("u [m/s]")

# Nested for loops
function main_loop(u::Array{Float64,2}, NT::Integer, dx::Float64, dy::Float64, c::Float64)
    for n in 0:NT ##loop across number of time steps
        un = copy(u)
        row, col = size(u)
        for j in 2:row-1
            for i in 2:col-1
                u[j,i] = un[j, i] - (c*dt/dx*(un[j,i] - un[j,i-1]))-(c*dt/dy*(un[j,i]-un[j-1,i]))
                u[1,:]   .= 1.0
                u[end,:] .= 1.0
                u[:,1]   .= 1.0
                u[:,end] .= 1.0
            end
        end
    end
end

@time main_loop(u, NT, dx, dy, c)

figure(2)
s2 = surf(x,y,u,cmap=ColorMap("coolwarm"))
xlabel("X [m]")
ylabel("Y [m]")
zlabel("u [m/s]")

################################################
# Array Operations using the same problem but instead of nested loops we use array operations
u = ones(Float64,NY,NX)
u[Int(0.5/dy):Int(1/dy),Int(0.5/dx):Int(1/dx)] .= 2.0 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
v[Int(0.5/dy):Int(1/dy),Int(0.5/dx):Int(1/dx)] .= 2.0 ##set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

function Linear_Conv2D() # This uses array operations. Is much faster
    for n in 1:NT #loop across number of time steps

        un = copy(u)
        vn = copy(v)

        u[2:end,2:end]=un[2:end,2:end]-(c*dt/dx*(un[2:end,2:end]-un[2:end,1:end-1]))-(c*dt/dy*(un[2:end,2:end]-un[1:end-1,2:end]))
        v[2:end,2:end]=vn[2:end,2:end]-(c*dt/dx*(vn[2:end,2:end]-vn[2:end,1:end-1]))-(c*dt/dy*(vn[2:end,2:end]-vn[1:end-1,2:end]))


        u[1,:]   .= 1.0
        u[end,:] .= 1.0
        u[:,1]   .= 1.0
        u[:,end] .= 1.0

        v[1,:]   .= 1.0
        v[end,:] .= 1.0
        v[:,1]   .= 1.0
        v[:,end] .= 1.0
    end
end

@time Linear_Conv2D()

figure(3)
s3 = surf(x,y,u,cmap=ColorMap("coolwarm"))
xlabel("X [m]")
ylabel("Y [m]")
zlabel("u [m/s]")

figure(4)
s4 = surf(x,y,v,cmap=ColorMap("coolwarm"))
xlabel("X [m]")
ylabel("Y [m]")
zlabel("v [m/s]")
