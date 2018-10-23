# 2D Nonlinear Convection (Advection) (Coupled) Equations
# Timestep will be discretized as a forward difference and both spatial steps will be discretized as backward differences

using PyPlot

# PHYSICAL PARAMTERS
Lx = 2.0
Ly = 2.0
# NUMERICAL PARAMETERS
NT = 80
NX = 101
NY = 101
c = 1.0
dx = Lx/(NX-1)
dy = Ly/(NY-1)
sigma = 0.2
dt = sigma*dx

x = range(0, length = NX, stop = Lx)
y = range(0, length = NY, stop = Ly)

u  = ones(Float64,NY,NX); #create a 1xn vector of 1's
v  = ones(Float64,NY,NX);
un = ones(Float64,NY,NX);
vn = ones(Float64,NY,NX);

#Assign initial conditions
s = Int(0.5/dy);
e = Int(1/dy);
u[s:e,s:e] .= 2.0 #set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
v[s:e,s:e] .= 2.0 #set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

figure(1)
ss1 = surf(x,y,u,cmap=ColorMap("coolwarm"))
xlabel("X [m]")
ylabel("Y [m]")
zlabel("u [m/s]")
figure(2)
ss2 = surf(x,y,v,cmap=ColorMap("coolwarm"))
xlabel("X [m]")
ylabel("Y [m]")
zlabel("v [m/s]")

function Nonlinear_Conv2D()
    for n in 1:NT #loop across number of time steps. This uses array operations (is fast)
        un = copy(u)
        vn = copy(v)

        u[2:end,2:end] = un[2:end,2:end]-(un[2:end,2:end].*(c*dt/dx*(un[2:end,2:end]-un[2:end,1:end-1])))-(vn[2:end,2:end].*(c*dt/dy*(un[2:end,2:end]-un[1:end-1,2:end])))
        v[2:end,2:end] = vn[2:end,2:end]-(un[2:end,2:end].*(c*dt/dx*(vn[2:end,2:end]-vn[2:end,1:end-1])))-(vn[2:end,2:end].*(c*dt/dy*(vn[2:end,2:end]-vn[1:end-1,2:end])))

        u[1,:]   .= 1.0;
        u[end,:] .= 1.0;
        u[:,1]   .= 1.0;
        u[:,end] .= 1.0;

        v[1,:]   .= 1.0;
        v[end,:] .= 1.0;
        v[:,1]   .= 1.0;
        v[:,end] .= 1.0;
    end
end

@time Nonlinear_Conv2D()

figure(3)
ss3 = surf(x,y,u,cmap=ColorMap("coolwarm"))
xlabel("X [m]")
ylabel("Y [m]")
zlabel("u [m/s]")
figure(4)
ss4 = surf(x,y,v,cmap=ColorMap("coolwarm"))
xlabel("X [m]")
ylabel("Y [m]")
zlabel("v [m/s]")
