# 2D Linear Diffusion Equation
# Used Central difference scheme

using PyPlot

# PHYSICAL PARAMTERS
Lx = 2.0
Ly = 2.0
# Variable declarations
NX = 101;
NY = 101;

nu = 0.2;
dx = 2/(NX-1)
dy = 2/(NY-1);
sigma = .25;
dt = sigma*dx*dy/nu;

x = range(0, length = NX, stop = Lx)
y = range(0, length = NY, stop = Ly)

u  = ones(NY,NX); #create a 1xn vector of 1's
un = ones(NY,NX);


#Assign initial conditions
s = Int(floor(0.5/dy))
e = Int(floor(1.0/dy))
u[s:e,s:e] .= 2.0; #set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

figure(1) # Its starts has a square wave
ss1 = surf(x,y,u)
xlim(0,2)
ylim(0,2)
zlim(1,2.5);

#Run through NT timesteps
function diffusion_2D(NT)
    for n in 1:NT
        un = copy(u) ;
        u[2:end-1,2:end-1] = un[2:end-1,2:end-1] + nu*dt/dx^2*(un[2:end-1,3:end]-2*un[2:end-1,2:end-1] + un[2:end-1,1:end-2]) +
                                                   nu*dt/dy^2*(un[3:end,2:end-1]-2*un[2:end-1,2:end-1] + un[1:end-2,2:end-1])
        u[1,:]  .= 1.0
        u[end,:].= 1.0
        u[:,1]  .= 1.0
        u[:,end].= 1.0
    end

    figure()
    surf(x,y,u,cmap=ColorMap("coolwarm"))
    xlim(0,2)
    ylim(0,2)
    zlim(1,2.5);
end

# Running the function for different time steps
diffusion_2D(10)
diffusion_2D(50)
diffusion_2D(100)
