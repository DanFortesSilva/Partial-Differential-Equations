# 1D Linear Diffusion Equation
# Used Central difference scheme

using PyPlot

function diffusion_1D(NX)
    # Nx = 41 ; # try changing this number from 41 to 81 and Run All ... what happens?
    dx = 2/(NX-1);
    NT = 20 ;   # NT is the number of timesteps we want to calculate
    nu = 0.3    # the value of viscosity
    sigma = 0.2 # sigma is a parameter, we'll learn more about it later
    dt = (sigma*(dx^2))/nu ; # dt is the amount of time each timestep covers (delta t)


    u = ones(Float64,NX);    # function ones()
    s = Int(floor(0.5/dx));
    e = Int(floor(1/dx));
    u[s:e] .= 2.0;
    un = ones(Float64,NX); # initialize a temporary array
    for n in 1:NT          # loop for values of n from 0 to nt, so it will run nt times
        un = copy(u)       # copy the existing values of u into un
        for i in 2:NX-1    # you can try commenting this line and...
            u[i] = un[i] + nu*dt/dx^2 * ( un[i+1] - 2*un[i] + un[i-1] );
        end
    end
    plot(range(0, length = NX ,stop = 2),u);
end

# For Plotting multiple solutions for time steps 41 to 81
for i in 41:81
    diffusion_1D(i)
end
