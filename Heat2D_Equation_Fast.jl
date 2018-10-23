# This script solves the 2D Laplace Equation using Finite Differences Method

using PyPlot

# PHYSICAL PARAMETERS
A  = 0.25  # Diffusion coefficient
Lx = 1.0   # Domain size x
Ly = 1.0   # Domain size y

# NUMERICAL PARAMETERS
NT = 5000     #Number of time steps
NX = 100        #Number of grid points in x
NY = 100         #Number of grid points in y

dx  = Lx/(NX-1)   #Grid step in x (space)
dx2 = dx^2
dy  = Ly/(NY-1)   #Grid step in y (space)
dy2 = dy^2
dt  = 0.00001     #Has to be small. preferably smaller than 0.00002

Thot  = 1000.0
Tcold =  295.0
T     = Tcold*ones(Float64,NX,NY) # Initial Temperature
RHS   = zeros(Float64,NX,NY)    # Right Hand Side

# Creating a meshgrid to plot solution u
x = range(0, length = NX, stop = Lx)
y = range(0, length = NY, stop = Ly)

function ini_cond(T::Array{Float64,2},dx::Float64,dy::Float64,NX::Integer,NY::Integer,Thot::Float64)
	for i = 1:NX
	    for j = 1:NY
		    if (((i*dx-0.5).^2+(j*dy-0.5).^2 <= 0.1 ) && ( (i*dx-0.5)^2+(j*dy-0.5)^2 >= 0.05 ) ) # These are two circle equations. Heat is added between these 2 circles
				T[i,j] = Thot # This adds heat as a circle. This adds heat as a ring.
		    end
	    end
        end
end

dx2_inv = 1/dx2;
dy2_inv = 1/dy2;

function Heat2D_Fast(T::Array{Float64,2},RHS::Array{Float64,2},dt::Float64,A::Float64,dx2_inv::Float64,dy2_inv::Float64,NT::Integer,NX::Integer,NY::Integer,Tcold::Float64)

    @inbounds for t=0:NT-1

		# FIXED BOUNDARY CONDITIONS AT EDGES (UNCOMENT BELOW IF NECESSERY)
		T[1,:]  .= Tcold # boundary condition T=295 at y=0
		T[NX,:] .= Tcold # boundary condition T=295 at y=1
		T[:,1]  .= Tcold # boundary condition T=295 at x=0
		T[:,NY] .= Tcold # boundary condition T=295 at x=1

        for  j=2:NY-1, i=2:NX-1

            RHS[i,j] =  dt*A*( (T[i-1,j]-2*T[i,j]+T[i+1,j])*dx2_inv  +
                               (T[i,j-1]-2*T[i,j]+T[i,j+1])*dy2_inv )
        end

        for j=2:NY-1, i=2:NX-1
            T[i,j] = T[i,j] + RHS[i,j]
        end

		if t%1000 == 0 # prints T value at every certain amount of iterations in the center to edge
			   println( "Iteration ",t, ",T value from center is:"," -- ",
			   round(T[50,50];digits=2)," -- ", round(T[50,60];digits=2)," -- ",
			   round(T[50,70];digits=2)," -- ", round(T[50,80];digits=2)," -- ",
			   round(T[50,90];digits=2)," -- ", round(T[50,100];digits=2) )
		end
    end

	figure() #2D plots
	surf(x,y,T, cmap=ColorMap("jet"))
    zlim(Tcold,Thot*0.5);
	xlabel("X")
	ylabel("Y")
	zlabel("T")

	figure() #Contour Plots
	contourf(x,y,T, cmap=ColorMap("jet"))
	xlabel("X")
	ylabel("Y")
	colorbar()

end

function main()
	 ini_cond(T,dx,dy,NX,NY,Thot)
	 Heat2D_Fast(T,RHS,dt,A,dx2_inv,dy2_inv,NT,NX,NY,Tcold) # call this function more times if you want to run the same simulation but for different time steps NT
end

@time main() # This runs the script
