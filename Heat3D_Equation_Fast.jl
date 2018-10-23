# This script solves the 3D Laplace Equation using Finite Differences Method

# PHYSICAL PARAMETERS
A  = 0.25  #Diffusion coefficient
Lx = 1.0   #Domain size x
Ly = 1.0   #Domain size y
Lz = 1.0   #Domain size z

# NUMERICAL PARAMETERS
NT = 1000        #Number of time steps
NX = 100         #Number of grid points in x
NY = 100         #Number of grid points in y
NZ = 100         #Number of grid points in z

dx  = Lx/(NX-1)   #Grid step in x (space)
dx2 = dx^2
dy  = Ly/(NY-1)   #Grid step in y (space)
dy2 = dy^2
dz  = Lz/(NZ-1)   #Grid step in y (space)
dz2 = dz^2

dt = 0.00002     # Has to be small

T   = zeros(Float64,NX,NY,NZ) # Temperature
RHS = zeros(Float64,NX,NY,NZ) # Right Hand Side

# Function for initial/boundary conditions. Definite the type of variables (ex Float64, Integer,etc) so it becomes is faster
function ini_cond(T::Array{Float64,3},dx::Float64,dy::Float64,dz::Float64,NX::Integer,NY::Integer,NZ::Integer)
	for i = 1:NX
	    for j = 1:NY
	    	for k = 1:NZ
		    if (((i*dx-0.5).^2+(j*dy-0.5).^2+(k*dz-0.5)^2 <= 0.1 ) && ( (i*dx-0.5)^2+(j*dy-0.5)^2+(k*dz-0.5)^2 >= 0.05 ) )# These are two circle equations. Heat is added between these 2 circles
				T[i,j,k] = 1.0 # This adds heat as a circle. This adds heat as a ring.
		    end
		end
	    end
        end
end

dx2_inv = 1/dx2;
dy2_inv = 1/dy2;
dz2_inv = 1/dz2;

function Heat3D_Fast(T::Array{Float64,3},RHS::Array{Float64,3},dt::Float64,A::Float64,dx2_inv::Float64,dy2_inv::Float64,dz2_inv::Float64,NX::Integer,NY::Integer,NZ::Integer)

    @inbounds for t=0:NT-1

		# FIXED BOUNDARY CONDITIONS AT EDGES (UNCOMENT BELOW IF NECESSERY)
		#T[1,:,:]  .= 295.0 # boundary condition T=295 at y=0
		#T[NX,:,:] .= 295.0 # boundary condition T=295 at y=1
		#T[:,1,:]  .= 295.0 # boundary condition T=295 at x=0
		#T[:,NY,:] .= 295.0 # boundary condition T=295 at x=1
		#T[:,:,1]  .= 295.0 # boundary condition T=295 at z=0
		#T[:,:,NZ] .= 295.0 # boundary condition T=295 at z=1

        for k=2:NZ-1, j=2:NY-1, i=2:NX-1

            RHS[i,j,k] =  dt*A*( (T[i-1,j,k]-2*T[i,j,k]+T[i+1,j,k])*dx2_inv  +
                                 (T[i,j-1,k]-2*T[i,j,k]+T[i,j+1,k])*dy2_inv  +
								 (T[i,j,k-1]-2*T[i,j,k]+T[i,j,k+1])*dz2_inv )
        end

        for k=2:NZ-1, j=2:NY-1, i=2:NX-1
            T[i,j,k] = T[i,j,k] + RHS[i,j,k]
        end

		if t%100 == 0 # prints T value at every certain amount of iterations in the center
			   println( "Iteration ",t, " value at center:"," -- ",
			   round(T[50,50,50];digits=4)," -- ", round(T[50,50,60];digits=4)," -- ",
			   round(T[50,50,70];digits=4)," -- ", round(T[50,50,80];digits=4)," -- ",
			   round(T[50,50,90];digits=4)," -- ", round(T[50,50,100];digits=4) )
		end
    end

end

function main()
	 ini_cond(T,dx,dy,dz,NX,NY,NZ)
	 Heat3D_Fast(T,RHS,dt,A,dx2_inv,dy2_inv,dz2_inv,NX,NY,NZ)
end

@time main() # This runs the script
