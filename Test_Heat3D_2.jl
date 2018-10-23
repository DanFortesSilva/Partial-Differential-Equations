
A  = 0.25  # Diffusion coefficient
Lx = 1.0   # Domain size x
Ly = 1.0   # Domain size y
Lz = 1.0   # Domain size z

NT = 100;  # Number of time steps
NX = 100;  # Number of grid points in x
NY = 100;  # Number of grid points in y
NZ = 100;  # Number of grid points in z

dx  = Lx/(NX-1) # Grid step in x (space)
dx2 = dx^2
dy  = Ly/(NY-1) # Grid step in y (space)
dy2 = dy^2
dz  = Lz/(NZ-1) # Grid step in z (space)
dz2 = dz^2

dt = 0.1

T   = zeros(NX,NY,NZ)
RHS = zeros(Float64,NX,NY,NZ);

# BOUNDARY CONDITIONS
#T[1,:,:]  .= 1;
#T[NX,:,:] .= 1;
#T[:,1,:]  .= 1;
#T[:,NY,:] .= 1;
#T[:,:,1]  .= 1;
#T[:,:,NZ] .= 0;


function test_main(T::Array{Float64,3},RHS::Array{Float64,3},dt::Float64,A::Float64,dx2::Float64,dy2::Float64,
dz2::Float64,NX::Integer,NY::Integer,NZ::Integer)

    dx2 = 1/dx2;
    dy2 = 1/dy2;
    dz2 = 1/dz2;

    @inbounds for t=0:NT-1

        for k=2:NZ-1, j=2:NY-1, i=2:NX-1

            RHS[i,j,k] =  dt*A*((T[i-1,j,k]-2*T[i,j,k]+T[i+1,j,k])*dx2  +
                                (T[i,j-1,k]-2*T[i,j,k]+T[i,j+1,k])*dy2  +
                                (T[i,j,k-1]-2*T[i,j,k]+T[i,j,k+1])*dz2 )
        end

        for k=2:NZ-1, j=2:NY-1, i=2:NX-1
            T[i,j,k] = T[i,j,k] + RHS[i,j,k]
        end
    end

end

test_main(T,RHS,dt,A,dx2,dy2,dz2,NX,NY,NZ)
@time test_main(T,RHS,dt,A,dx2,dy2,dz2,NX,NY,NZ)
