
A  = 0.25  #Diffusion coefficient
Lx = 1.0   #Domain size x
Ly = 1.0   #Domain size y

NT = 100;
NX = 100;
NY = 100;

dx  = Lx/(NX-1) #Grid step in x (space)
dx2 = dx^2
dy  = Ly/(NY-1) #Grid step in y (space)
dy2 = dy^2

dt = 0.1

T   = randn(NX,NY)
RHS = zeros(Float64,NX,NY);


function test_main(T::Array{Float64,2},RHS::Array{Float64,2},dt::Float64,A::Float64,dx2::Float64,dy2::Float64,NX::Integer,NY::Integer)

    dx2 = 1/dx2;
    dy2 = 1/dy2;

    @inbounds for t=0:NT-1

        for  j=2:NY-1, i=2:NX-1

            RHS[i,j] =  dt*A*((T[i-1,j]-2*T[i,j]+T[i+1,j])*dx2  +
                                (T[i,j-1]-2*T[i,j]+T[i,j+1])*dy2 )
        end

        for j=2:NY-1, i=2:NX-1
            T[i,j] = T[i,j] + RHS[i,j]
        end
    end

end

test_main(T,RHS,dt,A,dx2,dy2,NX,NY)
@time test_main(T,RHS,dt,A,dx2,dy2,NX,NY)
