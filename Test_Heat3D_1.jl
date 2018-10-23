NX = 100;
    NY = 20;
    NZ =20;
    T = randn(NX,NY,NZ);
    RHS = zeros(Float64,NX,NY,NZ);
    NT = 10;
    dx2 = 1/0.5;
    dy2 = 1/0.5;
    dz2 = 1/0.4;
    A = 1.1
    dt = 0.1

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
