

# PHYSICAL PARAMETERS
const A = 0.25    #Diffusion coefficient
const Lx = 1.0   #Domain size x
const Ly = 1.0   #Domain size y
const Lz = 1.0   #Domain size z

# NUMERICAL PARAMETERS
const NT = 100        #Number of time steps
const NX = 100        #Number of grid points in x
const NY = 100        #Number of grid points in y
const NZ = 100        #Number of grid points in z

const dx = Lx/(NX-1) #Grid step in x (space)
const dx2 = dx^2
const dy = Ly/(NY-1) #Grid step in y (space)
const dy2 = dy^2
const dz = Lz/(NZ-1) #Grid step in z (space)
const dz2 = dz^2

const dt = 0.00002

T = zeros(Float64,NX,NY,NZ)
RHS = zeros(Float64,NX,NY,NZ)

function ini_cond()
	for i = 1:NX
	    for j = 1:NY
	    	for k = 1:NZ
		    if (((i*dx-0.5).^2+(j*dy-0.5).^2+(k*dz-0.5)^2 <= 0.1 ) && ( (i*dx-0.5)^2+(j*dy-0.5)^2+(k*dz-0.5)^2 >= 0.05 ) )
				T[i,j,k] = 1
		    end
		end
	    end
        end
end

function Heat3D_Slow()
	 for n = 0:NT-1
	     RHS[2:NX-1,2:NY-1,2:NZ-1] = dt*A*( (T[1:NX-2,2:NY-1,2:NZ-1]-2*T[2:NX-1,2:NY-1,2:NZ-1]+T[3:NX,2:NY-1,2:NZ-1])/dx2  +
             			                (T[2:NX-1,1:NY-2,2:NZ-1]-2*T[2:NX-1,2:NY-1,2:NZ-1]+T[2:NX-1,3:NY,2:NZ-1])/dy2  +
                                                (T[2:NX-1,2:NY-1,1:NZ-2]-2*T[2:NX-1,2:NY-1,2:NZ-1]+T[2:NX-1,2:NY-1,3:NZ])/dz2 )

             T[2:NX-1,2:NY-1,2:NZ-1] = T[2:NX-1,2:NY-1,2:NZ-1] + RHS[2:NX-1,2:NY-1,2:NZ-1]

	     if n%10 == 0
	     	    println("Iteration ",n, "value at center:","--", T[50,50,50],"--", T[50,50,60],"--", T[50,50,70],"--", T[50,50,80],"--", T[50,50,90],"--", T[50,50,100])
	     end
	 end
end

function main()
	 ini_cond()
	 Heat3D_Slow()
end

main()
