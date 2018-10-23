# This script solves the 2D Laplace Equation using Finite Differences Method

using PyPlot # Calling python functions and plotting function matplotlib

# PHYSICAL PARAMETERS
A  = 0.25  #Diffusion coefficient
Lx = 1.0   #Domain size x
Ly = 1.0   #Domain size y

# NUMERICAL PARAMETERS
NT = 5000        #Number of time steps
NX = 100         #Number of grid points in x
NY = 100         #Number of grid points in y

dx = Lx/(NX-1)   #Grid step in x (space)
dx2 = dx^2
dy = Ly/(NY-1)   #Grid step in y (space)
dy2 = dy^2

dt = 0.00002     # Has to be small

T   = zeros(Float64,NX,NY) # Temperature
RHS = zeros(Float64,NX,NY) # Right Hand Side

function ini_cond()
	for i = 1:NX
	    for j = 1:NY
		    if (((i*dx-0.5).^2+(j*dy-0.5).^2 <= 0.1 ) && ( (i*dx-0.5)^2+(j*dy-0.5)^2 >= 0.05 ) )
				T[i,j] = 1.0 # This adds heat as a circle
		    end
	    end
        end
end

function Heat2D_Slow()
	 for t = 0:NT-1
	     RHS[2:NX-1,2:NY-1] = dt*A*( (T[1:NX-2,2:NY-1]-2*T[2:NX-1,2:NY-1]+T[3:NX,2:NY-1])/dx2  +
             			             (T[2:NX-1,1:NY-2]-2*T[2:NX-1,2:NY-1]+T[2:NX-1,3:NY])/dy2  )

         T[2:NX-1,2:NY-1] = T[2:NX-1,2:NY-1] + RHS[2:NX-1,2:NY-1]

	     if t%10 == 0 # prints T value at every certain amount of iterations in the center
	     	    println("Iteration ", t , " value at center:","--", T[50,50],"--", T[50,60],"--", T[50,70],"--", T[50,80],"--", T[50,90],"--", T[50,100])
	     end
	 end
end

function main()
	 ini_cond()
	 Heat2D_Slow()
end

@time main() # This runs the script

# Creating a meshgrid to plot solution T
x = range(0, length = NX, stop = Lx)
y = range(0, length = NY, stop = Ly)


figure(1)
HeatPlot = surf(x,y,T, cmap=ColorMap("jet"))
xlabel("X")
ylabel("Y")
zlabel("T")
colorbar(HeatPlot)

figure(2)
HeatContour = contourf(x,y,T, cmap=ColorMap("jet"))
xlabel("X")
ylabel("Y")
colorbar(HeatContour)
