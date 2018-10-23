using PyPlot

function buildUpB(rho, dt, dx, dy, u, v)
    b = zeros(size(u))
    b[2:end-1,2:end-1]=rho*(1/dt*((u[2:end-1,3:end]-u[2:end-1,1:end-2])/(2*dx)+(v[3:end,2:end-1]-v[1:end-2,2:end-1])/(2*dy))-
                      ((u[2:end-1,3:end]-u[2:end-1,1:end-2])/(2*dx)).^2-
                      2*((u[3:end,2:end-1]-u[1:end-2,2:end-1])/(2*dy).*(v[2:end-1,3:end]-v[2:end-1,1:end-2])/(2*dx))-
                      ((v[3:end,2:end-1]-v[1:end-2,2:end-1])/(2*dy)).^2)

    ####Periodic BC Pressure @ x = 2
    b[2:end-1,end]=rho*(1/dt*((u[2:end-1,1]-u[2:end-1,end-1])/(2*dx)+(v[3:end,end]-v[1:end-2,end])/(2*dy))-
        ((u[2:end-1,1]-u[2:end-1,end-1])/(2*dx)).^2-
    2*((u[3:end,end]-u[1:end-2,end])/(2*dy).*(v[2:end-1,1]-v[2:end-1,end-1])/(2*dx))-
    ((v[3:end,end]-v[1:end-2,end])/(2*dy)).^2)

    ####Periodic BC Pressure @ x = 0
    b[2:end-1,1]=rho*(1/dt*((u[2:end-1,2]-u[2:end-1,end])/(2*dx)+(v[3:end,1]-v[1:end-2,1])/(2*dy))-
    ((u[2:end-1,2]-u[2:end-1,end])/(2*dx)).^2-
    2*((u[3:end,1]-u[1:end-2,1])/(2*dy).*(v[2:end-1,2]-v[2:end-1,end])/(2*dx))-
                    ((v[3:end,1]-v[1:end-2,1])/(2*dy)).^2)

    return b
end
