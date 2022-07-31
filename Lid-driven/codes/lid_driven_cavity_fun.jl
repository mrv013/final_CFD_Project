#clearconsole()

using CPUTime
using Printf
using Plots
using FFTW
using PlotlyJS

font = Plots.font("Times New Roman", 18)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

#-----------------------------------------------------------------------------#
# Compute L-2 norm for a vector
#-----------------------------------------------------------------------------#
function compute_l2norm(nx, ny, r)
    norm = 0.0
    # println(residual)
    for j = 1:ny+1 for i = 1:nx+1
        norm = norm + r[i,j]^2
    end end
    # println(rms)
    norm = sqrt(norm/((nx+1)*(ny+1)))
    return norm
end

#-----------------------------------------------------------------------------#
# Fast poisson solver for homozeneous Dirichlet domain
   # f: source term of the Poisson equation
#∇^2(Ψ)= ω
#∇^2(u)= f
# get f in [nx+1,ny+1] and give Ψ in [2:nx,2:ny] without boundarys
#-----------------------------------------------------------------------------#

function fps_sine(nx,ny,dx,dy,f)

    data = Array{Complex{Float64}}(undef,nx-1,ny-1)
    data1 = Array{Complex{Float64}}(undef,nx-1,ny-1)
    e = Array{Complex{Float64}}(undef,nx-1,ny-1)

    u = Array{Complex{Float64}}(undef,nx-1,ny-1)

    for i = 1:nx-1
        for j = 1:ny-1
            data[i,j] = f[i+1,j+1]
        end
    end

    e = FFTW.r2r(data,FFTW.RODFT00)

    for i = 1:nx-1
        for j = 1:ny-1
            alpha = (2.0/(dx*dx))*(cos(pi*i/nx) - 1.0) +
                    (2.0/(dy*dy))*(cos(pi*j/ny) - 1.0)
            data1[i,j] = e[i,j]/alpha
        end
    end

    u = FFTW.r2r(data1,FFTW.RODFT00)/((2*nx)*(2*ny))

    return u
end


#-----------------------------------------------------------------------------#
# Compute L-2 norm for a vector
#-----------------------------------------------------------------------------#
#function compute_l2norm(nx, ny, r)


function bc2(nx,ny,dx,dy,w,s)
    # second order approximation
    # boundary condition for vorticity (Jensen) left and right
    u0=1
    for j = 1:ny+1
        w[1,j] = (3.5*s[1,j]-4.0*s[2,j]+0.5*s[3,j])/(dx*dx)
        w[nx+1,j]= (3.5*s[nx+1,j]-4.0*s[nx,j]+0.5*s[nx-1,j])/(dx*dx)
    end

    # boundary condition for vorticity (Jensen) bottom and top
    for i = 1:nx+1
        w[i,1] = (3.5*s[i,1]-4.0*s[i,2]+0.5*s[i,3])/(dy*dy)
        w[i,ny+1]= (3.5*s[i,ny+1]-4.0*s[i,ny]+0.5*s[i,ny-1])/(dy*dy) - u0*3.0/dy
    end
end

function bc_ψ(nx,ny,dx,dy,s)
    # second order approximation
    # boundary condition for ψ  left and right
    # boundary condition for ψ  bottom and top
    for i = 1:nx+1
        s[i,ny+1]= -dy/2
        s[i,1]= 0
    end
    for j = 1:ny+1
        s[nx+1,j]= 0
        s[1,j]= 0
    end
end

#-----------------------------------------------------------------------------#
# Compute numerical solution
#   - Time integration using Runge-Kutta third order
#   - 2nd-order finite difference discretization
#-----------------------------------------------------------------------------#
function numerical(nx,ny,nt,dx,dy,dt,re,wn,sn,rms)
    wn_all_time=Array{Float64}(undef, nx+1, ny+1,nt)
    sn_all_time=Array{Float64}(undef, nx+1, ny+1,nt)
    wt = Array{Float64}(undef, nx+1, ny+1) # temporary array during RK3 integration
    r = Array{Float64}(undef, nx+1, ny+1) # right hand side
    sp = Array{Float64}(undef, nx+1, ny+1) # old streamfunction

    """
    x_range=LinRange(0,L_x,nx+1)
    y_range=LinRange(0,L_y,ny+1)
    time_lapse=LinRange(0,,nt)
    """

    bc_ψ(nx,ny,dx,dy,sn)
    for k = 1:nt
        for i = 1:nx+1 for j = 1:ny+1
            sp[i,j] = sn[i,j]
        end end
        ##Runge-Kutta step1
        # Compute right-hand-side from vorticity  Runge-Kutta step1
        rhs(nx,ny,dx,dy,re,wn,sn,r)

        for i = 2:nx for j = 2:ny
            wt[i,j] = wn[i,j] + dt*r[i,j]# next step time ω : Runge-Kutta step1
        end end

        ##Runge-Kutta step2
        #bondry coundition :Runge-Kutta step2
        bc2(nx,ny,dx,dy,wt,sn)

        # compute streamfunction from vorticity:  Runge-Kutta step2
        sn[2:nx,2:ny] = fps_sine(nx,ny,dx,dy,-wt)
        #print(sn)
        # Compute right-hand-side from vorticity: Runge-Kutta step2
        rhs(nx,ny,dx,dy,re,wt,sn,r)

        for i = 2:nx for j = 2:ny
            wt[i,j] = 0.75*wn[i,j] + 0.25*wt[i,j] + 0.25*dt*r[i,j] #Runge-Kutta step2
        end end
        ##Runge-Kutta step3
        #bondry coundition :Runge-Kutta step3
        bc2(nx,ny,dx,dy,wt,sn)

        # compute streamfunction(Ψ) from vorticity (ω): Runge-Kutta step3
        sn[2:nx,2:ny] = fps_sine(nx,ny,dx,dy,-wt)

        # Compute right-hand-side from vorticity: Runge-Kutta step3
        rhs(nx,ny,dx,dy,re,wt,sn,r)

        for i = 2:nx for j = 2:ny
            wn[i,j] = (1.0/3.0)*wn[i,j] + (2.0/3.0)*wt[i,j] + (2.0/3.0)*dt*r[i,j]
        end end

        ##Runge-Kutta step4(next time step)
        #bondry coundition :Runge-Kutta step4(next time step)
        bc2(nx,ny,dx,dy,wn,sn)

        # compute streamfunction from vorticity
        sn[2:nx,2:ny] = fps_sine(nx,ny,dx,dy,-wn)

        sn_all_time[:,:,k]=sn
        wn_all_time[:,:,k]=wn

        rms[k] = 0.0 #  ?????
        for i = 1:nx+1 for j = 1:ny+1
            rms[k] = rms[k] + (sn[i,j] - sp[i,j])^2
        end end

        rms[k] = sqrt(rms[k]/((nx+1)*(ny+1)))
        println("time num=,$k,rms= $(rms[k])" )
    end
    return (end_sn=sn, end_wn=wn,
            sn_all_time=sn_all_time,
            wn_all_time=wn_all_time)
           #time_lapse=time_lapse, h_range=x_range, y_range=y_range)

end

#-----------------------------------------------------------------------------#
# Calculate right hand term of the inviscid Burgers equation
# nx,ny: total number of grids in x and y direction
# dx,dy: grid spacing in x and y direction
# re: Reynolds number of the flow
# w: vorticity field
# s: streamfunction
# r: right hand side of the Runge-Kutta scheme (-Jacobian+Laplacian terms)
# :r = -J(w,ψ) + ν ∇^2(w)
#-----------------------------------------------------------------------------#
function rhs(nx,ny,dx,dy,re,w,s,r)
    # Arakawa numerical scheme for Jacobian
    aa = 1.0/(re*dx*dx)
    bb = 1.0/(re*dy*dy)
    gg = 1.0/(4.0*dx*dy)
    hh = 1.0/3.0

    for i = 2:nx for j = 2:ny
        j1 = gg*((w[i+1,j]-w[i-1,j])*(s[i,j+1]-s[i,j-1]) -
                 (w[i,j+1]-w[i,j-1])*(s[i+1,j]-s[i-1,j]))

        j2 = gg*(w[i+1,j]*(s[i+1,j+1]-s[i+1,j-1]) -
                 w[i-1,j]*(s[i-1,j+1]-s[i-1,j-1]) -
                 w[i,j+1]*(s[i+1,j+1]-s[i-1,j+1]) +
                 w[i,j-1]*(s[i+1,j-1]-s[i-1,j-1]))

        j3 = gg*(w[i+1,j+1]*(s[i,j+1]-s[i+1,j]) -
                 w[i-1,j-1]*(s[i-1,j]-s[i,j-1]) -
            	 w[i-1,j+1]*(s[i,j+1]-s[i-1,j]) +
            	 w[i+1,j-1]*(s[i+1,j]-s[i,j-1]))

        jac = (j1+j2+j3)*hh

        #Central difference for Laplacian
        r[i,j] = -jac + aa*(w[i+1,j]-2.0*w[i,j]+w[i-1,j]) +
                        bb*(w[i,j+1]-2.0*w[i,j]+w[i,j-1])
        end end
end
#------------------------------------------------------------------------------#
function jac(nx,ny,dx,dy,w,s)
    # Arakawa numerical scheme for Jacobian
    jac=Array{Float64}(undef, nx+1, ny+1) # right hand side
    gg = 1.0/(4.0*dx*dy)
    hh = 1.0/3.0

    for i = 2:nx for j = 2:ny
        j1 = gg*((w[i+1,j]-w[i-1,j])*(s[i,j+1]-s[i,j-1]) -(-1)*
                 (w[i,j+1]-w[i,j-1])*(s[i+1,j]-s[i-1,j]))

        j2 = gg*(w[i+1,j]*(s[i+1,j+1]-s[i+1,j-1]) -
                 w[i-1,j]*(s[i-1,j+1]-s[i-1,j-1]) -(-1)*
                 w[i,j+1]*(s[i+1,j+1]-s[i-1,j+1]) +(-1)*
                 w[i,j-1]*(s[i+1,j-1]-s[i-1,j-1]))

        j3 = gg*(w[i+1,j+1]*(s[i,j+1]-s[i+1,j]) -
                 w[i-1,j-1]*(s[i-1,j]-s[i,j-1]) -(-1)*
            	 w[i-1,j+1]*(s[i,j+1]-s[i-1,j]) +(-1)*
            	 w[i+1,j-1]*(s[i+1,j]-s[i,j-1]))

        jac[i,j] = (j1+j2+j3)*hh

        end end
        return jac
end
#------------------------------------------------------------------------------#
function ditem_u_v(nx,ny,dx,dy,sn)
    # Arakawa numerical scheme for Jacobian
    u=Array{Float64}(undef, nx+1, ny+1) # right hand side
    v=Array{Float64}(undef, nx+1, ny+1) # right hand side
    u_P=Array{Float64}(undef, nx+3, ny+3) # right hand side
    v_P=Array{Float64}(undef, nx+3, ny+3) # right hand side
    # boundary condition in x dimention
    for j = 1:ny+1
        v[1,j] = 0;v[nx+1,j]=0
        u[1,j] = 0;u[nx+1,j]=0
    end
    # boundary condition in y dimention
    for i = 1:nx+1
        u[i,1] = 0;u[i,ny+1]=1
        v[i,1] = 0;v[i,ny+1]=0
    end

    for j = 1:ny+3
        v_P[1,j] = 0;v_P[nx+1,j]=0
        u_P[1,j] = 0;u_P[nx+1,j]=0
    end
    # boundary condition in y dimention
    for i = 1:nx+3
        u_P[i,1] = 0;u_P[i,ny+1]=1
        v_P[i,1] = 0;v_P[i,ny+1]=0
    end

    for i = 2:nx for j = 2:ny
        #x
        v[i,j] = (sn[i+1,j]+sn[i-1,j])/2/dx
        #y
        u[i,j] = (sn[i,j+1]+sn[i,j-1])/2/dy
    end end
    #0 1:ny+1 ny+2
    #1 2:ny+2 ny+3
    v_P[2:nx+2,2:ny+2]=v;u_P[2:nx+2,2:ny+2]=u
    return (u=u,v=v, u_P=u_P,v_P=v_P)
end
#------------------------------------------------------------------------------#
function ditem_P(nx,ny,dx,dy,v,u,re)
    # Arakawa numerical scheme for Jacobian
    #u=Array{Float64}(undef, nx+1, ny+1) # right hand side
    #v=Array{Float64}(undef, nx+1, ny+1) # right hand side
    P_sorce_term1=Array{Float64}(undef, nx+1, ny+1)
    P_sorce_term2=Array{Float64}(undef, nx+1, ny+1)
    P_sorce_term3=Array{Float64}(undef, nx+1, ny+1)
    P_sorce_term4=Array{Float64}(undef, nx+1, ny+1)
    P_sorce_term5=Array{Float64}(undef, nx+1, ny+1)
    P_sorce_term6=Array{Float64}(undef, nx+1, ny+1)
    P=Array{Float64}(undef, nx+1, ny+1)
    for i = 1:nx+1 for j = 1:ny+1
        P_sorce_term1[i,j] = 0.0 # initial condition
        P_sorce_term2[i,j]= 0.0
        P_sorce_term3[i,j]= 0.0
        P_sorce_term4[i,j] = 0.0
        P_sorce_term5[i,j] = 0.0
        P_sorce_term6[i,j] = 0.0
        P[i,j] = 00.0 # initial
    end end

    for i = 3:nx+1 for j = 3:ny+1
    P_sorce_term1[i-1,j-1]=u[i,j]*(u[i+1,j]-2*u[i,j]+u[i-1,j])/(dx)^2 + v[i,j]*(v[i,j+1]-2*v[i,j]+v[i,j-1])/(dy)^2

    P_sorce_term2[i-1,j-1]=(((u[i+1,j]-u[i-1,j])/(2*dx))^2+((v[i,j+1]-v[i,j-1])/(2*dy))^2)
    #println(i,j)
    P_sorce_term3[i-1,j-1]=(((u[i+1,j+1]-u[i+1,j-1])/(2*dy)-(u[i-1,j+1]-u[i-1,j-1])/(2*dy))/(2*dx)
                       +((v[i+1,j+1]-v[i+1,j-1])/(2*dy)-(v[i-1,j+1]-v[i-1,j-1])/(2*dy))/(2*dx))

    P_sorce_term5[i-1,j-1]=1/re*((u[i+2,j]-2*u[i+1,j]+u[i-1,j]-u[i-2,j])/dx^2
                           + (v[i,j+2]-2*v[i,j+1]+v[i,j-1]-v[i,j-2])/dy^2)

    P_sorce_term6[i-1,j-1]=1/re*(((u[i+1,j+1]-u[i-1,j+1])/(2*dx) - 2*(u[i+1,j]-u[i-1,j])/(2*dx)
                             +(u[i+1,j-1]-u[i-1,j-1])/(2*dx)/(dy)^2)
                             +((v[i+1,j+1]-v[i+1,j-1])/(2*dy) - 2*(v[i,j+1]-v[i,j-1])/(2*dy)
                             +(v[i-1,j+1]-v[i-1,j-1])/(2*dy)/(dx)^2 )  )
    end end


    P_sorce_term4[2:nx,2:ny]= jac(nx,ny,dx,dy,v[2:nx+2,2:ny+2],u[2:nx+2,2:ny+2])[2:nx,2:ny]
    #println("P_sorce_term",P_sorce_term)
    P_sorce_term_final= .-P_sorce_term1.-P_sorce_term2.-P_sorce_term3.-P_sorce_term4+1/re.*(P_sorce_term5.+P_sorce_term6)
    P[2:nx,2:ny]=fps_sine(nx,ny,dx,dy,P_sorce_term_final)

    #println("P",P)
    return (P=P)
end



index(ξ,Δξ)=floor(Int, (ξ+0.5Δξ)/Δξ)+1
print(index(0.6,0.1))
index_tim(ξ,Δξ)=floor(Int, (ξ+0.5Δξ)/Δξ)
slope(e,h)=@. log(e[2:end]/e[1:end-1])/log(h[2:end]/h[1:end-1])

#-----------------------------------------------------------------------------#

function plot_x_seris(ν, res_upwind, time_want, nRange1,var_name)
	#time_step_want=20
	show_x=res_upwind[1].h_range
	n_show_x=length(show_x)
	nr1=length(nRange1)
    plotc=Plots.scatter(xlabel="y",ylabel=var_name, title="solution for diffrent
    gride,grid line x=$ν,y=[0,1], time=$(time_want) s",legend= :bottomleft)
    for i in 1:nr1
	  time_step_want=index_tim.(time_want,res_upwind[i].time_lapse[1])
	  Plots.scatter!(plotc, res_upwind[i].h_range,
      res_upwind[i].all_time[:,time_step_want], label="n=$(nRange1[i]+1)" )# ** n=$(nRange1[i])
    end
	#plotc

    err=zeros(nr1-1,n)
	h3= @. 1.0/(nRange1) #[0.1, 0.05, 0.025, 0.0125]# **h3= @. 1.0/(nRange1-1)
	#τ=res_upwind[1].time_lapse[2]
	plo3=Plots.scatter(axis= :log,xlabel="\$ dy\$",ylabel="successive Error $var_name", title="solution convergence, ,grid line x=$ν,y=[0,1],dx=dy=h, time=$(time_want) s",legend= :bottomright)
    for i in 1:( nr1-1 )
		time_step_want=index_tim.(time_want,res_upwind[i].time_lapse[1])
        i1=index.(show_x,res_upwind[i].h_range[2])
	    i2=index.(show_x,res_upwind[i+1].h_range[2])
		for j in 2:(n_show_x-1)#**nRange1[1]
	        err[i,j]= abs.(res_upwind[i+1].all_time[i2[j],time_step_want]-res_upwind[i].all_time[i1[j],time_step_want])
		end
    end
	lable=["x=$ν,y=$ξ" for i in 1:1,ξ in show_x[2:end] ]
	successive_Error=Plots.scatter(plo3, h3[1:end-1], err[:,2:end], label=lable )

	plot_slope=Plots.plot(xlabel="\$ dy\$",ylabel="Slope $var_name ",title="solution Slope convergence, time=$(time_want) s",legend= :topleft)
    for i in 2:(n_show_x-1)#**length(nRange1)
    #for i in [3,5,7,9]#**length(nRange1)
		Plots.plot!(plot_slope,h3[1:end-2], slope(err[:,i],h3[1:end-1]),label="x=$ν,y=$(res_upwind[1].h_range[i])",lw=2)
    end


    return (solution=plotc, successive_Error=successive_Error, plot_slope=plot_slope)
end

function plot_x_seris_re(ν, res_upwind, time_want, Re_Range,var_name)
	#time_step_want=20
	show_x=res_upwind[1].h_range
	n_show_x=length(show_x)
	nr1=length(Re_Range)
    plotc=Plots.scatter(xlabel="y",ylabel=var_name, title="solution for diffrent
    Re,grid line x=$ν,y=[0,1], time=$(time_want) s",legend= :bottomleft)
    for i in 1:nr1
	  time_step_want=index_tim.(time_want,res_upwind[i].time_lapse[1])
	  Plots.scatter!(plotc, res_upwind[i].h_range,
      res_upwind[i].all_time[:,time_step_want], label="Re=$(Re_Range[i])" )# ** n=$(nRange1[i])
    end
	return (solution=plotc)
end

function plot_t_seris(time, ν_Range, res_upwind_t, x_want_t, x_line_want_L,method_name)
	dx_t=res_upwind_t[1].h_range[1]*2
	x_step_t=index_cv(x_want_t,dx_t)
	n_time=length(time)
	n_ν=length(ν_Range)
	err_t=zeros(n_ν-1,n_time)
	#h3= @. 1.0/(nRange1-1)
	plot_upwind_t=Plots.scatter(axis= :log,xlabel="dt",ylabel="successive Error",title="slution convergence $method_name, x=$(x_line_want_L),y=$(x_want_t) ",legend= :bottomright)
    println("m")
	for i in 1:( n_ν-1 )
        i1=index_tim.(time,ν_Range[i])
		println(i1)
	    i2=index_tim.(time,ν_Range[i+1])
		println(i2)
		for j in 1:n_time
	    err_t[i,j]= abs.(res_upwind_t[i+1].all_time[x_step_t,i2[j]]-res_upwind_t[i].all_time[x_step_t,i1[j]])
		end
    end
	lable_t=["t=$ξ   " for i in 1:1,ξ in time ]
	successive_Error=Plots.scatter(plot_upwind_t, ν_Range[1:end-1], err_t[:,2:end], label=lable_t )

	ploupwind_slope_t=Plots.plot(xlabel="dt",ylabel="Slope",title="slution convergence $method_name, x=$(x_line_want_L),y=$(x_want_t) ",legend= :bottomright)
    for i in 2:n_time
		Plots.plot!(ploupwind_slope_t,ν_Range[1:end-2], slope(err_t[:,i],ν_Range[1:end-1]),label="t=$(time[i])",lw=2)
    end
return (slope_t=ploupwind_slope_t , successive_Error=successive_Error)
end



function plot_t_seris(time, ν_Range, res_upwind_t, x_want_t, x_line_want_L,method_name)
	dx_t=res_upwind_t[1].h_range[1]*2
	x_step_t=index_cv(x_want_t,dx_t)
	n_time=length(time)
	n_ν=length(ν_Range)
	err_t=zeros(n_ν-1,n_time)
	#h3= @. 1.0/(nRange1-1)
	plot_upwind_t=Plots.scatter(axis= :log,xlabel="dt",ylabel="successive Error",title="slution convergence $method_name, x=$(x_line_want_L),y=$(x_want_t) ",legend= :bottomright)
    println("m")
	for i in 1:( n_ν-1 )
        i1=index_tim.(time,ν_Range[i])
		println(i1)
	    i2=index_tim.(time,ν_Range[i+1])
		println(i2)
		for j in 1:n_time
	    err_t[i,j]= abs.(res_upwind_t[i+1].all_time[x_step_t,i2[j]]-res_upwind_t[i].all_time[x_step_t,i1[j]])
		end
    end
	lable_t=["t=$ξ   " for i in 1:1,ξ in time ]
	successive_Error=Plots.scatter(plot_upwind_t, ν_Range[1:end-1], err_t[:,2:end], label=lable_t )

	ploupwind_slope_t=Plots.plot(xlabel="dt",ylabel="Slope",title="slution convergence $method_name, x=$(x_line_want_L),y=$(x_want_t) ",legend= :bottomright)
    for i in 2:n_time
		Plots.plot!(ploupwind_slope_t,ν_Range[1:end-2], slope(err_t[:,i],ν_Range[1:end-1]),label="t=$(time[i])",lw=2)
    end
return (slope_t=ploupwind_slope_t , successive_Error=successive_Error)
end

#
function generate_all_time_for_1_line(n,L_y,L_x,dt,end_time,re,w0,s0,rms,x_line_want_L)
    #L_y,L_x,dt,tf,re,wn,sn,rms


	nx =n; ny = n
    dx = (L_x)/nx ; dy = (L_y)/ny
	nt=Int64(end_time/dt)
    x_line_want=index(x_line_want_L,dx)
	w0 = Array{Float64}(undef, nx+1, ny+1)
	s0 = Array{Float64}(undef, nx+1, ny+1)
	for i = 1:nx+1 for j = 1:ny+1
        w0[i,j] = 0.0 # initial condition
        s0[i,j] = 0.0 # initial streamfunction
    end end

	#dt = 0.001;tf = 2.00;nt = Int64(tf/dt)
	x_range=LinRange(0,L_x,nx+1); y_range=LinRange(0,L_y,ny+1)
	time_lapse=LinRange(dt,end_time,nt)

	numerical_stor=numerical(nx,ny,nt,dx,dy,dt,re,w0,s0,rms)
	#print(numerical_stor.wn_all_time)
	wn_x_line=numerical_stor.wn_all_time[x_line_want,:,:]; sn_x_line=numerical_stor.sn_all_time[x_line_want,:,:]
	#wn_x_line[:,]=(wn[x_line_want,:])';sn_x_line[:,]=(sn[x_line_want,:])'
	return (all_time=sn_x_line , time_lapse=time_lapse', h_range=y_range)
end
