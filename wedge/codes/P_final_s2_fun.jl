### A Pluto.jl notebook ###
# v0.19.2

using Markdown
using InteractiveUtils

begin
using CPUTime
using Printf
using Plots
using PlotlyJS
end

##
function Initial_conditions(nx,gamma,qn,x,dx)
  # Sod's Riemann problem
  # Left side
  rhoL = 1.0
  uL = 1.0
  pL = 1.0
  # Right side
  rhoR = 0.125
  uR = 0.0
  pR = 0.1

  # nodal storage location (grid)
  #for i = 1:nx
  #	x[i] = -0.5*dx + dx*(i)
  #end


  xc = 0 # seperator location ## in bundery
  for i = 1:nx
	if (x[i] > xc)
		rho = rhoR
		u = uR
		p = pR
	else
		rho = rhoL
		u = uL
		p = pL
	end

	e = p/(rho*(gamma-1.0)) + 0.5*u*u

	#conservative variables
	qn[i,1] = rho
	qn[i,2] = rho*u
	qn[i,3] = rho*e
  end
  return qn
end
##
function fluxes(nx,gamma,q,f)
	for i = 1:nx+1
		p = (gamma-1.0)*(q[i,3]-0.5*q[i,2]*q[i,2]/q[i,1])
		f[i,1] = q[i,2]
		f[i,2] = q[i,2]*q[i,2]/q[i,1] + p
		f[i,3] = q[i,2]*q[i,3]/q[i,1] + p*q[i,2]/q[i,1]
	end
	return f
end

##
function roe(nx,gamma,uL,uR,f,fL,fR)
	dd = Array{Float64}(undef,3)
	dF = Array{Float64}(undef,3)
	V = Array{Float64}(undef,3)
	gm = gamma-1.0

	for i = 1:nx+1
		#Left and right states:
		rhLL = uL[i,1]
		uuLL = uL[i,2]/rhLL
		eeLL = uL[i,3]/rhLL
	    ppLL = gm*(eeLL*rhLL - 0.5*rhLL*(uuLL*uuLL))
	    hhLL = eeLL + ppLL/rhLL

		rhRR = uR[i,1]
		uuRR = uR[i,2]/rhRR
		eeRR = uR[i,3]/rhRR
	    ppRR = gm*(eeRR*rhRR - 0.5*rhRR*(uuRR*uuRR))
	    hhRR = eeRR + ppRR/rhRR

		alpha = 1.0/(sqrt(abs(rhLL)) + sqrt(abs(rhRR)))

		uu = (sqrt(abs(rhLL))*uuLL + sqrt(abs(rhRR))*uuRR)*alpha
		hh = (sqrt(abs(rhLL))*hhLL + sqrt(abs(rhRR))*hhRR)*alpha
		aa = sqrt(abs(gm*(hh-0.5*uu*uu)))

		D11 = abs(uu)
		D22 = abs(uu + aa)
		D33 = abs(uu - aa)

		beta = 0.5/(aa*aa)
		phi2 = 0.5*gm*uu*uu

		#Right eigenvector matrix
		R11, R21, R31 = 1.0, uu, phi2/gm
		R12, R22, R32 = beta, beta*(uu + aa), beta*(hh + uu*aa)
		R13, R23, R33 = beta, beta*(uu - aa), beta*(hh - uu*aa)

		#Left eigenvector matrix
		L11, L12, L13 = 1.0-phi2/(aa*aa), gm*uu/(aa*aa), -gm/(aa*aa)
		L21, L22, L23 = phi2 - uu*aa, aa - gm*uu, gm
		L31, L32, L33 = phi2 + uu*aa, -aa - gm*uu, gm

		for m = 1:3
			V[m] = 0.5*(uR[i,m]-uL[i,m])
		end

		dd[1] = D11*(L11*V[1] + L12*V[2] + L13*V[3])
		dd[2] = D22*(L21*V[1] + L22*V[2] + L23*V[3])
		dd[3] = D33*(L31*V[1] + L32*V[2] + L33*V[3])

		dF[1] = R11*dd[1] + R12*dd[2] + R13*dd[3]
		dF[2] = R21*dd[1] + R22*dd[2] + R23*dd[3]
		dF[3] = R31*dd[1] + R32*dd[2] + R33*dd[3]

		for m = 1:3
			f[i,m] = 0.5*(fR[i,m]+fL[i,m]) - dF[m]
		end
	end
	return f
end

##
function left_bc!()
        #K= 3*(n - 2)+1 .+(1:2)
	    #p=3*(n-3)
		rho = 0.125
        u = 0.0
        p = 0.1
	    gamma=1.4
	    e = p/(rho*(gamma-1.0)) + 0.5*u*u
        f1 = [rho rho*u rho*e]

	return f1
end

function right_bc!(tp)
        #K= 3*(n - 2)+1 .+(1:2)
	    #p=3*(n-3)
	  f1=0
	  if tp==1
		f1=[1 1 1]
	  end
	  if tp==0
		f1=[1 -1 1]
	  end

	return f1
end

#---------------------------------------------------------------------------#
#nonlinear weights for upwind direction
#---------------------------------------------------------------------------#
function wcL(v1,v2,v3,v4,v5)
    eps = 1.0e-6

    # smoothness indicators
    s1 = (13.0/12.0)*(v1-2.0*v2+v3)^2 + 0.25*(v1-4.0*v2+3.0*v3)^2
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)^2 + 0.25*(v2-v4)^2
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)^2 + 0.25*(3.0*v3-4.0*v4+v5)^2

    # computing nonlinear weights w1,w2,w3
    c1 = 1.0e-1/((eps+s1)^2)
    c2 = 6.0e-1/((eps+s2)^2)
    c3 = 3.0e-1/((eps+s3)^2)

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    # candiate stencils
    q1 = v1/3.0 - 7.0/6.0*v2 + 11.0/6.0*v3
    q2 =-v2/6.0 + 5.0/6.0*v3 + v4/3.0
    q3 = v3/3.0 + 5.0/6.0*v4 - v5/6.0

    # reconstructed value at interface
    f = (w1*q1 + w2*q2 + w3*q3)

    return f
end

#-----------------------------------------------------------------------------#
# WENO reconstruction for upwind direction (positive; left to right)
# u(i): solution values at finite difference grid nodes i = 1,...,N
# f(j): reconstructed values at nodes j = i-1/2; j = 1,...,N+1
#-----------------------------------------------------------------------------#
function wenoL(n,u,tp)
	f = Array{Float64}(undef,n+1,3)
        #Lb=left_bc!()
	    #Rb=right_bc!(tp)
	left_bc=b_c(tp[1])
	right_bc=b_c(tp[2])
	for m = 1:3
	    i = 0
	    v1 = u[i+3,m]*left_bc.coeff[m]  +left_bc.constant[m]#Lb[m]  #-3  u[i-2,m]
	    v2 = u[i+2,m]*left_bc.coeff[m]  +left_bc.constant[m]#Lb[m]  #-2 u[i-1,m]
	    v3 = u[i+1,m]*left_bc.coeff[m]  +left_bc.constant[m]#Lb[m]    #-1 u[i,m]
	    v4 = u[i+1,m] #1
	    v5 = u[i+2,m]
	    f[i+1,m] = wcL(v1,v2,v3,v4,v5) #f(1)  uR[0,m]

	    i = 1
	    v1 = u[i+1,m]*left_bc.coeff[m]  +left_bc.constant[m]#Lb[m] #2  u[i-2,m]
	    v2 = u[i,m]*left_bc.coeff[m]  +left_bc.constant[m]#Lb[m] #1  u[i-1,m]
	    v3 = u[i,m]   #1  u[i,m]
	    v4 = u[i+1,m] #2
	    v5 = u[i+2,m] #3
	    f[i+1,m] = wcL(v1,v2,v3,v4,v5)#fL(1+1/2)   uR[1,m]

	    i = 2
	    v1 = u[i-1,m]*left_bc.coeff[m]  +left_bc.constant[m]#Lb[m] #1  u[i-2,m]
	    v2 = u[i-1,m] #1  u[i-1,m]
	    v3 = u[i,m]   #2  u[i,m]
	    v4 = u[i+1,m] #3
	    v5 = u[i+2,m] #4
	    f[i+1,m] = wcL(v1,v2,v3,v4,v5)#f(3)  uL[2,m]

	    for i = 3:n-2
	        v1 = u[i-2,m]
	        v2 = u[i-1,m]
	        v3 = u[i,m]
	        v4 = u[i+1,m]
	        v5 = u[i+2,m]
	        f[i+1,m] = wcL(v1,v2,v3,v4,v5)
	    end

	    i = n-1
	    v1 = u[i-2,m]
	    v2 = u[i-1,m]
	    v3 = u[i,m]
	    v4 = u[i+1,m]
	    v5 = u[i+1,m]*right_bc.coeff[m]  +right_bc.constant[m]#u[i+1,m]* Rb[m]
	    f[i+1,m] = wcL(v1,v2,v3,v4,v5)

	    i = n
	    v1 = u[i-2,m]
	    v2 = u[i-1,m]
	    v3 = u[i,m]
	    v4 = u[i,m]*right_bc.coeff[m]  +right_bc.constant[m]#u[i,m]* Rb[m]
	    v5 = u[i-1,m]*right_bc.coeff[m]  +right_bc.constant[m]#u[i-1,m]* Rb[m]
	    f[i+1,m] = wcL(v1,v2,v3,v4,v5)
	end
	return f
end
#---------------------------------------------------------------------------#
#nonlinear weights for downwind direction
#---------------------------------------------------------------------------#
function wcR(v1,v2,v3,v4,v5)
    eps = 1.0e-6

    s1 = (13.0/12.0)*(v1-2.0*v2+v3)^2 + 0.25*(v1-4.0*v2+3.0*v3)^2
    s2 = (13.0/12.0)*(v2-2.0*v3+v4)^2 + 0.25*(v2-v4)^2
    s3 = (13.0/12.0)*(v3-2.0*v4+v5)^2 + 0.25*(3.0*v3-4.0*v4+v5)^2

    c1 = 3.0e-1/(eps+s1)^2
    c2 = 6.0e-1/(eps+s2)^2
    c3 = 1.0e-1/(eps+s3)^2

    w1 = c1/(c1+c2+c3)
    w2 = c2/(c1+c2+c3)
    w3 = c3/(c1+c2+c3)

    # candiate stencils
    q1 =-v1/6.0      + 5.0/6.0*v2 + v3/3.0
    q2 = v2/3.0      + 5.0/6.0*v3 - v4/6.0
    q3 = 11.0/6.0*v3 - 7.0/6.0*v4 + v5/3.0

    # reconstructed value at interface
    f = (w1*q1 + w2*q2 + w3*q3)
    return f
end
#-----------------------------------------------------------------------------#
# WENO reconstruction for downwind direction (negative; right to left)
# u(i): solution values at finite difference grid nodes i = 1,...,N+1
# f(j): reconstructed values at nodes j = i-1/2; j = 2,...,N+1
#-----------------------------------------------------------------------------#
function wenoR(n,u,tp)
	f = Array{Float64}(undef,n+1,3)
    #Lb=left_bc!()
	#Rb=right_bc!(tp)
	left_bc=b_c(tp[1])
	right_bc=b_c(tp[2])
	for m = 1:3
	    i = 1
	    v1 = u[i+1,m]*left_bc.coeff[m]  +left_bc.constant[m] #Lb[m]
	    v2 = u[i,m]*left_bc.coeff[m]    +left_bc.constant[m]#Lb[m]
	    v3 = u[i,m]
	    v4 = u[i+1,m]
	    v5 = u[i+2,m]
	    f[i,m] = wcR(v1,v2,v3,v4,v5)

	    i = 2
	    v1 = u[i-1,m]*left_bc.coeff[m]  +left_bc.constant[m]#Lb[m]
	    v2 = u[i-1,m]
	    v3 = u[i,m]
	    v4 = u[i+1,m]
	    v5 = u[i+2,m]
	    f[i,m] = wcR(v1,v2,v3,v4,v5)

	    for i = 3:n-2
	        v1 = u[i-2,m]
	        v2 = u[i-1,m]
	        v3 = u[i,m]
	        v4 = u[i+1,m]
	        v5 = u[i+2,m]
	        f[i,m] = wcR(v1,v2,v3,v4,v5)
	    end

	    i = n-1
	    v1 = u[i-2,m]
	    v2 = u[i-1,m]
	    v3 = u[i,m]
	    v4 = u[i+1,m]
	    v5 = u[i+1,m]*right_bc.coeff[m]  +right_bc.constant[m]#u[i+1,m]*Rb[m]
	    f[i,m] = wcR(v1,v2,v3,v4,v5)

	    i = n
	    v1 = u[i-2,m]
	    v2 = u[i-1,m]
	    v3 = u[i,m]
	    v4 = u[i,m]*right_bc.coeff[m]  +right_bc.constant[m]#u[i,m]*Rb[m]
	    v5 = u[i-1,m]*right_bc.coeff[m]  +right_bc.constant[m]#u[i-1,m]*Rb[m]
	    f[i,m] = wcR(v1,v2,v3,v4,v5)

	    i = n+1
	    v1 = u[i-2,m]
	    v2 = u[i-1,m]
	    v3 = u[i-1,m]*right_bc.coeff[m]  +right_bc.constant[m]#u[i-1,m]*Rb[m]
	    v4 = u[i-2,m]*right_bc.coeff[m]  +right_bc.constant[m]#u[i-2,m]*Rb[m]
	    v5 = u[i-3,m]*right_bc.coeff[m]  +right_bc.constant[m]#u[i-3,m]*Rb[m]
	    f[i,m] = wcR(v1,v2,v3,v4,v5)
	end
	return f
end


##
function rhs(nx,dx,gamma,q,b_tp)
    qL0 = Array{Float64}(undef,nx+1,3)
    qR0 = Array{Float64}(undef,nx+1,3)

	fL0 = Array{Float64}(undef,nx+1,3)
    fR0 = Array{Float64}(undef,nx+1,3)
	f0 = Array{Float64}(undef,nx+1,3)

    r=Array{Float64}(undef,nx,3)

	# WENO Reconstruction
	qL = wenoL(nx,q,b_tp) #left in face not valume
    qR = wenoR(nx,q,b_tp) #Right in face not valume

	# Computing fluxes
	fL=fluxes(nx,gamma,qL,fL0)
	fR=fluxes(nx,gamma,qR,fR0)

	# compute Riemann solver using Gudanov scheme
	f=roe(nx,gamma,qL,qR,f0,fL,fR)

	# RHS
	for i = 1:nx for m = 1:3
		r[i,m] = -(f[i+1,m] - f[i,m])/dx
	end end
	return r
end
##
function RK3(nx,gamma,qn,dx,dt,b_tp)
   qt = Array{Float64}(undef, nx,3) # temporary array during RK3 integration
   #r = Array{Float64}(undef, nx,3)
   #TVD RK3 for time integration
   r=rhs(nx,dx,gamma,qn,b_tp)
   for i = 1:nx for m = 1:3
	 qt[i,m] = qn[i,m] + dt*r[i,m]
   end end
   r=rhs(nx,dx,gamma,qt,b_tp)
   for i = 1:nx for m = 1:3
	qt[i,m] = 0.75*qn[i,m] + 0.25*qt[i,m] + 0.25*dt*r[i,m]
   end end
   r=rhs(nx,dx,gamma,qt,b_tp)
   for i = 1:nx for m = 1:3
	qt[i,m] = (1.0/3.0)*qn[i,m] + (2.0/3.0)*qt[i,m] + (2.0/3.0)*dt*r[i,m]

	if m==1 && qt[i,m]<0
		qt[i,m]=0
	end
   end end

   return qt
end
##
function callcu_1D_1time(nx,dx,dt,qn,b_tp)
    #x = Array{Float64}(undef, nx)
    qn = Array{Float64}(undef, nx,3) # numerical solsution at every time step
    #qt = Array{Float64}(undef, nx,3) # temporary array during RK3 integration
    #r = Array{Float64}(undef, nx,3)

    gamma = 1.4 # specific gas ratio

    #Initial_conditions(nx,gamma,qn,x,dx)

	qn=RK3(nx,gamma,qn,dx,dt,b_tp)

	# TVD RK3 for time integration

	return qn
end
#------------------------------------------------------------------------------#
function numerical(nx,ns,nt,dx,dt,q_stor,b_tp)
    x = Array{Float64}(undef, nx)
    qn0 = Array{Float64}(undef, nx,3) # numerical solsution at every time step
    #qt = Array{Float64}(undef, nx,3) # temporary array during RK3 integration
    #r = Array{Float64}(undef, nx,3)
    x
	for i = 1:nx
	x[i] = -0.5*dx + dx*(i)
    end
    gamma = 1.4 # specific gas ratio

    qn=Initial_conditions(nx,gamma,qn0,x,dx)
    print("Initial_conditions :$qn")
	ri = 1 # record index
    for i = 1:nx for m = 1:3
        q_stor[i,m,ri] = qn[i,m] # store solution at t=0
    end	end

	# TVD RK3 for time integration
    for n = 1:nt # time step
		println(n)
		qn=RK3(nx,gamma,qn,dx,dt,b_tp)
        println("after step $n:$qn")

        freq = Int64(nt/ns)
        if (mod(n,freq) == 0)
            ri = ri + 1
			for i = 1:nx for m = 1:3
            	q_stor[i,m,ri] = qn[i,m]
			end end
        end

    end
	return qn ,q
end

function index_x_yline(n_yline,dy,dx,alfa,y0,x0,ny_max,nx_max)

 #y0=ny0*dx
 #x0=nx0*dx
 ymin_element=(n_yline-1)*dy

 x_in_wedge_to_line=(ymin_element-y0)/tan(alfa/180*3.14) + x0

 nx_end_yline = floor(Int,x_in_wedge_to_line/dx) + 1

 if nx_end_yline>ny_max
 nx_end_yline=nx_max
 end

 return nx_end_yline
end

function index_y_xline(n_xline,dy,dx,alfa,y0,x0,ny_max,nx_max)#ny_max,nx_max

 #y0=ny0*dx
 #x0=nx0*dx
 xmax_element=(n_xline)*dx
 y_out_wedge_to_line = (xmax_element-x0)*tan(alfa/180*3.14)+y0
 ny_start_xline = floor(Int,y_out_wedge_to_line/dy)+1
 if ny_start_xline<=0
  ny_start_xline=1
 end

 return ny_start_xline
end

function Initial_conditions_2D_for_1line(nx,gamma,qn)
 # Sod's Riemann problem
 # Left side

 rho = 0.125
 u = 0.0
 v = 0.0
 p = 0.1
 e = p/(rho*(gamma-1.0)) + 0.5*u*u+ 0.5*v*v
 for i = 1:nx
  #conservative variables
  qn[i,1] = rho
  qn[i,2] = rho*u
  qn[i,3] = rho*v
  qn[i,4] = rho*e

 end
 return qn
end

function Initial_conditions_2D_for_all_line(dy,dx,alfa,y0,x0,ny,nx,gamma)
 D2_q=zeros(Float64,ny,nx,4)
 #print(D2_q)
 for j=1:ny
 nx_end_yline=index_x_yline(j,dy,dx,alfa,y0,x0,ny,nx)
 D2_q_1y = D2_q[j,1:nx_end_yline,:]
 D2_q_1y=Initial_conditions_2D_for_1line(nx_end_yline,gamma,D2_q_1y)
 D2_q[j,1:nx_end_yline,:]=D2_q_1y
 end

 for i=1:nx
  ny_start_xline=index_y_xline(i,dy,dx,alfa,y0,x0,ny,nx)
  num_y=ny-ny_start_xline+1
  D2_q_1x = D2_q[ny_start_xline:end,i,:]
  D2_q_1x=Initial_conditions_2D_for_1line(num_y,gamma,D2_q_1x)
  D2_q[ny_start_xline:end,i,:]=D2_q_1x
 end
 return D2_q
end

function b_c(tp)
	  if tp==1
		coeff=[1 1 1]
		constant=[0 0 0]
	  elseif tp==0
		coeff=[1 -1 1]
		constant=[0 0 0]
	  elseif tp==2
		rho = 1.0
        p = 1.0
	    gamma=1.4
		a=sqrt(gamma*p/rho)
		u = 1.2*a
	    e = p/(rho*(gamma-1.0)) + 0.5*u*u
		coeff=[0 0 0]
        constant = [rho rho*u rho*e]
	  end

	return (coeff=coeff,constant=constant)
end

function time_1step_all_xline(D2_q,dy,dx,dt,alfa,y0,x0,ny,nx,gamma)
 for i=1:nx
 #i=9
 ny_start_xline=index_y_xline(i,dy,dx,alfa,y0,x0,ny,nx)
 num_y=ny-ny_start_xline+1
 D2_q_1x_0 = D2_q[ny_start_xline:end,i,[1,3,4]] #[rho,rho*v,rho*e]
 	if (i*dx)<x0 && ny_start_xline==1
	   b_tp_end_y=[1,1] else b_tp_end_y=[0,1]
	end
 D2_q_1x_1 = RK3(num_y,gamma,D2_q_1x_0,dx,dt,b_tp_end_y)
 D2_q[ny_start_xline:end,i,[1,3,4]]=D2_q_1x_1
 end
	return D2_q
end

function time_1step_all_yline(D2_q,dy,dx,dt,alfa,y0,x0,ny,nx,gamma)
 for j=1:ny
 #j=3
 nx_end_yline=index_x_yline(j,dy,dx,alfa,y0,x0,ny,nx)
 D2_q_1y_0 = D2_q[j,1:nx_end_yline,[1,2,4]] #[rho,rho*u,rho*e]
	if nx_end_yline==nx
	   b_tp_end_x=[2,1] else b_tp_end_x=[2,0]
	end
 D2_q_1y_1 = RK3(nx_end_yline,gamma,D2_q_1y_0,dx,dt,b_tp_end_x)
 D2_q[j,1:nx_end_yline,[1,2,4]]=D2_q_1y_1
 end
	return D2_q
end


index_cv(ξ,Δξ)=floor(Int, (ξ+0.5*Δξ)/Δξ+0.5)
index_tim(ξ,Δξ)=floor(Int, (ξ+0.5*Δξ)/Δξ)
slope(e,h)=@. log(e[2:end]/e[1:end-1])/log(h[2:end]/h[1:end-1])

##
function plot_x_seris(y_want_t,ν, res_upwind, time_want, nRange1,var_name)
	#time_step_want=20
	show_x=y_want_t
	#res_upwind[1].h_range

	n_show_x=length(show_x)
	nr1=length(nRange1)
    plotc=Plots.scatter(xlabel="L",ylabel=var_name, title="solution $var_name for diffrent
    gride,grid line y=$ν,x=[0,1], time=$(time_want) s",legend= :bottomleft)
    for i in 1:nr1
	  time_step_want=index_tim.(time_want,res_upwind[i].time_lapse[1])
	  scatter!(plotc, res_upwind[i].h_range,
      res_upwind[i].all_time[:,time_step_want], label="n=$(nRange1[i])" )
	end
    #n=$(nRange1[i]+1)
    err=zeros(nr1-1,nRange1[1])
	h3= @. 1.0/(nRange1) #[0.1, 0.05, 0.025, 0.0125]# **h3= @. 1.0/(nRange1-1)
	#τ=res_upwind[1].time_lapse[2]
	plo3=Plots.scatter(axis= :log,xlabel="\$ h\$",ylabel="successive Error $var_name",
	 title="solution convergence, ,grid line y=$ν,x=[0,1],dx=dy=h, time=$(time_want) s",legend= :bottomright)

    for i in 1:( nr1-1 ) #all n in nRenge
		time_step_want=index_tim.(time_want,res_upwind[i].time_lapse[1])
		dx1=res_upwind[i].h_range[1]*2
        i1=index_cv.(show_x,dx1) # n  of  bace_x in dx-->nRenge1 i
		dx2=res_upwind[i+1].h_range[1]*2
	    i2=index_cv.(show_x,dx2) # n  of  bace_x in dx-->nRenge1 i+1
		for j in 1:(n_show_x)#**2:nRange1[1] # all point by smallest n
	    err[i,j]= abs.(res_upwind[i+1].all_time[i2[j],time_step_want]-
        res_upwind[i].all_time[i1[j],time_step_want])
		end
    end

	lable=["x=$ξ,y=$ν" for i in 1:1,ξ in show_x[1:end] ]#**2:end
	successive_Error=Plots.scatter(plo3, h3[1:end-1], err[:,2:end], label=lable )

	plot_slope=Plots.plot(xlabel="\$ h\$",ylabel="Slope $var_name ",title="solution Slope $var_name convergence, time=$(time_want) s",legend= :topleft)
    for i in 1:(n_show_x)#**length(nRange1)
    #for i in [3,5,7,9]#**length(nRange1)
		Plots.plot!(plot_slope,h3[1:end-2], slope(err[:,i],h3[1:end-1]),label="x=$ν,y=$(show_x[i])",lw=2)
    end
return (solution=plotc,successive_Error=successive_Error,plot_slope=plot_slope,err=err)
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


function generate_all_time_for_1_line(n,L_y,L_x,end_time,dt,alfa,y0,x0,gamma,x_line_want_L)
    #L_y,L_x,dt,tf,re,wn,sn,rms


	nx =n; ny = n
    dx = (L_x)/nx ; dy = (L_y)/ny
	nt=index_tim(end_time,dt)#Int64(end_time/dt)
    x_line_want=index_cv(x_line_want_L,dx)
    #println(x_line_want)
	D2_q_0=Initial_conditions_2D_for_all_line(dy,dx,alfa,y0,x0,ny,nx,gamma)#(dy,dx,alfa,y0,x0,ny,nx,gamma)
    D2_q=generate_all_time(D2_q_0,dy,dx,dt,alfa,y0,x0,ny,nx,nt,gamma)

	#dt = 0.001;tf = 2.00;nt = Int64(tf/dt)
	x_range=LinRange(0,L_x-dx,nx).+0.5*dx;y_range=LinRange(0,L_y-dy,ny).+0.5*dy
	#x_range=LinRange(0,L_x,nx+1); y_range=LinRange(0,L_y,ny+1)
	time_lapse=LinRange(dt,end_time,nt)

	#print(numerical_stor.wn_all_time)
	rho_x_line=(D2_q.stor[x_line_want,:,1,:]);
    rho_u_x_line=(D2_q.stor[x_line_want,:,2,:])
	rho_v_x_line=(D2_q.stor[x_line_want,:,3,:])
	rho_e_x_line=(D2_q.stor[x_line_want,:,4,:])
	#wn_x_line[:,]=(wn[x_line_want,:])';sn_x_line[:,]=(sn[x_line_want,:])'
return [(all_time=rho_x_line , time_lapse=time_lapse', h_range=y_range)
      ,(all_time=rho_u_x_line , time_lapse=time_lapse', h_range=y_range)
      ,(all_time=rho_v_x_line , time_lapse=time_lapse', h_range=y_range)
      ,(all_time=rho_e_x_line , time_lapse=time_lapse', h_range=y_range)]
end

function generate_all_time(D2_q_0,dy,dx,dt,alfa,y0,x0,ny,nx,nt,gamma)
 D2_q_stor=Array{Float64}(undef,ny,nx,4,nt)
 for t= 0:2:nt-2
 D2_qy1=time_1step_all_yline(D2_q_0,dy,dx,dt,alfa,y0,x0,ny,nx,gamma)
 println("time step: $(t+1),yline ")
 D2_q_stor[:,:,:,t+1]=D2_qy1
 D2_qx1=time_1step_all_xline(D2_qy1,dy,dx,dt,alfa,y0,x0,ny,nx,gamma)
 println("time step: $(t+2),xline ")
 D2_q_stor[:,:,:,t+2]=D2_qx1
 D2_q_0=D2_qx1
 #print(D2_q_0)
 end
 return (stor=D2_q_stor, endtime=D2_q_0)
end
