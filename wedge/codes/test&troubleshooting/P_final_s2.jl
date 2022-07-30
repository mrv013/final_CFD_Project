### A Pluto.jl notebook ###
# v0.19.2

using Markdown
using InteractiveUtils

# ╔═╡ 78610b9b-43c3-490d-b04e-a0a95b357d7b
begin
using CPUTime
using Printf
using Plots
end

# ╔═╡ ce187bce-2c98-4fad-8b84-15b249154a8e


# ╔═╡ 2c360e6d-ff1f-4653-8c9a-4fc650780926
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

xc = 0 # seperator location ## in bundery
for i = 1:nx
	if (x[i] >= xc)
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

# ╔═╡ 44d05ac3-0ae3-430f-a755-fac8aff0adef


# ╔═╡ f6b5833a-a39c-491a-9218-cd0c956c5d05
function fluxes(nx,gamma,q,f)
	for i = 1:nx+1
		p = (gamma-1.0)*(q[i,3]-0.5*q[i,2]*q[i,2]/q[i,1])
		f[i,1] = q[i,2]
		f[i,2] = q[i,2]*q[i,2]/q[i,1] + p
		f[i,3] = q[i,2]*q[i,3]/q[i,1] + p*q[i,2]/q[i,1]
	end
	return f
end


# ╔═╡ 5bc2f630-f76b-4e11-8901-6b6d6a1018a4
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


# ╔═╡ e95f8105-ab89-45a3-8324-e8bd838f654a
function left_bc!()
        #K= 3*(n - 2)+1 .+(1:2)
	    #p=3*(n-3)
		rho = 1.0
        u = 1.0
        p = 1.0
	    gamma=1.4
	    e = p/(rho*(gamma-1.0)) + 0.5*u*u
        f1 = [rho rho*u rho*e]

	return f1
end

# ╔═╡ 5b5726e6-d600-4f3b-af8e-a3c58c8cde9f
function b_c(tp)
	  if tp==1
		coeff=[1 1 1]
		constant=[0 0 0]
	  elseif tp==0
		coeff=[1 -1 1]
		constant=[0 0 0] 
	  elseif tp==2  
		rho = 1.0
        u = 1.0
        p = 1.0
	    gamma=1.4
	    e = p/(rho*(gamma-1.0)) + 0.5*u*u
		coeff=[0 0 0] 
        constant = [rho rho*u rho*e]
	  end
	   
	return (coeff=coeff,constant=constant)
end

# ╔═╡ 6e74ea04-2b8d-47a5-91e4-b437796303c9
b_c(1)

# ╔═╡ 08fc3d88-2d40-4b48-8434-b9fad1ab0b7c
left_bc!()[1,1]

# ╔═╡ ae7eff58-e358-40b0-ba31-e85507e641eb
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

# ╔═╡ 0f3a4c00-666f-46ac-8547-d728ef99a10e
tp4=[1,2]

# ╔═╡ 679bf7ab-5f35-4315-8dcc-0b52c3498f84

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


# ╔═╡ afd16e84-96b3-4715-a712-f0d8fc19b4c8

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


# ╔═╡ 5500059d-1aed-43c5-846f-fd56e8bef8c2

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


# ╔═╡ 98832477-0c6c-4770-8f45-74e244b7c42b

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



# ╔═╡ de66199c-155c-4e23-a2d6-ff665237037f
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


# ╔═╡ bc4a9a50-09c0-11ed-1051-0d7888c35c8a
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
end end
	return qt
end

# ╔═╡ 61aa2c26-0d88-4782-98c8-2acf857d8eb7
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


# ╔═╡ 73cf2abf-ed39-4ef0-9895-c736677a9be7
wcR(1,1,1,1,-1)#v1,v2,v3,v4,v5)

# ╔═╡ 08bca7eb-ae0e-4faa-9104-faa6fcef2cb4
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

# ╔═╡ 582ae9c6-c502-4430-b7cf-2d1f6ab782bf
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

# ╔═╡ a0f08d8f-9ef3-4c25-a970-bc973942caa2
function Initial_conditions_2D_for_1line(nx,gamma,qn)
 # Sod's Riemann problem
 # Left side
 rho = 1.0
 u = 1.0
 v = 1.0
 p = 1.0
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

# ╔═╡ 74c5b2d9-0e8c-4992-8f75-060bf5b41b75
function Initial_conditions_2D_for_all_line(dy,dx,alfa,y0,x0,ny,nx,gamma)
 D2_q=zeros(Float64,nx,nx,4)

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
	

# ╔═╡ 175c19f6-d11e-4c3c-adf3-b024846da01b
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

# ╔═╡ 0c85b3fa-1e00-4171-951a-24619a0bcc0b
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

# ╔═╡ 1a908140-b336-47db-9c54-6649103d1ec2
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

# ╔═╡ 942d57b9-0ade-4813-9298-292107d1215e
index_cv(ξ,Δξ)=floor(Int,(ξ+0.5*Δξ)/Δξ+0.5)

# ╔═╡ fbfa7298-6f8f-4011-87bb-03c32a6296f7
index_tim(ξ,Δξ)=floor(Int, (ξ+0.5*Δξ)/Δξ)

# ╔═╡ 56fcad51-0146-4a24-b8d6-231f8d72f8f5
slope(e,h)=@. log(e[2:end]/e[1:end-1])/log(h[2:end]/h[1:end-1])

# ╔═╡ 1b373135-757d-4286-add8-7ab0d6a1b3cf
function generate_all_time_for_1_line(n,L_y,L_x,end_time,dt,alfa,y0,x0,gamma,x_line_want_L)
    #L_y,L_x,dt,tf,re,wn,sn,rms


	nx =n; ny = n
    dx = (L_x)/nx ; dy = (L_y)/ny
	nt=Int64(end_time/dt)
    x_line_want=index_cv(x_line_want_L,dx)
    #println(x_line_want)
	D2_q_0=Initial_conditions_2D_for_all_line(dy,dx,alfa,y0,x0,ny,nx,gamma)#(dy,dx,alfa,y0,x0,ny,nx,gamma)
    D2_q=generate_all_time(D2_q_0,dy,dx,dt,alfa,y0,x0,ny,nx,nt,gamma)

	#dt = 0.001;tf = 2.00;nt = Int64(tf/dt)
	x_range=LinRange(0,L_x-dx,nx).+0.5*dx;y_range=LinRange(0,L_y-dy,ny).+0.5*dy
	#x_range=LinRange(0,L_x,nx+1); y_range=LinRange(0,L_y,ny+1)
	time_lapse=LinRange(dt,end_time,nt)

	#print(numerical_stor.wn_all_time)
	rho_x_line=(D2_q.stor[:,x_line_want,1,:]);
    rho_u_x_line=(D2_q.stor[:,x_line_want,2,:])
	rho_v_x_line=(D2_q.stor[:,x_line_want,3,:])
	rho_e_x_line=(D2_q.stor[:,x_line_want,4,:])
	#wn_x_line[:,]=(wn[x_line_want,:])';sn_x_line[:,]=(sn[x_line_want,:])'
	return [(all_time=rho_x_line , time_lapse=time_lapse', h_range=y_range),(all_time=rho_u_x_line , time_lapse=time_lapse', h_range=y_range),(all_time=rho_v_x_line , time_lapse=time_lapse', h_range=y_range),(all_time=rho_e_x_line , time_lapse=time_lapse', h_range=y_range)]
end


# ╔═╡ 7e178c53-14ab-4217-b551-cf0db73579ca
a[1]

# ╔═╡ 3bab86a7-37ff-4e98-b0b2-e4450a9ec126
function plot_t_seris(time, ν_Range, res_upwind_t, x_want_t, method_name)
	dx_t=res_upwind_t[1].h_range[1]*2
	x_step_t=index_cv(x_want_t,dx_t)
	n_time=length(time)
	n_ν=length(ν_Range)
	err_t=zeros(n_ν-1,n_time)
	#h3= @. 1.0/(nRange1-1)
	plot_upwind_t=scatter(axis= :log,xlabel="dt",ylabel="successive Error",title="slution convergence $method_name, x=$(x_want_t)*L,y=[0,1] ",legend= :bottomright)
    for i in 1:( n_ν-1 )
        i1=index_tim.(time,ν_Range[i])
	    i2=index_tim.(time,ν_Range[i+1])
		for j in 1:n_time
	    err_t[i,j]= abs.(res_upwind_t[i+1].all_time[x_step_t,i2[j]]-res_upwind_t[i].all_time[x_step_t,i1[j]])
		end
    end
	lable_t=["t=$ξ   " for i in 1:1,ξ in time ]
	successive_Error=scatter(plot_upwind_t, ν_Range[1:end-1], err_t[:,2:end], label=lable_t )

	ploupwind_slope_t=plot(xlabel="dt",ylabel="Slope",title="slution convergence $method_name, x=$((x_want_t))*L,y=[0,1] ",legend= :bottomright)
    for i in 2:n_time
		plot!(ploupwind_slope_t,ν_Range[1:end-2], slope(err_t[:,i],ν_Range[1:end-1]),label="t=$(time[i])",lw=2)
    end
	return (slope_t=ploupwind_slope_t , successive_Error=successive_Error)
end

# ╔═╡ ce673755-7701-4b33-b0fb-61ae7991731f
function plot_x_seris(ν, res_upwind, time_want, nRange1,var_name)
	#time_step_want=20
	show_x=res_upwind[1].h_range
	n_show_x=length(show_x)
	nr1=length(nRange1)
    plotc=Plots.scatter(xlabel="Y",ylabel=var_name, title="solution for diffrent
    gride,grid line x=$ν,y=[0,1], time=$(time_want) s",legend= :bottomleft)
    for i in 1:nr1
	  time_step_want=index_tim.(time_want,res_upwind[i].time_lapse[1])
	  scatter!(plotc, res_upwind[i].h_range, 
      res_upwind[i].all_time[:,time_step_want], label="n=$(nRange1[i])" )	
	end
#n=$(nRange1[i]+1)	
    err=zeros(nr1-1,nRange1[1])
	h3= @. 1.0/(nRange1) #[0.1, 0.05, 0.025, 0.0125]# **h3= @. 1.0/(nRange1-1)
	#τ=res_upwind[1].time_lapse[2]
	plo3=Plots.scatter(axis= :log,xlabel="\$ h\$",ylabel="successive Error $var_name", title="solution convergence, ,grid line x=$ν,y=[0,1],dx=dy=h, time=$(time_want) s",legend= :bottomright)
	
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
	
	lable=["x=$ν,y=$ξ" for i in 1:1,ξ in show_x[1:end] ]#**2:end
	successive_Error=Plots.scatter(plo3, h3[1:end-1], err[:,2:end], label=lable )

	plot_slope=Plots.plot(xlabel="\$ h\$",ylabel="Slope $var_name ",title="solution Slope convergence, time=$(time_want) s",legend= :topleft)
    for i in 1:(n_show_x)#**length(nRange1)
    #for i in [3,5,7,9]#**length(nRange1)
		Plots.plot!(plot_slope,h3[1:end-2], slope(err[:,i],h3[1:end-1]),label="x=$ν,y=$(show_x[i])",lw=2)
    end
    return (solution=plotc,successive_Error=successive_Error,plot_slope=plot_slope)
end

# ╔═╡ 89e5d667-a88f-4b58-a72b-f4409773f819
x_line_want_L=0.75

# ╔═╡ c5f47587-0c32-448a-afaf-4e649f161814
index_cv(0.745,0.01)

# ╔═╡ 37a4f19b-0dfc-4e72-b377-e9cf2ed994dd
x_want_t=0.75

# ╔═╡ 8ce3dc4a-a935-49d8-9c04-e94bb5aa71d1
time1=[0.25, 0.5, 0.75,1].*0.04

# ╔═╡ 874f362a-6d6e-4847-96bc-aae026037290
index_tim.([0.010000000000000002, 0.020000000000000004, 0.03, 0.04000000000000001, 0.05],0.01)

# ╔═╡ ff989c7c-3e4b-44b7-936e-6f2054ce7865
t_Range1=[1,0.5, 0.25, 0.125]*0.01

# ╔═╡ d94f8b07-54a7-4688-865c-4268888e1cb0
index_tim.(time1,t_Range1[4])

# ╔═╡ f5ef87de-dc97-4d84-8459-3c50091f0b47
CFL_t=(2.2*1.2).*t_Range1./0.01

# ╔═╡ e06079e2-39c5-441a-8542-a8836ca5972e


# ╔═╡ 1da4f341-b1d2-4c1e-99c2-b4482ce8ae9e
begin


n=10
nx =n;ny = n
x_left = 0.0 ;x_right = 1.0 ;y_bottom = 0.0 ;y_top = 1.0
L_x=x_right-x_left;L_y=y_top-y_bottom
dx = (L_x)/nx ; dy = (L_y)/ny
##dt = time element;tf = end time;nt =
dt = 0.01;end_time = 0.05;nt = Int64(end_time/dt)

##
re = 100.0
#x_range=LinRange(0,L_x,nx+1)
#y_range=LinRange(0,L_y,ny+1)
x_range=LinRange(0,L_x-dx,nx).+0.5*dx;
y_range=LinRange(0,L_y-dy,ny).+0.5*dy

alfa=35;y0=0;x0=0.4;
b_tp=[2,1]# type bc [L,R]


ns = nt
ds = end_time/ns
r = Array{Float64}(undef, nx,3)
q = zeros(Float64,nx,3,ns+1)
#numerical(nx,ns,nt,dx,dt,q)
# numerical(nx,ns,nt,dx,dt,q,1)
qn = Array{Float64}(undef, nx,3)
gamma = 1.4 # specific gas ratio
x = Array{Float64}(undef, nx)
for i = 1:nx
	x[i] = -0.5*dx + dx*(i)
end
	
q0=Initial_conditions(nx,gamma,qn,x,dx)

#q2=callcu_1D_1time(nx,dx,dt,q0,b_tp)

x = Array(0.5*dx:dx:1.0-0.5*dx)

#qn=Initial_conditions(nx,gamma,q,dx)
end

# ╔═╡ c02b74a9-f983-4c04-ac15-e52772bfa952
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


# ╔═╡ 2b6baee4-8246-4e5a-aeb1-2542b4adb0b2
left_bc=b_c(b_tp[1])

# ╔═╡ ec2e868a-9d90-456f-97a4-1b7e8a51c9e4
y_range2=[Array(0.5*dy:dy:1.0-0.5*dy)']

# ╔═╡ 302c3656-18ee-48c5-90bb-5c6becf99276
y_range3=LinRange(0,L_y-dy,ny).+0.5*dy

# ╔═╡ eede2ace-aa82-4fda-b0e0-cf8d6ba9f7a0
(rho_x_line=(all_time=rho_x_line , time_lapse=time_lapse', h_range=y_range),rho_u_x_line=(all_time=rho_u_x_line , time_lapse=time_lapse', h_range=y_range),rho_v_x_line=(all_time=rho_v_x_line , time_lapse=time_lapse', h_range=y_range),rho_e_x_line=(all_time=rho_e_x_line , time_lapse=time_lapse', h_range=y_range))

# ╔═╡ 1ec4b988-4336-4a2a-a816-c0b1cd1edebf
time_want=end_time

# ╔═╡ 720a9c0c-7426-40fd-872f-8b78f514c7ab
nRange1=[1,3,9,27].*(n)

# ╔═╡ 79d5bb94-9fc1-4de2-bd3e-e37aa72aae5c
(1)./nRange1

# ╔═╡ 2ec8e3aa-edbb-4666-a003-c7ef9d918c7f
CFL_x=(2.2*1.2).*0.001./((1)./nRange1)

# ╔═╡ 70a4ae48-16e6-4025-8ecf-7f8f8b0ef7ef
t_Range=[1,0.5, 0.25, 0.125]*dt

# ╔═╡ 0a62bca5-6baf-4247-a4b7-30105699d864
n_ν=length(t_Range)

# ╔═╡ 32d0cb63-9623-4ff1-8226-76d213b9fce9
time=[0.2, 0.4, 0.6, 0.8,1].*end_time

# ╔═╡ 05f57bcf-b9a6-49af-8844-56943d2ce285
plot_t_seris(time, t_Range, res_upwind_t, x_want_t,"method_name")

# ╔═╡ 93e251b4-b2a2-4479-809b-bedd4a638d43
s=(1,1)

# ╔═╡ d96fbabb-bdae-4c22-8b3b-18d146173834
v=([(s=s,b=dti) for dti in t_Range ])

# ╔═╡ c9e2a93b-8e91-47ac-a05a-a5f820da0156
tuple(v[2:end]...,5)

# ╔═╡ d70d2141-67bc-4cec-b387-bf43b89cfc97
for i=1:2
	b1=2
b2=tuple(b1,b2...)
end

# ╔═╡ 4c37c065-08f7-41a2-95c3-eb8360cc7101
res_upwind=[generate_all_time_for_1_line(ni,L_y,L_x,end_time,dt,alfa,y0,x0,gamma,x_line_want_L) for ni in nRange1 ]

# ╔═╡ 539e6fe0-afe1-4cf3-9a89-de551b58b76a
o=plot_x_seris(x_line_want_L, res_upwind, time_want, nRange1,"rho")

# ╔═╡ 172736e9-0870-4ac5-85dc-f83d985d60a0
"""
nRange=[1,2,4,8].*(n).+1
print(nRange1[1]+1)
res_upwind= [generate_all_time_for_1_line(ni,L_y,L_x,dt,end_time,re,wn,sn,rms,x_line_want_L) for ni in nRange1 ]
print(res_upwind[1].time_lapse[1])
print(res_upwind[1].all_time[:,1])
"""
#numerical(nx,ny,nt,dx,dy,dt,re,wn,sn,rms)

# ╔═╡ 5cb856fb-1aac-468d-b73f-b6d6a4dd0ecc


# ╔═╡ 7977faf1-740a-4b9d-9f71-5956b3332954


# ╔═╡ 152986b8-0a83-4a6b-9dbc-a560de9c2d89
x

# ╔═╡ 9ab3c5f0-71f2-4353-8f12-9acfa3f8b8a0
xg = Array(0.5*dx:dx:1.0-0.5*dx)

# ╔═╡ 817ae3fa-c2b0-43b2-b485-539d3fa189bc
A2=rhs(nx,dx,gamma,q0,b_tp)

# ╔═╡ 857ae133-12fe-4892-a9ad-26ab85b0aab6
RK3(nx,gamma,q0,dx,dt,b_tp)

# ╔═╡ 5ad64159-d980-4b0e-927d-44cca6c60b5c
tan(15/180*3.14)

# ╔═╡ f88c0bb2-ff66-4927-b689-a0a377346c6b
index_x_yline(2,0.1,0.1,15,0,0.4,10,10)

# ╔═╡ 0e6eddfe-1ab4-431d-bdc8-9067e7fe4a8e
index_y_xline(8,0.1,0.1,15,0.4,0,10,10)

# ╔═╡ acc0d931-110b-4efc-b814-1385a2e22124
D2_q=Initial_conditions_2D_for_all_line(dy,dx,alfa,y0,x0,ny,nx,gamma)

# ╔═╡ 2236232e-5d68-4b2c-ae4e-7fbf721d5602
"""begin
#D2_q= Array{Float64}(zeros,ny,nx,4)
D2_q=zeros(Float64,nx,nx,4)
for j=1:ny
D2_q_1y = D2_q[j,:,:]
	
nx_end_yline=index_x_yline(j,dy,dx,15,0,0.4,ny,nx)
	
D2_q_1y=Initial_conditions_2D_for_1line(nx_end_yline,gamma,D2_q_1y,x)
D2_q[j,:,:]=D2_q_1y 
end
end"""

# ╔═╡ 92518268-a8ef-4fba-900a-8b2b83d17104
# ╠═╡ disabled = true
#=╠═╡
D2_q3[:,:,1]=D2_q2[:,:,1]
D2_q3[:,:,2]=D2_q2[:,:,4]
D2_q3[:,:,3]=D2_q2[:,:,3]




  ╠═╡ =#

# ╔═╡ d359a1c2-2284-4785-ba46-18cc4b51a02b
#ny_start_xlin=index_y_xline(1,dy,dx,15,0,0.4,nx,nx)

# ╔═╡ 34d3b2e0-8fa9-4bf4-8847-ad6f4ae439c8
#num_y=ny-ny_start_xline+1

# ╔═╡ 5b5dc60d-1529-45ce-8b78-d0567ed498cb
#D2_q_1x = D2_q[ny_start_xline:end,10,:]

# ╔═╡ b1f783a9-e24e-44ba-b093-ccae670f0237
#Initial_conditions_2D_for_1line(num_y,gamma,D2_q_1x,x)

# ╔═╡ 24957064-d6ea-41f5-9d06-dbe917aa3caa
"""
begin
#D2_q= Array{Float64}(undef,ny,nx,3)
for i=1:nx

ny_start_xline=index_y_xline(i,dy,dx,15,0,0.4,ny,nx)
num_y=ny-ny_start_xline+1
D2_q_1x = D2_q[ny_start_xline:end,i,:]
D2_q_1x=Initial_conditions_2D_for_1line(num_y,gamma,D2_q_1x,x)
D2_q[ny_start_xline:end,i,:]=D2_q_1x
end
end
"""

# ╔═╡ 1e43fb15-5c63-4e7f-8b93-f403b22c5b39
D2_q[:,:,[1,3]]

# ╔═╡ e54ca7cd-0bc4-4a34-aecc-32f7f5cc23fc
D2_q

# ╔═╡ 96955631-b4c6-4770-9e09-bb1ded0631f4
#nx_end_yline=index_x_yline(1,dy,dx,alfa,y0,x0,ny,nx)

# ╔═╡ 9fc21ad4-ecd2-4477-9763-c89bc9b813a7
#D2_q_1y_0 = D2_q[1,1:nx_end_yline,[1,2,4]]

# ╔═╡ b2c42d36-04cc-42d7-a512-386a786a42f0
#time_1step_all_yline(D2_q,dy,dx,alfa,y0,x0,ny,nx,gamma)

# ╔═╡ 72f76f4f-d341-40bf-92b9-86f3afbb669a
t= [0:2:10]

# ╔═╡ a6daba88-eb8b-4b6a-a70d-37ca4570200c
generate_all_time(D2_q,dy,dx,dt,alfa,y0,x0,ny,nx,nt,gamma)

# ╔═╡ 87cb5bb1-fdbc-47d5-af80-60926405513a


# ╔═╡ bcf5fafa-eb57-4de3-b277-3021c941e50c

begin
#D2_q_stor=zeros(Float64,nx,nx,4)
D2_q_stor=Array{Float64}(undef,ny,nx,4,nt)
D2_q1=Initial_conditions_2D_for_all_line(dy,dx,alfa,y0,x0,ny,nx,gamma)
	
for t= 0:2:nt-2

D2_q1=time_1step_all_yline(D2_q1,dy,dx,dt,alfa,y0,x0,ny,nx,gamma)
D2_q_stor[:,:,:,t+1]=D2_q1
D2_q1=time_1step_all_xline(D2_q1,dy,dx,dt,alfa,y0,x0,ny,nx,gamma)
D2_q_stor[:,:,:,t+2]=D2_q1
end

end

# ╔═╡ ac605695-7383-48a1-bde5-b825a856d049
D2_q1

# ╔═╡ a3ebb31c-d426-4f5c-8965-1d18853dd2e6
D2_q_stor

# ╔═╡ c80016af-b49a-48c4-9d7e-3e4d551695cf
(D2_q_stor[3,:,2,:])'

# ╔═╡ 79a59406-fc78-4c51-b67a-8a45f89ae0dd
for t= 0:2:nt-2
	print(t)
end

# ╔═╡ f505589c-3e20-4869-9761-b291c01bce50
D2_q1

# ╔═╡ 9d8e4e9f-76c6-4e43-b63f-81979b2ab513

begin
#D2_q= Array{Float64}(zeros,ny,nx,4)
D2_q_test=zeros(Float64,nx,nx,4)#[rho,rho*u,rho*v,rho*e]
#for j=1:ny
j=3
nx_end_yline=index_x_yline(j,dy,dx,alfa,y0,x0,ny,nx)
D2_q_1y_0 = D2_q_test[j,1:nx_end_yline,[1,2,4]] #[rho,rho*u,rho*e]
	if nx_end_yline==nx 
	   b_tp_end_x=[2,1] else b_tp_end_x=[2,0]
	end
D2_q_1y_1 = RK3(nx_end_yline,gamma,D2_q_1y_0,dx,dt,b_tp_end_x)
D2_q_test[j,1:nx_end_yline,[1,2,4]]=D2_q_1y_1

	
#D2_q_1y=Initial_conditions_2D_for_1line(nx_end_yline,gamma,D2_q_1y,x)
 

 """
D2_q3[:,:,1]=D2_q2[:,:,1]
D2_q3[:,:,2]=-D2_q2[:,:,3]
D2_q3[:,:,3]=D2_q2[:,:,4]
 
D2_q3[:,:,1]=D2_q2[:,:,1]
D2_q3[:,:,3]=-D2_q2[:,:,2]
D2_q3[:,:,4]=D2_q2[:,:,3]
""" 
end


# ╔═╡ 2c08f977-678d-4a27-b8e8-ef5feb67f00c


# ╔═╡ f4ee884f-2ac9-447e-a842-e3b9059dd768
nt

# ╔═╡ 4bec6d83-6fc7-49f9-969d-0359149765f7
begin
#D2_q= Array{Float64}(undef,ny,nx,3)
#for i=1:nx
i=9
ny_start_xline=index_y_xline(i,dy,dx,alfa,y0,x0,ny,nx)
num_y=ny-ny_start_xline+1
D2_q_1x_0 = D2_q_test[ny_start_xline:end,i,[1,3,4]] #[rho,rho*v,rho*e]
 	if (i*dx)<x0 && ny_start_xline==1
	   b_tp_end_y=[0,1] else b_tp_end_y=[0,1]
	end
D2_q_1x_1 = RK3(num_y,gamma,D2_q_1x_0,dx,dt,b_tp_end_y)
D2_q_test[ny_start_xline:end,i,[1,3,4]]=D2_q_1x_1

end

# ╔═╡ ad912c9c-97f7-47b1-871d-dc634e17b4e0
x0

# ╔═╡ cadbed3b-c56b-46e9-9bfb-0059410d0a33
i*dx

# ╔═╡ b0d8d7c4-1d5c-419b-903b-2344cd2d46b6
D2_q_1x_1

# ╔═╡ 93f08a58-e64c-4adc-a42f-ffd9abc23f43
ny_start_xline

# ╔═╡ e43aaf77-96d9-4d1a-90b4-621ea6b6f810
D2_q_test

# ╔═╡ d13c587b-bedd-44ec-957c-f268038c74ec
 numerical(nx,ns,nt,dx,dt,q,1)

# ╔═╡ fddf6bb8-a0f6-4a5a-9ff1-cb06fd1f4425
 for n = 1:5
 print(n)
 end

# ╔═╡ 4c0184e3-6db1-42bb-b2da-64ad2fdb4f97
function internal_coeffs_upwind!(I, J, V, n, ν)
	for i= 2:(n-1)
		K=3*(i-2).+(1:3)# [3*n+1 ][3*o+1, 3*o+2, 3*o+3]
		I[K] .= i#[i, i, i]
		J[K] = [i-1, i,i+1]
		V[K] = [ν, 1-ν,0]
		#C[i] = c3*τ*T_inf/c1
	end
	#return I, J, V#, C
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CPUTime = "a9c8d775-2e2e-55fc-8582-045d282d599e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
CPUTime = "~1.0.0"
Plots = "~1.31.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CPUTime]]
git-tree-sha1 = "2dcc50ea6a0a1ef6440d6eecd0fe3813e5671f45"
uuid = "a9c8d775-2e2e-55fc-8582-045d282d599e"
version = "1.0.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "80ca332f6dcb2508adba68f22f551adb2d00a624"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.3"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "924cdca592bc16f14d2f7006754a621735280b74"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.1.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c5544d8abb854e306b7b2f799ab31cdba527ccae"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.0"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "ccd479984c7838684b3ac204b716c89955c76623"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "037a1ca47e8a5989cc07d19729567bb71bfabd0c"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.66.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "c8ab731c9127cd931c93221f65d6a1008dad7256"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.66.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ed47af35905b7cc8f1a522ca684b35a212269bd8"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.2.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "1a43be956d433b5d0321197150c2f94e16c0aaa0"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.16"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7c88f63f9f0eb5929f15695af9a4d7d3ed278a91"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.16"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "9f4f5a42de3300439cb8300236925670f844a555"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.1"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9888e59493658e476d3073f1ce24348bdc086660"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "5a1e85f3aed2e0d3d99a4068037c8582597b89cf"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.31.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "2690681814016887462cf5ac37102b51cd9ec781"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "22c5201127d7b243b9ee1de3b43c408879dff60f"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.3.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "23368a3313d12a2326ad0035f0db0c0966f438ef"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "66fe9eb253f910fe8cf161953880cfdaef01cdf0"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.0.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2c11d7290036fe7aac9038ff312d3b3a2a5bf89e"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.4.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "472d044a1c8df2b062b23f222573ad6837a615ba"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.19"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "ec47fb6069c57f1cee2f67541bf8f23415146de7"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.11"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═78610b9b-43c3-490d-b04e-a0a95b357d7b
# ╠═ce187bce-2c98-4fad-8b84-15b249154a8e
# ╟─c02b74a9-f983-4c04-ac15-e52772bfa952
# ╟─61aa2c26-0d88-4782-98c8-2acf857d8eb7
# ╟─2c360e6d-ff1f-4653-8c9a-4fc650780926
# ╠═44d05ac3-0ae3-430f-a755-fac8aff0adef
# ╟─bc4a9a50-09c0-11ed-1051-0d7888c35c8a
# ╟─de66199c-155c-4e23-a2d6-ff665237037f
# ╠═f6b5833a-a39c-491a-9218-cd0c956c5d05
# ╟─5bc2f630-f76b-4e11-8901-6b6d6a1018a4
# ╟─e95f8105-ab89-45a3-8324-e8bd838f654a
# ╟─5b5726e6-d600-4f3b-af8e-a3c58c8cde9f
# ╠═6e74ea04-2b8d-47a5-91e4-b437796303c9
# ╠═08fc3d88-2d40-4b48-8434-b9fad1ab0b7c
# ╠═ae7eff58-e358-40b0-ba31-e85507e641eb
# ╠═0f3a4c00-666f-46ac-8547-d728ef99a10e
# ╠═2b6baee4-8246-4e5a-aeb1-2542b4adb0b2
# ╟─afd16e84-96b3-4715-a712-f0d8fc19b4c8
# ╟─679bf7ab-5f35-4315-8dcc-0b52c3498f84
# ╟─98832477-0c6c-4770-8f45-74e244b7c42b
# ╠═73cf2abf-ed39-4ef0-9895-c736677a9be7
# ╟─5500059d-1aed-43c5-846f-fd56e8bef8c2
# ╟─08bca7eb-ae0e-4faa-9104-faa6fcef2cb4
# ╟─582ae9c6-c502-4430-b7cf-2d1f6ab782bf
# ╟─a0f08d8f-9ef3-4c25-a970-bc973942caa2
# ╟─74c5b2d9-0e8c-4992-8f75-060bf5b41b75
# ╟─175c19f6-d11e-4c3c-adf3-b024846da01b
# ╟─0c85b3fa-1e00-4171-951a-24619a0bcc0b
# ╟─1a908140-b336-47db-9c54-6649103d1ec2
# ╠═942d57b9-0ade-4813-9298-292107d1215e
# ╠═fbfa7298-6f8f-4011-87bb-03c32a6296f7
# ╠═56fcad51-0146-4a24-b8d6-231f8d72f8f5
# ╠═ec2e868a-9d90-456f-97a4-1b7e8a51c9e4
# ╠═302c3656-18ee-48c5-90bb-5c6becf99276
# ╟─1b373135-757d-4286-add8-7ab0d6a1b3cf
# ╠═eede2ace-aa82-4fda-b0e0-cf8d6ba9f7a0
# ╠═7e178c53-14ab-4217-b551-cf0db73579ca
# ╟─3bab86a7-37ff-4e98-b0b2-e4450a9ec126
# ╟─ce673755-7701-4b33-b0fb-61ae7991731f
# ╠═89e5d667-a88f-4b58-a72b-f4409773f819
# ╠═c5f47587-0c32-448a-afaf-4e649f161814
# ╠═1ec4b988-4336-4a2a-a816-c0b1cd1edebf
# ╠═720a9c0c-7426-40fd-872f-8b78f514c7ab
# ╠═79d5bb94-9fc1-4de2-bd3e-e37aa72aae5c
# ╠═2ec8e3aa-edbb-4666-a003-c7ef9d918c7f
# ╠═539e6fe0-afe1-4cf3-9a89-de551b58b76a
# ╠═37a4f19b-0dfc-4e72-b377-e9cf2ed994dd
# ╠═70a4ae48-16e6-4025-8ecf-7f8f8b0ef7ef
# ╠═0a62bca5-6baf-4247-a4b7-30105699d864
# ╠═8ce3dc4a-a935-49d8-9c04-e94bb5aa71d1
# ╠═874f362a-6d6e-4847-96bc-aae026037290
# ╠═d94f8b07-54a7-4688-865c-4268888e1cb0
# ╠═ff989c7c-3e4b-44b7-936e-6f2054ce7865
# ╠═f5ef87de-dc97-4d84-8459-3c50091f0b47
# ╠═e06079e2-39c5-441a-8542-a8836ca5972e
# ╠═32d0cb63-9623-4ff1-8226-76d213b9fce9
# ╠═05f57bcf-b9a6-49af-8844-56943d2ce285
# ╠═1da4f341-b1d2-4c1e-99c2-b4482ce8ae9e
# ╠═93e251b4-b2a2-4479-809b-bedd4a638d43
# ╠═d96fbabb-bdae-4c22-8b3b-18d146173834
# ╠═c9e2a93b-8e91-47ac-a05a-a5f820da0156
# ╠═d70d2141-67bc-4cec-b387-bf43b89cfc97
# ╠═4c37c065-08f7-41a2-95c3-eb8360cc7101
# ╠═172736e9-0870-4ac5-85dc-f83d985d60a0
# ╠═5cb856fb-1aac-468d-b73f-b6d6a4dd0ecc
# ╠═7977faf1-740a-4b9d-9f71-5956b3332954
# ╠═152986b8-0a83-4a6b-9dbc-a560de9c2d89
# ╠═9ab3c5f0-71f2-4353-8f12-9acfa3f8b8a0
# ╠═817ae3fa-c2b0-43b2-b485-539d3fa189bc
# ╠═857ae133-12fe-4892-a9ad-26ab85b0aab6
# ╠═5ad64159-d980-4b0e-927d-44cca6c60b5c
# ╠═f88c0bb2-ff66-4927-b689-a0a377346c6b
# ╟─0e6eddfe-1ab4-431d-bdc8-9067e7fe4a8e
# ╠═acc0d931-110b-4efc-b814-1385a2e22124
# ╠═2236232e-5d68-4b2c-ae4e-7fbf721d5602
# ╠═92518268-a8ef-4fba-900a-8b2b83d17104
# ╠═d359a1c2-2284-4785-ba46-18cc4b51a02b
# ╠═34d3b2e0-8fa9-4bf4-8847-ad6f4ae439c8
# ╠═5b5dc60d-1529-45ce-8b78-d0567ed498cb
# ╠═b1f783a9-e24e-44ba-b093-ccae670f0237
# ╠═24957064-d6ea-41f5-9d06-dbe917aa3caa
# ╠═1e43fb15-5c63-4e7f-8b93-f403b22c5b39
# ╠═e54ca7cd-0bc4-4a34-aecc-32f7f5cc23fc
# ╠═96955631-b4c6-4770-9e09-bb1ded0631f4
# ╠═9fc21ad4-ecd2-4477-9763-c89bc9b813a7
# ╠═b2c42d36-04cc-42d7-a512-386a786a42f0
# ╠═72f76f4f-d341-40bf-92b9-86f3afbb669a
# ╠═a6daba88-eb8b-4b6a-a70d-37ca4570200c
# ╠═87cb5bb1-fdbc-47d5-af80-60926405513a
# ╠═bcf5fafa-eb57-4de3-b277-3021c941e50c
# ╠═ac605695-7383-48a1-bde5-b825a856d049
# ╠═a3ebb31c-d426-4f5c-8965-1d18853dd2e6
# ╠═c80016af-b49a-48c4-9d7e-3e4d551695cf
# ╠═79a59406-fc78-4c51-b67a-8a45f89ae0dd
# ╠═f505589c-3e20-4869-9761-b291c01bce50
# ╟─9d8e4e9f-76c6-4e43-b63f-81979b2ab513
# ╠═2c08f977-678d-4a27-b8e8-ef5feb67f00c
# ╠═f4ee884f-2ac9-447e-a842-e3b9059dd768
# ╠═4bec6d83-6fc7-49f9-969d-0359149765f7
# ╠═ad912c9c-97f7-47b1-871d-dc634e17b4e0
# ╠═cadbed3b-c56b-46e9-9bfb-0059410d0a33
# ╠═b0d8d7c4-1d5c-419b-903b-2344cd2d46b6
# ╠═93f08a58-e64c-4adc-a42f-ffd9abc23f43
# ╠═e43aaf77-96d9-4d1a-90b4-621ea6b6f810
# ╠═d13c587b-bedd-44ec-957c-f268038c74ec
# ╟─fddf6bb8-a0f6-4a5a-9ff1-cb06fd1f4425
# ╠═4c0184e3-6db1-42bb-b2da-64ad2fdb4f97
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
