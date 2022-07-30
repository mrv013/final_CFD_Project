include("P_final_s2_fun.jl")

index_cv(ξ,Δξ)=floor(Int, (ξ+0.5*Δξ)/Δξ+0.5)
index_tim(ξ,Δξ)=floor(Int, (ξ+0.5*Δξ)/Δξ)
slope(e,h)=@. log(e[2:end]/e[1:end-1])/log(h[2:end]/h[1:end-1])

##
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

function plot_t_seris(time, ν_Range, res_upwind_t, x_want_t, nRange1,method_name)
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
return (all_time=rho_x_line , time_lapse=time_lapse', h_range=y_range)
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


##
n=200
nx =n;ny = n
x_left = 0.0 ;x_right = 1.0 ;y_bottom = 0.0 ;y_top = 1.0
L_x=x_right-x_left;L_y=y_top-y_bottom
dx = (L_x)/nx ; dy = (L_y)/ny
##dt = time element;tf = end time;nt =
dt = 0.0001;end_time = 1;nt = Int64(end_time/dt)
##
gamma=1.4
rho = 1.0
p = 1.0
gamma=1.4
a=sqrt(gamma*p/rho)
CFL=(1.2*a+a)*dt/dx
x_range=LinRange(0,L_x-dx,nx).+0.5*dx;
y_range=LinRange(0,L_y-dy,ny).+0.5*dy
##
alfa=15
y0=0;x0=0.4
##
b_tp=[2,1]# type bc [L,R]
##
D2_q_0=Initial_conditions_2D_for_all_line(dy,dx,alfa,y0,x0,ny,nx,gamma)#(dy,dx,alfa,y0,x0,ny,nx,gamma)
D2_q=generate_all_time(D2_q_0,dy,dx,dt,alfa,y0,x0,ny,nx,nt,gamma)

##
#print(p[:,:,1])

D2_q_0_contour_ρ = PlotlyJS.contour(x=x_range,y=y_range,z=D2_q_0[:,:,1])
layout_ρ0 = Layout(;xaxis_title="x",yaxis_title="y",font_size=14,title="ρ,Mach: α:$alfa,grid:$n*$n,dt:$dt,time:$end_time",titlefont_size=14)
PlotlyJS.plot(D2_q_0_contour_ρ,layout_ρ0)

using PyPlot;
close("all");pygui(true)
PyPlot.figure()
PyPlot.contour(x_range, y_range, D2_q.endtime[:,:,1],100)
PyPlot.xlabel("Y")
PyPlot.ylabel("X")
PyPlot.colorbar()
PyPlot.show()

D2_q_end_contour_1 = PlotlyJS.contour(x=x_range,y=y_range,z=D2_q.stor[:,:,1,8000])
layout_ρ = Layout(;xaxis_title="x",yaxis_title="y",font_size=14,title="ρ,Mach: α:$alfa,grid:$n*$n,dt:$dt,time:$end_time",titlefont_size=14)
PlotlyJS.plot(D2_q_end_contour_1,layout_ρ)



D2_q_end_contour_2 = PlotlyJS.contour(x=x_range,y=y_range,z=D2_q.endtime[:,:,2])
layout_ρu = Layout(;xaxis_title="x",yaxis_title="y",font_size=14,title="ρ*u,Mach: α:$alfa,grid:$n*$n,dt:$dt,time:$end_time",titlefont_size=14)
PlotlyJS.plot(D2_q_end_contour_2,layout_ρu)


D2_q_end_contour_3 = PlotlyJS.contour(x=x_range,y=y_range,z=D2_q.endtime[:,:,3])
layout_ρv = Layout(;xaxis_title="x",yaxis_title="y",font_size=14,title="ρ*v,Mach: α:$alfa,grid:$n*$n,dt:$dt,time:$end_time",titlefont_size=14)
PlotlyJS.plot(D2_q_end_contour_3,layout_ρv)

D2_q_end_contour_4 = PlotlyJS.contour(x=x_range,y=y_range,z=D2_q.endtime[:,:,4])
layout_e = Layout(;xaxis_title="x",yaxis_title="y",font_size=14,title="ρ*e,Mach: α:$alfa,grid:$n*$n,dt:$dt,time:$end_time",titlefont_size=14)
PlotlyJS.plot(D2_q_end_contour_4,layout_e)
#D2_q_for_1_line=generate_all_time_for_1_line(n,L_y,L_x,end_time,dt,alfa,y0,x0,gamma,x_line_want_L)


"""


#time error analysis
x_want_t=0.75
t_Range=[1,0.5, 0.25, 0.125]*dt
time=[0.2, 0.4, 0.6, 0.8,1].*end_time
res_solve_t=[generate_all_time_for_1_line(n,L_y,L_x,end_time,dti,alfa,y0,x0,gamma,x_line_want_L) for dti in t_Range ]
plot_t_seris(time, t_Range, res_solve_t, x_want_t, nRange1,"Y var name")


#grid error analysis
x_line_want_L=0.75
time_want=end_time
nRange1=[1,5,25,125].*(n)
res_solve_x=[generate_all_time_for_1_line(ni,L_y,L_x,end_time,dt,alfa,y0,x0,gamma,x_line_want_L) for ni in nRange1 ]
o=plot_x_seris(x_line_want_L, res_solve_x, time_want, nRange1,"rho")

"""

"""


Initial_conditions(nx,gamma,qn,x,dx)
D2_q= Array{Float64}(undef,ny,nx,4)
D2_q_1y = D2_q[j,:,1:3]
D2_q_1y=Initial_conditions(nx,gamma,D2_q_1y,x_range,dx)
D2_q[j,:,1:3]=D2_q_1y
#D2_qx = D2_q[:,i,1:2:4]
##
begin
nx = 10 # j
ny = 10 # i
dt = 0.001;tm = 0.005;

dx = 1.0/nx;
dy = 1.0/ny;
nt = Int64(tm/dt)
ns = nt
ds = tm/ns

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
r = Array{Float64}(undef, nx,3)
b_tp=1
#q2=callcu_1D_1time(nx,dx,dt,q0,b_tp)

x = Array(0.5*dx:dx:1.0-0.5*dx)

#qn=Initial_conditions(nx,gamma,q,dx)
end

D2_q= Array{Float64}(undef,nx,nx,3)

A2=rhs(nx,dx,gamma,q0,b_tp)
RK3(nx,gamma,q0,dx,dt,b_tp)
numerical(nx,ns,nt,dx,dt,q,1)
"""
