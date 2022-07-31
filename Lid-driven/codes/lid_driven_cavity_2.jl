include("lid_driven_cavity_fun.jl")


#---------------------------------------------------------------------------#
# main program
#---------------------------------------------------------------------------#
n=64
nx =n;ny = n
x_left = 0.0 ;x_right = 1.0 ;y_bottom = 0.0 ;y_top = 1.0
L_x=x_right-x_left;L_y=y_top-y_bottom
dx = (L_x)/nx ; dy = (L_y)/ny
##dt = time element;tf = end time;nt =
dt = 0.001;end_time = 10;nt = Int64(end_time/dt)
##
re =100
x_range=LinRange(0,L_x,nx+1)
y_range=LinRange(0,L_y,ny+1)


x = Array{Float64}(undef, nx+1)
y = Array{Float64}(undef, ny+1)
rms = Array{Float64}(undef, nt)

for i = 1:nx+1
    x[i] = dx*(i-1)
end
for i = 1:ny+1
    y[i] = dy*(i-1)
end

wn = Array{Float64}(undef, nx+1, ny+1)
sn = Array{Float64}(undef, nx+1, ny+1)

time = 0.0

for i = 1:nx+1 for j = 1:ny+1
    wn[i,j] = 0.0 # initial condition
    sn[i,j] = 0.0 # initial streamfunction
end end
test=fps_sine(nx,ny,dx,dy,-wn)

numerical_res=numerical(nx,ny,nt,dx,dy,dt,re,wn,sn,rms)

velocity=ditem_u_v(nx,ny,dx,dy,numerical_res.end_sn)

P=ditem_P(nx,ny,dx,dy,velocity.v_P,velocity.u_P,re)
#------------------------------------------------------------------------------#
#plot part : must run after total run by shift + Enter
#------------------------------------------------------------------------------#
#Pyplot test conture
using PyPlot;
close("all");pygui(true)
PyPlot.figure()
PyPlot.contour(x, y, numerical_res.end_wn,100)
PyPlot.xlabel("Y")
PyPlot.ylabel("X")
PyPlot.colorbar()
PyPlot.show()

stream_fun_contour = PlotlyJS.contour(x=x,y=y,z=numerical_res.end_sn)
layout_Ψ = Layout(;xaxis_title="Y",yaxis_title="X",font_size=14,title="Ψ,Re:$re,grid:$n*$n",titlefont_size=14)
contour_Ψ=PlotlyJS.plot(stream_fun_contour,layout_Ψ)

vorticity_contour = PlotlyJS.contour(x=x,y=y,z=numerical_res.end_wn)
layout_ω = Layout(;xaxis_title="Y",yaxis_title="X",font_size=14,title="ω,Re:$re,grid:$n*$n",titlefont_size=14)
contour_ω=PlotlyJS.plot(vorticity_contour,layout_ω)

#print(velocity.u)
u_contour = PlotlyJS.contour(x=x,y=y,z=(velocity.u))
layout_u = Layout(;xaxis_title="Y",yaxis_title="X",font_size=14,title="u,Re:$re,grid:$n*$n",titlefont_size=14)
contour_u=PlotlyJS.plot(u_contour,layout_u)
# v cantor
v_contour = PlotlyJS.contour(x=x,y=y,z=(velocity.v))
layout_v = Layout(;xaxis_title="Y",yaxis_title="X",font_size=14,title="v,Re:$re,grid:$n*$n",titlefont_size=14)
contour_v=PlotlyJS.plot(v_contour,layout_v)

P_contour = PlotlyJS.contour(x=x,y=y,z=(P))
layout_P = Layout(;xaxis_title="Y",yaxis_title="X",font_size=14,title="P,Re:$re,grid:$n*$n",titlefont_size=14)
PlotlyJS.plot(P_contour,layout_P)
