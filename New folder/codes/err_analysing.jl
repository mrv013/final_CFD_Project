include("lid_driven_cavity_fun.jl")


# main program
#---------------------------------------------------------------------------#
n=15
nx =n;ny = n
x_left = 0.0 ;x_right = 1.0 ;y_bottom = 0.0 ;y_top = 1.0
L_x=x_right-x_left;L_y=y_top-y_bottom
dx = (L_x)/nx ; dy = (L_y)/ny
##dt = time element;tf = end time;nt =
dt = 0.001;end_time = 2.00;nt = Int64(end_time/dt)
##
re = 100.0
x_range=LinRange(0,L_x,nx+1)
y_range=LinRange(0,L_y,ny+1)
rms = Array{Float64}(undef, nt)

wn = Array{Float64}(undef, nx+1, ny+1)
sn = Array{Float64}(undef, nx+1, ny+1)

time = 0.0

for i = 1:nx+1 for j = 1:ny+1
    wn[i,j] = 0.0 # initial condition
    sn[i,j] = 0.0 # initial streamfunction
end end
test=fps_sine(nx,ny,dx,dy,-wn)
x_line_want_L=0.3
print(index(0.5,0.1))
#a=generate_all_time_for_1_line(n,L_y,L_x,dt,end_time,re,wn,sn,rms,x_line_want_L)

nRange1=[1,2,4,8].*(n)
nRange=[1,2,4,8].*(n).+1
print(nRange1[1]+1)
time_want=2

"""
res_upwind= [generate_all_time_for_1_line(ni,L_y,L_x,dt,end_time,re,wn,sn,rms,x_line_want_L) for ni in nRange1 ]
print(res_upwind[1].time_lapse[1])
print(res_upwind[1].all_time[:,1])
#numerical(nx,ny,nt,dx,dy,dt,re,wn,sn,rms)

# ν, res_upwind, time_want, nRange1,var_name
o=plot_x_seris(0.7, res_upwind, time_want, nRange1,"Ψ")
#solution=plotc, successive_Error=successive_Error, plot_slope=plot_slope
o.solution
o.plot_slope
o.successive_Error
"""
Re_Range=[1,10,100,500]
res_Re= [generate_all_time_for_1_line(n,L_y,L_x,dt,end_time,Rei,wn,sn,rms,x_line_want_L) for Rei in Re_range]

Q=plot_x_seris_re(0.7, res_Re, time_want, Re_Range,"Ψ")
Q.solution

"""
n=64
ν=5
time=[0.1, 0.2, 0.3, 0.4, 0.5]
nRange1=[1,2,4,8].*(n-1).+1
res_upwind= [generate_all_time_upwind(ni,inational_value!(L,ni,2),L, a, ν,end_time,bc) for ni in nRange1 ]
res_upwind, time_want, nRange1,method_name
"""
