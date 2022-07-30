include("P_final_s2_fun.jl")


##
n=10
nx =n;ny = n
x_left = 0.0 ;x_right = 1.0 ;y_bottom = 0.0 ;y_top = 1.0
L_x=x_right-x_left;L_y=y_top-y_bottom
dx = (L_x)/nx ; dy = (L_y)/ny
##dt = time element;tf = end time;nt =
dt = 0.001000;end_time = 1;nt =index_tim(end_time,dt)# Int64(end_time/dt)
##
re = 10.0
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


#grid error analysis
x_line_want_L=0.55
time_want=0.3
y_want_t=[0.25,0.35,0.45,0.55]#,0.65,0.75,0.85,0.95]
nRange1=[1,3,9,27].*(n)

res_solve_x=[generate_all_time_for_1_line(ni,L_y,L_x,end_time,dt,alfa,y0,x0,gamma,x_line_want_L)[1] for ni in nRange1 ]

P_x_seris=plot_x_seris(y_want_t,x_line_want_L, res_solve_x, time_want, nRange1,"rho")
P_x_seris.solution
P_x_seris.plot_slope
P_x_seris.successive_Error



#time error analysis
"""
P_x_seris_ρ=plot_x_seris(time, t_Range, res_solve_x, y_want_t,x_line_want_L,"ρ")
P_x_seris_ρ.plot_slope
P_x_seris_ρ.successive_Error

P_x_seris_ρu=plot_x_seris(time, t_Range, res_solve_x[2], y_want_t,x_line_want_L,"ρu")
P_x_seris_ρu.plot_slope
P_x_seris_ρu.successive_Error

P_x_seris_ρv=plot_x_seris(time, t_Range, res_solve_x[3], y_want_t,x_line_want_L,"ρv")
P_x_seris_ρv.plot_slope
P_x_seris_ρv.successive_Error

P_x_seris_ρe=plot_x_seris(time, t_Range, res_solve_x[4], y_want_t,x_line_want_L,"ρe")
P_x_seris_ρe.plot_slope
P_x_seris_ρe.successive_Error
"""
