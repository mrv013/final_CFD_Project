include("P_final_s2_fun.jl")


##
n=100
nx =n;ny = n
x_left = 0.0 ;x_right = 1.0 ;y_bottom = 0.0 ;y_top = 1.0
L_x=x_right-x_left;L_y=y_top-y_bottom
dx = (L_x)/nx ; dy = (L_y)/ny
##dt = time element;tf = end time;nt =
dt = 0.001;end_time = 1;nt = Int64(end_time/dt)
##
re = 100.0
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

#time error analysis
x_line_want_L=0.305
y_want_t=0.745
t_Range=[1,0.5, 0.25, 0.125]*dt
time=[0.2, 0.4, 0.6, 0.8,1].*end_time
res_solve_t_ρ=[generate_all_time_for_1_line(n,L_y,L_x,end_time,dti,alfa,y0,x0,gamma,x_line_want_L)[1] for dti in t_Range]

P_t_seris_ρ=plot_t_seris(time, t_Range, res_solve_t_ρ, y_want_t,x_line_want_L,"ρ")
P_t_seris_ρ.slope_t
P_t_seris_ρ.successive_Error

res_solve_t_ρu=[generate_all_time_for_1_line(n,L_y,L_x,end_time,dti,alfa,y0,x0,gamma,x_line_want_L)[2] for dti in t_Range]

P_t_seris_ρu=plot_t_seris(time, t_Range, res_solve_t_ρu, y_want_t,x_line_want_L,"ρu")
P_t_seris_ρu.slope_t
P_t_seris_ρu.successive_Error


res_solve_t_ρv=[generate_all_time_for_1_line(n,L_y,L_x,end_time,dti,alfa,y0,x0,gamma,x_line_want_L)[3] for dti in t_Range]
P_t_seris_ρv=plot_t_seris(time, t_Range, res_solve_t_ρv, y_want_t,x_line_want_L,"ρv")
P_t_seris_ρv.slope_t
P_t_seris_ρv.successive_Error

res_solve_t_ρe=[generate_all_time_for_1_line(n,L_y,L_x,end_time,dti,alfa,y0,x0,gamma,x_line_want_L)[4] for dti in t_Range]
P_t_seris_ρe=plot_t_seris(time, t_Range, res_solve_t_ρe, y_want_t,x_line_want_L,"ρe")
P_t_seris_ρe.slope_t
P_t_seris_ρe.successive_Error

"""
#grid error analysis
x_line_want_L=0.75
time_want=end_time
nRange1=[1,5,25,125].*(n)
res_solve_x=[generate_all_time_for_1_line(ni,L_y,L_x,end_time,dt,alfa,y0,x0,gamma,x_line_want_L) for ni in nRange1 ]
P_x_seris=plot_x_seris(x_line_want_L, res_solve_x, time_want, nRange1,"rho")
P_x_seris.solution
P_x_seris.slope_t
P_x_seris.successive_Error
"""
