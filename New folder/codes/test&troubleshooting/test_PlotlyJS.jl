#using PlotlyJS

using PlotlyJS
"""
close("all");pygui(true)
PlotlyJS.plot(PlotlyJS.contour(
    z=[
        10      10.625      12.5       15.625     20
        5.625    6.25       8.125      11.25      15.625
        2.5      3.125      5.         8.125      12.5
        0.625    1.25       3.125      6.25       10.625
        0        0.625      2.5        5.625      10
    ]'
))

"""
using PyPlot
function test_PlotlyJs()
trace = PyPlot.contour(x, y, z)

#x = [-9, -6, -5 , -3, -1]
x =  [0, 2, 4, 6, 8]
    y = [0, 1, 3, 4, 5]
    z = [1  1  1  0  0
         1  1  1  0  0
         1  1  1  0  0
         0  0  0  0  0
         0  0  0  0  0]
    trace = PlotlyJS.contour(x=x, y=y, z=z)

    layout = Layout(title="Setting the X and Y Coordinates in a Contour Plot")
    PlotlyJS.plot(trace, layout)
trace = PyPlot.contour(x, y, z)
end
