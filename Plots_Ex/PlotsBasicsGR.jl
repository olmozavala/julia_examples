## http://docs.juliaplots.org/latest/basics/
# ENV["GRDISPLAY"] = "qt5"
using GR

histogram(randn(1000))

x = linspace(-2,2,21)
y = x
X, Y = meshgrid(x,y)