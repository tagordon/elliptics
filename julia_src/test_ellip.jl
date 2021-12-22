
using PyPlot
include("ellip.jl")
import Elliptic

function elc(x,kc,p)
  el1int(xi,kc) = 1/sqrt((1+xi^2)*(1+(kc*xi)^2))
  el2int(xi,kc) = 1/(1+xi^2)/sqrt((1+xi^2)*(1+(kc*xi)^2))
  el3int(xi,kc,p) = (1+xi^2)/(1+p*xi^2)/sqrt((1+xi^2)*(1+(kc*xi)^2))
  nxi = 10000; dxi = x/nxi
  xi = collect(0.5:1.0:nxi-0.5)*dxi
  return sum(dxi*el1int.(xi,kc)),sum(dxi*el2int.(xi,kc)),sum(dxi*el3int.(xi,kc,p))
end

ntest = 100
el1grid = zeros(5,ntest)
el2grid = zeros(5,ntest)
el3grid = zeros(5,ntest)
for i=1:ntest
  xtest = 10^(-3+6*rand()); kctest = 10^(-3+6*rand())
  p = rand()
  el1_test,el2_test,el3_test = elc(xtest,kctest,p)
  el1_an = el1(xtest,kctest)
  #el1_an = Elliptic.F(atan(xtest),1-kctest^2)
  el2_an = el2(xtest,kctest,1.0,0.0)
  el3_an = el3(xtest,kctest,p)
  #el3_an = Elliptic.Pi(1-p,atan(xtest),1-kctest^2)
  el1grid[:,i] = [xtest,kctest,el1_an,el1_test,el1_an/el1_test-1.0]
  el2grid[:,i] = [xtest,kctest,el2_an,el2_test,el2_an/el2_test-1.0]
  el3grid[:,i] = [xtest,kctest,el3_an,el3_test,el3_an/el3_test-1.0]
  #println(xtest," ",kctest," ",el1_an," ",el1_test," ",el1_an/el1_test-1.0," ",el2_an,
  #      " ",el2_test," ",el2_an/el2_test-1.0," ",el3_an," ",el3_test," ",el3_an/el3_test-1.0)
end

clf()
igood = abs.(el1grid[5,:]) .< 1e-8
loglog(el1grid[1,igood],el1grid[2,igood],".")
ibad = abs.(el1grid[5,:]) .>= 1e-8
loglog(el1grid[1,ibad],el1grid[2,ibad],".")
read(stdin,Char)
igood = abs.(el2grid[5,:]) .< 1e-8
clf()
loglog(el2grid[1,igood],el2grid[2,igood],".")
ibad = abs.(el2grid[5,:]) .>= 1e-8
loglog(el2grid[1,ibad],el2grid[2,ibad],".")
read(stdin,Char)
igood = abs.(el3grid[5,:]) .< 1e-8
clf()
loglog(el3grid[1,igood],el3grid[2,igood],".")
ibad = abs.(el3grid[5,:]) .>= 1e-8
loglog(el3grid[1,ibad],el2grid[2,ibad],".")
