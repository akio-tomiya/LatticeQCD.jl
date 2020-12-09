L = (4,4,4,4)
β = 6
NTRACE = 3
#gparam = Setup_Gauge_action(β)
gparam =  GaugeActionParam_standard(β,NTRACE)

BoundaryCondition=[1,1,1,-1]
Nwing = 1
initial="cold"
NC =3


hop= 0.141139#Hopping parameter
r= 1#Wilson term
eps= 1e-19
Dirac_operator= "Wilson"
MaxCGstep= 3000

fparam = FermiActionParam_Wilson(hop,r,eps,Dirac_operator,MaxCGstep)