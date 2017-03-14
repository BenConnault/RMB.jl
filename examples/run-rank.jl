reload("RMB")

module dev

using RMB
import RMB: inventorymodel, pidc1, pidc2, detensor, rank, pi!, rsm

dims=(5,3,10)
m2=inventorymodel(dims,5.0,20.0,[0.8,1.0,1.2])
pi!(m2,pidc1(m2))
u!(m2,[0.0,1.0])    #free stocking, expensive order
solvedp!(m2)


dd,dc,ds,dy,dr=m2.dims
dk,da=dd*dc*ds,dy*dr

ac=reshape([1.0*(ik==jk)*(ia==1) for ik=1:dk,jk=1:dk,ia=1:da],dk,dk*da)

function pi!(model::InventoryModel,eta::AbstractVector)
    pidc=RMB.z2q(eta)
    pi!(model,pidc)
end


eta=rand(dd*dc*(dd*dc-1))
j1=RMB.idrank(m2,ac,eta)


# pic1=rsm(dc,dc)
# pic2=rsm(dc,dc)
# pid1=rsm(dd,dd)
# pid2=rsm(dd,dd)

# function pi!(model::InventoryModel,eta::AbstractVector)
#     pic=eta[1]*pic1+(1-eta[1])*pic2
#     pid=eta[2]*pid1+(1-eta[2])*pid2
#     pi!(model,kron(pic,pid))
# end

# eta=[0.5,0.5]
# j2=RMB.idrank(m2,eta)


end #module


