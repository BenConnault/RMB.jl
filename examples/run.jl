module dev


using RMB
import RMB: inventorymodel, pidc1, pidc2, detensor, explicitinv

include("RMBplots")



# inventorymodel((dd,dc,ds),pricelevel,spotbuy,wholesaleprices)
# dims=(5,3,10)

# m1=inventorymodel(dims,5.0,20.0,[0.1,1.0,1.1])
# pi!(m1,pidc1(m1))
# # pid,pic=detensor(m1)   #for debugging
# u!(m1,[0,0.0])    #free stocking, free order
# solvedp!(m1)
# ka=rand(m1,200)
# tsplot(m1,ka)

# m2=inventorymodel(dims,5.0,20.0,[0.8,1.0,1.2])
# pi!(m2,pidc1(m2))
# u!(m2,[0.0,1.0])    #free stocking, expensive order
# solvedp!(m2)
# ka=rand(m2,200)
# tsplot(m2,ka)


# m3=inventorymodel(dims,8.0,20.0,[0.90,0.9,0.9])
# pi!(m3,pidc2(m3))
# u!(m3,[0.05,0.1])    #free stocking, free order
# solvedp!(m3)
# ka=rand(m3,300)
# tsplot(m3,ka)


end #module


