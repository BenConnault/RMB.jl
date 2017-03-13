using Plots
plotly()


function tsplot(model::RMB.InventoryModel,ka)
    T=size(ka,1)
    dd,dc,ds,dy,dr=model.dims
    ddd=[ind2sub((dd,dc,ds),ka[t,1])[1]-1 for t=1:T]
    cc=[model.wholesaleprices[ind2sub((dd,dc,ds),ka[t,1])[2]] for t=1:T]
    ss=[ind2sub((dd,dc,ds),ka[t,1])[3]-1 for t=1:T]
    p=plot(layout = (3,1),size=(1000,400))
    plot!(p[1],1:T,ss,label="stocks")
    plot!(p[2],1:T,cc,label="wholesale prices",l=:green)
    plot!(p[3],1:T,ddd,label="demand",l=:orange)
    display(p)
end