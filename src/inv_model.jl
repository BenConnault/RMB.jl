type InventoryModel <: RustModel
    dims::NTuple{5,Int}             #dd,dc,ds,dy,dr=dims
    pricelevel::Float64
    spotbuy::Float64
    wholesaleprices::Vector{Float64}
    uaa::Array{Float64,2}
    ubb::Array{Float64,1}

    ### RustModel interface
    beta::Float64
    u::Array{Float64,2}                 # u[k,a]: flow utilities                        
    v::Array{Float64,2}                 # v[k,a]: choice-specific discounted value
    V::Array{Float64,1}                 # V[k]: expected value
    pi::Array{Float64,3}                # pi[k,k',a]: condition state transitions 
    p::Array{Float64,2}                 # p[k,a]: conditional choice probabilities (CCP)
end


# dcs2k(dims::Tuple,id,ic,is)=sub2ind(dims,id,ic,is)
# k2dcs(dims::Tuple,ik)=ind2sub(dims,ik)


#id = 1+ demand
#ic = wholesaleprice category
#is = 1 + stock at the beginning of the period
#iy = 1 + sales
#ir = 1 + stock next perido (decision variable)
#iq (implicit) = orders = sales + stock next period - stock today = (iy-1) + (ir - 1) - (is - 1)
#I can "sell" above demand for 0 profit
#I can "sell" above stock by paying a high spotprice (should only happen with extremely high random utility shocks)
#I can destroy goods for free with iq=iy+ir-is <0 
uconstant(model::InventoryModel,id,ic,is,iy,ir)=
    model.pricelevel*(min(id-1,iy-1)-model.wholesaleprices[ic]*max(iy-1+ir-is,0)-model.spotbuy*max(iy-is,0))   
ualpha(model::InventoryModel,id,ic,is,iy,ir)=-model.pricelevel*(is-1)
ubeta(model::InventoryModel,id,ic,is,iy,ir)=-model.pricelevel*1(iy-1+ir-is>0)

function ua(model::InventoryModel,id,ic,is,iy,ir,itheta)
    itheta==1?ualpha(model::InventoryModel,id,ic,is,iy,ir):ubeta(model::InventoryModel,id,ic,is,iy,ir)
end

function inventorymodel(dims,pricelevel,spotbuy,wholesaleprices)
    (dd,dc,ds)=dims
    dtheta=2
    dy,dr=ds,ds
    dk=dd*dc*ds
    da=dy*dr
    model=InventoryModel(
        (dd,dc,ds,dy,dr),
        pricelevel,
        spotbuy,
        wholesaleprices,
        zeros(dk*da,dtheta),  #uaa
        zeros(dk*da),         #ubb
        0.95,
        zeros(dk,da),       #u
        zeros(dk,da),       #v
        zeros(dk),          #V
        zeros(dk,dk,da),    #pi
        zeros(dk,da)        #p
        )
    uaa=[ua(model,id,ic,is,iy,ir,itheta) for id=1:dd,ic=1:dc,is=1:ds,iy=1:dy,ir=1:dr,itheta=1:dtheta]
    model.uaa=reshape(uaa,dk*da,dtheta)
    ubb=[uconstant(model,id,ic,is,iy,ir) for id=1:dd,ic=1:dc,is=1:ds,iy=1:dy,ir=1:dr]
    model.ubb=vec(ubb)
    model
end


function u!(model::InventoryModel,theta)
    model.u[:]=model.uaa*theta+model.ubb
end

function pi!(model::InventoryModel,pidc::AbstractMatrix)
    dd,dc,ds,dy,dr=model.dims
    dk=dd*dc*ds
    da=dy*dr
    for ik=1:dk
        idc=ind2sub((dd*dc,ds),ik)[1]
        for jk=1:dk
            jdc,js=ind2sub((dd*dc,ds),jk)
            for ia=1:da
                iy,ir=ind2sub((dy,dr),ia)
                model.pi[ik,jk,ia]=pidc[idc,jdc]*1(ir==js)
            end
        end
    end
end

# inventorymodel((dd,dc,ds),pricelevel,spotbuy,wholesaleprices)
# im1(pricelevel)=inventorymodel((5,3,10),pricelevel,10.0,[0.9,1.0,1.1])


function pidc1(model::InventoryModel)
    dd,dc,ds,dy,dr=model.dims
    pid=Tridiagonal(fill(0.1,dd-1),fill(0.8,dd),fill(0.1,dd-1))
    @assert dc==3
    pid[1,2]=0.2
    pid[dd,dd-1]=0.2
    pic=[0.1 0.9 0; 0.05 0.9 0.05; 0 0.9 0.1]
    kron(pic,pid)    #reverse order is on purpose the get the right layout
end

function pidc2(model::InventoryModel)
    dd,dc,ds,dy,dr=model.dims
    @assert dd==5
    @assert dc==3
    pid=[
    0.95 0.0 0.0 0.0 0.05;
    0.0 0.0 0.0 0.0 1.0;
    0.0 0.0 0.0 0.0 1.0;
    0.0 0.0 0.0 0.0 1.0;
    0.05 0.0 0.0 0.0 0.95;
    ]
    pic=[0.1 0.9 0; 0.05 0.9 0.05; 0 0.9 0.1]
    kron(pic,pid)    #reverse order is on purpose the get the right layout
end


function detensor(model::InventoryModel)
    dd,dc,ds,dy,dr=model.dims
    pidc=reshape(model.pi[:,:,1],(dd,dc,ds,dd,dc,ds))
    pic=sum(pidc[1,:,1,:,:,:],[2,4])[:,1,:,1]
    pid=sum(pidc[:,1,1,:,:,:],[3,4])[:,:,1,1]
    pid,pic
end

function twostep(model::InventoryModel)
    "todo"
    # V=(I-model.beta*model.pi[:,:,1])\(model.u[:,1]-log(model.p[:,1]))
    # thetau2=(I-model.beta*model.pi[:,:,2])*V+log(model.p[:,2])
    # model.u2\thetau2
end

#################### Identification

#a=[ A_1  I-betaPi_1; ... ; A_da  I-betaPi_da ] 

# function idrank(model::InventoryModel,eta)
#     dk,da=size(model.u)
#     dtheta=size(model.uaa,2)
#     b=[1(ind2sub((dk,da),ika)[1]==jk)-model.beta*model.pi[ind2sub((dk,da),ika)[1],jk,ind2sub((dk,da),ika)[2]] for ika=1:dk*da,jk=1:dk]
#     a=hcat(sparse(model.uaa),sparse(b))
#     aa=full(a'*a)
#     ai=inv(aa)[1:dtheta,:]
#     function ff(eta1)
#         pi!(model,eta1)
#         piv=vcat([model.pi[:,:,ia]*model.V for ia=1:da]...)
#         model.beta*ai*(a'*piv)
#     end
#     jpi=Calculus.finite_difference_jacobian(ff,eta)
#     jpi
# end 


