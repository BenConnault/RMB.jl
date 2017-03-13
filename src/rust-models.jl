abstract RustModel
    # beta::Float64
    # u::Array{Float64,2}                 # u[k,a]: flow utilities                        
    # v::Array{Float64,2}                 # v[k,a]: choice-specific discounted value
    # V::Array{Float64,1}                 # V[k]: expected value
    # pi::Array{Float64,3}                # pi[k,k',a]: condition state transitions 
    # p::Array{Float64,2}                 # p[k,a]: conditional choice probabilities (CCP)



function solvedp!(model::RustModel)
    da=size(model.u)[2]
    
    function ff!(z,fv)
        fv[:]=-1.0
        for ia=1:da
            # fv[:]+=exp.(min(model.u[:,ia]+model.beta*model.pi[:,:,ia]*z-z,1))   
            fv[:]+=exp.(model.u[:,ia]+model.beta*model.pi[:,:,ia]*z-z)   
        end
    end

    function fj!(z,fm)
        fm[:]=0.0                           #inefficient reallocation of temp at each call
        for ia=1:da
            exponent=model.u[:,ia]+model.beta*model.pi[:,:,ia]*z-z
            # fm[:,:]+=diagm(1(exponent.<1).*exp.(exponent))*(model.beta*model.pi[:,:,ia]-eye(size(fm,1)))
            fm[:,:]+=diagm(exp.(exponent))*(model.beta*model.pi[:,:,ia]-eye(size(fm,1)))
        end
    end

    df = NLsolve.DifferentiableMultivariateFunction(ff!,fj!)

    res=NLsolve.nlsolve(df,model.V,autoscale=false) #will automatically start at previous value
    
    # println(res)
    model.V[:]=res.zero         
    for ia=1:da
        model.p[:,ia]=exp(model.u[:,ia]+model.beta*model.pi[:,:,ia]*model.V-model.V) 
    end
    println(checkdp(model), " DP solved?")
    # println(round(model.p,2))
end

function checkdp(model::RustModel)
    V1=\(eye(model.pi[:,:,1])-model.beta*model.pi[:,:,1],model.u[:,1]-log(model.p[:,1]))
    V2=\(eye(model.pi[:,:,2])-model.beta*model.pi[:,:,2],model.u[:,2]-log(model.p[:,2]))
    norm(V1-V2)
end

########## Identification

# "nonparametric" two-step projection is an implicit system of d_k(d_a+1) equations Phi((u,V),Pi)=A(Pi) (u,V) - B(P)
# to compute the Jacobian of eta ---> u(Pi(eta)) we use the implicit function theorem on the full X=(u,V)
# the implicit function theorem require inverting jx:=dPhi/dX. Fortunately this is closed form.
# A(Pi)=[A 0; I B] ---> A(Pi)^{-1} = [BC I-BCA;-C CA], C=(AB)^{-1}
# After three steps:
# (1) selecting the u block that we are interested in
# (2) using the fact the identifiying equation does not depend on Pi
# (3) developing dA(Pi)dPi using kronecker product formulas
# we find: du/dPi = beta (I - B(AB)^{-1}A) (I_{d_a} \otimes V' \otimes I_{d_k}) J_Pi
# where J_Pi = d (vec(Pi))/d alpha is (da dk dk x deta) is the Jacobian of Pi(eta) can be computed numerically
# the computational bottleneck is inverting (AB) which is (dk x dk)


function explicitinv(a,b)
    n=size(b,1)
    c=inv(a*b)
    [b*c eye(n)-b*c*a; -c c*a]
end



function rank(model::RustModel,a,eta)
    dk,da=size(model.u)
    b=[1(in2sub((dk,da),ika)[1]==jk)-model.beta*model.pi[in2sub((dk,da),ika)[1],jk,in2sub((dk,da),ika)[2]] for ika=1:dk*da,jk=1:dk]
    m=I-b*((a*b)\a)
    m*kron(eye(da),model.V',eye(dk))*jpi
end 





########## Simulation

function rand(model::RustModel,T)
    dk,da=size(model.p)
    ka=zeros(Int,T,2)
    ka[1,1]=1
    ka[1,2]=wsample(1:da,model.p[ka[1,1],:])
    for t=1:T-1
        ka[t+1,1]=wsample(1:dk,model.pi[ka[t,1],:,ka[t,2]])
        ka[t+1,2]=wsample(1:da,model.p[ka[t+1,1],:])
    end
    ka
end
