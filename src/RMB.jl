module RMB

importall NLsolve
import Base.rand
import Distributions: wsample


export solvedp!, u!, pi!

include("stochastic-matrices.jl")
include("rust-models.jl")
include("inv_model.jl")



end #module