module RMB

importall NLsolve
import Base.rand
import Distributions: wsample
import Calculus


export solvedp!, u!, pi!
export InventoryModel

include("stochastic-matrices.jl")
include("rust-models.jl")
include("inv_model.jl")



end #module