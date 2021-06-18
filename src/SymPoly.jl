using Symbolics, SymbolicUtils
using SymbolicUtils: istree, operation, arguments
using Symbolics: value, get_variables, solve_for, derivative
using SymbolicUtils.Rewriters
using SymbolicUtils.Code

include("poly.jl")
include("arith.jl")
include("roots.jl")
include("factors.jl")
include("tests.jl")