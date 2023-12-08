## Installation and Running the Oscar Implementation

Running the code and using multiple cores requires a dev version of Oscar.
If you use the project.toml (instructions below) this should get you the right 
version. This is subject to change soonish once the changes get approved and
merged into master.

```julia
julia> ] 
pkg> activate .
pkg> instantiate

julia> include("MultigradedImplicitization.jl")
```

Running script will run a parallel computation using n_cores (currently set to 8)
once computation is done, the script saves to file named kernel.json.
