
function Base.display(M::TensorMesh3D)
	
	println("$(M.dim)-dimensional tensor mesh of size $(M.n)")
	println("Number of cells:   $(M.nc)")
	println("Number of faces:   $(M.nf) = $(sum(M.nf))")
	println("Number of edges:   $(M.ne) = $(sum(M.ne))")
	println("Number of nodes:   $(prod(M.n.+1))")
	println("Coordinate origin: $(M.x0)")
	println("Domain size:       $(sum(M.h1))m x $(sum(M.h2))m x $(sum(M.h3))m")
	println("Minimum cell size: $(minimum(M.h1))m x $(minimum(M.h2))m x $(minimum(M.h3))m")
	println("Maximum cell size: $(maximum(M.h1))m x $(maximum(M.h2))m x $(maximum(M.h3))m")
	
end


function Base.display(M::RegularMesh)
	
	println("$(M.dim)-dimensional regular mesh of size $(M.n)")
	println("Number of cells:   $(M.nc)")
	println("Number of faces:   $(M.nf) = $(sum(M.nf))")
	println("Number of edges:   $(M.ne) = $(sum(M.ne))")
	println("Number of nodes:   $(prod(M.n.+1))")
	println("Coordinate origin: $(M.x0)")
	println("Domain size:       $(M.domain)")
	println("Cell size:         $(M.h)")
	
end
