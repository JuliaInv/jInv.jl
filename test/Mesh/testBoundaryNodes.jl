using jInv.Mesh
using Base.Test

# setup 3D Mesh
nc = [5, 7, 8]
x0 = [0, 0, 0]
domain = [x0[1], 1, x0[2], 1, x0[3], 1]
domain2D = [x0[1], 1, x0[2], 1]
h   = (domain[2:2:end]-domain[1:2:end])./nc
h1  = h[1]*ones(nc[1])
h2  = h[2]*ones(nc[2])
h3  = h[3]*ones(nc[3])

Mt = getTensorMesh3D(h1,h2,h3,x0)
Mr = getRegularMesh(domain,nc)

Mr2D = getRegularMesh(domain2D, [nc[1], nc[2]])

Meshes = (Mt, Mr, Mr2D);


for k=1:length(Meshes)


    M = Meshes[k]
    print("\ttesting boundary nodes for $(typeof(M))...")

    if M.dim==2

        X = getNodalGrid(M)[:,1];
        Y = getNodalGrid(M)[:,2];
        ib, iin = getBoundaryNodes(M);

        # Check if there are 0s in X and Y at boundary nodes
        @test length(X[ib])-length(find( x->(x < eps(Float32)), abs.(X[ib])))>0
        @test length(Y[ib])-length(find( x->(x < eps(Float32)), abs.(Y[ib])))>0

        # Check if there are 1s in X and Y at boundary nodes
        @test length(X[ib])-length(find( x->(x < eps(Float32)), abs.(X[ib]-1)))>0
        @test length(Y[ib])-length(find( x->(x < eps(Float32)), abs.(Y[ib]-1)))>0

        # Check that there aren't 0s in X, Y at inner nodes
        @test length(X[iin])-length(find(X[iin]))==0
        @test length(Y[iin])-length(find(Y[iin]))==0

        # Check that there aren't 1s in X, Y at inner nodes
        @test length(X[iin])-length(find(X[iin]-1))==0
        @test length(Y[iin])-length(find(Y[iin]-1))==0


        elseif M.dim==3

        X = getNodalGrid(M)[:,1];
        Y = getNodalGrid(M)[:,2];
        Z = getNodalGrid(M)[:,3];
        ib, iin = getBoundaryNodes(M);

        # Check if there are 0s in X, Y, and Z at boundary nodes
        @test length(X[ib])-length(find( x->(x < eps(Float32)), abs.(X[ib])))>0
        @test length(Y[ib])-length(find( x->(x < eps(Float32)), abs.(Y[ib])))>0
        @test length(Y[ib])-length(find( x->(x < eps(Float32)), abs.(Z[ib])))>0

        # Check if there are 1s in X, Y, and Z at boundary nodes
        @test length(X[ib])-length(find( x->(x < eps(Float32)), abs.(X[ib]-1)))>0
        @test length(Y[ib])-length(find( x->(x < eps(Float32)), abs.(Y[ib]-1)))>0
        @test length(Y[ib])-length(find( x->(x < eps(Float32)), abs.(Z[ib]-1)))>0

        # Check that there aren't 0s in X,Y,Z at inner nodes
        @test length(X[iin])-length(find(X[iin]))==0
        @test length(Y[iin])-length(find(Y[iin]))==0
        @test length(Z[iin])-length(find(Z[iin]))==0

        # Check that there aren't 1s in X,Y,Z at inner nodes
        @test length(X[iin])-length(find(X[iin]-1))==0
        @test length(Y[iin])-length(find(Y[iin]-1))==0
        @test length(Z[iin])-length(find(Z[iin]-1))==0

    end


    println("passed!\n")
end
