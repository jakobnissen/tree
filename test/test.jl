
push!(LOAD_PATH, "../src")

using Tree
using Base.Test

# Create a manual testing tree
root = Root("luca")
children = Dict{String, Child}()
children["a"] = Child("a", root, 1<<0)
children["b"] = Child("b", root, 1<<1)
children["c"] = Child("c", children["a"], 1<<2)
children["d"] = Child("d", children["b"], 1<<3)
children["e"] = Child("e", children["c"], 1<<4)
children["f"] = Child("f", children["d"], 1<<5)
children["g"] = Child("g", children["d"], 1<<6)
children["h"] = Child("h", children["f"], 1<<7)

@testset "Instantiation" begin
    # Test Root
    @test root.children == Child[children["a"], children["b"]]
    @test root.name == "luca"

    # Test Child
    @test children["e"].children == Child[]
    @test children["a"].children == Child[children["c"]]
    @test children["e"].name == "e"
    @test children["a"].parent == root
    @test children["c"].parent == children["a"]
    @test children["a"].length == 1
end

@testset "Parsing/Lexing Newick" begin
    newickstring = "(Bovine:0.69395,(Gibbon:0.36079,(:0.33636,(Gorilla:0.17147,(Chimp, Human:0.11927):0.08386):0.06124):0.15057):0.54939, Mouse:1.21460);"
    parsedstring = "(Bovine:0.69395,(Gibbon:0.36079,(3:0.33636,(Gorilla:0.17147,(Chimp:0.0,Human:0.11927)5:0.08386)4:0.06124)2:0.15057)1:0.54939,Mouse:1.2146)root;"
    luca = newick(newickstring)

    @test luca.name == "root"
    @test [child.name for child in luca.children] == ["Bovine", "1", "Mouse"]
    @test [child.length for child in luca.children] == [0.69395, 0.54939, 1.2146]
    @test length(descendantsof(luca.children[2])) == 8
    @test newick(luca) == parsedstring
end

@testset "Parsing/Lexing JSON" begin
    jsonstring = """{"name":"root","children":[{"length":0.69395,"name":"Bovine","parent":"root","children":[]},{"length":0.54939,"name":"1","parent":"root","children":[{"length":0.36079,"name":"Gibbon","parent":"1","children":[]},{"length":0.15057,"name":"2","parent":"1","children":[{"length":0.33636,"name":"3","parent":"2","children":[]},{"length":0.06124,"name":"4","parent":"2","children":[{"length":0.17147,"name":"Gorilla","parent":"4","children":[]},{"length":0.08386,"name":"5","parent":"4","children":[{"length":0.0,"name":"Chimp","parent":"5","children":[]},{"length":0.11927,"name":"Human","parent":"5","children":[]}]}]}]}]},{"length":1.2146,"name":"Mouse","parent":"root","children":[]}]}"""
    luca = json(jsonstring)

    @test luca.name == "root"
    @test [child.name for child in luca.children] == ["Bovine", "1", "Mouse"]
    @test [child.length for child in luca.children] == [0.69395, 0.54939, 1.2146]
    @test length(descendantsof(luca.children[2])) == 8
    @test json(luca) == jsonstring
end

@testset "namemapof" begin
    # On root
    namemap = namemapof(root)

    @test length(namemap) == 9
    @test namemap["a"] == root.children[1]

    vals = sort([i for i in values(namemap)], by=x -> x.name)[1:8]
    @test vals == sort([i for i in values(children)], by=x -> x.name)
    
    # On child
    namemap = namemapof(children["b"])
    @test length(namemap) == 5
    @test namemap["d"] == children["d"]
    
    # Temporarily add a new child with an existing name
    Child("d", children["c"], 9.2)
    
    namemap = namemapof(root)
    @test length(namemap) == 9
    
    # With unique
    @test_throws ErrorException namemap = namemapof(root, unique=true)
end

@testset "delete!" begin
    # Begin by deleting the extra child we added before
    @test length(descendantsof(root)) == 9
    
    delete!(children["c"].children[2])
    
    @test length(descendantsof(root)) == 8
    
    # Test that children are properly added and length preserved
    new_f_length = children["d"].length + children["f"].length
    new_g_length = children["d"].length + children["g"].length
    
    delete!(children["d"])
    
    @test children["b"].children == [children["f"], children["g"]]
    @test children["f"].length == new_f_length
    @test children["g"].length == new_g_length
end

@testset "insert" begin
    # Insert the "d" we deleted before
    newchild = insert("d", children["g"], children["g"].length - children["d"].length)
    
    @test children["g"].length == 1 << 6
    @test newchild.length == children["d"].length
    @test newchild.children == [children["g"]]
    children["d"] = newchild
end

@testset "move!" begin
    num_children = length(descendantsof(root))
    move!(children["f"], children["d"])
    children["f"].length -= 1 << 3
    
    @test length(descendantsof(root)) == num_children
    @test children["d"].children == [children["g"], children["f"]]
    @test descendantsof(children["d"]) == [children["g"], children["f"], children["h"]]
end

@testset "isleaf" begin
    @test isleaf(root) == false
    @test isleaf(children["a"]) == false
    @test isleaf(children["e"]) == true
end

@testset "isredundant" begin
    nonredundant = Root("1")
    child1 = Child("1", nonredundant, 1)
    child2 = Child("1", nonredundant, 2)
    child3 = Child("1", child2, 4)
    child4 = Child("1", child2, 8)
    
    @test isredundant(nonredundant) == false
    
    for child in descendantsof(nonredundant)
        @test isredundant(child) == false
    end
    
    @test isredundant(root) == true
    @test isredundant(children["d"]) == true
    @test isredundant(children["e"]) == false
end

@testset "ispolytomic" begin
    polytomic = Root("1")
    child1 = Child("1", polytomic, 1)
    child2 = Child("1", polytomic, 2)
    child3 = Child("1", polytomic, 4)
    child4 = Child("1", child1, 8)
    child5 = Child("1", child1, 16)
    child6 = Child("1", child2, 32)
    
    @test ispolytomic(polytomic) == true
    
    
    for child in descendantsof(polytomic)
        @test ispolytomic(child) == false
    end
    
    @test ispolytomic(root) == false
end

@testset "lineageof/childlineageof" begin
    @test lineageof(root) == Node[root]
    @test lineageof(children["h"]) == Node[children["h"], children["f"], children["d"], children["b"], root]
    
    @test_throws MethodError childlineageof(root)
    @test childlineageof(children["h"]) == Child[children["h"], children["f"], children["d"], children["b"]]
end

# Test descendantsof
@test descendantsof(root) == sort(collect(values(children)), by=x -> x.name)
@test descendantsof(children["e"]) == Child[]

# Test mrcaof
@test mrcaof([root]) == root
@test mrcaof([children["h"], children["h"]]) == children["h"]
@test mrcaof([children["a"], children["e"]]) == children["a"]
@test mrcaof([children["h"], children["g"]]) == children["d"]

# Test distance
@test distance(root, root) == 0.0
@test distance(children["b"], children["b"]) == 0.0
@test distance(root, children["b"]) == 2.0
@test distance(children["a"], children["b"]) == 3.0
@test distance(children["g"], children["h"]) == 224.0

@test distance(root, root, root) == 0.0
@test distance(children["b"], children["b"],  children["b"]) == 0.0
@test distance(root, children["b"], root) == 2.0
@test distance(children["a"], children["b"], root) == 3.0
@test distance(children["g"], children["h"], children["d"]) == 224.0

# Test simplify
root2 = deepcopy(root)
simplify!(root2)
@test [child.name for child in root2.children] == ["e", "d"]
@test length(descendantsof(root2)) == 4
@test length(descendantsof(root2.children[2])) == 2
@test distance(root2.children[1], root2.children[2].children[1]) == 191.0

# Test subtree
# Test furthestleaf - lots of methods!
# Test diameter
# Test shufflenames!
