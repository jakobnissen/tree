
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

# Test Newick parsing/lexing
newickstring = "(Bovine:0.69395,(Gibbon:0.36079,(:0.33636,(Gorilla:0.17147,(Chimp, Human:0.11927):0.08386):0.06124):0.15057):0.54939, Mouse:1.21460);"

parsedstring = "(Bovine:0.69395,(Gibbon:0.36079,(3:0.33636,(Gorilla:0.17147,(Chimp:0.0,Human:0.11927)5:0.08386)4:0.06124)2:0.15057)1:0.54939,Mouse:1.2146)root;"

luca = newick(newickstring)

@test luca.name == "root"
@test [child.name for child in luca.children] == ["Bovine", "1", "Mouse"]
@test [child.length for child in luca.children] == [0.69395, 0.54939, 1.2146]
@test newick(luca) == parsedstring

# Test JSON parsing/lexing
jsonstring = """{"name":"root","children":[{"length":0.69395,"name":"Bovine","parent":"root","children":[]},{"length":0.54939,"name":"1","parent":"root","children":[{"length":0.36079,"name":"Gibbon","parent":"1","children":[]},{"length":0.15057,"name":"2","parent":"1","children":[{"length":0.33636,"name":"3","parent":"2","children":[]},{"length":0.06124,"name":"4","parent":"2","children":[{"length":0.17147,"name":"Gorilla","parent":"4","children":[]},{"length":0.08386,"name":"5","parent":"4","children":[{"length":0.0,"name":"Chimp","parent":"5","children":[]},{"length":0.11927,"name":"Human","parent":"5","children":[]}]}]}]}]},{"length":1.2146,"name":"Mouse","parent":"root","children":[]}]}"""

luca = json(jsonstring)

@test luca.name == "root"
@test [child.name for child in luca.children] == ["Bovine", "1", "Mouse"]
@test [child.length for child in luca.children] == [0.69395, 0.54939, 1.2146]
@test json(luca) == jsonstring

# Test isleaf
@test isleaf(root) == false
@test isleaf(children["a"]) == false
@test isleaf(children["e"]) == true

# Test isredundant
nonredundant = Root("a")
x = Child("1", nonredundant, 1)
x = Child("1", nonredundant, 2)
Child("1", x, 4)
Child("1", x, 8)

@test isredundant(nonredundant) == false
@test isredundant(x) == false
@test isredundant(x.children[1]) == false

@test isredundant(root) == true
@test isredundant(children["d"]) == true
@test isredundant(children["e"]) == false

# Test ispolytomic
@test ispolytomic(luca) == true
@test ispolytomic(luca.children[1]) == false

@test ispolytomic(root) == false
@test ispolytomic(children["e"]) == false

# Test lineageof
@test lineageof(root) == Node[root]
@test lineageof(children["h"]) == Node[children["h"], children["f"], children["d"], children["b"], root]

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
