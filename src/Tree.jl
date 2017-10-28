
 module Tree

import JSON

abstract type Node end

mutable struct Child <: Node
    name::String
    parent::Node
    length::Float64
    children::Array{Child,1}

    function Child(name, parent, length)
        self = new(name, parent, length, Child[])
        push!(parent.children, self)
        return self
    end
end

mutable struct Root <: Node
    name::String    
    children::Array{Child,1}

    function Root(name)
        new(name, Child[])
    end
end

function _parsechildren(jsonparent, parent)
    for jsonchild in jsonparent["children"]
        if isa(jsonchild, String)
            jsonchild = JSON.parse(jsonchild)
        end
        child = Child(jsonchild["name"], parent, jsonchild["length"])
        
        _parsechildren(jsonchild, child)
    end
end

function json(string::T) where T <: AbstractString
    jsonroot::Dict{String, Any} = JSON.parse(string)
    root = Root(jsonroot["name"])
    parsechildren(jsonroot, root)
    return root
end

function _newickrecurser(nodestring, parent, nodearray)
    "Recursively instantiates children of a Newick string"
    
    name, branchlength, stringofchildren = _getnodeinfo(nodestring)
    
    if isempty(name)
        number = nodearray[1]
        nodearray[1] = number + 1
        name = string(number)
    end
    
    child = Child(name, parent, branchlength)
    
    for childstring in _newicksubstrings(stringofchildren)
        _newickrecurser(childstring, child, nodearray)
    end
end

function _newicksubstrings(string)
    "Given a newick tree representing zero or more children, properly splits them."
    
    if isempty(string)
        return SubString{String}[]
    end
    
    pars = 0
    left = 1
    substrings = SubString{String}[]
    
    for i in eachindex(string)
        if string[i] == '('
            pars += 1
        elseif string[i] == ')'
            pars -= 1
        elseif pars == 0 && string[i] == ','
            push!(substrings, strip(string[left: i - 1]))
            left = i + 1
        end
    end
    
    push!(substrings, strip(string[left: length(string)]))
    return substrings
end

function _getnodeinfo(nodestring)
    """Gets branchlength and name and the rest of the string for the top node
    of a Newick string"""
    
    rightparenpos = rsearch(nodestring, ')')
    
    colonpos = search(nodestring, ':', rightparenpos+1)
    
    if colonpos == 0
        name = String(nodestring[rightparenpos+1:end])
        branchlength = 0.0
    else
        name = String(nodestring[rightparenpos+1:colonpos-1])
        branchlength = parse(Float64, nodestring[colonpos+1:end])
    end
    
    return name, branchlength, nodestring[2:rightparenpos-1]
end

function newick(newickstring::T) where T<: AbstractString
    "Returns the root of a tree represented by a Newick string"
    
    nodearray = [1]
    
    newickstring = rstrip(newickstring, [';', '\n'])
    substring = SubString(newickstring, 1, length(newickstring))
    
    name, _, childrenstrings = _getnodeinfo(substring)
    
    if isempty(name)
        name = "root"
    end
    
    result = Root(name)
    
    for childstring in _newicksubstrings(childrenstrings)
        _newickrecurser(childstring, result, nodearray)
    end
    
    return result
end

function Base.show(io::IO, node::Root)
    print(io, "*" * node.name * "(" * join([c.name for c in node.children], ", ") * ")")
end

function Base.show(io::IO, node::Child)
    print(io, node.name * "(" * join([c.name for c in node.children], ", ") * ")")
end

function _jsonrecurse(node::Child)
    object = Dict{String, Any}([])
    object["name"] = node.name
    object["parent"] = node.parent.name
    object["length"] = node.length
    object["children"] = Any[]
    
    for child in node.children
        push!(object["children"], _jsonrecurse(child))
    end
    
    return object
end

function json(node::Child)
    object = Dict{String, Any}([])
    object["name"] = node.name
    object["parent"] = node.parent.name
    object["length"] = node.length
    object["children"] = Any[]
    
    for child in node.children
        push!(object["children"], _jsonrecurse(child))
    end
    
    return JSON.json(object)
end

function json(node::Root)
    object = Dict{String, Any}([])
    object["name"] = node.name
    object["children"] = Any[]
    
    for child in node.children
        push!(object["children"], _jsonrecurse(child))
    end
    
    return JSON.json(object)
end

function newick(node::Child, asroot=true)
    if isleaf(node)
        return node.name * ":" * string(node.length)
    else
        childstrings = [newick(c, false) for c in node.children]
        nodestring = "($(join(childstrings, ",")))$(node.name):$(node.length)"
        if asroot
            return nodestring * ";"
        else
            return nodestring
        end
    end
end

function newick(node::Root)
    childstrings = [newick(c, false) for c in node.children]
    return "($(join(childstrings, ",")))$(node.name);"
end

function isleaf(node::T) where T <: Node
    return isempty(node.children)
end

function lineageof(node::T) where T<:Node
    "Returns an array of ancestors from the node to the root."
    
    ancestors = Node[node]
    
    while isa(node, Child)
        node = node.parent
        push!(ancestors, node)
    end
    
    return ancestors
end

function descendantsof(node::T) where T<:Node
    "Returns an array of all children desceding from the node."
    
    result = copy(node.children)
    
    for child in result
        for grandchild in child.children
            push!(result, grandchild)
        end
    end
    return result
end

# Shortcutting saves 90% memory (15 KiB for 1000 node tree at LUCA)
# but basically no time for these two functions.

function isredundant(node::T) where T<:Node
    return length(node.children) == 1 || any(length(i.children) == 1 for i in descendantsof(node))
end

function ispolytomic(node::T) where T<:Node
    return length(node.children) > 2 || any(length(i.children) > 2 for i in descendantsof(node))
end

function mrcaof(nodes::Array{T, 1}) where T<:Node
    "Another implementation"
    
    lineages = [lineageof(node) for node in nodes]
    
    minlength = minimum([length(lineage) for lineage in lineages])
    generations = Array{Node}(length(nodes), minlength)
    
    for (i, lineage) in enumerate(lineages)
        generations[i,:] = lineage[end - minlength + 1 : end]
    end

    for generation in 1:minlength
        base = generations[1, generation]
        if all(i == base for i in generations[:, generation])
            return base
        end
    end
    return nothing
end

function distance(one::T1, two::T2, mrca::T3) where {T1<:Node, T2<:Node, T3<:Node}
    "Finds the vertical distance between two nodes"
    dist = 0.0
    
    for leaf in (one, two)
        head = leaf
        
        while head != mrca
            dist += head.length
            head = head.parent
        end
    end
    
    return dist
end

function distance(one::T1, two::T2) where {T1<:Node, T2<:Node}
    mrca = mrcaof([one, two])
    if mrca == nothing
        error("the two nodes do not share an ancestor")
    end
    
    return distance(one, two, mrca)
end

_sarrs(node::Child) = [node], [findfirst(node.parent.children, node)]
_sarrs(node::Root) = copy(node.children), [i for i in 1:length(node.children)]

function simplify!(node::T) where T<:Node
    "Removes as many nodes as possible without changing topology of tree."
    
    nodesleft, indexes = _sarrs(node)
    
    # Iterate over all nodes.
    while length(nodesleft) > 0
        child = pop!(nodesleft)
        childindex = pop!(indexes)
    
        if length(child.children) == 1
            parent = child.parent
            distance = child.length

            # If a string of single children are observed, go directly to last child
            while length(child.children) == 1
                child = child.children[1]
                distance += child.length
            end

            # Now connect the last child to the parent directly.
            child.parent = parent
            child.length = distance
            parent.children[childindex] = child
        end
        
        # Now there are guaranteed not 1 child, so add these.
        for (index, grandchild) in enumerate(child.children)
            push!(nodesleft, grandchild)
            push!(indexes, index)
        end
    end
end

function subtree(nodes::Array{T,1}) where T<:Node
    """Creates a new tree with the same topology, lengths and names, including
    only the nodes given and the ancestors necessary to keep topology same."""
    
    if isempty(nodes)
        error("no nodes given")
    end
    
    
    copyof = Dict{Node, Node}()
    lineage = reverse(lineageof(nodes[1]))
    
    # Create the root which is returned
    luca = Root(lineage[1].name)
    copyof[lineage[1]] = luca
    
    # For the first node, create all the rest of the lineage
    for template in lineage[2:end]
        copy = Child(template.name, copyof[template.parent], template.length)
        copyof[template] = copy
    end

    # For each subsequent node
    for node in nodes[2:end]
        lineage = Child[]
        
        # While the ancestors have never been seen, keep getting more distanct
        # ancestors
        while ~haskey(copyof, node)
            push!(lineage, node)
            node = node.parent
            
            # If hit a root, if it is the same root, stop, if it is another root,
            # raise an error as the nodes do not share a common ancestor
            if isa(node, Root)
                if ~haskey(copyof, node)
                    error("Not all nodes share an ancestor")
                else
                    break
                end
            end
        end
        
        # Now instantiate all these unseen ancestors
        for template in reverse(lineage)
            copy = Child(template.name, copyof[template.parent], template.length)
            copyof[template] = copy
        end
    end

    # Lastly, remove all the one-children nodes.
    simplify!(luca)
    return luca
end

function furthestleaf(node::Child, backto::T, among::Array{Child,1}) where T<: Node
    totarget = Dict{Node, Float64}()
    current::Child = node
    distance = 0.0
    
    if current != backto
        while current.parent != backto
            totarget[current] = distance
            distance += current.length
            current = current.parent
        end
    end
    
    totarget[current] = distance
    totarget[current.parent] = distance + current.length
    
    maxdistance, maxleaf = 0.0, node
    
    for otherleaf in among
        if otherleaf == node
            continue
        end
        
        current = otherleaf
        distance = 0.0
        encountereds = [current]
        traveleds = [0.0]
        
        while ~haskey(totarget, current.parent)
            distance += current.length
            current = current.parent
            push!(encountereds, current)
            push!(traveleds, distance)
        end
        
        distance += current.length
        for (encountered, traveled) in zip(encountereds, traveleds)
            totarget[encountered] = totarget[current.parent] + (distance - traveled)
        end

        distance += totarget[current.parent]
        if distance > maxdistance
            maxdistance = distance
            maxleaf = otherleaf
        end
    end
    
    return maxleaf, maxdistance
end

function furthestleaf(node::Root, backto::Root, among::Array{Child,1})
    if node != backto
        error("node is root, backto must be same root")
    end
    
    maxleaf = among[1]
    maxdistance = 0.0
    
    for leaf in among
        current::Child = leaf
        dist = current.length
        
        while isa(current.parent, Child)
            current = current.parent
            dist += current.length
        end
        
        if dist > maxdistance
            maxdistance = dist
            maxleaf = leaf
        end
    
    end
    
    return maxleaf, maxdistance
    
end

function furthestleaf(node::Child, backto::T) where T<: Node
    leaves = filter(isleaf, descendantsof(backto))
    return furthestleaf(node, backto, leaves)
end

function furthestleaf(node::Root, backto::Root)
    leaves = filter(isleaf, descendantsof(backto))
    return furthestleaf(node, backto, leaves)
end

function furthestleaf(node::Child)
    leaves = filter(isleaf, descendantsof(node))
    backto = lineageof(node)[end]
    return return furthestleaf(node, backto, leaves)
end

function furthestleaf(node::Root)
    leaves = filter(isleaf, descendantsof(node))
    return furthestleaf(node, node, leaves)
end

function furthestleaf(node::Child, among::Array{Child,1})
    leaves = copy(among)
    push!(node, leaves)
    backto = mrcaof(leaves)
    return return furthestleaf(node, backto, among)
end

function furthestleaf(node::Root, among::Array{Child,1})
    return furthestleaf(node, node, among)
end

function diameter(node::T) where T<:Node
    leaves = filter(isleaf, descendantsof(node))
    
    if length(leaves) < 2
        return 0.0
    end
    
    startleaf = leaves[1]
    
    leafone, _ = furthestleaf(startleaf, node, leaves)
    leaftwo, diameter = furthestleaf(leafone, node, leaves)
    
    return diameter, leafone, leaftwo
end

function shufflenames!(node::T) where T <: Node
    descendants = descendantsof(node)
    names = [child.name for child in descendants]
    shuffle!(names)

    for i in eachindex(descendants)
        descendants[i].name = names[i]
    end
    
    return nothing
end

export Node, Child, Root 
export newick, json, isleaf, lineageof, descendantsof, isredundant, ispolytomic
export simplify!, subtree, furthestleaf, diameter, shufflenames!, mrcaof, distance

 end # module Tree
