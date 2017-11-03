
module Tree

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

function Base.show(io::IO, node::Root)
    print(io, "*" * node.name * "(" * join([c.name for c in node.children], ", ") * ")")
end

function Base.show(io::IO, node::Child)
    print(io, node.name * "(" * join([c.name for c in node.children], ", ") * ")")
end

function Base.delete!(node::Child)
    "Deletes a node returning it, while connecting the parent and the children of it."
    
    length = node.length
    parent = node.parent
    
    index = findfirst(parent.children, node)
    deleteat!(parent.children, index)

    for child in node.children
        child.length += length
        child.parent = parent
        
        # Make sure the newly added children are in the same position in the
        # children array as the old child was.
        insert!(parent.children, index, child)
        index += 1
    end

    return node
end

function insert(name, child, length)
    "Inserts a node on a branch at 'length' distance ancestral from 'child'."
    
    if length > child.length
        error("length is longer than the branch leading to its child")
    end
    
    newchild = Child(name, child.parent, child.length - length)
    
    children = child.parent.children
    deleteat!(children, findfirst(children, child))
    
    child.parent = newchild
    child.length = length
    
    return newchild
end

function detach!(node::Child; asroot=false)
    """Deletes a subtree from a node, returning it.
    If asroot, convert the given node to a root."""
    
    deleteat!(parent.children, findfirst(parent.children, node))
    
    if asroot
        root = Root(node.name)
        root.children = copy(node.children)
        
        for child in root.children
            child.parent = root
        end
        
        return root
    else
        return node
    end
end

function namemapof(node::Node)
    "Returns a name => node dict of all nodes in subtree"
    
    namemap = Dict{String, Node}(node.name=>node)

    for node in descendantsof(node)
        if node.name in keys(namemap)
            error("name of nodes are not unique")
        end
        
        namemap[node.name] = node
    end
    
    return namemap
end

function isleaf(node::T) where T <: Node
    "Returns true if node has 0 children, else false"
    
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

function childlineageof(node::Child)
    "Returns an array of ancestors from the node to the oldest child."
    
    ancestors = Child[]

    while isa(node, Child)
        push!(ancestors, node)
        node = node.parent
    end
    
    return ancestors
end

function descendantsof(node::T) where T<:Node
    "Returns an array of all children descending from the node."
    
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
    "Returns true if nore or any descendant has one child, else false"
    
    return length(node.children) == 1 || any(length(i.children) == 1 for i in descendantsof(node))
end

function ispolytomic(node::T) where T<:Node
    "Returns true if nore or any descendant has more than two children, else false"
    return length(node.children) > 2 || any(length(i.children) > 2 for i in descendantsof(node))
end

function mrcaof(nodes::Array{T, 1}) where T<:Node
    """Finds the most recent common ancestor of two nodes.
    If the two nodes are not related, returns nothing"""
    
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
    "Finds the distance between two nodes"
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
    "Finds the distance between two nodes"
    mrca = mrcaof([one, two])
    if mrca == nothing
        error("the two nodes do not share an ancestor")
    end
    
    return distance(one, two, mrca)
end

function simplify!(node::T) where T<:Node
    "Deletes all redundant nodes in subtree"
    
    while length(node.children) == 1
        delete!(node.children[1])
    end
    
    for child in descendantsof(node)
        if length(child.children) == 1
            delete!(child)
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

function _recursor(array, node, length)
    newlength = node.length + length
    push!(array, (node, newlength))
    for child in node.children
        _recursor(array, child, newlength)
    end
end

function _ndchild(node::Child, ancestors::Array{Child, 1}, withroot)
    pairs = Tuple{Node, Float64}[]
    
    dist = 0.0
    previousancestor = node
    for i in eachindex(ancestors)
        ancestor = ancestors[i]
        push!(pairs, (ancestor, dist))
        
        if i == 1
            for child in ancestor.children
                _recursor(pairs, child, dist)
            end
        else
            for child in filter(c->c != previousancestor, ancestor.children)
                _recursor(pairs, child, dist)
            end
        end
        
        dist += ancestor.length
        previousancestor = ancestor
    end
    
    # Add the root is asked to
    if withroot  
        root::Root = ancestors[end].parent
        push!(pairs, (root, dist))

        for child in filter(c->c != ancestors[end], root.children)
            _recursor(pairs, child, dist)
        end
    end
    
    return pairs
end

nodedistances(node::Child) = _ndchild(node, lineageof(node, true), true)

function nodedistances(node::Child, backto::Child)
    
    ancestors = lineageof(node, true)
    backtoindex = findfirst(ancestors, backto)
    if backtoindex == 0
        error("backto must be an ancestor of node")
    end
    return _ndchild(node, ancestors[1:backtoindex], false)
end

function nodedistances(node::Child, backto::Root)
    ancestors = lineageof(node, true)
    if backto != ancestors[end].parent
        error("backto must be an ancestor of node")
    end
    return _ndchild(node, ancestors, true)
end

function nodedistances(node::Root)
    pairs = Tuple{Node, Float64}[(node, 0.0)]
    
    for child in node.children
        _recursor(pairs, child, 0)
    end
    
    return pairs
end

function nodedistances(node::Root, backto::Root)
    if backto != node
        error("backto must be an ancestor of node")
    end
    return nodedistances(node)
end

function _furthestleaf(node, backto)
    leafdists::Array{Tuple{Child, Float64}} = filter(x->isleaf(x[1]), nodedistances(node, backto))
    maxdist = -1.0
    result = leafdists[1][1]
    for (leaf, dist) in leafdists
        if dist > maxdist
            result = leaf
        end
    end
    return result
end

function diameter(node::T) where T <: Node
    leafone = _furthestleaf(node, node)
    leaftwo = _furthestleaf(leafone, node)
    
    return distance(leafone, leaftwo)
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

function reroot!(node::Child)
    # Rerooting is done by reversing the lineage.
    lineage = lineageof(node)
    
    # Since we reverse, lengths are now attached to the other end of the branch
    lengths = [n.length for n in lineage[1:end-1]]
    insert!(lengths, 1, pop!(lengths))
    for i in eachindex(lengths)
        lineage[i].length = lengths[i]
    end
    
    # Instead of making new child and root, we swap properties of last and first node
    child, root = lineage[1], lineage[end]
    child.name, root.name = root.name, child.name 
    child.children, root.children = root.children, child.children
    
    # When iterating, compensate for the fact that the two above "traded places"
    lineage[1], lineage[end] = lineage[end], lineage[1]
    
    # Update parents and children for all
    for i in eachindex(lineage)
        
        # Update parent for all but the new root and delete its previous child
        if i == 2
            # This is special since its previous child is the new root.
            lineage[i].parent = lineage[1]
            deleteat!(lineage[i].children, findfirst(lineage[i].children, lineage[end]))
        elseif i > 2
            lineage[i].parent = lineage[i-1]
            deleteat!(lineage[i].children, findfirst(lineage[i].children, lineage[i-1]))
        end
        
        # For all except old root, add the old parent as child
        if i < length(lineage)
            push!(lineage[i].children, lineage[i+1])
        end
        
    end
    
    return root
end

include("parsetree.jl")

export Node, Child, Root
export newick, json
export delete!, insert, detach!, namemapof, lineageof, childlineageof
export isleaf, descendantsof, isredundant, ispolytomic, distance, nodedistances
export simplify!, subtree, diameter, shufflenames!, mrcaof, reroot!

end # module Tree


