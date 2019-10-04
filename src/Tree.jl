module Tree

#=
TODO: A function that calculates if two trees are identical

=#

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

Child(name::String, parent::Node) = Child(name, parent, 0.0)

mutable struct Root <: Node
    name::String
    children::Array{Child,1}
    function Root(name)
        new(name, Child[])
    end
end

function Base.getproperty(x::Child, sym::Symbol)
    if sym === :parent
        return getfield(x, :parent)::Union{Child, Root}
    else
        return getfield(x, sym)
    end
end

function Base.show(io::IO, node::Root)
    print(io, "*$(node.name)($(join([c.name for c in node.children], ", ")))")
end

function Base.show(io::IO, node::Child)
    print(io, "$(node.name)($(join([c.name for c in node.children], ", ")))")
end

function Base.deepcopy_internal(x::Child, dict::IdDict)
    if haskey(dict, x)
        return dict[x]
    elseif isempty(dict)
        y = Child(x.name, Root("x"), x.length)
        y.parent = x.parent
    else
        y = Child(x.name, dict[x.parent], x.length)
    end

    dict[x] = y

    for child in x.children
        Base.deepcopy_internal(child, dict)
    end

    return y
end

function Base.deepcopy_internal(x::Root, dict::IdDict)
    y = Root(x.name)
    dict[x] = y
    for child in x.children
        Base.deepcopy_internal(child, dict)
    end
    return y
end

function deepcopyasroot(x::Child)
    dict = IdDict()
    y = Root(x.name)
    dict[x] = y
    for child in x.children
        Base.deepcopy_internal(child, dict)
    end

    return y
end

deepcopyasroot(x::Root) = deepcopy(x)

function removechild!(child::Child)
    deleteat!(child.parent, findfirst(isequal(child), child.parent.children))
end

function Base.delete!(node::Child)
    "Deletes a node returning it, while connecting the parent and the children of it."

    length = node.length
    parent = node.parent

    index = findfirst(isequal(node), parent.children)
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

function insert(name, node::Child, length)
    "Inserts a node on a branch at 'length' distance ancestral from 'child'."

    if length > node.length
        error("length is longer than the branch leading to its child")
    end
    newchild = Child(name, node.parent, node.length - length)
    push!(newchild.children, node)
    removechild!(node)
    node.parent = newchild
    node.length = length

    return newchild
end

function detach!(node::Child; asroot=false)
    """Deletes a subtree from a node, returning it.
    If asroot, convert the given node to a root."""

    removechild!(node)
    asroot || return node
    root = Root(node.name)
    root.children = copy(node.children)
    for child in root.children
        child.parent = root
    end
    return root
end

function move!(node::Child, parent::Node)
    "Moves a child to a new parent"
    removechild!(node)
    push!(parent.children, node)
    node.parent = parent
end

function namemapof(node::T; unique=true) where T <: Node
    "Returns a name => node dict of all nodes in subtree"

    namemap = Dict{String, Node}(node.name=>node)
    for node in DepthsFirst(node)
        if unique && haskey(namemap, node.name)
            error("node.name is present more than once in tree")
        end
        namemap[node.name] = node
    end
    return namemap
end

abstract type AbstractLineageIterator{T} end
Base.IteratorEltype(::Type{<:AbstractLineageIterator}) = Base.HasEltype()
Base.IteratorSize(::Type{<:AbstractLineageIterator}) = Base.SizeUnknown()

struct ChildLineageIterator <: AbstractLineageIterator{Child}
    x::Child
end

Base.iterate(x::ChildLineageIterator, state=x.x) = (state isa Root) ? nothing : (state, state.parent)
Base.eltype(::Type{ChildLineageIterator}) = Child

struct LineageIterator{T} <: AbstractLineageIterator{T}
    x::Union{T}
end
LineageIterator(x::T) where {T<:Node} = LineageIterator{T}(x)

function Base.iterate(x::LineageIterator, state=x.x)
    state === nothing && return nothing
    return (state isa Root) ? (state, nothing) : (state, state.parent)
end
Base.eltype(::Type{LineageIterator{Child}}) = Union{Root, Child}
Base.eltype(::Type{LineageIterator{Root}}) = Root

lineageof(node::Node) = LineageIterator(node)
childlineageof(node::Child) = ChildLineageIterator(node)

rootof(node::Root) = node
function rootof(node::Child)::Root
    for i in lineageof(node)
        isa(i, Root) && return i
    end
end

struct DepthsFirst{T}
    x::T
    remaining::Vector{Union{Child, T}}
    DepthsFirst(x::T) where {T<:Node} = new{T}(x, Union{Child, T}[x])
end

function Base.iterate(x::DepthsFirst, state=x)
    isempty(x.remaining) && return nothing
    node = pop!(x.remaining)
    isleaf(node) || append!(x.remaining, node.children)
    return (node, x)
end

Base.IteratorEltype(::Type{<:DepthsFirst}) = Base.HasEltype()
Base.IteratorSize(::Type{<:DepthsFirst}) = Base.SizeUnknown()
Base.eltype(::Type{DepthsFirst{T}}) where T = Union{T, Child}

struct BreadthFirst{T}
    start::T
    currentlevel::Vector{Union{Child, T}}
    nextlevel::Vector{Union{Child, T}}
    BreadthFirst(x::T) where {T<:Node} = new{T}(x, Union{Child, T}[x], Union{Child, T}[])
end

function Base.iterate(x::BreadthFirst, state=x)
    if isempty(x.currentlevel)
        copy!(x.currentlevel, x.nextlevel)
        empty!(x.nextlevel)
    end
    isempty(x.currentlevel) && return nothing
    node = pop!(x.currentlevel)
    isleaf(node) || append!(x.nextlevel, node.children)
    return (node, x)
end

Base.IteratorEltype(::Type{<:BreadthFirst}) = Base.HasEltype()
Base.IteratorSize(::Type{<:BreadthFirst}) = Base.SizeUnknown()
Base.eltype(::Type{BreadthFirst{T}}) where T = Union{T, Child}

subtreeof(node::Node) = DepthsFirst(node)
descendantsof(node::Node) = iterate(DepthsFirst(node))[2]
isleaf(node::Node) = isempty(node.children)
isredundant(node::Node) = any(length(i.children) == 1 for i in subtreeof(node))
ispolytomic(node::Node) = any(length(i.children) > 2 for i in subtreeof(node))

function mrcaof(nodes::Array{T, 1}) where T<:Node
    lineages = [collect(lineageof(node)) for node in nodes]
    minlength = minimum(map(length, lineages))
    generations = Array{Node}(undef, (length(nodes), minlength))

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

mrcaof(one::Child, two::Root) = rootof(one) === two ? two : nothing
mrcaof(one::Root, two::Child) = mrcaof(two, one)
function mrcaof(one::Child, two::Child)
    lin1, lin2 = collect(lineageof(one)), collect(lineageof(two))
    len1, len2 = length(lin1), length(lin2)
    smallen, largelen = minmax(len1, len2)
    if len1 != len2
        smallest, largest = ifelse(len1 < len2, (lin1, lin2), (lin2, lin1))
        unsafe_copyto!(largest, 1, largest, largelen-smallen+1, smallen)
    end
    @inbounds for i in 1:smallen
        if lin1[i] === lin2[i]
            return lin1[i]
        end
    end
    return nothing
end

function distance(one::Node, two::Node, mrca::Node)
    "Finds the distance between two nodes"
    dist = 0.0

    for leaf in (one, two)
        leaf isa Root && continue
        leafchild::Child = leaf
        for child in ChildLineageIterator(leafchild)
            child === mrca && break
            dist += child.length
        end
    end
    return dist
end

function distance(one::T1, two::T2) where {T1<:Node, T2<:Node}
    "Finds the distance between two nodes"
    mrca = mrcaof(one, two)
    if mrca === nothing
        error("the two nodes do not share an ancestor")
    end
    return distance(one, two, mrca)
end

function simplify!(node::T) where T<:Node
    "Deletes all redundant nodes in subtree"
    for child in subtreeof(node)
        length(child.children) == 1 && delete!(child)
    end
end

function subtree(nodes::Array{T,1}) where T<:Node
    """Creates a new tree with the same topology, lengths and names, including
    only the nodes given and the ancestors necessary to keep topology same."""

    if isempty(nodes)
        error("no nodes given")
    end

    copyof = Dict{Node, Node}()
    lineage = reverse!(collect(lineageof(nodes[1])))

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
        while !haskey(copyof, node)
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

nodedistances(node::Child) = _ndchild(node, collect(childlineageof(node)))

function nodedistances(node::Child, backto::Child)
    ancestors = collect(childlineageof(node))
    backtoindex = findfirst(isequal(backto), ancestors)
    if backtoindex === nothing
        error("backto must be an ancestor of node")
    end
    return _ndchild(node, ancestors[1:backtoindex], false)
end

function nodedistances(node::Child, backto::Root)
    ancestors = collect(childlineageof(node))
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

    return distance(leafone, leaftwo), leafone, leaftwo
end

function shufflenames!(node::T) where T <: Node
    descendants = descendantsof(node)
    names = [child.name for child in descendants]
    shuffle!(names)

    @inbounds for i in eachindex(descendants)
        descendants[i].name = names[i]
    end

    return nothing
end

reroot!(root::Root) = root
function reroot!(node::Child)
    lineage = collect(childlineageof(node))
    lengths = [node.length for node in lineage]
    root::Root = lineage[end].parent

    # Switch properties of node and root
    root.name, node.name = node.name, root.name
    oldrootchildren = copy(root.children)
    root.children = node.children
    push!(root.children, node.parent)
    for child in root.children
        child.parent = root
    end

    node.parent = lineage[end]
    node.length = lengths[end]
    deleteat!(oldrootchildren, findfirst(isequal(lineage[end]), oldrootchildren))
    node.children = oldrootchildren

    for child in node.children
        child.parent = node
    end

    for childindex in 2:length(lineage)
        child = lineage[childindex]
        # Set length
        child.length = lengths[childindex - 1]

        # Set new parent
        if childindex == 2
            child.parent = root
        else
            child.parent = lineage[childindex - 1]
        end

        # Identify new child
        if childindex == length(lineage)
            newchild = node
        else
            newchild = lineage[childindex + 1]
        end

        # Identify old child to be removed and replace old with new
        oldchild = lineage[childindex - 1]
        @assert oldchild in child.children
        i = findfirst(isequal(oldchild), child.children)
        child.children[i] = newchild
    end

    return root
end

function reroot!(node::Child, length::Number)
    root = rootof(node)
    newnode = insert(root.name, node, length)
    newroot = reroot!(newnode)

    if length(newnode.children) == 1
        delete!(newnode)
    end

    return newroot
end

function midpointof(node::T) where T <: Node
    "Finds the midpoint - i.e. halfway through the diameter"

    diameter, leafone, leaftwo = diameter(node)
    radius = diameter / 2
    mrca = mrcaof(leafone, leaftwo)

    if distance(leafone, mrca, mrca) >= radius
        longestleaf = leafone
    else
        longestleaf = leaftwo
    end

    traveled = 0.0
    for child in childlineageof(longestleaf)
        traveled += child.length

        if traveled > radius
            return child, radius + child.length - traveled
        end
    end

    error("Should not reach this line")
end

function reroot_midpoint!(node::Root)
    "Reroots the tree at the midpoint."

    child, distance = midpointof(node)
    newnode = insert(node.name, child, distance)
    newroot = reroot!(newnode)
    delete!(newnode)
    return newroot
end

struct CompactTree
    rootname::String
    names::Vector{String}
    parentids::Vector{Int}
    lengths::Vector{Float64}
end

function pack(node::Node)
    nodemap = Dict{Node,Int}(node => 0)
    names = Vector{String}()
    parents = Vector{Int}()
    lengths = Vector{Float64}()
    for (n, child) in enumerate(descendantsof(node))
        push!(names, child.name)
        push!(lengths, child.length)
        parentid = nodemap[child.parent]
        push!(parents, parentid)
        nodemap[child] = n
    end
    return CompactTree(node.name, names, parents, lengths)
end

function unpack(c::CompactTree)
    root = Root(c.rootname)
    nodemap = Dict{Int, Node}(0 => root)
    for i in eachindex(c.names)
        parent = nodemap[c.parentids[i]]
        child = Child(c.names[i], parent, c.lengths[i])
        nodemap[i] = child
    end
    reverse!(root.children)
    for child in descendantsof(root)
        isleaf(child) || reverse!(child.children)
    end
    return root
end

include("parsetree.jl")

export Node, Child, Root,
    # Instantiations
    newick,  json, deepcopyasroot, pack, unpack,
    # Single-node mutations
    delete!, insert, detach!, move!,
    # Transversal
    lineageof, childlineageof, rootof, descendantsof, subtreeof, mrcaof,
    # Boolean
    isredundant, ispolytomic, isleaf,
    # Metrics
    distance, nodedistances, diameter, midpointof,
    # Whole-tree rearrangements
    simplify!, shufflenames!, reroot!, reroot_midpoint!
    # Methods returning new trees
    subtree,
    # Misc
    namemapof,

end # module Tree
