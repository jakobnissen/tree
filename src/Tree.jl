module Tree

#=
TODO: Re-implement nodistances and diameter. It's WRONG and slow.


A function that calculates if two trees are identical
With and without considering node labels
With and without considering node lengths
    (can be implemented by operating on copied 0-length tree)
=#

struct Pairs{T}
    x::Vector{T}
end

function Base.iterate(x::Pairs, state::Tuple{Int,Int}=(1,2))
    ai, bi = state
    len = length(x.x)
    if bi > len
        ai += 1
        bi = ai + 1
    end
    bi > len && return nothing
    val = @inbounds (x.x[ai], x.x[bi])
    return val, (ai, bi+1)
end

Base.IteratorSize(::Type{Pairs}) = Base.HasLength()
Base.IteratorEltype(::Type{Pairs}) = Base.HasEltype()
Base.eltype(::Type{Pairs{T}}) where T = Tuple{T, T}
Base.length(x::Pairs) = (length(x.x) * (length(x.x)-1)) >>> 1

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

# This funcion only exists because mutually recursive types cannot exist
# in Julia as of 20191004, so I need some way to tell the compiler that
# child.parent returns Union{Child, Root}, not Node
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
    deleteat!(child.parent.children, findfirst(isequal(child), child.parent.children))
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
    for child in descendantsof(node)
        if unique && haskey(namemap, child.name)
            error("$(node.name) is present more than once in tree")
        end
        namemap[child.name] = child
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

siblingsof(node::Child)::Vector{Child} = node.parent.children

struct DepthsFirst{T}
    remaining::Vector{Union{Child, T}}
end

DepthsFirst(x::T) where {T<:Node} = DepthsFirst{T}(Union{Child, T}[x])

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
    currentlevel::Vector{Union{Child, T}}
    nextlevel::Vector{Union{Child, T}}
end

BreadthFirst(x::T) where {T<:Node} = BreadthFirst{T}(Union{Child, T}[x], Union{Child, T}[])

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

struct GenerationIterator{T}
    currentlevel::Vector{Union{Child, T}}
end

function GenerationIterator(x::T) where {T<:Node}
    return GenerationIterator{T}(Union{Child, T}[x])
end

function Base.iterate(x::GenerationIterator, state=x)
    isempty(x.currentlevel) && return nothing
    res = copy(x.currentlevel)
    empty!(x.currentlevel)
    for i in res
        append!(x.currentlevel, i.children)
    end
    return (res, x)
end

Base.IteratorEltype(::Type{<:GenerationIterator}) = Base.HasEltype()
Base.IteratorSize(::Type{<:GenerationIterator}) = Base.SizeUnknown()
Base.eltype(::Type{GenerationIterator{T}}) where T = Vector{Union{T, Child}}

generationsof(x::Node) = GenerationIterator(x)

nextgenerations(node::Node) = GenerationIterator{Child}(copy(node.children))

subtreeof(node::Node) = DepthsFirst(node)
descendantsof(node::Node) = DepthsFirst{Child}(copy(node.children))

isleaf(node::Node) = isempty(node.children)
leavesof(node::Child) = (child for child in subtreeof(node) if isleaf(child))
leavesof(node::Root) = (child for child in descendantsof(node) if isleaf(child))

function isredundant(node::Node)
    # Use descendantsof not subtree of because type stability
    length(node.children) == 1 && return true
    return any(length(i.children) == 1 for i in descendantsof(node))
end

function ispolytomic(node::Node)
    # Use descendantsof not subtree of because type stability
    length(node.children) > 2 && return true
    return any(length(i.children) > 2 for i in descendantsof(node))
end

function mrcaof(nodes::Vector{Root})
    isempty(nodes) && error("empty vector")
    root = first(nodes)
    length(nodes) == 1 && return root
    return any(nodes[i] !== root for i in 2:length(nodes)) ? nothing : root
end


function mrcaof(nodes::Vector{Node})
    length(nodes) == 1 && return first(nodes)
    roots = Vector{Root}()
    children = Vector{Children}()
    for node in nodes
        (node isa Root) ? push!(roots, node) : push!(children, node)
    end
    if length(roots) == 0
        return mrcaof(children)
    end
    root = first(roots)
    if length(roots) > 1 && any(roots[i] !== root for i in 2:length(roots))
        return nothing
    end
    if all(rootof(child) === root for child in children)
        return root
    else
        return nothing
    end
end

function mrcaof(nodes::Vector{Child})
    isempty(nodes) && error("empty vector")
    length(nodes) == 1 && return @inbounds nodes[1]
    length(nodes) == 2 && return mrcaof((@inbounds nodes[1]), (@inbounds nodes[2]))
    lineages = [collect(lineageof(node)) for node in nodes]
    minlength = minimum(map(length, lineages))
    generations = Array{Node}(undef, (length(nodes), minlength))

    for (i, lineage) in enumerate(lineages)
        generations[i,:] = lineage[end - minlength + 1 : end]
    end

    for generation in 1:minlength
        base = generations[1, generation]
        if all(i === base for i in generations[:, generation])
            return base
        end
    end
    return nothing
end

mrcaof(one::Node, two::Root) = rootof(one) === two ? two : nothing
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

treelength(node::Node) = sum(c.length for c in descendantsof(node))

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
    "Deletes all redundant nodes in subtree, never deletes self"
    for child in descendantsof(node)
        length(child.children) == 1 && delete!(child)
    end
end

function subtree(nodes::Vector{<:Node})
    if isempty(nodes)
        error("no nodes given")
    elseif length(nodes) == 1
        return Root(first(nodes.name))
    end

    mrca = mrcaof(nodes)
    if mrca === nothing
        throw(ArgumentError("Not all nodes are related"))
    end
    root = Root(mrca.name)
    copyof = Dict{Node, Node}()
    copyof[mrca] = root

    # For the first node, create all the rest of the lineage
    for template in reverse!(collect(childlineageof(first(nodes))))
        copyof[template] = Child(template.name, copyof[template.parent], template.length)
    end

    # For each subsequent node
    lineage = Child[]
    for node in nodes[2:end]
        empty!(lineage)
        # While the ancestors have never been seen, keep getting more distanct
        # ancestors
        for template in childlineageof(node)
            if haskey(copyof, template)
                break
            end
            push!(lineage, template)
        end

        # Now instantiate all these unseen ancestors
        for template in reverse!(lineage)
            copyof[template] = Child(template.name, copyof[template.parent], template.length)
        end
    end

    simplify!(root)
    return root
end

# distance
# nodedistances

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
    if backto !== node
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

function diameter(node::Node)
    leafone = _furthestleaf(node, node)
    leaftwo = _furthestleaf(leafone, node)

    return distance(leafone, leaftwo), leafone, leaftwo
end


@inline npairs(n::Int) = (n * (n-1)) >>> 1

function _update(node, cache)
    @inbounds begin
    nchildren = length(node.children)
    if nchildren == 0
        return (1, 0.0, 0.0)
    elseif nchildren == 1
        child = node.children[1]
        value = cache[child]
        return (value[1], value[2], value[3] + child.length)
    # These next two do the same thing, but the first one is faster for nchildren == 2
    elseif nchildren == 2
        return _update2(node, cache)
    else
        return _updaten(node, cache, nchildren)
    end
    end # inbounds
end

function _update2(node, cache)
    # See comments in _updaten
    @inbounds begin
    child1, child2 = node.children[1], node.children[2]
    L1, L2 = child1.length, child2.length
    v1, v2 = cache[child1], cache[child2]
    N1, N2 = v1[1], v2[1]
    I1, I2 = v1[2], v2[2]
    D1, D2 = v1[3], v2[3]

    N = N1 + N2
    I = (npairs(N1)*I1 + npairs(N2)*I2 + N1*N2*(D1+D2+L1+L2)) / (npairs(N1) + npairs(N2) + N1*N2)
    D = (N1*(D1+L1) + N2*(D2+L2)) / N
    end # inbounds
    return (N, I, D)
end

function _updaten(node, cache, nchildren)
    # We cache this to avoid looking up in the dictionary tonnes of times
    @inbounds begin
    childdata = Vector{Tuple{Int, Float64, Float64, Float64}}(undef, nchildren)
    for i in 1:nchildren
        child = node.children[i]
        value = cache[child]
        childdata[i] = (value[1], value[2], value[3], child.length)
    end

    # Here, we calculate the mean internal path distance to all the leaves
    sumd = 0.0
    npaths = 0
    for i in 1:nchildren
        # First add the paths within each child's leaves
        np = npairs(childdata[i][1])
        npaths += np
        sumd += np * childdata[i][2]
        # Then add the paths between the children
        for j in (i+1):nchildren
            np = childdata[i][1] * childdata[j][1]
            npaths += np
            sumd += np * (childdata[i][3] + childdata[i][4] + childdata[j][3] + childdata[j][4])
        end
    end
    internal = sumd / npaths

    # Calculate mean distances
    npaths = 0
    sumd = 0.0
    for i in 1:nchildren
        np = childdata[i][1]
        npaths += np
        sumd += np * (childdata[i][3] + childdata[i][4])
    end
    meandist = sumd / npaths
    end # inbounds
    return (npaths, internal, meandist)
end

function meandistance(node::Node)
    # We keep track of 3 properties for each node:
    # 1) Number of leaves descendants, 2) Mean internal distance between leaves
    # 3) Mean distance from leaf to the node
    cache = Dict{Child, Tuple{Int, Float64, Float64}}()
    generations = reverse!(collect(nextgenerations(node)))
    # We then iteratively merge siblings together WITHOUT needing to visit
    # all the sibling's descendants.
    for generation in generations
        for child in generation
            cache[child] = _update(child, cache)
        end
    end
    return _update(node, cache)[2]
end

function shufflenames!(node::Node)
    descendants = descendantsof(node)
    names = [child.name for child in descendants]
    shuffle!(names)

    @inbounds for i in eachindex(descendants)
        descendants[i].name = names[i]
    end

    return nothing
end

#= This function does two things:
1) Updated leafdistances for a node using the leafdistances of its children
2) Returns an array of one array per child, each containing the leaves of that child,
and the distance to the node. The distance between the leaves are then the sum of all leaves
across two of these inner arrays
=#
function _update_dist_parent!(node::Node, leafdistances::Dict{Node, Vector{Tuple{Int, Float64}}},
                              indexof::Dict{Child,Int})
    node_leafdistances = valtype(leafdistances)()
    dists_by_child = valtype(leafdistances)[]
    for child in node.children
        child_distances = valtype(leafdistances)()
        push!(dists_by_child, child_distances)

        # If child is leaf, it's not in leafdistances, but it leaves are just itself
        if isleaf(child)
            val = (indexof[child], child.length)
            push!(node_leafdistances, val)
            push!(child_distances, val)
        # Else, for each of its leaves, we add the length of the child
        else
            for (leafindex, dist) in leafdistances[child]
                val = (leafindex, dist + child.length)
                push!(node_leafdistances, val)
                push!(child_distances, val)
            end
            pop!(leafdistances, child)
        end
    end
    leafdistances[node] = node_leafdistances
    return dists_by_child
end

function _write_distances!(dists::Matrix{Float64}, dists_by_child::Vector{Vector{Tuple{Int, Float64}}})
    @inbounds for (vec1, vec2) in Pairs(dists_by_child)
        for (index1, dist1) in vec1, (index2, dist2) in vec2
            # This is distance from leaf A to node + distance from leaf B to node
            dist = dist1 + dist2
            dists[index1, index2] = dist
            dists[index2, index1] = dist
        end
    end
end

#= Algorithm: Start from youngest leaves, moving up to the MRCA. Along the way,
keep track of the distance to the node you're at from all its leaf descendants.
Every time you reach a new node, that node is the MRCA of leaf descendants of two different
children, and you already got the distance from those leaves to the node, so just add them.
=#
function distance_matrix(nodes::Vector{Child})
    all(isleaf, nodes) || throw(ArgumentError("All nodes must be leaves"))
    # This dict is node => [(index, dist_from_leaf_to_node) ..] for all node's leaves
    leafdistances = Dict{Node, Vector{Tuple{Int, Float64}}}()
    mrca = mrcaof(nodes)
    generations = reverse!(collect(generationsof(mrca)))
    result = zeros(length(nodes), length(nodes))
    indexof = Dict(i=>n for (n,i) in enumerate(nodes))
    # Start at youngest generation, skipping first that only has leaves
    for generation in generations[2:end]
        for node in generation
            if !isleaf(node)
                dists_by_child = _update_dist_parent!(node, leafdistances, indexof)
                _write_distances!(result, dists_by_child)
            end
        end
    end
    return result
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

function midpointof(node::Node)
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
    DepthsFirst, BreadthFirst, generationsof, nextgenerations, leavesof,
    # Boolean
    isredundant, ispolytomic, isleaf,
    # Metrics
    distance, nodedistances, diameter, midpointof, treelength, meandistance,
    distance_matrix,
    # Whole-tree rearrangements
    simplify!, shufflenames!, reroot!, reroot_midpoint!,
    # Methods returning new trees
    subtree,
    # Misc
    namemapof

end # module Tree
