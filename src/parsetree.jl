
import JSON

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
    _parsechildren(jsonroot, root)
    return root
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

    # Must work with "bear:6.80041"

    rightparenpos = findlast(')', nodestring)

    if rightparenpos === nothing
        stringofchildren = ""
        namestart = 1
    else
        stringofchildren = nodestring[2:rightparenpos-1]
        namestart = rightparenpos + 1
    end

    colonpos = findnext(':', nodestring, namestart)
    if colonpos === nothing
        branchlength = 0.0
        nameend = length(nodestring)
    else
        branchlength = Float64(parse(BigFloat, nodestring[colonpos+1:end]))
        nameend = colonpos - 1
    end

    name = String(nodestring[namestart:nameend])
    name, branchlength, stringofchildren
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
