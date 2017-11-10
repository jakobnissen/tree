###  Getting or setting simple node/tree properties
__rename node__

`node.name = "newname"`

__change branch length__

Let `node` be the node at the end of the branch relative to the root:

`node.length = 45.11`

__check whether a node is a leaf__

`isleaf(node)`

__check whether the node is the root__

`isa(node, Root)`

__check whether a subtree has polytomic nodes__

A polytomic node is a node with more than two children. All fully resolved phylogenetic trees are neither polytomic nor redundant.

`ispolytomic(node)`

__check whether a subtree has redundant nodes__

A redundant node is a node with exactly one child. All fully resolved phylogenetic trees are neither polytomic nor redundant.

`isredundant(node)`

### Functions adding, removing or copying nodes
__creating a new tree__

`root = Root("newrootname")`

__adding a new node__

Let `node` be the parent of the node you wish to add, and `branchlength` the length of the edge from `node` to the new node you are creating. 

`childnode = Child("newnodename", node, branchlength)`

__inserting a node on a branch__

Let `child` be the child of the branch you want to insert a new node. Let `length` be the distance from `child` to the newly inserted node.

`insert("newnodename", child, length)`

__deleting a node from a tree__

This removes references to the node from its parent and children, and connects the children of the node to the parent of the node, perserving the total distance. The new children of the parental node are 
inserted in the same position in the array of children as the original node was. 

`delete!(node)`

__removing a node and all its descendants from a tree, returning it__

`subtree = detach!(node)`

If you want to returned node to be an independent tree with its own root:

`newtree = detach!(node, asroot=true)`

__copy a subtree__

`subtreecopy = deepcopy(node)`

Note: If the ancestral node is a Child, its parent will be the original node's parent. The newly created copy will not feature among the parent's children.

If you want to copy a subtree without any weird parent/child linkage to the old tree, the ancestral node must be a root. To this end, you can use another function that functions like `deepcopy`, except 
that it always returns a root.

`subtreecopy = deepcopyasroot(node)`

__extract a minimal subtree with the given nodes__

This function copies the given nodes, as well as the minimal number of ancestral nodes needed to preserve the topology. The ancestor to all nodes is returned.

`nodes = [node1, node2, node3]`

`mysubtree = subtree(nodes)`

### Parsing/lexing:
__reading a Newick string__

`root = newick(string)`

__converting a subtree to Newick format__

`string = newick(node)`

__reading a JSON string__

`root = json(string)`

__converting a subtree to JSON format__

`string = json(node)`

### Functions returning nodes in tree
__getting all the ancestors of a node__

This returns an array of nodes from the given node back to its root.

`ancestor_array = lineageof(node)`

If you only want its ancestral children and not the ancestral root for type-stability reasons:

`childancestor_array = childlineageof(node)`

__getting a dictionary mapping node names to nodes in subtree__

`nodemap = nodemapof(node)`

__getting the root node__

`root = lineageof(node)[end]`

__getting the children of a node__

`child_array = node.children`

__getting the parent of a node__

`parent = node.parent`

__getting all descendants of a node__

This returns a 1D array of all descendants of the node, not including the node, depth-first.

`child_array = descendantsof(node)`

__getting all leaf descendants of a node__

`child_array = filter(isleaf, descendantsof(node))`

__getting the siblings of a node__

`child_array = [child for child in node.parent.children if child != node]`

__getting the most recent common ancestor of several nodes__

`nodes = [node1, node2, node3]`

`mrca = mrcaof(nodes)`

### Tree measurements
__getting the distance between two nodes__

`distance(node1, node2)`

This function can be significantly sped up by providing the most recent common ancestor to the nodes, if it is already computed. Note that this assumes the given MRCA to be correct.

`distance(node1, node2, mrca)`

__getting distances to all nodes in tree__

This returns a 1D array of (Node, Float64) tuples of all nodes in the tree - *including the descendant of ancestral nodes* - and their distance to the given node.

`nodedistances(node)`

If you only want to include all the descendants of a particular ancestor in the array, you can add that as an argument. It is assumed that `ancestor` is ancestral to `node`.

`nodedistances(node, ancestor)`

__getting tree diameter (longest path in tree), and/or the two leaves at endpoints__

`treediameter, node1, node2 = diameter(node)`

__getting the subtree length (sum of branch lengths)__

`treelength = sum(node.length for node in descendantsof(node))`

### Tree mutation
__rerooting the tree__

Rerooting is done by converting an existing node to the root and changing the tree accordingly. If you want to insert a new root, you must first use `insert` followed by `reroot!`.

Let `node` be the child you want to convert to the new root.

`reroot!(node)`

__shuffle the names of all nodes in a subtree__

The given node is not included, only the descendants.

`shufflenames!(node)`

__make tree nonredundant__

This function deletes all redundant nodes in the subtree.

`simplify!(node)`

### Current operations are not implemented
* Parsing/writing PhyloXML
* Check whether two trees are exactly equal (names, lengths, topology). Implement a hash function.
* Make an efficient iterator of nodes, breath-first
* Make an efficient iterator of nodes, depth-first
* Find the midpoint outgroup
