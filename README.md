# Tree

A Julia library for creating and manipulating trees

### Design

The tree is implemented as a hierarchy of nodes which refer to their parent and children. All operations done on subtrees takes the ancestral node as argument. Thus a "tree" is really just a reference to 
the ancestral node. As such, there is no tree object. Instead the properties of trees emerge from the node objects.

In each tree, one node is special, namely the root node. This is the only one with no parent. All other nodes have a parent, and a length associated with them. This represents the distance from the node 
to its parent.

Unrooted trees and dendrograms are special cases of trees - namely the ones where the location of the root node and the lengths, respectively, are unimportant.

### Tutorial

Trees can be instantiated manually. This is not particularly easy, but it serves for demonstration purposes. As the only node without a parent, the root must be instantiated first by providing a name for 
it:

`root = Root("sample_098")`

After this, Child nodes can be instantiated. These require a name, a parent, and a distance to parent.

`zebrafish = Child("Danio rerio", root, 7.711)`

There is no fundamental difference between an internal node and a leaf. Let's instantiate a primate ancestor and its descendants, humans and chimps, and also some rodents:

```
primate = Child("primate ancestor", root, 3.829)
Child("Homo sapiens", primate, 4.332)
Child("Pan troglodytes", primate, 4.165)
rodent = Child("rodent ancestor", root, 2.991)

Child("Mus musculus", rodent, 4.003)
Child("Rattus rattus", rodent, 9.195)
```

The tree now looks something like this:

                                     /------------------------- Rattus rattus
                /------rodent ancestor
               /                     \------------- Mus musculus
              /
    sample_098--------------------- Danio rerio
              \
               \                          /----- Pan troglodytes
                \-------- primate ancestor
                                          \------ Homo sapiens
                                
And can be interpreted either as rooted or unrooted.

This tree is not fully resolved. The common ancestor of mouse and primates is not included. To create it, we need to insert it at an already existing branch, namely the one ancestral to the primate 
ancestor.

`mammal = insert("mammal ancestor", primate, 1.4)`

And move the rodent ancestor to the mammal node.

`move!(rodent, mammal)`

The resulting tree looks something like this:

               /--------------------- Danio rerio
              /
    sample_098                                          /------------------------- Rattus rattus
              \                  /-------rodent ancestor
               \                /                       \------------- Mus musculus
                \mammal ancestor 
                                \                 /Pan troglodytes
                                 \primate ancestor
                                                  \Homo sapiens

In reality, you would probably want to instantiate trees from a Newick or JSON string. Here's how to load the same tree:

    newickstring = "(Danio rerio:7.711,((Homo sapiens:4.332,Pan troglodytes:4.165)primate ancestor:1.4,(Mus musculus:4.003,Rattus rattus:9.195)rodent ancestor:2.991)mammal ancestor:2.429)sample_098;"
    root = newick(newickstring)

For tree manipulation functions, check the associated documentation.
