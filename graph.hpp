// graph.hpp

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   Graph utilities
*
*/


#ifndef GRAPH_HPP_73854  // random number
#define GRAPH_HPP_73854


#include "common.hpp"
using namespace Common_sp;



namespace Common_sp
{



struct DiGraph : Root
// Directed graph
// Bi-directional access
// n = number of nodes
// m = number of arcs
{
  struct Arc; 
  struct Node;
  

  typedef  unordered_map <const Node*, const Node*>  Node2Node;  
    // !nullptr


  struct Node : VirtNamed, DisjointCluster
  {
    friend DiGraph;
    
    const DiGraph* graph {nullptr};
      // nullptr <=> *this is detach()'ed
  private:
    List<Node*>::iterator graphIt;
      // In graph->nodes
  public:
    List<Arc*> arcs [2 /*bool: out*/];

    Node* scc {nullptr};  
      // Strongly-connected component
      // = root of SCC-subtree of DFS tree
      // scc->orderDfs <= orderDfs 
    size_t orderDfs {0};
      // Order by depth-first search
      // 0 <=> not visited by DFS
  private:
    // Auxiliary      
    bool inStack {false};  
      // => orderDfs
  public:
    bool reachable {false};


    explicit Node (DiGraph &graph_arg)
			{ attach (graph_arg); }
    Node (const Node &other)
      : VirtNamed (other)
      , DisjointCluster (other)
      , scc      (other. scc)
      , orderDfs (other. orderDfs)
      , inStack  (other. inStack)
      {}
      // To be followed by: attach()
    Node* copy () const override
      { return new Node (*this); } 
   ~Node ();
      // Remove this from graph 
      // Time: O(m) for all nodes
    void qc () const override;
    void saveText (ostream& os) const override;
      // Invokes: getName(), saveContent()


  protected:
    virtual void saveContent (ostream &/*os*/) const 
      {}
  public:
    string getName () const override
      // Return: !empty()
	    { ostringstream oss;
		  	oss << this;
			  return oss. str ();
		  }

    Node* getDisjointCluster ()
      { return static_cast <Node*> (DisjointCluster::getDisjointCluster ()); }
  	void attach (DiGraph &graph_arg);
      // Requires: !graph; no Arc's
      // Invokes: graph_arg.nodes.push_back(this)
      // Time: O(1)
		virtual string getHumanName () const
		  { return getName (); }
		virtual string getLeafName () const
		  { return getName (); }
    bool isIncident (const Node* n,
                     bool out) const;
      // Return: n is among arcs[out]->node[out]
    bool isIncidentExcept (const Node* n,
                           bool out) const;
    size_t getDegree () const
      { return   arcs [false]. size () 
      	       + arcs [true].  size (); 
      }
    VectorPtr<Node> getNeighborhood (bool out) const;
    VectorPtr<Node> getNeighborhood () const
      { return getNeighborhood (false) << getNeighborhood (true); }
      // Return: may have duplicates
    VectorPtr<Node> getChildren () const
      { return getNeighborhood (false); }
    void deleteNeighborhood (bool out);
    void getDependents (bool out,
                        unordered_set<const Node*> &dependents) const
      { dependents. insert (this);
        for (const Arc* arc : arcs [out])
          arc->node [out] -> getDependents (out, dependents);
      }
    void setReachable (bool reachItself);
      // From parent to children
      // Output: reachable = true
      // Invokes: setReachable(true)
  private:
    Node* setScc (size_t &visitedNum,
                  stack<Node*,vector<Node*>> &sccStack);
      // If node n is reachable from this then
      //   the SCC of n is reachable from this and the SCC is a subtree of DFS tree
      // Output: scc, orderDfs in nodes reachable from this
      // Tarjan's alogorithm
      // Return: n s.t. n->inStack and there is a path from this to n
      // Requires: n->inStack <=> there is a path from n to this
      // Time: O(n + m) for all nodes
    virtual void contractContent (const Node* /*from*/) {}
      // Required time: O(1)
  public:      
    void contract (Node* from);
      // Update: this: No parallel arcs
      // Invokes: contractContent(from), Arc::contractContent(), delete from
      // Requires: from != this
      //           No parallel arcs
      // Time: O(n + m log n) for all nodes
    void isolate ();
      // Make degree = 0
    void detach ();
      // Output: graph = nullptr
      // Requires: No Arc's
      // Invokes: list::erase()
    Node* copyGraph (bool out,
                     Node2Node &node2node) const;
  };


  struct Arc : Root
  {
    friend DiGraph;
    friend Node;
   
    array <Node*, 2 /*bool: out*/> node;
      // !nullptr
  private:
    array <List<Arc*>::iterator, 2 /*bool: out*/> arcsIt;
      // in node[b]->arcs[!b]
  public:


    Arc (Node* start,
         Node* end)
      { node [false] = nullptr;
      	node [true]  = nullptr;
      	attach (start, end);
      }
  private:
  	void attach (Node* start,
                 Node* end);
      // Adds *this to the graph
      //   node[i]->arcs[!i].push_back(this)
      // Requires: !node[i]
      // Time: O(1)
  public:
  	Arc (const Arc& other)
  	  : Root (other)
  	  { node [false] = nullptr;
      	node [true]  = nullptr;
      }
      // To be followed by: attach()
    Arc* copy () const override
      { return new Arc (*this); } 
   ~Arc ();
      // Remove this from node->graph
      // Time: O(1)


    virtual void saveContent (ostream &/*os*/) const 
      {}
  private:
    virtual void contractContent (const Arc* /*from*/) 
      {}
      // Required time: O(1)
  public:
    void setNode (Node* newNode,
                  bool out);
      // Update: node[out] = newNode
      // Preserves the ordering of node[!out]->arcs[out]
    bool selfLoop () const
      { return node [false] == node [true]; }
  };


  List<Node*> nodes;
    // size() == n


  DiGraph () = default;
  typedef  unordered_map <const Node* /*old*/, Node* /*new*/>  Old2new;
  DiGraph (const DiGraph &other)
    { Old2new old2new (other. nodes. size ());
    	init (other, old2new);
    }
  DiGraph (const DiGraph &other,
           Old2new &old2new)
    { init (other, old2new); }
private:
  void init (const DiGraph &other,
             Old2new &old2new);
    // Output: old2new
public:
  DiGraph* copy () const override
    { return new DiGraph (*this); } 
 ~DiGraph ()
    { deleteNodes (); }
  void qc () const override;
  void saveText (ostream &os) const override;
  bool empty () const final
    { return nodes. empty (); }
  void clear () final
    { return deleteNodes (); }


  void deleteNodes ();
    // Invokes: Node::delete
    // Time: O(n + m)
  void connectedComponents ();
    // Output: Node::getConnectedComponent()  
    // Invokes: DisjointCluster::init()
  void scc (); 
    // Output: Node::{scc,orderDfs}
    // Invokes: Node::setScc()
    // Time: O(n + m)
  void contractScc ();
    // Output: DAG
    // Requires: After scc()
    // Invokes: Node::contract()
    // Time: O(n + m log n)
  void clearReachable ()
    { for (Node* node : nodes)
        node->reachable = false;
    }
  VectorPtr<Node> getEnds (bool out) const;
    // Return: distinct, !nullptr
    // Input: out: false - roots
    //             true  - leaves
  const Node* getRoot (bool out) const
		{ const VectorPtr<Node> ends (getEnds (out));
			if (ends. size () == 1)
			  return * ends. begin ();
			return nullptr;
		}
  static Node2Node reverse (const Node2Node& old2new);
  void borrowArcs (const Node2Node &other2this,
                   bool parallelAllowed);
    // Input: other2this: other node to *this node
    // Invokes: new Arc
    // Time: O(|other2this| log|other2this| outdegree_max (parallelAllowed ? 1 : outdegree'_max)),
    //          where outdegree_max  = max(outdegree(other2this.keys()),
    //                outdegree'_max = max(outdegree(other2this.values())
};



struct Tree : DiGraph
// m = n - 1
// Parent <=> out = true
{
	struct TreeNode : DiGraph::Node
	{
	  friend Tree;
	  bool frequentChild {false};
	    // For a directed tree
	  size_t frequentDegree {0};
	    // For an undirected tree
	  size_t leaves {0};

		TreeNode (Tree &tree,
		          TreeNode* parent)
			: DiGraph::Node (tree)
			{ setParent (parent); }
		  // Input: parent_arg: may be nullptr
		void qc () const override;
   	void saveText (ostream &os) const override;
   	  // Invokes: getName(), saveContent(), getSaveSubtreeP()

		string getHumanName () const final
		  { return getLcaName (); }

    virtual bool getSaveSubtreeP () const 
      { return true; }
	  virtual double getParentDistance () const
	    { return -1.0; }
	    // Return: -1 || >= 0
	  virtual string getNewickName (bool /*minimal*/) const
	    { return getName (); }
		static string name2newick (const string &s)
		  { return to_url (s); }
	private:
	  void printNewick_ (ostream &os,
	                     bool internalNames,
	                     bool minimalLeafName) const;
	    // Input: os.setprecision
	    // Invokes: getParentDistance(), getNewickName(), name2newick()
  public:
	  const Tree& getTree () const
  	  { return * static_cast <const Tree*> (graph); }
		bool isLeaf () const
		  { return arcs [false]. empty (); }
		virtual bool isLeafType () const
		  { return isLeaf (); }
		virtual bool isInteriorType () const
		  { return false; }
		const TreeNode* getParent () const
			{ return arcs [true]. empty () ? nullptr : static_cast <TreeNode*> (arcs [true]. front () -> node [true]); }
		  // Return: nullptr <=> root
		const TreeNode* getAncestor (size_t height) const;
		  // Return: !nullptr
		  // getAncestor(0) = this
		void setParent (TreeNode* newParent);
		  // Update: *newParent; makes *this the last child of *newParent
		  //         getTree()->root if !newParent
    void printAncestors (const TreeNode* end) const;
    struct TipName : Root
    { string name; 
      size_t depth {0}; 
      TipName () = default;
      TipName (const string &name_arg,
               size_t depth_arg)
        : name (name_arg)
        , depth (depth_arg)
        {}
      void saveText (ostream &os) const override
        { os << name << '\t' << depth; }
    };
    TipName getTipName () const;
      // Return: identification of *this by a tip
		size_t getTopologicalDepth () const
		  { if (const TreeNode* parent_ = getParent ())
		  		return parent_->getTopologicalDepth () + 1;
		  	return 0;
		  }
		struct NodeDist
		{ const TreeNode* node;
		  double dist;
		  static bool distLess (const NodeDist &x,
		                        const NodeDist &y) 
		    { return x. dist < y. dist; }
		  bool operator< (const NodeDist &other) const
		    { return node < other. node; }
		  bool operator== (const NodeDist &other) const
		    { return node == other. node; }
		};
		void getSubtreeHeights (Vector<NodeDist> &nodeHeights) const;
		  // Append: nodeHeights: interior nodes
		  // Invokes: getParentDistance()
    void getLeafDepths (Vector<NodeDist> &leafDepths) const
      { getLeafDepths_ (leafDepths, true); }
  private:
    void getLeafDepths_ (Vector<NodeDist> &leafDepths,
                         bool first) const;
		  // Append: leafDepths: leaves
		  // Invokes: getParentDistance()
	public:
		size_t getHeight () const;
		  // Return: 0 <=> isLeaf()
		size_t getInteriorHeight () const;
		double getDistanceHeight () const;
		  // Invokes: getParentDistance()
    void getBifurcatingInteriorBranching (size_t &bifurcatingInteriorNodes,
                                          size_t &branches) const;
      // Update: bifurcatingInteriorNodes, branches
      //         branches >= bifurcatingInteriorNodes
      //         branches <= 2 bifurcatingInteriorNodes
	  double getRootDistance () const
		  { if (const TreeNode* parent_ = getParent ())
		  		return parent_->getRootDistance () + getParentDistance ();
		  	return 0;
		  }
		bool descendantOf (const TreeNode* ancestor) const
		  { if (! ancestor)
		  	  return true;
		  	if (this == ancestor)
		  		return true;
		  	if (const TreeNode* parent_ = getParent ())
		  		return parent_->descendantOf (ancestor);
		  	return false;
		  }
		const TreeNode* getPrevAncestor (const TreeNode* ancestor) const
		  { if (! ancestor)
		      return nullptr;
		    const TreeNode* parent_ = getParent ();
		    if (parent_ == ancestor)
		      return this;
		    if (! parent_)
		      return nullptr;
		    return parent_->getPrevAncestor (ancestor);
		  }
		double getPathLength (const TreeNode* ancestor) const
		  { if (this == ancestor)
		  		return 0.0;
		  	if (const TreeNode* parent_ = getParent ())
		  		return getParentDistance () + parent_->getPathLength (ancestor);
		  	if (! ancestor)
		  	  return 0.0;
		  	return numeric_limits<double>::quiet_NaN ();
		  }
		size_t getSubtreeSize (bool countLeaves) const;
		  // Return: number of Arc's in the subtree
		  //         0 <= isLeaf()
		  // Input: !countLeaves => leaf arcs are ignored
    double getSubtreeLength () const;
		  // Return: 0 <= isLeaf()
		// Output: leaves in the subtree
		void subtreeSize2leaves ();
		  // Output: leaves: number of Arc's in the subtree
		void setLeaves ();
    void setLeaves (size_t leaves_arg);
    //
		size_t getLeavesSize () const;
		void children2frequentChild (double rareProb);
		  // Input: leaves
		  // Output: TreeNode::frequentChild
		  // Invokes: isInteriorType()
    void getLeaves (VectorPtr<TreeNode> &leafVec) const;
      // Update: leafVec
		const TreeNode* getClosestLeaf (size_t &leafDepth) const;
		  // Return: !nullptr
		  // Output: leafDepth; nullptr <=> Return = this
    const TreeNode* getOtherChild (const TreeNode* child) const;
      // Return: May be nullptr; != child
      // Requires: getChildren().size() <= 2
    const TreeNode* getDifferentChild (const TreeNode* child) const;
      // Return: !nullptr; != child
    const TreeNode* getFirstDecendant () const;
    const TreeNode* getLastDecendant () const;
    string getLcaName () const;
	  void childrenUp ();
	    // Children->setParent(getParent())
	    // Post-condition: arcs[false].empty()
	  TreeNode* isTransient () const
	    { return arcs [false]. size () == 1 
	    	        ? static_cast <TreeNode*> (arcs [false]. front () -> node [false]) 
	    	        : nullptr; 
	    }
	    // Return: Single child of *this
	  void detachChildrenUp ();
	    // Invokes: childrenUp(), detach()
	  TreeNode* isolateTransient ()
			{ TreeNode* transient = isTransient ();
				if (transient)
					detachChildrenUp ();
			  return transient;
			}
	  bool deleteTransient ()
	    { if (! isolateTransient ())
		 			return false;
		    delete this;
        return true;
	    }
	  void deleteSubtree ();
	    // Does not delete *this
	    // Postcondition: isLeaf()
	  const TreeNode* makeRoot ();
	    // Redirect TreeArc's so that this = getTree()->root
	    // Return: old getTree()->root, !nullptr
    void getArea (uint radius,
                  VectorPtr<TreeNode> &area,
                  VectorPtr<TreeNode> &boundary) const
      { getArea_ (radius, nullptr, area, boundary); }
      // Update: area: connected TreeNode's with one root, distinct
      //         boundary: distinct; degree = 1 in the subgraph
      //         area.contains(boundary)
  private:
    void getArea_ (uint radius,
                   const TreeNode* prev,
                   VectorPtr<TreeNode> &area,
                   VectorPtr<TreeNode> &boundary) const;
      // Update: area, boundary
      //         area.contains(boundary)
  public:
    void getDistanceArea (double radius,
                          VectorPtr<TreeNode> &area,
                          VectorPtr<TreeNode> &boundary) const
      { getDistanceArea_ (radius, nullptr, area, boundary); }
      // Update: area: connected TreeNode's with one root, distinct
      //         boundary: distinct; degree = 1 in the subgraph
      //         area.contains(boundary)
  private:
    void getDistanceArea_ (double radius,
                           const TreeNode* prev,
                           VectorPtr<TreeNode> &area,
                           VectorPtr<TreeNode> &boundary) const;
  public:
    void getClosestLeaves (size_t neighbors_max,
                           Vector<NodeDist> &neighbors) const
      { getClosestLeaves_ (nullptr, 0.0, neighbors_max, neighbors); }
      // Output: neighbors: sorted by distLess
  private:
    void getClosestLeaves_ (const Tree::TreeNode* prev,
                            double distance,
                            size_t neighbors_max,
                            Vector<NodeDist> &neighbors) const;
      // Update: neighbors: sorted by distLess
  public:
    void getSubtreeArea (const VectorPtr<Tree::TreeNode> &possibleBoundary,
	                       VectorPtr<Tree::TreeNode> &area,
                         VectorPtr<Tree::TreeNode> &boundary) const;
      // Output: area: connected TreeNode's with one root, distinct
      //         boundary: distinct; degree = 1 in the subgraph
      //         area.contains(boundary)
    template <typename StrictlyLess>
    	void sort (const StrictlyLess &strictlyLess)
  			{ VectorPtr<DiGraph::Node> children (getChildren ());
  				for (const DiGraph::Node* child : children)
  				  const_static_cast <TreeNode*> (child) -> sort (strictlyLess);
  				Common_sp::sort (children, strictlyLess);
  				// To reorder arcs[false]
  				for (const DiGraph::Node* child : children)
  				{ const TreeNode* s = static_cast <const TreeNode*> (child);
  				  var_cast (s) -> setParent (var_cast (s->getParent ()));  
  				}
  			}
	};
	const TreeNode* root {nullptr};
	  // nullptr <=> nodes.empty()
  static const char objNameSeparator {':'};


  Tree () = default;
  void qc () const override;
	void saveText (ostream &os) const override
	  { if (root)
	  	  root->saveText (os);
      os << endl;
    }


  void printNewick (ostream &os,
                    bool internalNames,
                    bool minimalLeafName) const
    { root->printNewick_ (os, internalNames, minimalLeafName);
    	os << ';' << endl;
    }
    // Input: internalNames <=> print name at each internal node
  void printAsn (ostream &os) const;	
    // http://www.ncbi.nlm.nih.gov/tools/treeviewer/biotreecontainer/
  void printArcLengths (ostream &os) const;
    // Output: os: <printArcLengthsColumns()>
    // Requires: getParentDistance() > 0 for all nodes except root
  static string printArcLengthsColumns ()
    { return "<node name> <arc length> <depth length> <log(<parent arc length>/<arc length>)"; }
  double getLength () const
    { return root->getSubtreeLength (); }
  double getAveArcLength () const;
  struct Patristic
  { const TreeNode* leaf1 {nullptr};
    const TreeNode* leaf2 {nullptr};  
      // != nullptr
      // leaf1->getName() < leaf2->getName()
    double distance {0.0};
    Patristic (const TreeNode* leaf1_arg, 
               const TreeNode* leaf2_arg,
               double distance_arg);        
    Patristic () = default;
  };
  Vector<Patristic> getLeafDistances () const;
  void setLeaves ()
    { if (root)
        var_cast (root) -> setLeaves (); 
    }
  size_t size (bool countLeaves) const
    { return nodes. size () <= 1
               ? countLeaves
               : 1 + root->getSubtreeSize (countLeaves); 
    }
  size_t countInteriorNodes () const;
    // Input: TreeNode::isInteriorType()
  bool isStar () const
    { return countInteriorNodes () == 1; }
  size_t getInteriorHeight () const
    { if (root && root->isInteriorType ())
        return root->getInteriorHeight ();
      return 0;
    }
  static size_t radius2boundarySize (uint radius) 
    { return radius 
               ? 3 * powInt (2, radius - 1) 
               : 1; 
    }
    // Requires: binary tree
  double getBifurcatingInteriorBranching () const;
    // For unrooted tree
    // Return: if !root->isInteriotType() then -1 else between 1 and 2
    // # Bifurcating interior nodes = [1 1]' [[branching 0]' [1 1]']^depth [1 0] = \sum_{i=0}^depth branching^i = branching^{depth+1} - 1
    //   Vector meaning: [open_nodes closed_nodes]
    // # Leaves = # Bifurcating interior nodes + 1 = branching^{depth+1}
  size_t countInteriorUndirectedArcs () const;
    // Arc is interior <=> arc's nodes are interior
    // Return: <= countInteriorNodes()
    // Invokes: countInteriorNodes()
    
  struct LcaBuffer
  { VectorPtr<TreeNode> vec1;
  	VectorPtr<TreeNode> vec2;
  	void clear ()  
  	  { vec1. clear (); vec2. clear (); }
  };
  static const TreeNode* getLca (const TreeNode* n1,
                                 const TreeNode* n2,
                                 LcaBuffer &buf);
    // Return: nullptr <=> !n1 || !n2
  static const TreeNode* getLca (const VectorPtr<TreeNode> &nodeVec,   
	                               Tree::LcaBuffer &buf);
    // Return: nullptr <= nodeVec.empty()
    // Input: nodeVec: may be nullptr
  static Set<const TreeNode*> getParents (const VectorPtr<TreeNode> &nodeVec, 
	                                        Tree::LcaBuffer &buf);
    // Return: !nullptr, !contains(getLca(nodeVec)), contains(nodeVec)
    // Invokes: getLca(nodeVec)
  static VectorPtr<TreeNode>& getPath (const TreeNode* n1,
								                       const TreeNode* n2,
								                       const TreeNode* ca,
									                     const TreeNode* &lca,
								                       LcaBuffer &buf);
    // Return: reference to buf (buf.vec1 or buf.vec2)
    // Input: ca: may be nullptr
    // Output: path: sequential arcs on the path from n1 to n2 or reversed, distinct, !nullptr
    //         lca: !nullptr
    // Requires: ca is the common ancestor of n1 and n2
  VectorPtr<TreeNode> leaves2lcas (const Set<const TreeNode*> &leaves) const;
    // Opposite to TreeNode::getLeaves()
    // Return: size() = 1 => Return.front().getLeaves().contains(leaves)
    // Invokes: setLeaves()
    // Time: n log(n)
  void setFrequentChild (double rareProb);
    // Input: 0 <= rareProb < 0.5
    // Output: TreeNode::frequentChild: statistically consistent estimate
    // Invokes: setLeaves(), children2frequentChild()
  void setFrequentDegree (double rareProb);
    // Input: 0 <= rareProb < 0.3
    // Output: TreeNode::frequentDegree: statistically consistent estimate
    // Invokes: setLeaves(), isLeafType()
  void setRoot ();
    // Output: root
  size_t deleteTransients ();
    // Return: # TreeNode's delete'd
  virtual void deleteLeaf (TreeNode* leaf,
                           bool /*deleteTransientAncestor*/) 
    { delete leaf; }
  size_t restrictLeaves (const StringVector &leafNames,
                         bool deleteTransientAncestor);
    // Return: # leaves delete'd
    // Input: leafNames: sort()'ed
    // Invokes: isLeafType(), deleteLeaf()

  template <typename StrictlyLess>
    void sort (const StrictlyLess &strictlyLess)
      { if (root)
      	  var_cast (root) -> sort (strictlyLess); 
      }
private:
	static bool strictlyLess_std (const DiGraph::Node* a,
	                              const DiGraph::Node* b);

public:
  void sort ()
    { setLeaves ();
      sort (strictlyLess_std); 
    }
};



struct TopologicalSort : Root
// Usage: 
//   while (Node* n = ts.getFront ())  ...
// Time: O(graph.n + graph.m)
{
	DiGraph &graph;
	const bool out;
private:
	VectorPtr<DiGraph::Node> order;
	size_t index {0};
	  // <= order.size()
public:
	
	
	TopologicalSort (DiGraph &graph_arg,
	                 bool out_arg);
	  
	  
	const DiGraph::Node* getFront ();
		// delete Arc's of graph
		// graph is DAG <=> all Arc's are delete'd, all Node's are returned
};



}



#endif

