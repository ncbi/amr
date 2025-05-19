// graph.cpp

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


#undef NDEBUG

#include "graph.hpp"

#include "common.inc"



namespace Common_sp
{
 


//////////////////////////////// DiGraph /////////////////////////////////////

// DiGraph::Node

void DiGraph::Node::attach (DiGraph &graph_arg)
{
  ASSERT (! graph);
#ifndef NDEBUG
	for (const bool b : {false, true})
	  ASSERT (arcs [b]. empty ());
#endif
  
  graph = & graph_arg;
  graph_arg. nodes << this;
  graphIt = graph_arg. nodes. end (); 
  graphIt--;
}

 

DiGraph::Node::~Node ()
{
	for (const bool b : {false, true})
    deleteNeighborhood (b);
  if (graph)
    const_cast <DiGraph*> (graph) -> nodes. erase (graphIt);
}

 

void DiGraph::Node::qc () const
{
  if (! qc_on)
    return;
  Root::qc ();
    
  QC_IMPLY (graph, *graphIt == this);
  if (! graph)
  	for (const bool b : {false, true})
  	  { QC_ASSERT (arcs [b]. empty ()); }
}



void DiGraph::Node::saveText (ostream &os) const
{
  os << getName () << "  ";
  saveContent (os);
  
  if (orderDfs)
    os << "  DFS_order = " << orderDfs;
  if (scc)
    os << "  SCC: " << scc->getName ();

  os << endl;  
	for (const bool b : {false, true})
  {
    os << "  " << (b ? "Out" : "In") << ":" << endl;
    for (const Arc* arc : arcs [b])
    {
      os << "    " << arc->node [b] -> getName ();
      os << ": ";
      arc->saveContent (os);
      os << endl;
    }
  }
}

 

bool DiGraph::Node::isIncident (const DiGraph::Node* n,
                                bool out) const
{
	ASSERT (n);
	ASSERT (graph);
	ASSERT (n->graph == graph);

  for (const Arc* arc : arcs [out])
    if (arc->node [out] == n)
    	return true;
  return false;
}



bool DiGraph::Node::isIncidentExcept (const DiGraph::Node* n,
	                                    bool out) const
{
	ASSERT (n);
	ASSERT (graph);
	ASSERT (n->graph == graph);

  for (const Arc* arc : arcs [out])
    if (arc->node [out] != n)
    	return true;
  return false;
}



VectorPtr<DiGraph::Node> DiGraph::Node::getNeighborhood (bool out) const
{
	VectorPtr<Node> s;  s. reserve (arcs [out]. size ());
  for (const Arc* arc : arcs [out])
    s << arc->node [out];
	return s;
}

 

void DiGraph::Node::deleteNeighborhood (bool out)
{
  for (Iter <List<Arc*> > iter (arcs [out]); iter. next (); )
    delete *iter;  // changes arcs[i]
}



void DiGraph::Node::setReachable (bool reachItself)
{
  if (reachable)
    return;
  if (reachItself)
    reachable = true;
	for (Arc* arc : arcs [false])
		static_cast <Node*> (arc->node [false]) -> setReachable (true);
}



DiGraph::Node* DiGraph::Node::setScc (size_t &visitedNum,
                                      stack<DiGraph::Node*,vector<DiGraph::Node*>> &sccStack)
// Tarjan's alogorithm:
/*
  Input: Graph G = (V, E), Start node v0
  
  index = 0                       // DFS node number counter 
  S = empty                       // An empty stack of nodes
  tarjan(v0)                      // Start a DFS at the start node
  
  procedure tarjan(v)
    v.index = index               // Set the depth index for v
    v.lowlink = index
    index = index + 1
    S.push(v)                     // Push v on the stack
    forall (v, v') in E do        // Consider successors of v 
      if (v'.index is undefined)  // Was successor v' visited? 
        tarjan(v')                // Recurse
        v.lowlink = min(v.lowlink, v'.lowlink)
      elseif (v' in S)            // Is v' on the stack?
        v.lowlink = min(v.lowlink, v'.lowlink)
    if (v.lowlink == v.index)     // Is v the root of an SCC?
      print "SCC:"
      repeat
        v' = S.pop
        print v'
      until (v' == v)
*/
{
  if (orderDfs)
  {
    if (inStack)
      return scc;
    else
      return nullptr;
  }
  ASSERT (! inStack);
  
  visitedNum++;
  orderDfs = visitedNum;

  sccStack. push (this);
  inStack = true;

  scc = this;

  for (Arc* arc : arcs [true])
    if (Node* lowNode = arc->node [true] -> setScc (visitedNum, sccStack))
      if (lowNode->orderDfs < scc->orderDfs)
        scc = lowNode;

  ASSERT (scc->orderDfs <= orderDfs);  
  if (scc->orderDfs < orderDfs)  
  {
    ASSERT (scc != this);
    return scc;
  }
  
  // this is the root of a SCC
  ASSERT (scc->orderDfs == orderDfs);  
  ASSERT (scc == this);
  for (;;)
  {
    Node* n = sccStack. top ();
    sccStack. pop ();
    n->inStack = false;
    if (n == this)
      break;
    n->scc = this;
  }

  return nullptr;
}

 

void DiGraph::Node::contract (Node* from)
{
  ASSERT (from);
  ASSERT (this != from);
  
  contractContent (from);

	for (const bool b : {false, true})
  {  
    unordered_map <Node*, Arc*> m;  m. rehash (arcs [b]. size ());
    for (Arc* arc : arcs [b])
      m [arc->node [b]] = arc;

    for (Iter <List<Arc*> > iter (from->arcs [b]); iter. next (); )
    {
      Arc* arc = *iter;
      if (m. find (arc->node [b]) == m. end ())
      {
        iter. erase ();
        arc->node [! b] = this;
        arcs [b]. push_back (arc);
        arc->arcsIt [! b] = arcs [b]. end ();
        arc->arcsIt [! b] --;
      }
      else
        m [arc->node [b]] -> contractContent (arc);
    }
  }
  
  delete from;
}



void DiGraph::Node::isolate ()
{ 
  for (const bool b : {false, true})
  	while (! arcs [b]. empty ())
  	{ 
  		Arc* a = arcs [b]. front ();
  		delete a;
    }
}



void DiGraph::Node::detach ()
{
  ASSERT (graph);
#ifndef NDEBUG
	for (const bool b : {false, true})
	  ASSERT (arcs [b]. empty ());
#endif
  var_cast (graph) -> nodes. erase (graphIt);  
  graphIt = var_cast (graph) -> nodes. end (); 
  graph = nullptr;
}



DiGraph::Node* DiGraph::Node::copyGraph (bool out,
                                         Node2Node &node2node) const
{
  {
    const Node* n = nullptr;
    if (find (node2node, this, n))
      return var_cast (n);
  }
  Node* n = copy ();
  node2node [this] = n;
  for (const Arc* arc : arcs [out])
  {
    const Node* child = arc->node [out];
    Node* child1 = child->copyGraph (out, node2node);
    new Arc (child1, n);
  }
  return n;
}


 

// DiGraph::Arc

void DiGraph::Arc::attach (Node* start,
                           Node* end) 
{      
  ASSERT (start);
  ASSERT (end);
  ASSERT (start->graph);
  ASSERT (start->graph == end->graph);
  ASSERT (! node [false]);
  ASSERT (! node [true]);

  node [false] = start;
  node [true]  = end;

	for (const bool b : {false, true})
  {
  	List<Arc*>& arcs = node [b] -> arcs [! b];
    arcs. push_back (this);
    arcsIt [b] = arcs. end ();
    arcsIt [b] --;
  }
}

 

DiGraph::Arc::~Arc ()
{
	for (const bool b : {false, true})
    node [b] -> arcs [! b]. erase (arcsIt [b]);
}

 

void DiGraph::Arc::setNode (Node* newNode,
                            bool out)
{
  ASSERT (newNode);
  ASSERT (node [out]);
  ASSERT (node [out] -> graph == newNode->graph);
  
  if (node [out] == newNode)
    return;
  
  node [out] -> arcs [! out]. erase (arcsIt [out]);
  node [out] = newNode;
  newNode->arcs [! out]. push_back (this);
  arcsIt [out] = newNode->arcs [! out]. end ();
  arcsIt [out] --;
}




// DiGraph

void DiGraph::init (const DiGraph &other,
	                  Old2new &old2new)
{
  ASSERT (old2new. empty ());
  ASSERT (old2new. bucket_count () == other. nodes. size ());

  for (const Node* node : other. nodes)
  {
  	Node* n = node->copy ();
  	n->attach (*this);
  	old2new [node] = n;
  } 
  ASSERT (old2new. size () == other. nodes. size ());


  for (const Node* node : other. nodes)
  {
  	Node* n = old2new [node];
  	ASSERT (n);
  	
  	EXEC_ASSERT (n->parentDC = old2new [static_cast <Node*> (n->parentDC)]);
  	
  	if (n->scc)
  	{
	  	n->scc = old2new [n->scc];
	  	ASSERT (n->scc);
	  }
  	
	  for (const Arc* arc : node->arcs [true])
	  {
	  	Arc* a = arc->copy ();
	  	a->attach (n, old2new [arc->node [true]]);
	  }
  }   
}



void DiGraph::qc () const
{
  if (! qc_on)
    return;

  Set<const Node*> nodes_;
  Set<const Arc*> arcs_ [2];
  Set<string> names;
  size_t arcs = 0;
  for (const Node* node : nodes)
  {
    QC_ASSERT (node);
    QC_ASSERT (node->graph == this);
    nodes_ << node;
    try { node->qc (); }
      catch (const exception &e)
        { throw runtime_error ("In node " + node->getHumanName () + ": " + e. what ()); }
  	for (const bool b : {false, true})
      for (const Arc* arc : node->arcs [b])
      {
        QC_ASSERT (arc);
        arcs_ [b] << arc;
        arcs++;
        if (b)
        {
    	  	Unverbose unv;
          arc->qc ();
        }
      }
    if (names. contains (node->getName ()))
      throw runtime_error ("Duplicate name: " + node->getName ());
    names << node->getName ();
  }
  QC_ASSERT (nodes. size () == nodes_. size ());
	for (const bool b : {false, true})
    QC_ASSERT (2 * arcs_ [b]. size () == arcs);
  QC_ASSERT (arcs_ [false] == arcs_ [true]);
}

 

void DiGraph::saveText (ostream &os) const
{
  for (const Node* node : nodes)
  {
    node->saveText (os);
    os << endl;
  }
}

 

void DiGraph::deleteNodes ()
{
  while (! nodes. empty ())
  {
    Node* n = nodes. front ();
    if (! n)
      errorExit ("DiGraph::Node is nullptr");
    delete n;
  }
}



void DiGraph::connectedComponents ()
{
  for (Node* node : nodes)
	  node->DisjointCluster::init ();
  for (Node* node : nodes)
  	for (const bool b : {false, true})
      for (Arc* arc : node->arcs [b])
	    	arc->node [b] -> DisjointCluster::merge (* arc->node [! b]);
}



void DiGraph::scc ()
{
  size_t visitedNum = 0;
  vector<Node*> vec (nodes. size ());
  vec. clear ();
  stack<Node*,vector<Node*>> sccStack (vec);
  for (Node* node : nodes)
  {
    if (! node->orderDfs)
      node->setScc (visitedNum, sccStack);
    ASSERT (sccStack. empty ());
  }
}

 

void DiGraph::contractScc ()
{
  for (Iter <List<Node*> > iter (nodes); iter. next (); )
  {
    Node* n = *iter;
    ASSERT (n->scc);
    if (n->scc != n)
      n->scc->contract (n);
  }
}



VectorPtr<DiGraph::Node> DiGraph::getEnds (bool out) const
{
  VectorPtr<Node> s;
  for (const Node* node : nodes)
    if (node->arcs [out]. empty ())
      s << node;
  return s;
}



DiGraph::Node2Node DiGraph::reverse (const Node2Node& old2new)
{
  Node2Node new2old (old2new. size ());
  for (const auto& it : old2new)
    new2old [it. second] = it. first;    
  return new2old;
}



void DiGraph::borrowArcs (const Node2Node &other2this,
                          bool parallelAllowed)
{
#ifndef NDEBUG
  const DiGraph* otherGraph = nullptr;
#endif
  for (const auto& it : other2this)
  {
    const Node* other = it. first;
    const Node* from  = it. second;
    ASSERT (other);
    ASSERT (from);
    ASSERT (other->graph);
    ASSERT (from ->graph);
  #ifndef NDEBUG
    if (otherGraph)
      { ASSERT (otherGraph == other->graph); }
    else
      otherGraph = other->graph;
  #endif
    ASSERT (this == from->graph);
    ASSERT (otherGraph != this);
    const VectorPtr<Node> otherNeighborhood (other->getNeighborhood (true));
    for (const Node* otherNeighbor : otherNeighborhood)
      if (const Node* to = findPtr (other2this, otherNeighbor))
        if (parallelAllowed || ! from->isIncident (to, true))
          new Arc ( var_cast (from)
                  , var_cast (to)
                  );
  }
}




////////////////////////////////////////// Tree ////////////////////////////////////////////

// Tree::TreeNode

void Tree::TreeNode::qc () const
{
  if (! qc_on)
    return;
  DiGraph::Node::qc ();
    
  if (graph)
  {
    QC_ASSERT (! (isLeafType () && isInteriorType ()));
    QC_IMPLY (isLeaf (), isLeafType ());
    QC_IMPLY (isInteriorType () && getParent (), getParent () -> isInteriorType ());
  }
  if (contains (getName (), objNameSeparator))
    throw runtime_error (strQuote (getName ()) + " contains separator " + strQuote (string (1, objNameSeparator)));
}
  


void Tree::TreeNode::saveText (ostream &os) const
{
  os << getName () << ": ";
  saveContent (os);

  
  bool saveSubtreeP = false;
	for (const Arc* arc : arcs [false])
	  if (static_cast <const TreeNode*> (arc->node [false]) -> getSaveSubtreeP ())
	  {
	    saveSubtreeP = true;
	    break;
	  }

  if (saveSubtreeP)
  {
  	Offset ofs;
  	for (const Arc* arc : arcs [false])
  	{ 
  		Offset::newLn (os);
  	  static_cast <TreeNode*> (arc->node [false]) -> saveText (os);
  	}
  }
  else
    if (! isLeaf ())
      os << " (" << getLeavesSize () << ")";  // Time ??
}



void Tree::TreeNode::printNewick_ (ostream &os,
	                                 bool internalNames,
	                                 bool minimalLeafName) const
{
  // Cf. saveText() ??

	if (isLeaf ())
		os << name2newick (getNewickName (minimalLeafName));
	else
	{
	  string internalName;
		if (internalNames)
			internalName = name2newick (getNewickName (minimalLeafName));
	  ASSERT (! arcs [false]. empty ());
	  
	#if 0
	  bool done = false;
	  if (   internalName. empty () 
	      && arcs [false]. size () == 1
	     )
		{
		  const TreeNode* n = static_cast <TreeNode*> (arcs [false]. front () -> node [false]);
		  ASSERT (n);
    	const double dist = n->getParentDistance ();
    	if (dist == dist && ! dist)
    	{
    	  n->printNewick_ (os, internalNames, minimalLeafName);
    	  done = true;
    	}
		}			
		if (! done)
  #endif
		{	  
  		os << "(";
  		bool first = true;
  		for (const Arc* arc : arcs [false])
  		{
  			const TreeNode* n = static_cast <TreeNode*> (arc->node [false]);
  			if (! first)
  			  os << ",";
  			n->printNewick_ (os, internalNames, minimalLeafName);
  			first = false;
  		}	  
  		os << ")" << internalName;
    }
	}
	
	const double dist = getParentDistance ();
	if (dist == dist && dist != -1.0)
		os << ":" << fixed << dist;
}



#if 0
string Tree::TreeNode::name2newick (const string &s) 
{
  string s1 (s);
  replace (s1, "\"\' ():;,[]<>=", '_');
  return s1;
//return s1. substr (0, 100 /*50*/);  // PAR
    // Newick allowes the name length >= 50 characters, but http://www.trex.uqam.ca/ breaks with 50
}
#endif



const Tree::TreeNode* Tree::TreeNode::getAncestor (size_t height) const
{
  const TreeNode* n = this;
  FOR (size_t, i, height)
  	if (n)
      n = n->getParent ();
    else
      break;
  if (n)
    return n;
  return getTree(). root;  
}



void Tree::TreeNode::setParent (TreeNode* newParent)
{ 
	ASSERT (newParent != this);

#if 0	
  // Breaks sort()'ing
	if (newParent && getParent () == newParent)
	  return;
#endif
	
	if (! arcs [true]. empty ())
	{
		Arc* a = arcs [true]. front ();
		delete a;
  }
  ASSERT (arcs [true]. empty ());

	if (newParent)
	{
	  new Arc (this, newParent);
	  ASSERT (getParent () == newParent);
	}
	else
		var_cast (getTree ()). root = this;
}



void Tree::TreeNode::printAncestors (const TreeNode* end) const
{
  const TreeNode* node = this;
  for (;;)
  {
    cout << ' ' << node;
    if (! node || node == end)
      break;
    node = node->getParent ();
  }
  cout << endl;
}



Tree::TreeNode::TipName Tree::TreeNode::getTipName () const
{ 
  if (isLeaf ()) 
    return TipName (getName (), 0);
    
  TipName tn_best;
	for (const DiGraph::Arc* arc : arcs [false])
	{
	  const TipName tn = static_cast <const TreeNode*> (arc->node [false]) -> getTipName ();
	  if (tn_best. name. empty () || tn_best. name > tn. name)
	    tn_best = tn;
	}
	tn_best.depth ++;

	return tn_best;
}



void Tree::TreeNode::getSubtreeHeights (Vector<Tree::TreeNode::NodeDist> &nodeHeights) const
{
  if (isLeaf ())
    return;

  double height = 0.0;
	for (const DiGraph::Arc* arc : arcs [false])
	{
	  const TreeNode* node = static_cast <const TreeNode*> (arc->node [false]);
	  const double dist = node->getParentDistance ();
	  ASSERT (dist >= 0.0);
	  node->getSubtreeHeights (nodeHeights);
	  maximize (height, dist + (node->isLeaf () ? 0.0 : nodeHeights. back (). dist));
  }
  nodeHeights << NodeDist {this, height};
}



void Tree::TreeNode::getLeafDepths_ (Vector<Tree::TreeNode::NodeDist> &leafDepths,
                                     bool first) const
{
  const size_t start = leafDepths. size ();
  
  if (isLeaf ())
    leafDepths << NodeDist {this, 0.0};
  else
  	for (const DiGraph::Arc* arc : arcs [false])
  	{
  	  const TreeNode* node = static_cast <const TreeNode*> (arc->node [false]);
  	  node->getLeafDepths_ (leafDepths, false);
    }
    
  if (! first)
  {
    const double parentDistance = getParentDistance ();
    FFOR_START (size_t, i, start, leafDepths. size ())
      leafDepths [i]. dist += parentDistance;
  }
}



size_t Tree::TreeNode::getHeight () const
{
  size_t n = 0;
	for (const Arc* arc : arcs [false])
	  maximize (n, 1 + static_cast <const TreeNode*> (arc->node [false]) -> getHeight ());
	return n;
}



size_t Tree::TreeNode::getInteriorHeight () const
{
  ASSERT (isInteriorType ());
  
  size_t n = 0;
	for (const Arc* arc : arcs [false])
	{
	  const TreeNode* node = static_cast <const TreeNode*> (arc->node [false]);
	  if (node->isInteriorType ())
	    maximize (n, 1 + node->getInteriorHeight ());
	}
	return n;
}



double Tree::TreeNode::getDistanceHeight () const
{
  double n = 0.0;
	for (const Arc* arc : arcs [false])
	{
		const TreeNode* child = static_cast <const TreeNode*> (arc->node [false]);
		const double parentDistance = child->getParentDistance ();
		ASSERT (parentDistance >= 0.0);
	  maximize (n, parentDistance + child->getDistanceHeight ());
	}
	return n;
}



void Tree::TreeNode::getBifurcatingInteriorBranching (size_t &bifurcatingInteriorNodes,
                                                      size_t &branches) const
{
  ASSERT (isInteriorType ());
  ASSERT (/*a*/ arcs [false]. size () >= 2);  // transient nodes are prohibited

  // Make *this bifurcating
  branches                 += arcs [false]. size () - 2;  // b
  bifurcatingInteriorNodes += arcs [false]. size () - 2;  // n

  size_t interiorArcs = 0;  // i
	for (const Arc* arc : arcs [false])
	{
	  const TreeNode* node = static_cast <const TreeNode*> (arc->node [false]);
	  if (node->isInteriorType ())
	  {
	    interiorArcs++;
	    node->getBifurcatingInteriorBranching (bifurcatingInteriorNodes, branches);
	  }
	}  
	
	branches += interiorArcs;
	if (interiorArcs)
	  bifurcatingInteriorNodes++;
	ASSERT (branches >= bifurcatingInteriorNodes);

  ASSERT (branches <= 2 * bifurcatingInteriorNodes);  
/*Proof:
  a >= 2.
  i <= a.
  For each parent node:
  b = a - 2 + i.
  n = a - 2 + (bool) i.
  x := 2 n - b = a - 2 + 2 (bool) i - i = (a - i) + 2 ((bool) i - 1).
  if i > 0 then x = a - i >= 0.
  if i = 0 then x = a - 2 >= 0.
*/
}



size_t Tree::TreeNode::getSubtreeSize (bool countLeaves) const
{
	size_t n = 0;
	for (const Arc* arc : arcs [false])
	{
	  const TreeNode* child = static_cast <TreeNode*> (arc->node [false]);
    if (! countLeaves && child->arcs [false]. empty ())
      continue;
	  n += 1 + child->getSubtreeSize (countLeaves);
	}
	return n;
}



void Tree::TreeNode::subtreeSize2leaves () 
{
	leaves = 0;
	for (const Arc* arc : arcs [false])
	{
	  const TreeNode* child = static_cast <TreeNode*> (arc->node [false]);
	  var_cast (child) -> subtreeSize2leaves ();
	  leaves += 1 + child->leaves;
	}
}



double Tree::TreeNode::getSubtreeLength () const
{
	double len = 0.0;
	for (const Arc* arc : arcs [false])
	{
	  const TreeNode* node = static_cast <TreeNode*> (arc->node [false]);
	  len += node->getParentDistance () + node->getSubtreeLength ();
	}
	return len;
}



void Tree::TreeNode::setLeaves () 
{
	leaves = 0;
	for (const Arc* arc : arcs [false])
	{
	  TreeNode* child = static_cast <TreeNode*> (arc->node [false]);
	  child->setLeaves ();
	  ASSERT (child->leaves);
	  leaves += child->leaves;
	}
	if (! leaves)
	  leaves = 1;
}



void Tree::TreeNode::setLeaves (size_t leaves_arg) 
{
	leaves = leaves_arg;
	for (const Arc* arc : arcs [false])
	{
	  TreeNode* child = static_cast <TreeNode*> (arc->node [false]);
	  child->setLeaves (leaves_arg);
	}
}



size_t Tree::TreeNode::getLeavesSize () const
{
	size_t n = 0;
	for (const Arc* arc : arcs [false])
	  n += static_cast <TreeNode*> (arc->node [false]) -> getLeavesSize ();
	return max<size_t> (n, 1);
}



namespace
{
  struct ChildFreq
  {
    Tree::TreeNode* node;
    size_t freq;
    ChildFreq (Tree::TreeNode* node_arg,
               size_t freq_arg)
      : node (node_arg)
      , freq (freq_arg)
      { ASSERT (node);
        ASSERT (freq);
      }
    static bool strictlyLess (const ChildFreq& a,
                              const ChildFreq& b)
      { LESS_PART (b, a, freq);
        LESS_PART (a, b, node);
        return false;
      }
  };
}



void Tree::TreeNode::children2frequentChild (double rareProb)
{
  ASSERT (rareProb >= 0.0);
  ASSERT (rareProb < 0.5);  // => in a bifurcating node at least one child is frequentChild
  ASSERT (frequentChild);  
  
  Vector<ChildFreq> childFreqs;
  {  
    const VectorPtr<DiGraph::Node> children (getChildren ());
    childFreqs. reserve (children. size ());
    for (const DiGraph::Node* child : children)
    {
      const TreeNode* node = static_cast <const TreeNode*> (child);
      ASSERT (node->leaves);
      childFreqs << ChildFreq (var_cast (node), node->leaves); 
    }
  }
  Common_sp::sort (childFreqs, ChildFreq::strictlyLess);  // Bifurcatization
  
  size_t sum = 0;
#ifndef NDEBUG
  size_t freq_prev = numeric_limits<size_t>::max ();
#endif
  for (ChildFreq& cf : childFreqs)
  {
  #ifndef NDEBUG
    ASSERT (freq_prev >= cf. freq);
    freq_prev = cf. freq;
  #endif
    sum += cf. freq;
    ASSERT (sum);
    if ((double) cf. freq / (double) sum < rareProb)
      continue;
    cf. node->frequentChild = true;      
    cf. node->children2frequentChild (rareProb);    
  }
}
  


void Tree::TreeNode::getLeaves (VectorPtr<TreeNode> &leafVec) const
{
  if (arcs [false]. empty ())
    leafVec << this;
  else
  	for (const Arc* arc : arcs [false])
  	  static_cast <TreeNode*> (arc->node [false]) -> getLeaves (leafVec);
}



const Tree::TreeNode* Tree::TreeNode::getClosestLeaf (size_t &leafDepth) const
{
	const TreeNode* leaf = nullptr;
  leafDepth = SIZE_MAX;
	for (const Arc* arc : arcs [false])
	{
		size_t depth1;
		const TreeNode* leaf1 = static_cast <TreeNode*> (arc->node [false]) -> getClosestLeaf (depth1);
		if (minimize (leafDepth, depth1))
			leaf = leaf1;
	}
	
	if (leaf)
		leafDepth++;
	else
	{
		leaf = this;
		leafDepth = 0;
	}

	ASSERT (! leafDepth == (leaf == this));
		
	return leaf;
}



const Tree::TreeNode* Tree::TreeNode::getOtherChild (const TreeNode* child) const
{
  ASSERT (child);
  ASSERT (child->getParent () == this);

  const TreeNode* otherChild = nullptr;
	for (const Arc* arc : arcs [false])
	{
		const TreeNode* n = static_cast <TreeNode*> (arc->node [false]);
	  if (n != child)
	  {
	  	ASSERT (! otherChild);
	  	otherChild = n;
	  }
	}
  return otherChild;
}



const Tree::TreeNode* Tree::TreeNode::getDifferentChild (const TreeNode* child) const
{
  ASSERT (child);
  ASSERT (child->getParent () == this);

	for (const Arc* arc : arcs [false])
	{
		const TreeNode* n = static_cast <TreeNode*> (arc->node [false]);
	  if (n != child)
	    return n;
	}
  throw runtime_error ("getDifferentChild(): Transient parent");
}



const Tree::TreeNode* Tree::TreeNode::getFirstDecendant () const
{
  const TreeNode* n = this;
  while (! n->isLeaf ())
    n = static_cast <TreeNode*> (n->arcs [false]. front () -> node [false]);
  return n;
}



const Tree::TreeNode* Tree::TreeNode::getLastDecendant () const
{
  const TreeNode* n = this;
  while (! n->isLeaf ())
    n = static_cast <TreeNode*> (n->arcs [false]. back () -> node [false]);
  return n;
}



string Tree::TreeNode::getLcaName () const
{ 
  const TreeNode* left  = getFirstDecendant ();
	const TreeNode* right = getLastDecendant ();
	string name;
	if (left == right)
		name = left->getLeafName ();
	else
	  name = left->getLeafName () + objNameSeparator + right->getLeafName (); 
	return name;
}



void Tree::TreeNode::childrenUp ()
{
  const VectorPtr<DiGraph::Node> children (getChildren ());
	for (const DiGraph::Node* node : children)
	{	
		TreeNode* n = const_static_cast <TreeNode*> (node);
		n->setParent (var_cast (getParent ()));  
	}
	ASSERT (arcs [false]. empty ());
}



void Tree::TreeNode::detachChildrenUp ()
{ 
  childrenUp ();
	if (! arcs [true]. empty ())
	{ 
		Arc* a = arcs [true]. front ();
		delete a;
  }  
  detach ();
}



void Tree::TreeNode::deleteSubtree ()
{
	for (const DiGraph::Arc* arc : arcs [false])
		static_cast <TreeNode*> (arc->node [false]) -> deleteSubtree ();
	while (! arcs [false]. empty ())
	{
		TreeNode* n = static_cast <TreeNode*> (arcs [false]. front() -> node [false]);
		delete n;
	}
}



const Tree::TreeNode* Tree::TreeNode::makeRoot ()
{
	const TreeNode* root_old = getTree (). root;
	ASSERT (root_old);
	
	if (this == root_old)
	  return this;
	
	TreeNode* parent_new = nullptr;
	TreeNode* node = this;
	while (node)
	{
		TreeNode* parent_old = var_cast (node->getParent ());
	  node->setParent (parent_new);
	  parent_new = node;
	  node = parent_old;
	}
	ASSERT (arcs [true]. empty ());
	ASSERT (this == getTree (). root);
	ASSERT (root_old != getTree (). root);
	
	return root_old;
}



void Tree::TreeNode::getArea_ (uint radius,
                               const Tree::TreeNode* prev,
                               VectorPtr<Tree::TreeNode> &area,
                               VectorPtr<Tree::TreeNode> &boundary) const
{
  area << this;

  size_t degree = (size_t) (prev ? 1 : 0);  // Degree of *this in area
  if (radius)
  {
    const TreeNode* parent_ = getParent ();
    if (parent_ && parent_ != prev)
    {
      parent_->getArea_ (radius - 1, this, area, boundary);
      degree++;
    }
    for (const Arc* arc : arcs [false])
    {
      const TreeNode* child = static_cast <const TreeNode*> (arc->node [false]);
      if (child != prev)
      {
        child->getArea_ (radius - 1, this, area, boundary);
        degree++;
      }
    }
  }
  
  if (degree <= 1)
    boundary << this;
}



void Tree::TreeNode::getDistanceArea_ (double radius,
                                       const Tree::TreeNode* prev,
                                       VectorPtr<Tree::TreeNode> &area,
                                       VectorPtr<Tree::TreeNode> &boundary) const
{
  if (radius < 0.0)
    return;
  
  area << this;

  size_t degree = (size_t) (prev ? 1 : 0);  // Degree of *this in area

  const TreeNode* parent_ = getParent ();
  if (parent_ && parent_ != prev)
  {
    parent_->getDistanceArea_ (radius - getParentDistance (), this, area, boundary);
    degree++;
  }

  for (const Arc* arc : arcs [false])
  {
    const TreeNode* child = static_cast <const TreeNode*> (arc->node [false]);
    if (child != prev)
    {
      child->getDistanceArea_ (radius - child->getParentDistance (), this, area, boundary);
      degree++;
    }
  }
  
  if (degree <= 1)
    boundary << this;
}



void Tree::TreeNode::getClosestLeaves_ (const Tree::TreeNode* prev,
                                        double distance,
                                        size_t neighbors_max,
                                        Vector<NodeDist> &neighbors) const
{
  ASSERT (neighbors_max);
  
  if (   neighbors. size () == neighbors_max
      && neighbors. back (). dist <= distance
     )
    return;

  if (prev && isLeaf ())
  {
    neighbors << NodeDist {this, distance};
    neighbors. sort (NodeDist::distLess);
    if (neighbors. size () > neighbors_max)
      neighbors. pop ();
  }
  ASSERT (neighbors. size () <= neighbors_max);
      
  const TreeNode* parent_ = getParent ();
  if (parent_ && parent_ != prev)
    parent_->getClosestLeaves_ (this, distance + getParentDistance (), neighbors_max, neighbors);
  for (const Arc* arc : arcs [false])
  {
    const TreeNode* child = static_cast <const TreeNode*> (arc->node [false]);
    if (child != prev)
      child->getClosestLeaves_ (this, distance + child->getParentDistance (), neighbors_max, neighbors);
  }
}



void Tree::TreeNode::getSubtreeArea (const VectorPtr<Tree::TreeNode> &possibleBoundary,
	                                   VectorPtr<Tree::TreeNode> &area,
                                     VectorPtr<Tree::TreeNode> &boundary) const
{
	area << this;
	if (possibleBoundary. contains (this))
		boundary << this;
  else
  	if (arcs [false]. empty ())
  		boundary << this;
  	else
		  for (const Arc* arc : arcs [false])
		    static_cast <const TreeNode*> (arc->node [false]) -> getSubtreeArea (possibleBoundary, area, boundary);
}




// Tree

void Tree::qc () const
{
  if (! qc_on)
    return;
	DiGraph::qc ();

  Set<string> names;
  bool transient = false;
	for (const DiGraph::Node* node : nodes)
	{
	  QC_ASSERT (node->arcs [true]. size () <= 1);
	#if 0
	  {
	    cout << "Multiple parents of " << node->getHumanName () << endl;
	    const VectorPtr<DiGraph::Node> neighbors (node->getNeighborhood (true));
	    QC_ASSERT (neighbors. size () >= 2);
	    for (const DiGraph::Node* neighbor : neighbors)
	      cout << " " << neighbor->getHumanName ();
	    cout << endl << endl;
	    if (verbose ())
	      saveText (cout);
	    ERROR;
	  }
	#endif
	  const TreeNode* n = static_cast <const TreeNode*> (node);
	  if (n->isLeaf ())
	  {
	    const string newickName (TreeNode::name2newick (n->getNewickName (true)));
      if (names. contains (newickName))
        throw runtime_error ("Duplicate name: " + newickName);
      names << newickName;
    }
    if (n->isTransient ())
      transient = true;
	}
	QC_ASSERT (! root == nodes. empty ());
	QC_IMPLY (root, getRoot (true) == root);
	QC_IMPLY (root, (nodes. size () > 1) == ! root->isLeaf ());
	
	QC_IMPLY (! transient, nodes. size () <= 2 * root->getLeavesSize () - 1);
}



namespace 
{
  typedef  map<size_t, string> AsnFeatures;


  void printAsnFeatures (ostream &os,
                         const AsnFeatures &features)
  {
    if (features. empty ())
      return;
    os << ",\n\
      features {\n";
    bool first = true;
    for (const auto& it : features)
    {
      if (! first)
        os << ',';
      os << "\
        {\n\
          featureid " << it. first << ",\n\
          value \"" << it. second << "\"\n\
        }";
      first = false;
    }
    os << "\n\
      }";
  }
}



void Tree::printAsn (ostream &os) const
{
  unordered_map<const DiGraph::Node*, size_t/*index*/> node2index;  node2index. rehash (nodes. size ());
  size_t index = 0;
  for (const DiGraph::Node* n : nodes)
  {
    node2index [n] = index;
    index++;
  }    
  ASSERT (node2index. size () == nodes. size ());

  os << "BioTreeContainer ::= {fdict {\n\
    {\n\
      id 0,\n\
      name \"label\"\n\
    },\n\
    {\n\
      id 1,\n\
      name \"dist\"\n\
    },\n\
    {\n\
      id 2,\n\
      name \"info\"\n\
    }\n\
  },\n\
  nodes {";

  bool first = true;
  for (const DiGraph::Node* n : nodes)
  {
    const TreeNode* node = static_cast <const TreeNode*> (n);
    const TreeNode* parent = node->getParent ();
    if (! first)
      os << ',';
    os << "\n\
    {\n\
      id " << node2index [n];
    if (parent)
      os << ",\n\
      parent " << node2index [parent];
    string label;
    string info;
    if (node->isLeaf ())
    {
      label = n->getName ();
      const size_t pos = label. find (' ');
      if (pos != string::npos)
      {
        info = label. substr (pos + 1);
        label. erase (pos);
      }
    }
    AsnFeatures features;
    if (! label. empty ())
      features [0] = label;
    if (parent)
    {
      ostringstream oss;
      oss << node->getParentDistance ();
      features [1] = oss. str ();
    }
    if (! info. empty ())
      features [2] = info;
    printAsnFeatures (os, features);
    os << "\n\
    }";
    first = false;
  }
  os << "\n\
  }\n\
}\n\
";
}



void Tree::printArcLengths (ostream &os) const
{
  for (const DiGraph::Node* n : nodes)
  {
    const TreeNode* node = static_cast <const TreeNode*> (n);
    if (n == root)
      continue;
  	const double dist = node->getParentDistance ();
  	if (dist == dist && dist != -1.0 && dist > 0.0)
  	{
  		os << node->getLcaName () << " " << dist << " " << node->getRootDistance ();
  		double distPar = 0.0;
      if (const TreeNode* parent = node->getParent ())
        if (parent != root)
        {
        	distPar = parent->getParentDistance ();
        	if (distPar == distPar && distPar != -1.0)
        	{
        	  ASSERT (distPar > 0.0);
        	  os << " " << log (distPar) - log (dist);
        	}
        }
      if (! distPar)
     	  os << " ?";
      os << endl;
    }
  }  
}



double Tree::getAveArcLength () const
{
  double len = 0.0;
  size_t n = 0;
  for (const DiGraph::Node* node : nodes)
  {
    if (node == root)
      continue;
    const double arcLen = static_cast <const TreeNode*> (node) -> getParentDistance ();
    ASSERT (arcLen >= 0.0);
    len += arcLen;
    n++;
  }
  return len / (double) n;
}



Tree::Patristic::Patristic (const TreeNode* leaf1_arg, 
                            const TreeNode* leaf2_arg,
                            double distance_arg)
: leaf1 (leaf1_arg)
, leaf2 (leaf2_arg)
, distance (distance_arg)
{
  ASSERT (leaf1);
  ASSERT (leaf2);
  ASSERT (distance == distance);  // != NaN
  ASSERT (leaf1->graph);
  ASSERT (leaf1->graph == leaf2->graph);
  ASSERT (leaf1->getName () != leaf2->getName ());
  if (leaf1->getName () > leaf2->getName ())
    swap (leaf1, leaf2);
}



namespace 
{

typedef  map <const Tree::TreeNode*, double>  Leaf2dist;
 


Vector<Tree::Patristic> node2leafDistances (const Tree::TreeNode* node,
                                            Leaf2dist &leaf2dist) 
// Output: leaf2dist
{
  ASSERT (node);
  ASSERT (leaf2dist. empty ());
  
  Vector<Tree::Patristic> res;
  if (node->isLeaf ())
    leaf2dist [node] = 0.0;
  else
		for (const Tree::Arc* arc : node->arcs [false])
		{
			const Tree::TreeNode* n = static_cast <Tree::TreeNode*> (arc->node [false]);
			Leaf2dist nodeLeaf2dist;
			res << std::move (node2leafDistances (n, nodeLeaf2dist));
			ASSERT (! nodeLeaf2dist. empty ());
			const double dist = n->getParentDistance ();
			ASSERT (dist == dist);  // != NaN
			for (auto& it : nodeLeaf2dist)
			  it. second += dist;
		  for (const auto& it1 : leaf2dist)
			  for (const auto& it2 : nodeLeaf2dist)
			    res << Tree::Patristic (it1. first, it2. first, it1. second + it2. second);
			for (const auto& it : nodeLeaf2dist)
			  leaf2dist [it. first] = it. second;
	  }
    
  return res;
}

}



Vector<Tree::Patristic> Tree::getLeafDistances () const
{
  Leaf2dist leaf2dist;
  return node2leafDistances (root, leaf2dist);
}



size_t Tree::countInteriorNodes () const
{ 
  size_t n = 0;
  for (const DiGraph::Node* node : nodes)
    if (static_cast <const TreeNode*> (node) -> isInteriorType ())
      n++;
  return n;
} 



double Tree::getBifurcatingInteriorBranching () const
{ 
  if (! (root && root->isInteriorType ()))
    return -1;
    
  size_t bifurcatingInteriorNodes = 0;
  size_t branches = 0;
  root->getBifurcatingInteriorBranching (bifurcatingInteriorNodes, branches);
  const double branching = (double) branches / (double) bifurcatingInteriorNodes;
  ASSERT (branching >= 1.0);
  ASSERT (branching <= 2.0);
  return branching;
}



size_t Tree::countInteriorUndirectedArcs () const
{ 
  size_t n = countInteriorNodes ();
  if (root->isInteriorType ())
    n--;

  const VectorPtr<DiGraph::Node> children (root->getChildren ());
  switch (children. size ())
  {
    case 1: if (static_cast <const TreeNode*> (children [0]) -> isInteriorType ())
              n--;
            break;
    case 2: // root is transient in the undirected tree
      {
        size_t m = 0;
        for (const DiGraph::Node* child : children)
          if (static_cast <const TreeNode*> (child) -> isInteriorType ())
            m++;
        if (m)
          n--;
      }
  }

  return  n;
}



namespace 
{

bool getParentsOrTarget (const Tree::TreeNode* from,
			                   const Tree::TreeNode* target,
			                   VectorPtr<Tree::TreeNode> &parents) 
// Return: true <=> from->descendantOf(target)
{
	ASSERT (from)
	ASSERT (target)
	ASSERT (from->graph);
  ASSERT (from->graph == target->graph);
  ASSERT (parents. empty ());
  
	while (from)
	{
		if (from == target)
			return true;
		parents << from;
		from = from->getParent ();
	}
	
	return false;
}

}



const Tree::TreeNode* Tree::getLca (const TreeNode* n1,
	                                  const TreeNode* n2,
                                    LcaBuffer &buf) 
{
  IMPLY (n1 && n2, n1->graph && n1->graph == n2->graph);
  
  if (   ! n1 
  	  || ! n2
  	 )
  	return nullptr;
  	
  if (n1 == n2)
    return n1;
	
  buf. clear ();
  auto& vec1 = buf. vec1;
  auto& vec2 = buf. vec2;

	if (getParentsOrTarget (n1, n2, vec1))
		return n2;
	ASSERT (! vec1. empty ());

	if (getParentsOrTarget (n2, n1, vec2))
		return n1;
	ASSERT (! vec2. empty ());
	
	const TreeNode* m = nullptr;
	size_t i1 = vec1. size () - 1;
	size_t i2 = vec2. size () - 1;
	ASSERT (vec1 [i1] == vec2 [i2]);
	while (vec1 [i1] == vec2 [i2])
	{
		m = vec1 [i1];
		ASSERT (i1);
		ASSERT (i2);
		i1--;
		i2--;
	}
	ASSERT (m);
	return m;
}



const Tree::TreeNode* Tree::getLca (const VectorPtr<TreeNode> &nodeVec,   
	                                  Tree::LcaBuffer &buf) 
{
  if (nodeVec. empty ())
    return nullptr;
    
	const TreeNode* n = nodeVec [0];
	FFOR_START (size_t, i, 1, nodeVec. size ())
    n = getLca (n, nodeVec [i], buf);
	return n;
}



Set<const Tree::TreeNode*> Tree::getParents (const VectorPtr<TreeNode> &nodeVec,   
	                                           Tree::LcaBuffer &buf) 
{
	Set<const TreeNode*> s;  
  const TreeNode* lca = getLca (nodeVec, buf);
  for (const TreeNode* n : nodeVec)
	  while (n != lca)
	  {
	  	ASSERT (n);
		  s << n;
		  n = n->getParent ();
		}

	return s;
}



VectorPtr<Tree::TreeNode>& Tree::getPath (const TreeNode* n1,
											                    const TreeNode* n2,
											                    const TreeNode* ca,
												                  const TreeNode* &lca,
											                    LcaBuffer &buf)
{ 
  ASSERT (n1);
  ASSERT (n2);
  ASSERT (n1->graph);
  ASSERT (n1->graph == n2->graph);
  IMPLY (ca, n1->graph == ca->graph);
  
  buf. clear ();
  auto& vec1 = buf. vec1;
  auto& vec2 = buf. vec2;

  if (n1 == n2)
  {
    lca = n1;
    return vec1;
  }
  
  const TreeNode* n1_init = n1;
  
	while (n1 != ca && n1 != n2)
	{
	  ASSERT (n1);
		vec1 << n1;
		n1 = n1->getParent ();
	}
	if (n1 == n2)
	{
    lca = n2;
		return vec1;
  }

  
	while (n2 != ca && n2 != n1_init)
	{
	  ASSERT (n2);
		vec2 << n2;
		n2 = n2->getParent ();
	}
	if (n2 == n1_init)
	{
    lca = n1_init;
		return vec2;
  }

	ASSERT (! vec1. empty ());
	ASSERT (! vec2. empty ());

  lca = ca;  
	size_t i1 = vec1. size () - 1;
	size_t i2 = vec2. size () - 1;
	while (vec1 [i1] == vec2 [i2])
	{
	  lca = vec1 [i1];
    vec1. pop_back ();
    vec2. pop_back ();
		ASSERT (i1);
		ASSERT (i2);
		i1--;
		i2--;
	}
	if (vec1. size () < vec2. size ())
	{
	  CONST_ITER_REV (VectorPtr<TreeNode>, it, vec1)
	    vec2 << *it;
	  return vec2;
  }
	else
	{
	  CONST_ITER_REV (VectorPtr<TreeNode>, it, vec2)
	    vec1 << *it;
	  return vec1;
  }
}



VectorPtr<Tree::TreeNode> Tree::leaves2lcas (const Set<const TreeNode*> &leaves) const
{
  ASSERT (root);
  
  VectorPtr<TreeNode> lcas;
  
  if (leaves. empty ())
    return lcas;
    
  if (qc_on)
    for (const TreeNode* leaf : leaves)
    {
      QC_ASSERT (leaf);
      QC_ASSERT (leaf->isLeaf ());
      QC_ASSERT (& leaf->getTree () == this);
    }

  // TreeNode::leaves
  var_cast (root) -> setLeaves ();

  Vector<Set<const TreeNode*>> leaves2nodes;  
    // Index = TreeNode::leaves - 1
  leaves2nodes. resize (root->leaves);
  for (const Node* node : nodes)
  {
    const TreeNode* treeNode = static_cast <const TreeNode*> (node);
    ASSERT (treeNode);
    ASSERT (treeNode->leaves);
    leaves2nodes [treeNode->leaves - 1] << treeNode;
  }    
    
  
  // leaves2nodes: erase()
  // Bottom-up
  for (Set<const TreeNode*>& treeNodes : leaves2nodes)
  {
    VectorPtr<TreeNode> bads;  bads. reserve (treeNodes. size ());
    for (const TreeNode* parent : treeNodes)
    {
      ASSERT (parent);
      
      bool good = true;
      const VectorPtr<Node> children (parent->getChildren ());
        // Current lcas of subtrees
      if (children. empty ())
      {
        ASSERT (parent->leaves == 1);
        ASSERT (parent->isLeaf ());
        if (! leaves. contains (parent))
          good = false;
      }
      else
        for (const Node* child_ : children)
        {
          const TreeNode* child = static_cast <const TreeNode*> (child_);
          ASSERT (child);
          ASSERT (child->leaves < parent->leaves);
          ASSERT (child->leaves);
          if (! leaves2nodes [child->leaves - 1]. contains (child))
          {
            good = false;
            break;
          }
        }
        
      if (good)
        for (const Node* child_ : children)
        {
          const TreeNode* child = static_cast <const TreeNode*> (child_);
          EXEC_ASSERT (leaves2nodes [child->leaves - 1]. erase (child) == 1);
        }
      else
        bads << parent;
    }
    for (const TreeNode* bad : bads)
      EXEC_ASSERT (treeNodes. erase (bad) == 1);
  }
  
  
  // lcas
  ASSERT (lcas. empty ());
  for (const Set<const TreeNode*>& treeNodes : leaves2nodes)
    for (const TreeNode* node : treeNodes)
      lcas << node;
  ASSERT (! lcas. empty ());
  
  return lcas;
}



void Tree::setFrequentChild (double rareProb)
{ 
  if (! root)
    return;
    
  setLeaves ();
    
  for (DiGraph::Node* node : nodes)
    static_cast <TreeNode*> (node) -> frequentChild = false;
  var_cast (root) -> frequentChild = true;
  var_cast (root) -> children2frequentChild (rareProb);
  
  VectorPtr<TreeNode> transients;  transients. reserve (nodes. size ());
  for (const DiGraph::Node* node : nodes)
  {
    const TreeNode* tn = static_cast <const TreeNode*> (node);
    if (! tn->frequentChild)
      continue;
    const VectorPtr<DiGraph::Node> children (tn->getChildren ());
    size_t freqChildren = 0;
    for (const DiGraph::Node* child : children)
      if (static_cast <const TreeNode*> (child) -> frequentChild)
        freqChildren++;
    if (freqChildren == 1)  // transient, not leaf
      transients << tn;
  }  
  for (const TreeNode* tn : transients)
    var_cast (tn) -> frequentChild = false;
}



void Tree::setFrequentDegree (double rareProb)
{ 
  ASSERT (rareProb >= 0.0);
  ASSERT (rareProb < 0.3);
  
  if (! root)
    return;
    
  setLeaves ();
    
  const size_t allLeaves = root->leaves;
#if 0
  Binomial bin;
  bin. setParam ((int) allLeaves, rareProb);
  const double pValue = 0.01;  // PAR
#endif
  for (DiGraph::Node* node : nodes)
  {
    size_t degree = 0;
    {
      size_t sum = 0;
      const VectorPtr<DiGraph::Node> children (node->getChildren ());
      for (const DiGraph::Node* child : children)
      {
        const TreeNode* tn = static_cast <const TreeNode*> (child);
        const size_t leaves = tn->leaves /*getLeavesSize ()*/; 
        ASSERT (leaves);
        if ((double) leaves / (double) allLeaves >= rareProb)
      //if (1 - bin. cdf ((int) leaves - 1) <= pValue)
          degree++;
        sum += leaves;
      }
      if (node != root)
      {
        ASSERT (allLeaves > sum);
        const size_t leaves = allLeaves - sum;
        if ((double) leaves / (double) allLeaves >= rareProb)
      //if (1 - bin. cdf ((int) leaves - 1) <= pValue)
          degree++;
      }
    }
    static_cast <TreeNode*> (node) -> frequentDegree = degree;
  }

  // Leaves
  for (const DiGraph::Node* node : nodes)
  {
    const TreeNode* tn = static_cast <const TreeNode*> (node);
    if (! tn->isLeafType ())
      continue;
    ASSERT (tn->frequentDegree == 1);
    const TreeNode* parent = tn;
    while (parent && parent->frequentDegree == 1)
      parent = parent->getParent ();
    if (parent && parent->frequentDegree < 3)  // 0 or 2
      var_cast (tn) -> frequentDegree = 0;
  }  
}



void Tree::setRoot ()
{
  root = nullptr;
  for (const DiGraph::Node* node : nodes)
    if (! static_cast <const TreeNode*> (node) -> getParent ())
    {
      ASSERT (! root);
      const_static_cast <TreeNode*> (node) -> setParent (nullptr);
    }
  IMPLY (! nodes. empty (), root);
}



size_t Tree::deleteTransients ()
{
	size_t n = 0;
 	for (List<DiGraph::Node*>::const_iterator it = nodes. begin (); 
 		   it != nodes. end ();
 		  )
 	{
 		TreeNode* node = static_cast <TreeNode*> (*it);
 		it++;
 		if (node->deleteTransient ())
      n++;
  }  
  return n;
}



size_t Tree::restrictLeaves (const StringVector &leafNames,
                             bool deleteTransientAncestor)
{
  size_t n = 0;
 	Vector<DiGraph::Node*> nodeVec;  nodeVec. reserve (nodes. size ());
 	insertAll (nodeVec, nodes);
 	for (DiGraph::Node* node_ : nodeVec)  
 	{
 	  const TreeNode* node = static_cast <const TreeNode*> (node_);
    if (node->isLeafType ())
      if (! leafNames. containsFast (node->getName ()))
      {
        deleteLeaf (var_cast (node), deleteTransientAncestor);
        n++;
      }
  }
  return n;
}



bool Tree::strictlyLess_std (const DiGraph::Node* a,
	                           const DiGraph::Node* b)
{
	ASSERT (a);
	ASSERT (b);
	ASSERT (a->graph);
	ASSERT (a->graph == b->graph);
	const TreeNode* a_ = static_cast <const TreeNode*> (a);
	const TreeNode* b_ = static_cast <const TreeNode*> (b);

  ASSERT (a_->leaves);
  ASSERT (b_->leaves);
	LESS_PART (*b_, *a_, leaves);
  LESS_PART (*a_, *b_, getFirstDecendant  () -> getName ());

  return false;
}




// TopologicalSort

TopologicalSort::TopologicalSort (DiGraph &graph_arg,
                                  bool out_arg)
: graph (graph_arg)
, out (out_arg)
{ 
	order. reserve (graph. nodes. size ());
	for (DiGraph::Node* node : graph. nodes)
	{
	  for (Iter<List<DiGraph::Arc*>> iter (node->arcs [! out]); iter. next (); )
	    if ((*iter) -> selfLoop ())
	      delete *iter;  		
    if (node->arcs [! out]. empty ())
    	order << node;
  }
}



const DiGraph::Node* TopologicalSort::getFront ()
{	
	ASSERT (index <= order. size ());	
	if (index == order. size ())
		return nullptr;
		
	const DiGraph::Node* n = order [index];
	ASSERT (n);
	ASSERT (n->arcs [! out]. empty ());
	VectorPtr<DiGraph::Node> add;
  for (const DiGraph::Arc* arc : n->arcs [out])
  {
  	const DiGraph::Node* next = arc->node [out];
  	ASSERT (next);
  	ASSERT (! next->arcs [! out]. empty ());
  	if (! next->isIncidentExcept (n, ! out))
  		if (! add. contains (next))
  		  add << next;
  }
  order << add;
  add. clear ();
	var_cast (n) -> isolate ();
	index++;
	
	return n;
}


}


