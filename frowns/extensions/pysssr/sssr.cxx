// sssr source code and helpers
#include "sssr.h"
#include <queue>
#include "assert.h"

//#define DEBUG

// Courtesy functions for debugging
void dumpAtoms(AtomList &x) {
	AtomList::iterator ix=x.begin();
	for(; ix != x.end(); ix++) {
		(*ix)->dump();
	}
}

void dumpAtoms(AtomVector &x) {
	AtomVector::iterator ix=x.begin();
	for(; ix != x.end(); ix++) {
		(*ix)->dump();
		
	}
}

Ring getRing(Atom *atom, AtomVector &atoms);

RingSet& ringDetection(Molecule *mol, int RSIZE) {
	AtomVector &fullSet = mol->atoms;
	AtomList trimset;
	AtomList temp1;
	AtomList temp2;
	AtomList NodesN2, startNodes;
	RingSet *ringSet = new RingSet(mol);
	

	AtomList::iterator atomiter;
	AtomVector::iterator atomviter;

	for (atomviter = fullSet.begin(); atomviter != fullSet.end(); atomviter++) {
		temp1.push_back(*atomviter);
#ifdef DEBUG
		(*atomviter)->dump();
#endif
	}

	AtomList &current = temp1;
	AtomList &next = temp2;
	int minimumDegree = 10000;
	Atom *init = NULL;

	int index = 0;
	int sizeBeforeLast = current.size();
	int lastsize = current.size();

	while (current.size()) {
		index++;
		next.clear();
		NodesN2.clear();
		startNodes.clear();

		// remove all the degree 0 nodes from the trimset
		minimumDegree = 10000;
		init = NULL;

		for (atomiter = current.begin(); atomiter != current.end(); atomiter++) {
			Atom *atom = *atomiter;
			int degree = atom->degree();
			// store the atom of minimum degree
			if (degree && degree < minimumDegree) {
				minimumDegree = atom->degree();
				init = atom;
			}

			if (degree == 0 ) {
				trimset.push_back(atom);
			} else if(degree == 2) {
				NodesN2.push_back(atom);
				next.push_back(atom);
			} else {
				next.push_back(atom);
			}
		}

		if (!init) {
#ifdef DEBUG
			std::cout << "warning init is NULL" << std::endl;
#endif
			break;
		}

#ifdef DEBUG
		std::cout << "Trial " << index << " minimumDegree " << minimumDegree << std::endl;
		std::cout << "\tcurrent size " << current.size() << std::endl;
		std::cout << "\ttrimset size " << trimset.size() << std::endl;
		std::cout << "\tMinimum degree is " << minimumDegree << std::endl;
#endif

		if (minimumDegree == 1) {
#ifdef DEBUG
			std::cout << "\tremoving mindegree 1 atom ";
			init->dump();
			std::cout << std::endl << "\t";
#endif
			init->trim(fullSet);
#ifdef DEBUG
			init->dump();
			std::cout << std::endl;
#endif

			trimset.push_back(init);
			next.remove(init);
#ifdef DEBUG
			std::cout << "\tnnext size " << next.size() << std::endl;
#endif

		} else if (minimumDegree == 2) {
			for (atomiter = NodesN2.begin(); atomiter != NodesN2.end(); atomiter++) {
				Ring ring = getRing((*atomiter), fullSet);

				if (ring.count()) {
#ifdef DEBUG
				  std::cout << "\t====>found ring of size " << ring.count() << std::endl;
#endif
					int start = ringSet->add(ring);

					if(start) {
						startNodes.push_back((*atomiter));
					}
#ifdef DEBUG
					else {
					  std::cout << "\t ===> ring already found\n";
					}
#endif
				}
			}

			if (startNodes.size() == 0) throw SSSRError();
			//			assert(startNodes.size() > 0);
			for (atomiter = startNodes.begin(); atomiter != startNodes.end(); atomiter++) {
				(*atomiter)->trim(fullSet);
			}

		} else if ( minimumDegree >= 3 ) {

			Ring ring = getRing(init, fullSet);
			if (ring.count() > 0) {
				ringSet->add(ring);
				checkEdges(init, ring, fullSet);
			} else {
				next.remove(init);
			}
		}

		// swap the fullsets
#ifdef DEBUG
		std::cout << "\tnext size " << next.size() << std::endl;
#endif
		if (&current == &temp1) {
		  current = temp2;
		  next = temp1;
		} else {
		  current = temp1;
		  next = temp2;
		}

		if (lastsize == current.size()) {
		  // Error, we didn't remove an atom!!!!
		  // this will cause an infinite loop
#ifdef DEBUG	      
		  std::cout << "******************************************\nFailed to reduce ringsize\n";
#endif
		}

		if (sizeBeforeLast == current.size() && index > 3) {
#ifdef DEBUG	      
		  std::cout << "******************************************\nFailed to reduce ringsize in two iterations, breaking!\n";
		  //		  break;
#endif
		}
		sizeBeforeLast = lastsize;
		lastsize = current.size();	
	}

#ifdef DEBUG
	std::cout << "Found " << ringSet->size() << " rings" << std::endl;
	ringSet->dump();


	if (RSIZE != -1 && ringSet->size() != RSIZE) {
		std::cout << "Failed:: Got " << ringSet->size() << " wanted " << RSIZE << std::endl;
		return *ringSet;
	}
#endif
	// this is a test of memory management of course
	//delete &(ringSet.to_atom_lists());

	return *ringSet;
}

class Bond{
	// This helps us keep track of bonds in the BFS
public:
	int source;   // Where we are coming from
	int next;     // Where we are going to

	Bond(int s, int n) {
		source=s; 
		next=n;
	}
};


typedef std::queue<Bond* > BondQue;

Ring getRing(Atom *atom, AtomVector &atoms) {
	Ring ringSet;
#ifdef DEBUG 
	std::cout << "Getting a ring from atom " << atom->index << std::endl;
#endif
	// I'm assuming that a bitset ends up initialized to empty....
	Ring *paths = new Ring[atoms.size()];

	BondQue queue;
	IntList::iterator iiter;

	// load the queue up with the initial breadth first elements
	int source = atom->index;
	int node = -1;

	int qsize = 0;
	for(iiter = atom->neighbors.begin(); iiter != atom->neighbors.end(); iiter++) {
		node = *iiter;
		Bond *bond = new Bond(source, node);
		queue.push(bond);
		qsize++;

		paths[node][source] = 1;
		paths[node][node] = 1;
	}

	
	while(!queue.empty()) {
		Bond *b = queue.front();
		source = b->source;
		node = b->next;
		delete b;  // Don't need the bond anymore
		queue.pop();
		qsize--;


		// This is the bfs traversal here
		Atom &n = *(atoms[node]);
		Ring &path = paths[node];

#ifdef DEBUG 
		std::cout << "queue size is " << qsize << std::endl;	
		std::cout << source << "->" << node << std::endl;
		std::cout << "path[" << node << "] is " << path << std::endl;
		// advance through the neighbors
		std::cout << "advancing to ";
		for(iiter = n.neighbors.begin(); iiter != n.neighbors.end(); iiter++) {
			std::cout << *iiter << ", ";
		}

		std::cout << std::endl;
#endif
		for(iiter = n.neighbors.begin(); iiter != n.neighbors.end(); iiter++) {
			Atom &m = *(atoms[*iiter]);
			

			if(m.index != source) {
#ifdef DEBUG
				std::cout << " traversing to atom " << *iiter << std::endl;
#endif
				Ring &mpath = paths[m.index];
				//{ Collision occurs then attached node m's path is
		                //  no longer empty }
				if (paths[m.index].count()) {
					Ring intersection = path & mpath;
#ifdef DEBUG
					std::cout << path.count() << " mpath size " << mpath.count() << " intersection size " << \
					  intersection.count() << "\n";

					std::cout << "\tcollision!" << std::endl;
					std::cout << "\t path=" << path << std::endl;
					std::cout << "\tmpath=" << mpath << std::endl;
					std::cout << "\t   &  " << intersection << std::endl;
#endif
					int intersectionSize = intersection.count();

					if (intersectionSize == 1) {
						Ring result = path | mpath;
						result[m.index] = 1;
#ifdef DEBUG
						std::cout << "==>returning result " << result << std::endl;
#endif
						// Memory management
						delete [] paths;
						while (!queue.empty()) {
							b = queue.front();
							delete b;
							queue.pop();
						}
						return result;
					}
			
				} else {
					// {Update the path m}
					paths[m.index] = path | mpath;
					paths[m.index][m.index] = 1;
#ifdef DEBUG
					std::cout << "\tupdating path " << m.index << std::endl;
					std::cout << "\t" << paths[m.index] << std::endl;	
					std::cout << "\tadding " << m.index << "->" << node << " to queue" << std::endl;
#endif
					Bond *addbond = new Bond(node, m.index);
					queue.push(addbond);					
					qsize++;
				}
			}
		}
	}

	// If we got here all bonds have been consumed;
	// free up memory
	delete[] paths;
	return ringSet;
}


void checkEdges(Atom *init, Ring &ring, AtomVector &fullSet) {
	std::list<Bond *> edges;

	std::list<int > atom_ids;
	std::list<int >::iterator iiter, niter;

	if (ring == 0) throw SSSRError();
	//	assert(ring != 0);
	// Get all the atoms in the ring
	for(unsigned int i=0; i < fullSet.size(); i++) {
	  if (ring[i]) {
	    atom_ids.push_back(i);
	  }
	}

	// Get all the edges in the ring by looking for two atoms
	// that are connected.
	// Note that this assumes the ring is the smallest ring and
	// doesn't have smaller branches between atoms
	for (iiter = atom_ids.begin(); iiter != atom_ids.end(); iiter++) {
		Atom &atom = (*fullSet[*iiter]);
		int &from_atom = *iiter;
		for (niter = atom.neighbors.begin(); niter != atom.neighbors.end(); niter++) {
			int &to_atom = *niter;

			if (ring[to_atom] && from_atom < to_atom) {
				Bond *b = new Bond(from_atom, to_atom);
				edges.push_back(b);
			}
		}
	}

	if(edges.size() == 0) throw SSSRError();
				//	assert(edges.size());
	// Cycle through all edges in the ring and break them!
	int r1, r2, temp;
	int min=0;
	Bond *minEdge = NULL;
	
	while (edges.size()) {
		Bond *edge = edges.front();
		edges.pop_front();
		Atom &atom1 = (*fullSet[edge->source]);
		Atom &atom2 = (*fullSet[edge->next]);

		atom1.break_bond(atom2);
		r1 = getRing(&atom1, fullSet).count();
		r2 = getRing(&atom2, fullSet).count();

		if (r2 > r1) {
			temp = r2;
		} else {
			temp = r1;
		}

		if( min == 0) {
			min = temp;
			minEdge = edge;
		} else if (temp < min) {
			delete minEdge;
			minEdge = edge;
		}

		atom1.add_bond(atom2);
	}

	
	Atom &atom1 = (*fullSet[minEdge->source]);
	Atom &atom2 = (*fullSet[minEdge->next]);

	atom1.break_bond(atom2);
	delete minEdge;
	
}


