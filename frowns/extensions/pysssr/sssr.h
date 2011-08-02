# ifndef __SSSR_BPK__
# define __SSSR_BPK__

#include <iostream>
#include <bitset>
#include <list>
#include <vector>
#include "assert.h"

#define MAXATOMS 512
typedef std::bitset<MAXATOMS > Ring; // we'll be in trouble if we have more than 500 bits...
typedef std::list<int > IntList;

class SSSRError {
};

class Atom {
	// An "atom": this holds the neighbor indices in the molecule
	// and has the degree as well
public:
	int index;
	IntList neighbors;
	IntList mark;

	Atom(int i) {
		index = i;
	}

	void finalize() {
		// finalize the atom (by keeping a copy of the neighbors list
		// the neighbors list will be modified at some point
		mark = neighbors;
	}

	void trim(std::vector<Atom *> atoms) {
		std::list<int >::iterator i;
		for (i=neighbors.begin(); i != neighbors.end(); i++) {
			atoms[*i]->neighbors.remove(index);
		}

		neighbors.clear();
	}

	void break_bond(Atom &other) {
		other.neighbors.remove(index);
		neighbors.remove(other.index);
	}

	void add_bond(Atom &other) {
		neighbors.push_back(other.index);
		other.neighbors.push_back(index);
	}

	inline int degree() {
		return neighbors.size();
	}

	void dump() {
		std::cout << "Atom(index=" << index << ", neighbors=[";
		std::list<int >::iterator i;
		for (i=neighbors.begin(); i != neighbors.end(); i++) {
			std::cout << *i << ", ";

		}
		std::cout << "])" << std::endl;
	}
};

typedef std::list<Atom* > AtomList;
typedef std::vector<Atom* > AtomVector;

class Molecule {
	// A molecule is a collection of atoms
	// the atoms know which atoms are linked to each other
public:
	AtomVector atoms;
	int numberOfAtoms;
	Molecule(int numAtoms) {
		numberOfAtoms = numAtoms;

		if (numAtoms > MAXATOMS) {
			// throw an error here, don't know what yet
		}

		atoms = AtomVector(numAtoms);
	
		for(unsigned int i = 0; i < atoms.size(); i++) {
			atoms[i] = new Atom(i);
		}
	}

	~Molecule() {
		// Memory management
		for(unsigned int i = 0; i < atoms.size(); i++ ) {
			delete atoms[i];
		}
	}

	void finalize() {
		// finalize the atoms (this keeps their original structure
		// intact)
		for (unsigned int i = 0; i< atoms.size(); i++) {
			atoms[i]->finalize();
		}
	}

	void add_bond(int atom1, int atom2) {
		// connect two atoms
		atoms[atom1]->add_bond(*(atoms[atom2]));
	}

};

class RingSet {
public:
	std::list<Ring > rings;
	Molecule *mol;

	RingSet(Molecule *molecule) {
		mol = molecule;
	}

	int size() {
		return rings.size();
	}

	int add(Ring r) {
	  if (r == 0) {
	    throw SSSRError();
	  }
	  //		assert( r != 0 ); // no empty rings!
		std::list<Ring >::iterator i;
		for (i = rings.begin(); i != rings.end(); i++) {
			if(r == *i) {
				return 0;
			}
		}

		rings.push_back(r);
		return 1;
	}

	void dump() {
		std::list<Ring >::iterator i;
		for (i = rings.begin(); i != rings.end(); i++) {
			std::cout << "\t" << *i << std::endl;
		}
	}

	std::list<IntList* > &to_atom_lists() {
		// convert the rings sets into a list of atom ids
		std::list<IntList* > *atomList = new std::list<IntList* >;
		
		std::list<Ring >::iterator i;
		int count=0;
		for (i = rings.begin(); i != rings.end(); i++, count++) {
			std::cout << "Ring " << count << std::endl;
			Ring ring = *i;
			IntList *list = new IntList;
			
			std::cout << "\t";
			for(int index=0; index < mol->numberOfAtoms; index++) {
				if(ring[index]) {
					list->push_back(index);
					std::cout << index << " ";
				}

			}
			std::cout << std::endl;
			atomList->push_back(list);
		}

		
		return *atomList;
	}


};
void dumpAtoms(AtomList &x);
void dumpAtoms(AtomVector &x);
RingSet& ringDetection(Molecule *mol, int RSIZE=-1);
void checkEdges(Atom *init, Ring &ring, AtomVector &fullSet);

# endif

