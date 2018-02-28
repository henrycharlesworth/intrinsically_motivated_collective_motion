#ifndef BIRDCLASSES_H_INCLUDED
#define BIRDCLASSES_H_INCLUDED

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>

using namespace std;

//DEFINE SOME GLOBAL PARAMETERS OF THE SIMULATION
const int nBirds = 50;
const int nF = 4; //number of future timesteps to consider.
const int nSteps = 1000; //number of simulation steps
const double v0 = 10.0, dt = 1, birdRad = 1, dv = 2; // v0 is default speed. Other values are v0+dv, v0-dv. dt is timestep, birdRad is radius of each bird. Basically we use these as our units of time/distance

int numStates(int numMoves, int nFut) {
	int num = 1; int multiplier = numMoves;
	for (int i = 1; i<nFut; i++) {
		num += multiplier;
		multiplier *= numMoves;
	}
	return num;
}

//Node class for a linked list (used in constructing the visual states).
class Node {
public:
	int identifier; // 0 is left side of interval, 1 is right side
	double value; //angular value.
	Node *next; //pointer to the next interval.
	void display(Node *start);
};

//printout linked list if necessary (mainly for debugging purposes).
void Node::display(Node *start) {
	if (start != 0) {
		double inter = start->value;
		cout << inter << " ";
		display(start->next);
	}
}

//bird class.
class Bird {
	double currX, currY;
	double currV;
	double updatedX, updatedY, updatedV;
	double currOrientation;
	Node *visualState;
public:
	Bird() {
		currOrientation = 0.0; currX = 0.0; currY = 0.0, currV = v0;
		visualState = new Node;
		visualState->value = 0.0;
		visualState->next = new Node;
		visualState->next->value = 0.0;
		visualState->next->next = 0;
	}
	Bird(double x, double y, double o, double v) {
		currX = x; currY = y; currOrientation = o, currV = v;
		visualState = new Node;
		visualState->value = 0.0;
		visualState->next = new Node;
		visualState->next->value = 0.0;
		visualState->next->next = 0;
	}
	void set_position(double x, double y) {
		currX = x; currY = y;
	}
	void set_V(double v) {
		updatedV = v;
	}
	double get_xPos() {
		return currX;
	}
	double get_yPos() {
		return currY;
	}
	double get_orientation() {
		return currOrientation;
	}
	double get_V() {
		return currV;
	}
	//return pointer to first node.
	Node* get_visualState() {
		return visualState;
	}
	void update_Pos(double x, double y) {
		updatedX = x;
		updatedY = y;
	}
	//run this after all birds have updated positions:
	void finishUpdate() {
		currX = updatedX;
		currY = updatedY;
		currV = updatedV;
	}
	void update_Orientation(double o) {
		currOrientation += o;
	}

	//add the interval defined by [l r] to the visual state.
	void addInterval(double l, double r) {
		int placed = 0; double cL = 0.0; double cR = 0.0;
		if (visualState->value == 0.0 && visualState->next->value == 0.0) { //then this is first interval to place.
			visualState->value = l;
			visualState->next->value = r;
			placed = 1;
			return;
		}
		Node *curr_L = visualState;
		Node *prev_L = visualState;
		while (placed == 0) {
			cL = curr_L->value;
			cR = curr_L->next->value;
			if (l<cL && r<cL) { //add new interval before this one.
				Node *newRoot = new Node;
				newRoot->value = l;
				newRoot->identifier = 0;
				newRoot->next = new Node;
				newRoot->next->value = r;
				newRoot->next->next = curr_L;
				if (curr_L == visualState) {
					visualState = newRoot;
				}
				else {
					prev_L->next->next = newRoot;
				}
				placed = 1;
			}
			else if (l <= cL && r >= cR) {
				curr_L->value = l;
				curr_L->next->value = r;
				//technically should check if next values in linked list are also < r. but resolved by just not letting sensorFrac > 1.
				placed = 1;
			}
			else if (l <= cL && r <= cR) {
				curr_L->value = l;
				placed = 1;
			}
			else if (l >= cL && r <= cR) {
				placed = 1; //dont need to do anything.
			}
			else if (l >= cL && l <= cR && r >= cR) {
				curr_L->next->value = r;
				placed = 1;
			}

			if (l > cR && r > cR) {
				if (curr_L->next->next != 0) {
					prev_L = curr_L;
					curr_L = curr_L->next->next;
				}
				else {
					Node *newEndL = new Node;
					newEndL->value = l;
					newEndL->identifier = 0;
					newEndL->next = new Node;
					newEndL->next->value = r;
					newEndL->next->identifier = 1;
					newEndL->next->next = 0;
					curr_L->next->next = newEndL;
					placed = 1;
				}
			}
		}
	}
	//remove any remaining overlaps.
	void cleanUp(Node *start) {
		Node *NP, *NNP;
		NP = start->next->next;
		if (NP == 0) return;
		NNP = start->next->next->next->next;
		double cR = start->next->value, nL = start->next->next->value, nR = start->next->next->next->value;

		if (nL < cR) {
			if (nR > cR) {
				start->next->value = nR;
				cR = nR;
			}

			start->next->next = NNP;
			delete NP->next;
			delete NP;
		}
		if (NNP != 0) {
			nL = NNP->value;
		}
		if (NNP != 0 && nL < cR) {//go again
			cleanUp(start);
		}
		else {
			if (NNP != 0) cleanUp(start->next->next);
		}
	}
	//reset the visual state.
	void resetVisualState() {
		Node *cNode = visualState;
		Node *nNode = visualState->next;
		while (nNode != 0) {
			delete cNode;
			cNode = nNode;
			nNode = nNode->next;
		}
		delete cNode;
		delete nNode;
		visualState = new Node;
		visualState->identifier = 0;
		visualState->value = 0.0;
		visualState->next = new Node;
		visualState->next->identifier = 1;
		visualState->next->value = 0.0;
		visualState->next->next = 0;
		return;
	}
};

#endif