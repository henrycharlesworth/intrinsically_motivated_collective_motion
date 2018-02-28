//This was one of my first C++ projects so the code is a bit of a mess.
// Written by Henry Charlesworth, University of Warwick, 2018.

#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <ctime>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <random>
#include <cstdlib>

#include "birdClasses.h"

using namespace std;

//nBirds, nSteps, nF, v, dt, birdRad defined in "birdClasses.h"

//define some other parameters which can be accessed globally.
const int nSensors = 40;
const int nMovesRO = 2; //reorientation moves.
const int nMovesSC = 3; //speed change moves
const int nMovesTot = nMovesRO+nMovesSC; //no. possible moves at each step.
double dTheta = 15 * M_PI / 180.0; //angle that birds can change their orientation by in a timestep.
double movesRO[nMovesRO] = { -dTheta, dTheta }; //possible moves.
double movesSC[nMovesSC] = { v0 - dv, v0, v0 + dv };
double decNoise = 0.00; //noise added to reorientiation decisions.
double speedNoise = 0.00; //noise added to speed

//these are used by the getVisualState(...) function
double sensorFrac[nSensors];
double sensorRef[nSensors];
double sensorRange = 2 * M_PI / ((double)nSensors);
double sensorTol = 0.5;

int equilibration = 0;

int nps = numStates(nMovesTot, nF); //given nF how many future states could there possibly be along a given branch.
vector< vector<int> > possibleStates(nps, vector<int>(nSensors, 0)); //stores possible future states.
int nUnique = 0;

//stores the modelled future positions of the other birds.
vector< vector< vector<double> > > futurePositions(nBirds, vector< vector<double> >(nF + 1, vector<double>(4, 0.0)));

//variables to record positions and orientations.
double xPositions[nSteps][nBirds], yPositions[nSteps][nBirds], orientations[nSteps][nBirds], speeds[nSteps][nBirds];

//set up random number generator.
default_random_engine generator(time(NULL));
uniform_real_distribution<double> uniformReal(0.0, 1.0);
normal_distribution<double> normalDist(0.0, 1.0);

//function prototypes
bool checkCollision(int i, int nFut, Bird *birds, double xi, double yi, double oi, double vi);
void getVisualState(Bird *birdList, int nFut, int i, double cX, double cY, double cAng);
void updateTree(double exploreX, double exploreY, double exploreO, double exploreV, Bird *bird, int bn, int nFut);

int main(int, const char * const[])
{
	//initialize - give angular range of each sensor (so sensor[0] accounts for angle 0 -> sensorRange, sensor[1] from sensorRange->2*sensorRange etc)
	sensorRef[0] = sensorRange;
	for (int u = 1; u<nSensors; u++) sensorRef[u] = sensorRef[u - 1] + sensorRange;

	//size of the box that birds are initialized in.	
	double initBoxX = 100, initBoxY = 100;

	//set up output to save data.
	ostringstream namestring;
	namestring << "outputN=" << nBirds << ".txt";
	string str1 = namestring.str();
	ofstream output(str1.c_str());

	//initialize birds in a box randomly, all with the similar orientations
	Bird birdList[nBirds];
	for (int i = 0; i<nBirds; i++) {
		birdList[i].set_position(initBoxX*uniformReal(generator), initBoxY*uniformReal(generator));
		birdList[i].update_Orientation(normalDist(generator)*dTheta);
	}


	int uniqueVisStates[nMovesTot]; //store the number of unique visual states for each available move.
	double cX, cY, fX, fY, exploreX, exploreY, exploreO, exploreV;

	//main time step loop
	for (int ts = 0; ts<nSteps; ts++) {

		//save current positions
		for (int i = 0; i<nBirds; i++) {
			xPositions[ts][i] = birdList[i].get_xPos();
			yPositions[ts][i] = birdList[i].get_yPos();
			orientations[ts][i] = birdList[i].get_orientation();
			speeds[ts][i] = birdList[i].get_V();
		}

		//fill in current positions
		for (int bn = 0; bn < nBirds; bn++) {
			futurePositions[bn][0][0] = birdList[bn].get_xPos();
			futurePositions[bn][0][1] = birdList[bn].get_yPos();
			futurePositions[bn][0][2] = birdList[bn].get_orientation();
			futurePositions[bn][0][3] = birdList[bn].get_V();
		}

		//loop over birds to choose how they each choose to update their orientation.
		for (int bn = 0; bn<nBirds; bn++) {

			//fill in first move they all make assuming that every other bird moves in a straight line at v0.
			for (int bn2 = 0; bn2 < nBirds; bn2++) {
				if (bn2 == bn) continue;
				double dAngle = 0.0;
				double dSpeed = 0.0;
				double updAngle = futurePositions[bn2][0][2] + dAngle;
				double updSpeed = v0;
				futurePositions[bn2][1][0] = futurePositions[bn2][0][0] + updSpeed*dt*cos(updAngle);
				futurePositions[bn2][1][1] = futurePositions[bn2][0][1] + updSpeed*dt*sin(updAngle);
				futurePositions[bn2][1][2] = futurePositions[bn2][0][2] + dAngle;
				futurePositions[bn2][1][3] = updSpeed;
			}

			//we now loop over each of the possible moves that the bird can make in the present.
			for (int l = 0; l<nMovesTot; l++) {
				uniqueVisStates[l] = 0; //clear
			}
			for (int mn = 0; mn<nMovesTot; mn++) {
				for (int l = 0; l<nps; l++) {
					for (int d = 0; d < nSensors; d++) {
						possibleStates[l][d] = 0; //clear
					}
				}
				nUnique = 0;
				if (mn < nMovesRO) {
					exploreO = birdList[bn].get_orientation() + movesRO[mn];
					exploreV = v0;
				}
				else {
					exploreO = birdList[bn].get_orientation();
					exploreV = movesSC[mn - nMovesRO];
				}
				exploreX = birdList[bn].get_xPos() + cos(exploreO)*exploreV*dt;
				exploreY = birdList[bn].get_yPos() + sin(exploreO)*exploreV*dt;
				//exploreX, exploreY, exploreO and exploreV are this birds position, orientation and speed at the next timestep if it were to choose this move.
				//we pass this into updateTree which then goes on to recursively explore a tree of future states from this next position. Note that nUnique (global) gets updated along the way.
				updateTree(exploreX, exploreY, exploreO, exploreV, &birdList[0], bn, 0);
				uniqueVisStates[mn] = nUnique; //number of unique visual states down this future branch.

			}
			//now we evaluate which presently available move leads to the largest number of unique possible visual states.
			int maxInd = 0, maxVal = uniqueVisStates[0];
			for (int h = 1; h<nMovesTot; h++) {
				if (uniqueVisStates[h] > maxVal) {
					maxInd = h; maxVal = uniqueVisStates[h];
				}
				else if (uniqueVisStates[h] == maxVal && h<nMovesRO) {
					if (abs(movesRO[h])<abs(movesRO[maxInd])) { //by default do the thing that reorients you less (if same number of unique visual states).
						maxInd = h;
					}
					else { //randomly choose whether to do this one or not.
						double randNum = uniformReal(generator);
						if (randNum < 0.5) {
							maxInd = h;
						}
					}
				}
			}
			//now we actually apply this update.
			if (maxInd < nMovesRO) {
				double cv = v0 + speedNoise*normalDist(generator);
				birdList[bn].update_Orientation(movesRO[maxInd] + decNoise*normalDist(generator));
				birdList[bn].update_Pos(birdList[bn].get_xPos() + cos(birdList[bn].get_orientation())*cv*dt, birdList[bn].get_yPos() + sin(birdList[bn].get_orientation())*cv*dt);
				birdList[bn].set_V(v0);
			}
			else {
				double cv = movesSC[maxInd - nMovesRO] + speedNoise*normalDist(generator);
				birdList[bn].set_V(cv);
				birdList[bn].update_Pos(birdList[bn].get_xPos() + cos(birdList[bn].get_orientation())*cv*dt, birdList[bn].get_yPos() + sin(birdList[bn].get_orientation())*cv*dt);
			}
			
		}

		for (int bn = 0; bn<nBirds; bn++) birdList[bn].finishUpdate();
		

		cout << "Timestep: " << ts << " completed!\n";
	} //end main timestep loop

	//OUTPUT DATA INTO A TEXT FILE.
	for (int ts = 0; ts<(nSteps - 1); ts++) {
		for (int bn = 0; bn<nBirds; bn++) {
			output << xPositions[ts][bn] << " " << yPositions[ts][bn] << " " << orientations[ts][bn] << " " << speeds[ts][bn] << "\n";
		}
	}

	return 0;
}

//check whether a bird has collided during its move from its current position to its new position (xi,yi,oi,vi).
//assume that it and other birds move straight during the timestep (oi and vi are updated at the start of the timestep).
bool checkCollision(int i, int nFut, Bird *birds, double xi, double yi, double oi, double vi) {
	double vrx, vry, xip, yip, xjp, yjp, rx, ry, rxp, ryp, m, c, srx, sry, val;
	int cond = 1; int index, counti = 0;
	for (int index = 0; index<nBirds; index++) {
		if (index == i) continue;
		double xj = futurePositions[index][nFut][0];
		double yj = futurePositions[index][nFut][1];
		double oj = futurePositions[index][nFut][2];
		double vj = futurePositions[index][nFut][3];
		//if they are now overlapping, then they have collided
		if ((xi - xj)*(xi - xj) + (yi - yj)*(yi - yj) < 4 * birdRad*birdRad) {
			return 1;
		}
		else { //check whether trajectories have crossed.
			   //get positions at previous timestep.
			xip = xi - vi*dt*cos(oi); yip = yi - vi*dt*sin(oi);
			xjp = xj - vj*dt*cos(oj); yjp = yj - vj*dt*sin(oj);
			//if they were overlapping at the start, actually we return 0. Only time this will happen is when a collision
			//happens by chance and two birds start a timestep overlapping (further tree branches where this happens will be cut off
			//as a coll is detected at the previous timestep). Need a way (other than noise) for birds to move away.
			if ((xip - xjp)*(xip - xjp) + (yip - yjp)*(yip - yjp) < 4 * birdRad*birdRad) {
				return 0;
			}
			//move into frame of bird i
			vrx = vj*cos(oj) - vi*cos(oi); vry = vj*sin(oj) - vi*sin(oi); //relative velocity.
			rx = xj - xi; ry = yj - yi; rxp = xjp - xip; ryp = yjp - yip; //relative positions.
			m = (ry - ryp) / (rx - rxp); c = ry - m*rx;
			srx = -c / (m + (1 / m));
			sry = -srx / m;
			if ((rx<rxp && srx > rx && srx < rxp) || (rx > rxp && srx < rx && srx > rxp)) {
				//min dist lies between timestep start/endpoints
				if (srx*srx + sry*sry < 4 * birdRad*birdRad) {
					return 1;
				}
			}
		}
		counti++;
		if (counti == nBirds) break;
	}
	return 0;
}

//finds the visual state of bird i based on its current "exploring position" (cX, cY, cAng) and the predicted positions of other birds at timestep nFut.
//visual state is defined by discretizing the bird's field of view into nSensors (relative to current orientation) and creating a vector of
//0s and 1s depending on whether each sensor is < half (or sensorFrac) filled or not.
void getVisualState(Bird *birdList, int nFut, int i, double cX, double cY, double cAng) {
	 
	double relX, relY, relDist, dAng, s, dTheta, ang1, ang2;
	//clear current visual state.
	birdList[i].resetVisualState();
	//clear sensorFrac vector
	for (int w = 0; w < nSensors; w++) sensorFrac[w] = 0.0;

	for (int j = 0; j<nBirds; j++) {
		if (i == j) continue;
		relX = futurePositions[j][nFut][0] - cX;
		relY = futurePositions[j][nFut][1] - cY;
		relDist = sqrt(relX*relX + relY*relY);
		dAng = acos((cos(cAng)*relX + sin(cAng)*relY) / relDist);
		dTheta = atan(birdRad / relDist);
		if (isnan(dAng) || isnan(dTheta)) { //very rarely we get a weird result (probably if they're very close to overlapping). This happens so infrequently it's easiest to ignore.
			continue;
		}
		s = cos(cAng)*relY - sin(cAng)*relX;
		if (s<0) dAng = 2 * M_PI - dAng;
		ang1 = dAng - dTheta; ang2 = dAng + dTheta;
		if (ang1 < 0) {
			birdList[i].addInterval(0, ang2);
			birdList[i].addInterval(2 * M_PI + ang1, 2 * M_PI);
		}
		else if (ang2 > 2 * M_PI) {
			birdList[i].addInterval(0, fmod(ang2, 2 * M_PI));
			birdList[i].addInterval(ang1, 2 * M_PI);
		}
		else {
			birdList[i].addInterval(ang1, ang2);
		}
	}
	Node *sI = birdList[i].get_visualState();
	birdList[i].cleanUp(sI);
	//now the angular interval due to each other bird has been added and merged together into a set of non-overlapping intervals.
	//no discretization until this point - now we look at how much each discrete sensor is filled.
	int ind1, ind2;
	for (int k = 0; k<nSensors; k++) sensorFrac[k] = 0.0; //initialize.
	while (sI->next->next != 0) {
		ang1 = sI->value; ang2 = sI->next->value;
		ind1 = floor(ang1 / sensorRange); ind2 = floor(ang2 / sensorRange);
		if (ind2 == nSensors) ind2--; //this happens if ang2 = 2pi (which can happen a lot).
		if (ind1 == ind2) {
			sensorFrac[ind1] += (ang2 - ang1) / sensorRange;
		}
		else if (ind2 - ind1 == 1) {
			sensorFrac[ind1] += (sensorRef[ind1] - ang1) / sensorRange;
			sensorFrac[ind2] += (ang2 - sensorRef[ind1]) / sensorRange;
		}
		else {
			sensorFrac[ind1] += (sensorRef[ind1] - ang1) / sensorRange;
			sensorFrac[ind2] += (ang2 - sensorRef[ind2 - 1]) / sensorRange;
			for (int y = ind1 + 1; y<ind2; y++) sensorFrac[y] = 1.0;
		}
		sI = sI->next->next;
	}
	//do final interval separately.
	ang1 = sI->value; ang2 = sI->next->value;
	ind1 = floor(ang1 / sensorRange); ind2 = floor(ang2 / sensorRange);
	if (ind2 == nSensors) ind2--; //this happens if ang2 = 2pi (which can happen a lot).
	if (ind1 == ind2) {
		sensorFrac[ind1] += (ang2 - ang1) / sensorRange;
	}
	else if (ind2 - ind1 == 1) {
		sensorFrac[ind1] += (sensorRef[ind1] - ang1) / sensorRange;
		sensorFrac[ind2] += (ang2 - sensorRef[ind1]) / sensorRange;
	}
	else {
		sensorFrac[ind1] += (sensorRef[ind1] - ang1) / sensorRange;
		sensorFrac[ind2] += (ang2 - sensorRef[ind2 - 1]) / sensorRange;
		for (int y = ind1 + 1; y<ind2; y++) sensorFrac[y] = 1.0;
	}
	//vector sensorFrac contains the visual state.
	int val;
	//now add the state to list if it is unique.
	if (nUnique == 0) {
		for (int y = 0; y < nSensors; y++) {
			if (sensorFrac[y] > sensorTol) {
				possibleStates[0][y] = 1;
			}
			else {
				possibleStates[0][y] = 0;
			}
		}
		nUnique++;
	}
	else { //check if state is unique
		int isUnique = 1; //assume it is
		int maxCount = 0, cCount = 0;
		for (int y = (nUnique - 1); y > -1; y--) {
			cCount = 0;
			for (int p = 0; p < nSensors; p++) {
				val = (sensorFrac[p] > sensorTol) ? 1 : 0;
				if (val == possibleStates[y][p]) cCount++;
			}
			if (cCount > maxCount) maxCount = cCount;
			if (maxCount == nSensors) {
				isUnique = 0;
				break;
			}
		}
			
		if (isUnique == 1) {
			for (int v = 0; v < nSensors; v++) {
				if (sensorFrac[v] > sensorTol) {
					possibleStates[nUnique][v] = 1;
				}
				else {
					possibleStates[nUnique][v] = 0;
				}
			}
			nUnique++;
		}
	}

	return;
}

//this function gets called recursively to explore the tree of possible futures.
void updateTree(double exploreX, double exploreY, double exploreO, double exploreV, Bird *bird, int bn, int nFut) {
	double o, x, y, v;
	if (checkCollision(bn, nFut + 1, bird, exploreX, exploreY, exploreO, exploreV)) return;
	getVisualState(bird, nFut + 1, bn, exploreX, exploreY, exploreO);

	if (nFut == (nF - 1)) return;

	//fill in current positions of this exploring bird (for other birds to make decision w/ heuristic later on.
	futurePositions[bn][nFut + 1][0] = exploreX;
	futurePositions[bn][nFut + 1][1] = exploreY;
	futurePositions[bn][nFut + 1][2] = exploreO;
	futurePositions[bn][nFut + 1][3] = exploreV;
	//use this to calculate the next positions (if necessary)
	if (nFut < (nF - 1)) {
		for (int bn2 = 0; bn2 < nBirds; bn2++) {
			if (bn2 == bn) continue;
			double dAngle = 0.0;
			double dSpeed = 0.0;
			double newSpeed = v0;
			double updAngle = futurePositions[bn2][nFut + 1][2] + dAngle;
			futurePositions[bn2][nFut + 2][0] = futurePositions[bn2][nFut + 1][0] + newSpeed*dt*cos(updAngle);
			futurePositions[bn2][nFut + 2][1] = futurePositions[bn2][nFut + 1][1] + newSpeed*dt*sin(updAngle);
			futurePositions[bn2][nFut + 2][2] = futurePositions[bn2][nFut + 1][2] + dAngle;
			futurePositions[bn2][nFut + 2][3] = newSpeed;
		}
	}

	//if we haven't yet explored out to the final number of future steps nF we call this function again for each
	//possible move we can make.
	if (nFut < (nF - 1)) {
		for (int m = 0; m<nMovesTot; m++) {
			if (m < nMovesRO) {
				o = exploreO + movesRO[m];
				v = v0;
			}
			else {
				o = exploreO;
				v = movesSC[m - nMovesRO];
			}
			x = exploreX + cos(o)*v*dt;
			y = exploreY + sin(o)*v*dt;
			updateTree(x, y, o, v, bird, bn, nFut + 1);
		}
	}
	else {
		return;
	}
}