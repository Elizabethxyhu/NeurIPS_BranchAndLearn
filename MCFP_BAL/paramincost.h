#include <iostream>
#include <list>
#include <algorithm>
#include <vector>
#include <iomanip>
#include <string>
#include <sstream>
#include <queue>
#include <map>
#include <limits>
#include <cmath>

#define DBL_PRECISION 1E-5

// A parameterized edge is represented as a struct
// Fields:
//      destination   -  denote the ending node of an edge. For example, 'v' in u-->v
//      capacity      -  the maximum capacity of an edge
//      residualFlow  -  the residual amount of flow that can flow through the edge
//      cost          -  a tuple (w, b) to represent the linear function w\alpha + b
//      counterEdge   -  a pointer to the counter edge in residual graph for performance optimization
struct ParaEdge{
    int destination;
    int capacity;
    int residualFlow;
    std::pair<double,double> cost;
    double realCost; 
    ParaEdge* counterEdge;
};

// A graph is represented as a struct
// Fields:
//      numVertices - denotes the number of vertices in the graph
//      adj         - Adjacency list : Collection of unordered lists one for each vertex
struct ParaGraph {
    int numVertices;
    vector<ParaEdge*> *adj;
};

struct ParaValue {
    std::pair<double, double> distance; // linear function for distance 
    int parentVertex; // vertices number 
    ParaEdge* parentEdge; // parent edge  
}; 

// piecewise data structure of table for BF 
// upper break point as the key 
typedef std::map<double, ParaValue> ParaCell; 
typedef std::map<double, double> piecewiseConst; 

// Generates a new edge (allocating space dynamically) and returns a pointed to the edge
ParaEdge* genEdge(int destination, int capacity, int residualFlow, double realCost, double w, double b){
    ParaEdge* e1 = new ParaEdge;
    e1->destination = destination;
    e1->capacity = capacity;
    e1->residualFlow = residualFlow;
    e1->cost = std::make_pair(w, b);
    e1->realCost = realCost; 
    return e1;
}

// Parameterized version of BellmanFord algorithm to find function to map alpha to the shortest path 
void ParaBF(ParaGraph resGraph, int source, int sink, ParaCell resultTable[], double lowerBound, double upperBound){
    // Initialize variables that will be needed 
    int numVertices = resGraph.numVertices;
    vector<ParaEdge*> *adj = resGraph.adj; 
    // paraParentVertex is a piecewise integer (constant) function: interval -> int
    // paraParentEdge is a piecewise pointer function: interval -> pointer  
    // distance is a piecewise linear function: interval -> linear function or tuple (w, b)
    // Initialize parameterized table 
    // all cells of the paraTable are empty, except for the source 
    for(int i=0; i < numVertices; i++){
        ParaValue temp = {std::make_pair(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()), -1, NULL}; 
        resultTable[i][upperBound] = temp; 
    }
    resultTable[source][upperBound].distance = std::make_pair(0,0);
    // do the parameterize BellmanFord algorithm 
    for(int i=0; i<numVertices-1; i++){
        // loop on all edges 
        // create a new paraTable and initialize it to 
        ParaCell paraTable[numVertices]; 
        for(int k=0; k < numVertices; k++){
            paraTable[k] = resultTable[k];
        }

        for(int u=0; u<numVertices; u++){
            for(int e=0; e<adj[u].size(); e++) {
                if(adj[u][e]->residualFlow > 0) {
                    // if(false) {
                    if( lowerBound != std::numeric_limits<double>::infinity() && upperBound != std::numeric_limits<double>::infinity() && abs(lowerBound-upperBound) < DBL_PRECISION ) {
                        // if lowerBound and upper Bound are too close, compute as determinisitc BF algorithm 
                        // skip if the souce vertex is not reachable now 
                        if( resultTable[u][upperBound].distance.second == std::numeric_limits<double>::infinity()) 
                            continue; 
                        double alpha = lowerBound;
                        int v = adj[u][e]->destination;
                        double w = adj[u][e]->cost.first * alpha + adj[u][e]->cost.second;
                        if( (paraTable[v][upperBound].distance.second - resultTable[u][upperBound].distance.second - w) > DBL_PRECISION ) {
                            // std::cout << "\t\tFind an edge (" << u << "," << v << "), " << "distance[v]: " << paraTable[v][upperBound].distance.second << ",distance[u]: " << resultTable[u][upperBound].distance.second  << ", edge cost: " << w << std::endl;
                            paraTable[v][upperBound].distance.second = resultTable[u][upperBound].distance.second + w; 
                            paraTable[v][upperBound].parentVertex = u; 
                            paraTable[v][upperBound].parentEdge = adj[u][e]; 
                        }
                    }
                    else {
                        int v = adj[u][e]->destination;
                        std::pair<double,double> w = adj[u][e]->cost;
                        // compute the sum distance[u]+w 
                        std::map<double, ParaValue> tempCell; 
                        for(std::map<double, ParaValue>::iterator iter=resultTable[u].begin(); iter != resultTable[u].end(); iter++) {
                            ParaValue temp = {std::make_pair(w.first + (*iter).second.distance.first, w.second + (*iter).second.distance.second), u, adj[u][e]}; 
                            tempCell[(*iter).first] = temp; 
                        }

                        // skip if the vertex is not reachable now 
                        if( (*tempCell.begin()).second.distance.first == std::numeric_limits<double>::infinity() )
                            continue;

                        // compute the InfComb of distance[v] and (distance[u]+w) 
                        std::map<double, ParaValue>::iterator alphaIter = paraTable[v].begin(); 
                        std::map<double, ParaValue>::iterator betaIter = tempCell.begin();
                        std::map<double, ParaValue> resultCell; 
                        double lowerPt = lowerBound;

                        // // output tempCell 
                        // cout << std::endl << "==========================================================================" << std::endl;
                        // cout << "\tFind an edge (" << u << "," << adj[u][e]->destination << ") with nonzero residual flow, cost function: " << adj[u][e]->cost.first << " * alpha + " << adj[u][e]->cost.second << std::endl;
                        // lowerPt = lowerBound;
                        // std::cout << "\t\tcurrent paraTable[" << v << "]: " << std::endl;
                        // for(std::map<double, ParaValue>::iterator iter = paraTable[v].begin(); iter != paraTable[v].end(); iter++ ) { 
                        //     std::cout << "\t\t\t(" << lowerPt << "," << (*iter).first << "], distFunc: " << (*iter).second.distance.first << " * alpha + " << (*iter).second.distance.second;
                        //     std::cout << ", ParentVertex: " << (*iter).second.parentVertex << std::endl; 
                        //     lowerPt = (*iter).first;
                        // }

                        // lowerPt = lowerBound;
                        // std::cout << "\t\tcurrent paraTable[" << u << "]: " << std::endl;
                        // for(std::map<double, ParaValue>::iterator iter = resultTable[u].begin(); iter != resultTable[u].end(); iter++ ) { 
                        //     std::cout << "\t\t\t(" << lowerPt << "," << (*iter).first << "], distFunc: " << (*iter).second.distance.first << " * alpha + " << (*iter).second.distance.second;
                        //     std::cout << ", ParentVertex: " << (*iter).second.parentVertex << std::endl; 
                        //     lowerPt = (*iter).first;
                        // }

                        // lowerPt = lowerBound;
                        // std::cout << "\t\tparaTable[" << u << "] + cost of (" << u << "," << v << "): " << std::endl;
                        // for(std::map<double, ParaValue>::iterator iter = tempCell.begin(); iter != tempCell.end(); iter++ ) { 
                        //     std::cout << "\t\t\t(" << lowerPt << "," << (*iter).first << "], distFunc: " << (*iter).second.distance.first << " * alpha + " << (*iter).second.distance.second;
                        //     std::cout << ", ParentVertex: " << (*iter).second.parentVertex << std::endl; 
                        //     lowerPt = (*iter).first;
                        // }

                        lowerPt = lowerBound;
                        while(lowerPt < upperBound && (alphaIter !=  paraTable[v].end() || betaIter != tempCell.end()) ) {
                            bool updateAlpha = true; 
                            if(betaIter == tempCell.end() ) {
                                resultCell[(*alphaIter).first] = (*alphaIter).second; // copy ParaCell to resultCell 
                                alphaIter++; // update the alphaIter 
                            }
                            else if(alphaIter == paraTable[v].end()){
                                resultCell[(*betaIter).first] = (*betaIter).second; // copy tempCell to result Cell 
                                betaIter++; // update the betaIter 
                            } else {
                                // std::cout << "\t\tMerge linear functions" << std::endl;
                                // compute the new break point 
                                double breakPtA = (*alphaIter).first; 
                                double breakPtB = (*betaIter).first;
                                // compute the intersection point of two linear functions 
                                std::pair<double, double> linFuncA = (*alphaIter).second.distance; 
                                std::pair<double, double> linFuncB = (*betaIter).second.distance; 
                                double interPt = ( linFuncB.second - linFuncA.second ) / ( linFuncA.first - linFuncB.first ); 
                                // std::cout << "\t\tbreakPtA: " << breakPtA << ", breakPtB: " << breakPtB << ", interPt: " << interPt << std::endl;
                                // std::cout << "\t\tlinFuncA: " << linFuncA.first << " * alpha + " << linFuncA.second << std::endl;
                                // std::cout << "\t\tlinFuncB: " << linFuncB.first << " * alpha + " << linFuncB.second << std::endl;
                                // std::cout << "\t\tlowerPt: " << lowerPt << ", interPt: " << interPt << ", breakPtA: " << breakPtA << ", breakPtB: " << breakPtB << std::endl;
                                // std::cout << "\t\tlowerPt < interPt: " << (lowerPt < interPt) << ", interPt < breakPtA:" << (interPt < breakPtA) << ", interPt < breakPtB:" << (interPt < breakPtB) << std::endl;
                                // std::cout << "\t\tinterPt - lowerPt: " << (interPt - lowerPt) << ", breakPtA - interPt:" << (breakPtA - interPt) << ", breakPtB - interPt:" << (breakPtB - interPt) << std::endl;
                                // std::cout << "\t\tlinFuncA.first == linFuncB.first: " << (linFuncA.first == linFuncB.first) << "(linFuncA.first - linFuncB.first): " << (linFuncA.first - linFuncB.first) << std::endl;
                                // std::cout << "\t\t(abs(linFuncA.first - linFuncB.first) < DBL_PRECISION): " << (abs(linFuncA.first - linFuncB.first) < DBL_PRECISION) << std::endl;
                                // note that lowerPt < breakPtA and lowerPt < breakPtB 
                                // double lowerTransiPt = min(breakPtA, breakPtB); 
                                if(linFuncA.first == std::numeric_limits<double>::infinity() || linFuncA.second == std::numeric_limits<double>::infinity()) {
                                    // handle the initial case where the node is first reached 
                                    // replace the linear function 
                                    // std::cout << "\t\thandle the initial case where the node is first reached" << std::endl; 
                                    ParaValue temp; 
                                    temp.distance = (*betaIter).second.distance; 
                                    temp.parentVertex = (*betaIter).second.parentVertex; 
                                    temp.parentEdge = adj[u][e]; 
                                    resultCell[min(breakPtA, breakPtB)] = temp; 
                                    alphaIter++; 
                                    betaIter++;
                                }
                                else if( abs(linFuncA.first - linFuncB.first) < DBL_PRECISION ) {
                                    // two linear function are parallel, and there exists no intersection point 
                                    // std::cout << "\t\ttwo linear function are parallel, and there exists no intersection point" << std::endl;
                                    // std::cout << "\t\tlinFuncA.second: " << linFuncA.second <<  ", linFuncB.second: " << linFuncB.second << std::endl;
                                    if((linFuncA.second < linFuncB.second) || abs(linFuncA.second - linFuncB.second) < DBL_PRECISION ) { 
                                        resultCell[min(breakPtA, breakPtB)] = (*alphaIter).second; 
                                    } else {
                                        ParaValue temp; 
                                        temp.distance = (*betaIter).second.distance; 
                                        temp.parentVertex = u; 
                                        temp.parentEdge = adj[u][e]; 
                                        // std::cout << temp.distance.first << " * alpha + " << temp.distance.second << std::endl;
                                        resultCell[min(breakPtA, breakPtB)] = temp; 
                                    }
                                    // update the lower Pt 
                                    if(breakPtA < breakPtB || breakPtA == breakPtB || abs(breakPtA - breakPtB) < DBL_PRECISION) {
                                        lowerPt = breakPtA; 
                                        alphaIter++;  
                                    }
                                    if(breakPtA > breakPtB || breakPtA == breakPtB || abs(breakPtA - breakPtB) < DBL_PRECISION) {
                                        lowerPt = breakPtB; 
                                        betaIter++; 
                                    }
                                }
                                else if( interPt - lowerPt > DBL_PRECISION  && breakPtA - interPt > DBL_PRECISION  && breakPtB - interPt > DBL_PRECISION ) {
                                    // std::cout << "\t\tintersection point is within the range" << std::endl; 
                                    // intersection point is within the range 
                                    // insert two linear function into the resultCell 
                                    ParaValue temp; 
                                    temp.distance = (*betaIter).second.distance; 
                                    temp.parentVertex = u; 
                                    temp.parentEdge = adj[u][e]; 
                                    if(linFuncA.first < linFuncB.first) {
                                        resultCell[interPt] = temp; 
                                        resultCell[min(breakPtA, breakPtB)] = (*alphaIter).second; 
                                    } else {
                                        resultCell[interPt] = (*alphaIter).second; 
                                        resultCell[min(breakPtA, breakPtB)] = temp; 
                                    }

                                    // update the lower Pt 
                                    if(breakPtA < breakPtB || breakPtA == breakPtB || abs(breakPtA - breakPtB) < DBL_PRECISION) {
                                        lowerPt = breakPtA; 
                                        alphaIter++;  
                                    }
                                    if(breakPtA > breakPtB || breakPtA == breakPtB || abs(breakPtA - breakPtB) < DBL_PRECISION) {
                                        lowerPt = breakPtB; 
                                        betaIter++; 
                                    }
                                } else {
                                    // std::cout << "\t\tthere is an intersection point but it is out of range " << std::endl;
                                    // there is an intersection point but it is out of range 
                                    // use the function value at lowerPt to determine the relative magnitude of two linear functions 
                                    if(interPt < lowerPt || abs(interPt - lowerPt) < DBL_PRECISION) {
                                    // if(interPt < lowerPt ) {
                                        // std::cout << "intersection point is on the left " << std::endl;
                                        // intersection point is on the left 
                                        if(linFuncA.first < linFuncB.first) 
                                            resultCell[min(breakPtA, breakPtB)] = (*alphaIter).second; 
                                        else {
                                            ParaValue temp; 
                                            temp.distance = (*betaIter).second.distance; 
                                            temp.parentVertex = u; 
                                            temp.parentEdge = adj[u][e]; 
                                            resultCell[min(breakPtA, breakPtB)] = temp; 
                                        }
                                    }
                                    else {
                                        // std::cout << "intersection point is on the right " << std::endl;
                                        // intersection point is on the right 
                                        if(linFuncA.first < linFuncB.first) {
                                            ParaValue temp; 
                                            temp.distance = (*betaIter).second.distance; 
                                            temp.parentVertex = u; 
                                            temp.parentEdge = adj[u][e]; 
                                            resultCell[min(breakPtA, breakPtB)] = temp; 
                                        }
                                        else
                                            resultCell[min(breakPtA, breakPtB)] = (*alphaIter).second; 
                                    }

                                    // udpate the pointer 
                                    if(breakPtA < breakPtB || breakPtA == breakPtB || abs(breakPtA - breakPtB) < DBL_PRECISION) {
                                        lowerPt = breakPtA; 
                                        alphaIter++; 
                                    } 

                                    if(breakPtA > breakPtB || breakPtA == breakPtB || abs(breakPtA - breakPtB) < DBL_PRECISION) {
                                        lowerPt = breakPtB;
                                        betaIter++; 
                                    }
                                }
                            }
                        }
                        
                        // std::cout << std::endl << "\t\tAfter resultCell: " << std::endl;
                        // lowerPt = lowerBound;
                        // for(std::map<double, ParaValue>::iterator iter = resultCell.begin(); iter != resultCell.end(); iter++ ) { 
                        //     std::cout << "\t\t\t(" << lowerPt << "," << (*iter).first << "], distFunc: " << (*iter).second.distance.first << " * alpha + " << (*iter).second.distance.second;
                        //     std::cout << ", ParentVertex: " << (*iter).second.parentVertex << std::endl;
                        //     lowerPt = (*iter).first;
                        // }

                        // merge the consecutive linear functions 
                        lowerPt = lowerBound;
                        std::vector<double> removedPts; 
                        double preSlopt = 0.0;
                        double preIntercept = 0.0;
                        int preParent = -1; 
                        for(std::map<double, ParaValue>::iterator iter = resultCell.begin(); iter != resultCell.end(); iter++ ) { 
                            if(abs(preSlopt - (*iter).second.distance.first) < DBL_PRECISION && (preIntercept - (*iter).second.distance.second) < DBL_PRECISION && preParent == (*iter).second.parentVertex )
                                removedPts.push_back(lowerPt);
                            preSlopt = (*iter).second.distance.first;
                            preIntercept = (*iter).second.distance.second;
                            preParent = (*iter).second.parentVertex; 
                            lowerPt = (*iter).first;
                        }
                        for(std::vector<double>::iterator iter = removedPts.begin(); iter != removedPts.end(); iter++)
                            resultCell.erase(*iter);

                        // assign the resultCell 
                        paraTable[v] = resultCell; 

                        // lowerPt = lowerBound;
                        // std::cout << std::endl << "\t\tAfter paraTable[" << v << "]: " << std::endl;
                        // for(std::map<double, ParaValue>::iterator iter = paraTable[v].begin(); iter != paraTable[v].end(); iter++ ) { 
                        //     std::cout << "\t\t\t(" << lowerPt << "," << (*iter).first << "], distFunc: " << (*iter).second.distance.first << " * alpha + " << (*iter).second.distance.second;
                        //     std::cout << ", ParentVertex: " << (*iter).second.parentVertex << std::endl;
                        //     lowerPt = (*iter).first;
                        // }
                    }   
                }
            }
        }

        // std::cout << std::endl << "================================= Iteration " << i << " =================================" << std::endl;
        // for(int k=0; k < numVertices; k++) {
        //     std::cout << "paraTable[" << k << "]: " << std::endl;
        //     double lowerPt = lowerBound;
        //     for(std::map<double, ParaValue>::iterator iter = paraTable[k].begin(); iter != paraTable[k].end(); iter++ ) { 
        //         std::cout << "\t\t(" << lowerPt << "," << (*iter).first << "), distFunc: " << (*iter).second.distance.first << " * alpha + " << (*iter).second.distance.second;
        //         std::cout << ", parentVertex: " << (*iter).second.parentVertex; 
        //         std::cout << ", (parentEdge == NULL): " << ((*iter).second.parentEdge == NULL) << std::endl;
        //         lowerPt=(*iter).first; 
        //     }
        // }
        
        // copy the paraTable to resultTable
        for(int k=0; k < numVertices; k++)
            resultTable[k] = paraTable[k]; 
        
    }

    for(int i=0; i < resGraph.numVertices; i++) {
        ParaCell tempCell; 
        double currentAlpha = -std::numeric_limits<double>::infinity(); 
        ParaValue currentValue = (*resultTable[i].begin()).second;
        for(std::map<double, ParaValue>::iterator iter = resultTable[i].begin(); iter != resultTable[i].end(); iter++ ) { 
            if(currentValue.parentVertex != (*iter).second.parentVertex) {
              tempCell[currentAlpha] = currentValue;
              currentValue = (*iter).second; 
            }
            currentAlpha = (*iter).first; 
        }
        tempCell[currentAlpha] = currentValue;
        resultTable[i] = tempCell;
    }

    // // output the parameterized table 
    // for(int i=0; i < resGraph.numVertices; i++) {
    //   std::cout << "resultTable[" << i << "]: " << std::endl;
    //   double lowerPt = lowerBound;
    //   for(std::map<double, ParaValue>::iterator iter = resultTable[i].begin(); iter != resultTable[i].end(); iter++ ) { 
    //     std::cout << "\t\t(" << lowerPt << "," << (*iter).first;
    //     if(abs((*iter).first - upperBound) < DBL_PRECISION)
    //         std::cout << "), distFunc: ";
    //     else
    //         std::cout << "], distFunc: ";
    //     std::cout << (*iter).second.distance.first << " * alpha + " << (*iter).second.distance.second;
    //     std::cout << ", parentVertex: " << (*iter).second.parentVertex; 
    //     std::cout << ", (parentEdge == NULL): " << ((*iter).second.parentEdge == NULL) << std::endl;
    //     lowerPt = (*iter).first; 
    //   }
    // }

    // cin.ignore(); 

}

piecewiseConst paraMinCost(ParaGraph resGraph, int source, int sink, int requiredFlow, double lowerBound, double upperBound){
    piecewiseConst returnFunc; 

    // run the parameterized bellmanford algorithm 
    ParaCell paraTable[resGraph.numVertices]; 
    ParaBF(resGraph, source, sink, paraTable, lowerBound, upperBound); 

    // check whether the sink is reachable 
    if(paraTable[sink].begin()->second.parentVertex == -1) {
        // construct and return the piecewise constant function 
        returnFunc[upperBound] = std::numeric_limits<double>::infinity(); 
        return returnFunc; 
    }

    // search over the parameterized table 
    ParaCell::iterator iter[resGraph.numVertices]; 
    for(int i = 0; i < resGraph.numVertices; i++)
      iter[i] = paraTable[i].begin(); 
    double lowerPt = lowerBound; 
    double currentAlpha = 0.0; 
    while(true) {
        double currentCost=0;
        int currentVertex = sink; 
        std::set<int> path; 
        // std::cout << "path: "; 
        do {
            if(path.count(currentVertex)){
                std::cout << std::endl << "contains cycle, exit!" << std::endl;
                exit(0);
            }
            // std::cout << currentVertex << " <- ";
            path.insert(currentVertex); 
            currentVertex = (*iter[currentVertex]).second.parentVertex;
        } while(currentVertex != source);
        // std::cout << currentVertex << std::endl;

        currentAlpha = upperBound; 
        for(int i = 0; i < resGraph.numVertices; i++) {
            if( path.count(i) != 0 )
                currentAlpha = min(currentAlpha, iter[i]->first); 
        }
        // std::cout << "currentAlpha: " << currentAlpha << std::endl;
        if(currentAlpha - lowerPt > DBL_PRECISION){
            // compute the min flow 
            int pathFlow = requiredFlow;
            currentVertex = sink;
            while(currentVertex != source) {
                // std::cout << "currentVertex: " << currentVertex << ", (*iter[currentVertex]).second.parentEdge == NULL: " << ((*iter[currentVertex]).second.parentEdge == NULL) << std::endl;
                // std::cout << "residualFlow for edge (" << (*iter[currentVertex]).second.parentVertex << "," << currentVertex << ") is " << (*iter[currentVertex]).second.parentEdge->residualFlow << std::endl;
                pathFlow = min(pathFlow, (*iter[currentVertex]).second.parentEdge->residualFlow);
                currentVertex = (*iter[currentVertex]).second.parentVertex;
            }
            // std::cout << "pathFlow: " << pathFlow << std::endl;
        
            // update the residual graph 
            currentVertex = sink;
            while(currentVertex != source) {
                ParaEdge *te1 = (*iter[currentVertex]).second.parentEdge;
                ParaEdge *te2 = te1->counterEdge;
                te1->residualFlow -= pathFlow;
                te2->residualFlow += pathFlow;
                currentCost += double(pathFlow)*(te1->realCost);
                currentVertex = (*iter[currentVertex]).second.parentVertex;
            }
            // std::cout << "currentCost: " << currentCost << std::endl;

            // output the current exploring path 
            // currentVertex = sink;
            // std::cout << "range: (" << lowerPt << "," << currentAlpha << "), cost: " << currentCost << ", pathFlow: " << pathFlow << ", path: ";
            // while(currentVertex != 0) {
            //     std::cout << currentVertex << " <- ";
            //     currentVertex = (*iter[currentVertex]).second.parentVertex;
            // }
            // std::cout << currentVertex << std::endl;

            // recursively call paraMinCost and add currentCost and tempFunc 
            // todo: handle when lowerPt is too close to currentAlpha
            // std::cout << "requiredFlow-pathFlow: " << requiredFlow-pathFlow << ", abs(lowerPt-currentAlpha): "<< abs(lowerPt-currentAlpha) << std::endl;
            if(requiredFlow-pathFlow > 0) {
                piecewiseConst tempFunc = paraMinCost(resGraph, source, sink, requiredFlow-pathFlow, lowerPt, currentAlpha);
                for(piecewiseConst::iterator iter = tempFunc.begin(); iter != tempFunc.end(); iter++ )
                    returnFunc[iter->first] = iter->second + currentCost;     
            } else {
                returnFunc[currentAlpha] = currentCost;
            }

            currentVertex = sink;
            // std::cout << "return from recursive call, range: (" << lowerPt << "," << currentAlpha << "), cost: " << currentCost << ", pathFlow: " << pathFlow << ", path: ";
            // while(currentVertex != 0) {
            //     std::cout << currentVertex << " <- ";
            //     currentVertex = (*iter[currentVertex]).second.parentVertex;
            // }
            // std::cout << currentVertex << std::endl;

            // restore the resGraph 
            currentVertex = sink;
            while(currentVertex != source) {
                ParaEdge *te1 = (*iter[currentVertex]).second.parentEdge;
                ParaEdge *te2 = te1->counterEdge;
                te1->residualFlow += pathFlow;
                te2->residualFlow -= pathFlow;
                currentVertex = (*iter[currentVertex]).second.parentVertex;
            }

            // // output the returnFunc 
            // std::cout << "In paraMinCost requiredFlow: " << requiredFlow << ", lowerBound: " << lowerBound << ", upperBound: " << upperBound << std::endl;
            // double tempPt = lowerBound; 
            // std::cout << "returnFunc: " << std::endl;
            // for(piecewiseConst::iterator it = returnFunc.begin(); it != returnFunc.end(); it++) {
            //     std::cout << "\t(" << tempPt << "," << it->first << "): " << it->second << std::endl;
            //     tempPt = it->first; 
            // }
            // cin.ignore();

            // std::cout << "currentAlpha: " << currentAlpha << ", upperBound: " << upperBound << std::endl;
            if( abs(currentAlpha - upperBound) < DBL_PRECISION || currentAlpha == std::numeric_limits<double>::infinity())
                break;
        }

        // std::cout << "currentAlpha: " << currentAlpha << std::endl;
        for(int i = 0; i < resGraph.numVertices; i++) {
            while(currentAlpha >= iter[i]->first)
                (iter[i])++; 
        }
        lowerPt = currentAlpha;
        // std::cout << "udpate lowerPt to " << lowerPt << std::endl;
    }

    // std::cout << "return from paraMinCost, source: " << source << ", sink: " << sink << ", requiredFlow: " << requiredFlow << ", lowerBound: " << lowerBound << ", upperBound: " << upperBound << std::endl;
    return returnFunc;
}

