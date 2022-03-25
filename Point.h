//
//  point.h
//  Cary Code Testing
//
//  Created by Andrew Froning on 3/19/20.
//  Copyright (c) 2020 Andrew Froning. All rights reserved.
//


#ifndef __Cary_Code_Testing__point__
#define __Cary_Code_Testing__point__

#include <stdio.h>

//#ifndef CaryCodeTesting_Point_h
//#define CaryCodeTesting_Point_h
/////////////////////
class Point
{
private:
    float x, y;//coordinate
    
    int rCellId;// The cell that its root belongs to
    int s;	// Negative size of the cluster (if the point is a root)
    int fiberNum;//The fiber number the disk is on
    
    short rId;//position in vector subCell[rCellId][rId] we dont need this since it always be the first element in the array
    short dx, dy;//displacement to root
    short diskType; //Whether it is perc. cluster (1), perc cluster on perc fiber(2), or neither(3)
    
public:
    Point();
    Point(float, float, int, int);
        
    void setId(int,int);	//Set up the Root: rCellId and rId
    void setDis(int, int);
    void setS(int);
    void setFiberNum(int);
    void setDiskType(int);
    
    int getRCellId();
    int getRootId();
    
    float getCoorX();
    float getCoorY();
    
    int getDisX();
    int getDisY();
    int getS();
    int getFiberNum();
    int getDiskType();
    
    void print();
};
///////////////////////

//#endif

#endif /* defined(__Cary_Code_Testing__point__) */
