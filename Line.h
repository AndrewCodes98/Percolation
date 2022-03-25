//
//  Line.h
//  Code Playground
//
//  Created by Andrew Froning on 4/28/20.
//  Copyright (c) 2020 Andrew Froning. All rights reserved.
//

#ifndef __Code_Playground__Line__
#define __Code_Playground__Line__

#include <stdio.h>

#endif /* defined(__Code_Playground__Line__) */

class Line
{
private:
    float x, y;//coordinate
    float length, theta, leftX, rightX, leftY, rightY;
    
    int rCellId;// The cell that its root belongs to
    int count;
    int s;    // Negative size of the cluster (if the point is a root)

    short dx, dy;//displacement to root
    short rId;//position in vector subCell[rCellId][rId] we dont need this since it always be the first element in the array
    
    unsigned short fiberType; //Whether it is perc. cluster (1), or not(0)

public:
    //INITIALIZERS
    Line();
    Line(float, float, float, float);
    
    //SETTERS
    void setId(int,int);	//Set up the Root: rCellId and rId
    void setDis(int, int);
    void setS(int);
    
    void incrementCount();
    void setEndPoints();
    
    void setFiberType(int);
    
    //GETTERS
    int getRCellId();
    int getRootId();
    
    float getCoorX();
    float getCoorY();
    float getTheta();
    float getThetaDegrees();
    float getLeftX();
    float getRightX();
    float getLeftY();
    float getRightY();
    
    int getCount();
    int getDisX();
    int getDisY();
    int getS();
    int getFiberType();
    
    void print();
};
