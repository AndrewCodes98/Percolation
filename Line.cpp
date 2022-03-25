//
//  Line.cpp
//  Code Playground
//
//  Created by Andrew Froning on 4/28/20.
//  Copyright (c) 2020 Andrew Froning. All rights reserved.
//

#include "Line.h"

#include <cstdio>
#include <cmath>
#include <vector>

#define PI 3.1415926535897

//Initializers
Line::Line(){
    x = y = -1;
    length = theta = -1;
}
Line::Line(float inputX, float inputY, float inputTheta, float inputLength){
    x = inputX;
    y = inputY;
    length = inputLength;
    theta = inputTheta;
    
    setEndPoints();
    s = count = fiberType = -1;
    dx = dy = 0;
}

//SETTERS
void Line::setId(int cell, int id){   //SUBCELL AND POS IN SUBCELL VECTOR
    if (cell >= 0 && id >=0){
        rCellId = cell;
        rId = id;
    }
    else
        printf("Not a valid fiber ID");
}
void Line::setDis(int disX, int disY){
    dx = disX; dy = disY;
}
void Line::setS(int size){
    s=size;
}

void Line::incrementCount(){
    count++;
}
void Line::setEndPoints(){
    leftX = x - 0.5 * length * cos(theta);
    leftY = y - 0.5 * length * sin(theta);
    rightX = x + 0.5 * length * cos(theta);
    rightY = y + 0.5 * length * sin(theta);
}

void Line::setFiberType(int t){
    fiberType = t;
}

//GETTERS
int Line::getRCellId(){
    return rCellId;
}
int Line::getRootId(){
    return rId;
}

float Line::getCoorX(){
    return x;
}
float Line::getCoorY(){
    return y;
}
float Line::getTheta(){
    return theta;
}
float Line::getThetaDegrees(){
    return theta * 180 / PI;
}
float Line::getLeftX(){
    return leftX;
}
float Line::getLeftY(){
    return leftY;
}
float Line::getRightX(){
    return rightX;
}
float Line::getRightY(){
    return rightY;
}

int Line::getDisX(){
    return dx;
}
int Line::getDisY(){
    return dy;
}
int Line::getS(){
    return s;
}
int Line::getCount(){
    return count;
}
int Line::getFiberType(){
    return fiberType;
}

void Line::print(){
    printf("%i\t\t%i\t%f\t%f\t%i\t%i\t%i\t%i\n", rCellId, rId, x, y, s, dx, dy, fiberType);
}
