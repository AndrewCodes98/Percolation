//
//  point.cpp
//  Cary Code Testing
//
//  Created by Andrew Froning on 3/19/20.
//  Copyright (c) 2020 Andrew Froning. All rights reserved.
//

#include "Point.h"
#include <cstdio>
#include <cmath>
#include <sstream>

#define d 1

//INITIALIZERS
Point::Point(){
    rCellId = rId = -1;
    x = y = -1;
    dx = dy = 0;
    s = 0;
    fiberNum = diskType = -1;
}

Point::Point(float InputX, float InputY, int subCell, int vecPos){
    rCellId = subCell;
    rId = vecPos;
    x = InputX;
    y = InputY;
    dx = dy = 0;
    s = -1;
    fiberNum = diskType = -1;
}

//SETTERS
void Point::setId(int cell, int id){//SUBCELL AND POS IN SUBCELL VECTOR
    if (cell >= 0 && id >=0){
        rCellId = cell;
        rId = id;
    }
    else
        printf("Not a valid disk ID");
}

void Point::setDis(int disX, int disY){
    dx = disX; dy = disY;
}

void Point::setS(int size){
    s=size;
}
void Point::setFiberNum(int fnum){
    fiberNum = fnum;
}
void Point::setDiskType(int t){
    diskType = t;
}

//GETTERS
int Point::getRCellId(){//SUBCELL
    return rCellId;
}

int Point::getRootId(){//POS IN SUBCELL VECTO{
    return rId;
}

float Point::getCoorX(){
    return x;
}

float Point::getCoorY(){
    return y;
}

int Point::getDisX(){
    return dx;
}

int Point::getDisY(){
    return dy;
}

int Point::getS(){
    return s;
}

int Point::getFiberNum(){
    return fiberNum;
}

int Point::getDiskType(){
    return diskType;
}
void Point::print(){
    printf("%i\t\t%i\t%f\t%f\t%i\t%i\t%i\t%i\t\t%i\n", rCellId, rId, x, y, s, dx, dy, diskType, fiberNum);
}
