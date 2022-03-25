//  Percolation Code V4.0
//  3/16/21

#include <unistd.h>
#include <cstdio>
#include <iomanip>
#include <ctime>
#include <vector>	// Changing to <thrust/vector> in the future will allow for CUDA vectors
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <sstream>

#include "Point.h"
#include "Line.h"

//FUCTION DECLARATIONS
const int findL();
const int calcLengthOfFiberSubCells();
void clusterPercentages(int totalDisks, int totalFibers, int trialNum);
void clusterPercentageData();
void diskChronoDataForImport(std::vector<Line> vec);
void lineChronoDataForImport(std::vector<Line> vec);
void importLines();
void fiberChronoDataForPythonCode(std::vector<Line> vec);
void plotFiberCoor();
void plotDiskCoor();
void plotStdDeviation(float SD);
void printToMasterLog(float time, float etaDisk, float etaTubule);
void printPercThresholds(int n[], float SD);
float calcEtaDisk(int n[]);
float calcEtaFiber(int n[]);
float calcStdDeviation(int n[]);
void importPoints();
Point findRoot(Point &p, int pCell, int pCellId);
void nineCells(int cId, int *array);
void bondvector(int j, int &dx, int &dy);
//float distance(Point A, Point B);
bool detectWrapping(int j, Point a, Point b);
bool detectWrapping(int j, Line a, Line b);
void setFibers();
bool onSegment(float x1, float y1, float x2, float y2, float x3, float y3);
int orientation(float x1, float y1, float x2, float y2, float x3, float y3);
bool intersect(Line A, Line B);
Line findRoot(Line &p, int pCell, int pCellId);
void nineFiberCells(int cId, int *array);
int getPointOnFiber(float &x, float &y, int count);
int populatedFibers();
void findClustersDisks(int rCell, int rId);
void findClustersFibers(int rCell, int rId);
void printAllDisks();
void printAllFibers();
int percolateCircles();
int percolateFibers();

//CONSTANTS TO LEAVE ALONE
#define d 1    //disk diameter
#define PI 3.1415926
#define ETA 1.12808737
#define fiberETA 5.6372858
float fiber_density = 10; //Use fiberETA here only for non-percolating fibers. 10 is here to give bad data if accidentally used.

#define usingDisks      true    //true = disk percolation
#define usingFibers     true    //true = using fibers
#define percFibers      true   //true = percolating fibers, false = non-percolating (pre-determined)

#define importData   false   //go to "importData()", change "fileName" to desired file. (percFiber must = true)
#define File         false
#define DEBUGGING    false
#define printClusters   true
//#define readInputs   true    //Inputs for L and fiber_length

#define NUM_RUN 1000

//Switch between either of the two lines below; to run a specific case or to use findL() for doing multiple cases.
const int L = 128;
//const int L = findL();
const int N = L * L;

#ifdef usingFibers

float fiber_length = 1.0; //will change if findL() is executed for const int L
float dtheta = 90 * PI / 180; //Maximum angle (degrees) that fiber makes w/ horizontal (Theta in radians)
const int L_FIBERS = calcLengthOfFiberSubCells();
const int N_FIBERS = L_FIBERS * L_FIBERS;

std::vector<Line> fiberSubCells[512*512]; //Can't make dynamic so let's just stick with array, doesn't affect speed

//int num_fibers = int(ETA * N * 2); //eta_c = 1.7938
int num_fibers = int(fiber_density * N / fiber_length / fiber_length);

#endif //usingFibers

std::vector<float> coorVec; //Only if importing points
std::vector<Point> subCells[512*512];
Line singleDimensionFiberArr[512*512*20];
std::vector<float> clusPercentages;
float clusPercentageArr[6][NUM_RUN]; // 3 cluster percentages and respective sigma values


/*
 /////////////////////////////////////////////////////////////
 -Steps to Plot (Manually):
 1) Run 4 system sizes of lengths L
 2) Copy and Paste 4 files into Conv. Code Directory
 3) Run Python Code
 4) Manually Type eta value and stick length into text file for gnuplot
 5) Repeat Steps 1-4
 6) Run gnuplot via terminal
 -Steps to Plot (Automatically):
 1) Allow L to be an input of Perc. Code
 2) Create new Master code that can run perc. code at various L sizes
 a) Run 4 L sizes to produce 4 files
 b) Have Conv. Code use those 4 files produced
 c) print conv. code eta val into a text file with stick size
 V4.0
 ------------------------

 -unnecessary stuff in lineChronoDataForImport()?
 
 3/19/21
 -Update file writing to append rather than rewrite.
 -Standarize file writing for lines to use two points rather than point and theta. Include appending 0 or 1 for clustering in Python code.
 -
 
 /////////////////////////////////////////////////////////////
 */


const int findL(){
    // Need to read value from text file, eliminate first line, then print remaining back
    /* Text file should be format:
     L  fiber_length
     32,1.0
     64,2.5
     etc.
     */
    
    //Go to the directory containing L_and_fibers_to_run.dat
        
//    chdir("/Users/AndrewFroning/Library/Developer/Xcode/DerivedData/Perc_4_Clustering_with_Old-djzvrccruekkbqbkrozmkgedxkco/Build/Products/Debug");
    chdir("/Users/AndrewFroning/Desktop/PercOutput");

    //This just helps when I move things around to make sure I am in the right directory.
//    char buffer[200];
//    getcwd(buffer, sizeof(buffer));
//    printf("Current directory: %s\n", buffer);
    
    std::string fileName = "L_and_fibers_to_run.dat";
    
    std::vector <std::string> lineVec;
    std::ifstream read(fileName);
    std::ofstream write;
    write << std::setprecision(8);
    
    int len = 0;
    float fiberLen = 0;
    std::string line, num;
    
    if(!read.is_open()){
        printf("L_and_fibers_to_run.dat file not found.\n");
        exit(0);
    }else{
        
        std::string line;
        getline(read, line);
        
        std::stringstream ss(line);
        
        //System size L
        if(getline(ss, num, ',')){
            len = stoi(num);
            
            //fiber length
            getline(ss,num, ',');
            fiberLen = stod(num);
            fiber_length = stod(num);
            
        }else{
            printf("L_and_fibers_to_run.dat is empty.\n");
            exit(0);
        }
        
        while(getline(read, line))
            lineVec.push_back(line);
        
        write.open(fileName);
        for(auto const& value: lineVec)
            write << value << std::endl;
        write.close();
    }
    
    read.close();

//    chdir("/Users/AndrewFroning/Desktop/Perc Code4");

    return len;
}
const int calcLengthOfFiberSubCells(){
    float l; //num of fiberSubCells along an axis
    if(modf(L / fiber_length, &l) != 0.0)
        l++;
    return l;
}

void clusterPercentages(long int totalDisks, long int totalFibers, int trialNum){
    /*
     State type of cluster percentage in clusPercentages:
     0 = Perc Disk on Perc Fiber vs Perc Disk
     1 = SD of 0
     2 = Perc Disk vs Disk          (val0 + val1) vs
     3 = SD of 2
     4 = Perc Fiber vs Fiber
     5 = SD of 4
     */
    float per0, SD0, per1, SD1, per2, SD2;
    per0 = SD0 = per1 = SD1 = per2 = SD2 = 0;
    
    for(int i = 0; i < N; i++)
        for(int j = 0; j < subCells[i].size(); j++){
            
            if(subCells[i][j].getDiskType()==2)
                per0++;//Count disks on Perc cluster and perc fibers
            else if(subCells[i][j].getDiskType()==1)
                per1++;//Count disks on perc cluster
        }
    
    for(int i = 0; i < N_FIBERS; i++)
        for(int j = 0; j < fiberSubCells[i].size(); j++)
            if(fiberSubCells[i][j].getFiberType()==1)
                per2++;//Count perc Fibers
    
    float percDiskTotal = per0 + per1;
    
    //Convert values to percentages
    per0 /= percDiskTotal;
    per1 = percDiskTotal / totalDisks;
    per2 /= totalFibers;
    
    per0 *= 100;
    per1 *= 100;
    per2 *= 100;
    
    //For multiple trials
    clusPercentageArr[0][trialNum] = per0;
    clusPercentageArr[2][trialNum] = per1;
    clusPercentageArr[4][trialNum] = per2;
}
void clusterPercentageData(){
    
    std::vector <std::string> lineVec;
    std::string fileName = "clusterPercentageData.dat";
    std::ifstream read(fileName);
    std::ofstream write;
    write << std::setprecision(8);
    
    //Format
    // L    fiber_length    per0    SD0    per1    SD1    per2    SD2
    
    //We want to find the standard deviations of the percentages
    const int len = 6;
    float arr[len] = {0}; //first row averages, second row StDev
    
    //Avergaes
    for(int i = 0; i < NUM_RUN; i++){
        arr[0] += clusPercentageArr[0][i];
        arr[2] += clusPercentageArr[2][i];
        arr[4] += clusPercentageArr[4][i];
    }

    arr[0] /= NUM_RUN;
    arr[2] /= NUM_RUN;
    arr[4] /= NUM_RUN;
    
    //StDev
    for(int i = 0; i < NUM_RUN; i++){
        arr[1] += pow(arr[0] - clusPercentageArr[0][i], 2);
        arr[3] += pow(arr[2] - clusPercentageArr[2][i], 2);
        arr[5] += pow(arr[4] - clusPercentageArr[4][i], 2);
    }
    
    arr[1] /= (NUM_RUN - 1);
    arr[3] /= (NUM_RUN - 1);
    arr[5] /= (NUM_RUN - 1);
    
    arr[1] = sqrt(arr[1]);
    arr[3] = sqrt(arr[3]);
    arr[5] = sqrt(arr[5]);

    
    //Note: clusPercentageArr has [1][], [3][], and [5][] memory locations unused since stand devitating after all runs
    
    //Print to File
    if(read.is_open()){
        std::string line;
        while(getline(read, line))
            lineVec.push_back(line);
        read.close();
    }else{
        printf("clusterPercentageData.dat file not found. Creating one.\n");
        
        write << "Columns 0 & 1 are system size and fiber length.\n"
            << "0 = Perc Disk on Perc Fiber vs Perc Disk\n"
            << "1 = SD of 0\n"
            << "2 = Perc Disk vs Disk\n"
            << "3 = SD of 2\n"
            << "4 = Perc Fiber vs Fiber\n"
            << "5 = SD of 4\n"
            << "L\tfibLen\ttrials\tclus0 \t\tclus1 \t\tclus2 \t\tclus3 \t\tclus4 \t\tclus5\n";
    }
    write.open(fileName);
    for(auto const& value: lineVec)
        write << value << std::endl;
    
    write << L <<",  "<< fiber_length <<",   "<< NUM_RUN <<"\t"<< arr[0] <<",\t"<< arr[1]
                <<",\t"<< arr[2] <<",\t"<< arr[3] <<",\t"<< arr[4] <<",\t"<< arr[5] << "\n";

    write.close();
}
void diskChronoDataForImport(std::vector<Point> vec){
    //Debugging Code
    
    //We run a normal case and use this function to create a file of fibers' x, y, and theta values.
    //We then use that file in the importLines() code to recreate the original simulation.
    
    std::string fileName = "DisksChronoForImport.dat";
    
    std::vector <std::string> lineVec;
    std::ofstream write;
    write.open(fileName);
    write << std::setprecision(8);
    
    if(write.is_open()){
        write << "\n\n" << L << "," << L << "\n";
        for(int i = 0; i < vec.size(); i++)
            write << vec[i].getCoorX() << "," << vec[i].getCoorY() << "\n";
    }
    write.close();
}
void lineChronoDataForImport(std::vector<Line> vec){
    //Debugging Code
    
    //We run a normal case and use this function to create a file of fibers' x, y, and theta values.
    //We then use that file in the importLines() code to recreate the original simulation.
    
    std::string fileName = "LinesChronoForImport.dat";
    
    std::vector <std::string> lineVec;
    std::ofstream write;
    write.open(fileName);
    write << std::setprecision(8);
    
    if(write.is_open()){
        write << "\n\n\n";
        for(int i = 0; i < vec.size(); i++)
            write << vec[i].getCoorX() << "," << vec[i].getCoorY() << "," << vec[i].getTheta() << "\n";
    }
    write.close();
}
void importLines(){
    
    float coor;
    std::string line, num;
    
    std::string fileName = "LinesChronoForImport.dat";
    
    std::vector <std::string> pointVec;
    std::ifstream lines(fileName);
    std::ofstream write;
    
    if(!lines.is_open()){
        printf("LinesChronoForImport.dat file not found.\n");
        write.open(fileName);
        write.close();
    }
    else{
        getline(lines,line);
        getline(lines,line);
        getline(lines,line);
        while(getline(lines, line)){//get each line of text
            //Each line with be "x,y,theta(radians)" from "lineDataForImport()"
            
            std::stringstream ss(line);
            
            while(getline(ss, num, ',')){ //Stores x and y and theta values in a vector
                coor = stod(num);
                coorVec.push_back(coor);
            }
        }
        lines.close();
    }
}
void fiberChronoDataForPythonCode(std::vector<Line> vec){
    // Debugging Code
    
    //This function reads and writes to a file contianing left and right points of fibers
    
    std::string fileName = "FiberChronoDataForPythonCode.dat";
    
    std::vector <std::string> lineVec;
    std::ofstream write;
    write.open(fileName);
    write << std::setprecision(8);
    
    
    if(write.is_open()){
        write << "\n\n\n";
        for(int i = 0; i < vec.size(); i++){
            write << vec[i].getLeftX() << "," << vec[i].getLeftY() << ","
            << vec[i].getRightX() << "," << vec[i].getRightY() << ",0\n";
        }
    }
    write.close();
}

void plotFiberCoor(){
    //This function reads and writes to a file contianing left and right points of fibers
    // as well as circle coordinates if used
    
    std::string fileName = "Lines.dat";
    std::ofstream write;
    write << std::setprecision(8);
    write.open(fileName);
    
    write << "This file made from the Percolation Code stores the start and end points of microtubules.\n"
    << "We implement wrapping, truncating\n"
    << "The format is (xleft, yleft, xright, yright).\n";
    if(percFibers){
        for(int i = 0; i < N_FIBERS; i++){
            for(int j = 0; j < fiberSubCells[i].size(); j++){
                write << fiberSubCells[i][j].getLeftX() << "," << fiberSubCells[i][j].getLeftY()
                << "," << fiberSubCells[i][j].getRightX() << "," << fiberSubCells[i][j].getRightY()
                << "," << fiberSubCells[i][j].getFiberType() << "\n";
            }
        }
    }else{
        for(int i = 0; i < num_fibers; i++){
            write << fiberSubCells[0][i].getLeftX() << "," << fiberSubCells[0][i].getLeftY()
            << "," << fiberSubCells[0][i].getRightX() << "," << fiberSubCells[0][i].getRightY() << "\n";
        }
    }
    
    write.close();
}
void plotDiskCoor(){
    std::string fileName = "Points.dat";
    std::ofstream write;
    write << std::setprecision(8);
    write.open(fileName);
    
    write << "This file stores the positions of points on microtubules.\n"
    << "The next row has the dimensions of the window.  After that is coordinate pairs.\n"
    << L << "," << L << "\n";
    
    for(int i = 0; i < N; i++)
        for(int j = 0; j < subCells[i].size(); j++)
            write << subCells[i][j].getCoorX() << "," << subCells[i][j].getCoorY() << "," << subCells[i][j].getDiskType() << "\n";
    
    write.close();
}
void plotStdDeviation(float SD){
    // This function writes the fiber length, system size and the standard deviation (SD) of the percolating disks
    // to a file so we may observe fluctuations in the disk SD as a function of fiber length.
    
    std::vector <std::string> lineVector;
    std::string fileName = "DiskSDvsFiberLength.dat";
    std::ifstream read(fileName);
    std::ofstream write;
    write << std::setprecision(8);
    
    if(read.is_open()){
        std::string line;
        while(getline(read, line))
            lineVector.push_back(line);
        
        write.open(fileName);
        for(auto const& value: lineVector)
            write << value << std::endl;
        
        write << fiber_length << "\t\t" << L << "\t" << SD << "\n";
    }else{
        printf("DiskSDvsFiberLength.dat does not exist. Creating one.\n");
        write.open(fileName);
        write << "fiber length\tL\tSD\n";
        write << fiber_length << "\t\t" << L << "\t" << SD << "\n";
    }
    write.close();
}
void printToMasterLog(float time, float etaDisk, float etaTubule){
    //This function records all the valuable data of all run on CaryCode
    
    std::vector <std::string> etaFileVec;
    std::string fileName = "MasterLog.txt";
    std::ifstream read(fileName);
    std::ofstream write;
    write << std::setprecision(8);
    
    if(read.is_open()){
        std::string line;
        while(getline(read, line))
            etaFileVec.push_back(line);
        
        write.open(fileName);
        for(auto const& value: etaFileVec)
            write << value << std::endl;
        
        write << L <<"\t"<< NUM_RUN <<"\t"<< time <<"\t"<< etaDisk <<"\t"<< fiber_length << "\t\t" << etaTubule << std::endl;
    }else{
        printf("File does not exist. Creating one.\n");
        write.open(fileName);
        write << "L\tTrials\tTime\tDiskEta\tDisk SD\tLineEta\tLine SD\n";
        write << L <<"\t"<< NUM_RUN <<"\t"<< time <<"\t"<< etaDisk <<"\t"<< fiber_length << "\t\t" << etaTubule <<"\t" << std::endl;
    }
    write.close();
}
void printPercThresholds(long int n[], float SD){
    //This function prints a list of the number of disks placed to achieve percolation
    
    //    FILE *fCountData;
    std::string fileName = "Population_Data_" + std::to_string(L);
    
    if(usingFibers){
        std::string append = "_" + std::to_string(fiber_length) + "_" + std::to_string(time(0)) +".txt";
        fileName += append;
    }
    else{
        std::string append = "_" + std::to_string(time(0)) + ".txt"; //Appends length and Unix time stamp
        fileName += append;
    }

    //Write file
    std::ofstream write;
    write.open(fileName);
    
    if(write.is_open()){
        write << "Percolation threshold data using ";
        if(usingFibers){
            write << "fiber length: " << fiber_length << "\n";
            if(!percFibers)
                write << "\nFiber density: " << fiber_density << ". ";
            else
                write << "Percolating Fibers. ";
        }
        else
            write << "Percolating Disks. ";
        
        write << "Total number of cells: " << N << "\n\n";
        
        write << "Run\tNum\t\n";
        write << "---\t---\t\n";
        
        for (int j = 0; j < NUM_RUN; j++)
            write << j << "\t" << n[j] << "\n";
    }
    write.close();
}

float calcEtaDisk(long int n[]){ //type = 1 for disks, 2 for sticks
    //Calculate the eta value by using the average of all the disks or sticks used.
    
    long long int tally = 0LL;
    for(int i = 0; i < NUM_RUN; i++)
        tally += n[i];
    
    return (float)tally / NUM_RUN * PI / 4 / N;//N = L * L, r = 0.5 * 0.5
}
float calcEtaFiber(long int n[]){ //type = 1 for disks, 2 for sticks
    //Calculate the eta value by using the average of all the disks or sticks used.
    
    long long int tally = 0LL;
    for(int i = 0; i < NUM_RUN; i++)
        tally += n[i];
    
    //num_sticks * tub_length * tub_length = 5.637
    return (float)tally / NUM_RUN / N * fiber_length * fiber_length;
}
float calcStdDeviation(long int n[]){
    //Calculate the Standard Deviation of the objects required to percolate
    
    long long int total = 0LL;
    float sd = 0;
    
    for(int i = 0; i < NUM_RUN; i++)
        total += n[i];
    
    float mean = total / NUM_RUN;

    for(int i = 0; i < NUM_RUN; i++)
        sd += pow(n[i] - mean, 2);
    
    return sqrt(sd / NUM_RUN);
}
void importPoints(){
    
    float coor;
    std::string line, num;
    std::string fileName = "Points.dat";
    
    std::vector <std::string> pointVec;
    std::ifstream points(fileName);
    std::ofstream write;
    
    if(!points.is_open()){
        printf("File not found.\n");
        write.open(fileName);
        write.close();
    }
    else{
        getline(points, line);//first two lines are text
        getline(points, line);
        getline(points, line);//length L of data
        
        while(getline(points, line)){//get each line of text
            std::stringstream ss(line);
            
            while(getline(ss, num, ',')){ //Stores x and y values in a vector
                coor = stod(num);
                coorVec.push_back(coor);
            }
        }
        points.close();
    }
}
Point findRoot(Point &p, int pCell, int pCellId){
    //This function takes a point (disk) and finds it's root
    
    // p is disk who's root we do not know
    // pCell subCell disk p was found in
    // pCellId is root vector disk p was found in
    
    //We check if the disk's
//    p.print();
    
    int pRcell = p.getRCellId();
    int pRid = p.getRootId();
    int disX = p.getDisX();
    int disY = p.getDisY();
    
    if (subCells[pRcell][pRid].getCoorX() == p.getCoorX() && subCells[pRcell][pRid].getCoorY() == p.getCoorY())
        return p;
    else{
        Point a = subCells[pRcell][pRid];//let a is the old root
        Point b = findRoot(a, pRcell, pRid);//recursive, findroot of the old root
        
        //update point a's size, displacement, rootCellId, rootID
        p.setS(b.getS());
        p.setDis(a.getDisX() + disX, a.getDisY() + disY);
        p.setId(b.getRCellId(), b.getRootId());
        subCells[pCell][pCellId] = p;
        
        return b;
    }
}
void nineCells(int cId, int *array){
    //input a cellID, output an array that stores the neighbor 9 cells
    
    array[0] = (cId + N - L - 1) % N;
    array[1] = (cId + N - L) % N;
    array[2] = (cId + N - L + 1) % N;
    array[3] = (cId + N - 1) % N;
    array[4] = cId % N;
    array[5] = (cId + 1) % N;
    array[6] = (cId + L - 1) % N;
    array[7] = (cId + L) % N;
    array[8] = (cId + L + 1) % N;
    if (cId % L == 0){
        array[0] = (cId + N - 1) % N;
        array[3] = (cId + N + L -1) % N;
        array[6] = (cId + N + 2*L - 1) % N;
    }
    if ((cId + 1) % L == 0){
        array[2] = (cId + N - 2*L + 1) % N;
        array[5] = (cId + N - L + 1)% N;
        array[8] = (cId + N + 1 )% N;
    }
}
void bondvector(int j, int &dx, int &dy){
    //This function establishes a corresponding vector to the neighboring subcells
    //input a neighbor id (0->9), output the bondvector
    
    switch (j) {
        case 0:
            dx = dy = d;
            break;
        case 1:
            dx = 0;
            dy = d;
            break;
        case 2:
            dx = -d;
            dy = d;
            break;
        case 3:
            dx = d;
            dy = 0;
            break;
        case 4:
            dx = dy = 0;
            break;
        case 5:
            dx = -d;
            dy = 0;
            break;
        case 6:
            dx = d;
            dy = -d;
            break;
        case 7:
            dx = 0;
            dy = -d;
            break;
        case 8:
            dx = -d;
            dy = -d;
            break;
    }
}
//float distance(Point A, Point B){
//    //Calculate the distance between two points
//
//    float x = A.getCoorX(), y = A.getCoorY();
//    float xB = B.getCoorX(), yB = B.getCoorY();
//    float delx, dely;
//
//    if (std::abs(x - xB) < (L - std::abs(x - xB)))
//        delx = (x - xB)*(x - xB);
//    else
//        delx = (L - std::abs(x - xB))*(L - std::abs(x - xB));
//    if (std::abs(y - yB) < (L - std::abs(y - yB)))
//        dely = (y - yB)*(y - yB);
//    else
//        dely = (L - std::abs(y - yB))*(L - std::abs(y - yB));
//    return std::sqrt(delx + dely);
//}
float distance(Point &A, Point &B){
    //Calculate the distance between two points
    
    float dx = std::abs(A.getCoorX()-B.getCoorX());
    float dy = std::abs(A.getCoorY()-B.getCoorY());

    float delx, dely;
    
    if (dx < (L - dx))
        delx = dx * dx;
    else
        delx = (L - dx)*(L - dx);
    if (delx > d)
        return 1.1;//if x length is too long just end function
    
    if (dy < (L - dy))
        dely = dy * dy;
    else
        dely = (L - dy)*(L - dy);
    return std::sqrt(delx + dely);
}

bool detectWrapping(int j, Point a, Point b){
    //Detects whether wrapping has occurred or not
    
    int oldDx = a.getDisX(), oldDy = a.getDisY();//the path through it own root
    int newDx, newDy;
    int disx, disy;
    
    bondvector(j, disx, disy);
    newDx = b.getDisX() - disx;
    newDy = b.getDisY() - disy;
    return ((std::abs(oldDx-newDx))==L||(std::abs(oldDy-newDy)==L));
}
bool detectWrapping(int j, Line a, Line b, int disx, int disy){
    //Detects whether wrapping has occurred or not
    
    int oldDx = a.getDisX(), oldDy = a.getDisY();//the path through it own root
//    int disx, disy;
    
//    bondvector(j, disx, disy);
    
//    if(j < 9)
//        bondvector(j, disX, disY);
//    else{
//        disX = bondVecArrDisX[j-9];
//        disY = bondVecArrDisY[j-9];
//    }
    
    int newDx = b.getDisX() - disx;
    int newDy = b.getDisY() - disy;
//    printf("%i, %i\n",std::abs(oldDx-newDx),std::abs(oldDy-newDy));
    return ((std::abs(oldDx-newDx))==L_FIBERS||(std::abs(oldDy-newDy)==L_FIBERS));
    //BEWARE OF L_fibers
}

void setFibers(){
    /*
     Purpose: Initialize the x, y, and theta values for all the fibers that will be created before percolation process.
     
     (1) Create random number generator
     (2) For all fibers to be created (num_fibers), create values for x, y, and theta from the x-axis
     (3) Create a fiber (Line) and push onto the first subCell in fiberSubCells
     --> this last step allows us to randomly select a fiber from the 1st slot in the vector. No need for
     assigning fibers to appropriate subCells since fibers are made on idea percolation is occurring.
     */
    
    std::random_device randDev;
    std::seed_seq seed{randDev()};
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> dis(0,L);
    
    for(int i = 0; i < num_fibers; i++){
        float x, y, theta;
        x = dis(gen);
        y = dis(gen);
        theta = dtheta * (dis(gen) / L * 2 - 1);
        
        Line A(x, y, theta, fiber_length);
        fiberSubCells[0].push_back(A);
        //        printf("i: %i\n", i);
        //        printf("fiber: %f\n", fiberSubCells[0][i].getCoorX());
        
    }
}
bool onSegment(float x1, float y1, float x2, float y2, float x3, float y3){
    if(x2 <= fmax(x1, x3) && x2 >= fmin(x1, x3) &&
       x2 <= fmax(x1, x3) && x2 >= fmin(x1, x3))
        return true;
    return false;
}
int orientation(float x1, float y1, float x2, float y2, float x3, float y3){
    //Use points 1, 2, 3 to find CW/CCW or linear orientation
    float val = (y2 - y1)*(x3 - x2) - (x2 - x1)*(y3 - y2);
    
//    if(val == 0) return 0;//collinear
    
    return (val > 0)? 1: 2; //Clock and counterclock wise
}
bool intersect(Line A, Line B){
    //Returns true if two lines intersect, false if not.
    float Ax1 = A.getLeftX();
    float Ay1 = A.getLeftY();
    float Ax2 = A.getRightX();
    float Ay2 = A.getRightY();
    float Bx1 = B.getLeftX();
    float By1 = B.getLeftY();
    float Bx2 = B.getRightX();
    float By2 = B.getRightY();
    
    int o1 = orientation(Ax1, Ay1, Ax2, Ay2, Bx1, By1);
    int o2 = orientation(Ax1, Ay1, Ax2, Ay2, Bx2, By2);
    int o3, o4;
    //General Intersection case
    if(o1 != o2){
        o3 = orientation(Bx1, By1, Bx2, By2, Ax1, Ay1);
        o4 = orientation(Bx1, By1, Bx2, By2, Ax2, Ay2);
        if(o3 != o4)
            return true;
    }
    /*
     //after 5E6 sticks placed, onSegment was never called. If one had been called, it would be 2E-5 % change. Our bounds can take that just fine. So we can ignore onSegment and dave cpu time by using if/else statement for orientation()
     //Colinear cases
     if(o1 == 0 && onSegment(Ax1, Ay1, Bx1, By1, Ax2, Ay2))
     return true;
     if(o2 == 0 && onSegment(Ax1, Ay1, Bx2, By2, Ax2, Ay2))
     return true;
     if(o1 != o2){
     if(o3 == 0 && onSegment(Bx1, By1, Ax1, Ay1, Bx2, By2))
     return true;
     if(o4 == 0 && onSegment(Bx1, By1, Ax2, Ay2, Bx2, By2))
     return true;
     }//*/
    
    //FYI fibers rarely encounter intersections
    
    //For edge pieces: execute same code again but with new coordinates
    bool changed = false;
//    printf("edge?\n");

    if(abs(B.getCoorX() - A.getCoorX()) > L - 2 * fiber_length){
        if(A.getCoorX() > B.getCoorX()){
            Bx1 += L;
            Bx2 += L;
        }else{
            Ax1 += L;
            Ax2 += L;
        }
        changed = true;
    }
    
    if(abs(B.getCoorY() - A.getCoorY()) > L - 2 * fiber_length){
        if(A.getCoorY() > B.getCoorY()){
            By1 += L;
            By2 += L;
        }else{
            Ay1 += L;
            Ay2 += L;
        }
        changed = true;
    }
    
    if(changed){
        o1 = orientation(Ax1, Ay1, Ax2, Ay2, Bx1, By1);
        o2 = orientation(Ax1, Ay1, Ax2, Ay2, Bx2, By2);
        if(o1 != o2){
            o3 = orientation(Bx1, By1, Bx2, By2, Ax1, Ay1);
            o4 = orientation(Bx1, By1, Bx2, By2, Ax2, Ay2);
            if(o3 != o4)
                return true;
        }
        /*
         if(o1 == 0 && onSegment(Ax1, Ay1, Bx1, By1, Ax2, Ay2))
         return true;
         if(o2 == 0 && onSegment(Ax1, Ay1, Bx2, By2, Ax2, Ay2))
         return true;
         if(o3 == 0 && onSegment(Bx1, By1, Ax1, Ay1, Bx2, By2))
         return true;
         if(o4 == 0 && onSegment(Bx1, By1, Ax2, Ay2, Bx2, By2))
         return true;
         //*/
    }
//    printf("no intersect\n");
    
    return false;
}
Line findRoot(Line &p, int pCell, int pCellId){
    //This function takes a point (a disk) and finds it's root
    
    int pRcell = p.getRCellId();
    int pRid = p.getRootId();
    int disX = p.getDisX();
    int disY = p.getDisY();
    
    if (fiberSubCells[pRcell][pRid].getCoorX() == p.getCoorX() && fiberSubCells[pRcell][pRid].getCoorY() == p.getCoorY())
        return p;
    else{
        Line a = fiberSubCells[pRcell][pRid];//let a is the old root
        Line b = findRoot(a, pRcell, pRid);//recursive, findroot of the old root
        
        //update point a's size, displacement, rootCellId, rootID
        p.setS(b.getS());
        p.setDis(a.getDisX() + disX, a.getDisY() + disY);
        p.setId(b.getRCellId(), b.getRootId());
        fiberSubCells[pCell][pCellId] = p;
        
        return b;
    }
}

void nineFiberCells(int cId, int *array){
    //input a cellID, output an array that stores the neighbor 9 stick cells
    array[0] = (cId + N_FIBERS - L_FIBERS - 1) % N_FIBERS;
    array[1] = (cId + N_FIBERS - L_FIBERS) % N_FIBERS;
    array[2] = (cId + N_FIBERS - L_FIBERS + 1) % N_FIBERS;
    array[3] = (cId + N_FIBERS - 1) % N_FIBERS;
    array[4] = cId % N_FIBERS;
    array[5] = (cId + 1) % N_FIBERS;
    array[6] = (cId + L_FIBERS - 1) % N_FIBERS;
    array[7] = (cId + L_FIBERS) % N_FIBERS;
    array[8] = (cId + L_FIBERS + 1) % N_FIBERS;
    if (cId % L_FIBERS == 0){
        array[0] = (cId + N_FIBERS - 1) % N_FIBERS;
        array[3] = (cId + N_FIBERS + L_FIBERS -1) % N_FIBERS;
        array[6] = (cId + N_FIBERS + 2 * L_FIBERS - 1) % N_FIBERS;
    }
    if ((cId + 1) % L_FIBERS == 0){
        array[2] = (cId + N_FIBERS - 2 * L_FIBERS + 1) % N_FIBERS;
        array[5] = (cId + N_FIBERS - L_FIBERS + 1)% N_FIBERS;
        array[8] = (cId + N_FIBERS + 1 )% N_FIBERS;
    }
}
int getPointOnFiber(float &x, float &y, int count){
    
    int t = x / L * num_fibers; //Represents the t-numbered fiber to be selected
    float s = (y / L - 0.5) * fiber_length; //pick a random distance from center of line
    
    if(percFibers){
        x = singleDimensionFiberArr[t].getCoorX() + s * cos(singleDimensionFiberArr[t].getTheta());
        y = singleDimensionFiberArr[t].getCoorY() + s * sin(singleDimensionFiberArr[t].getTheta());
        singleDimensionFiberArr[t].incrementCount();
    }else{
        //When not use fiber percolation, but already created num_fibers and am plotting disks on to them randomly.
        x = fiberSubCells[0][t].getCoorX() + s * cos(fiberSubCells[0][t].getTheta());
        y = fiberSubCells[0][t].getCoorY() + s * sin(fiberSubCells[0][t].getTheta());
        fiberSubCells[0][t].incrementCount();
    }
    
    //Make sure 0 < x, y < L
    if(x < 0)
        x += L;
    else if(x >= L)
        x -= L;
    if(y < 0)
        y += L;
    else if(y >= L)
        y -= L;
    return t;
}
int populatedFibers(){
    //For each trial, store the number of tubules containing at least one disk
    
    int tally = 0;
    if(percFibers){
        for(int i = 0; i < N_FIBERS; i++)
            for(int j = 0; j < fiberSubCells[i].size(); j++)
                if(fiberSubCells[i][j].getCount() > 0)
                    //                    printf("count: %i\n", fiberSubCells[i][j].getCount());
                    tally++;
    }
    else
        for(int i = 0; i < num_fibers; i++)
            if(fiberSubCells[0][i].getCount() > 0)
                tally++;
    return tally;
}

void findClustersDisks(int rCell, int rId){

    //Just Disks
    for(int i = 0; i < N; i++)
        for(int j = 0; j < subCells[i].size(); j++){
            findRoot(subCells[i][j], i, j);
            
            if(subCells[i][j].getRCellId() == rCell && subCells[i][j].getRootId() == rId)
                subCells[i][j].setDiskType(1);//Percolating cluster
            else
                subCells[i][j].setDiskType(0);//All other clusters
            
            //Disks and Fibers
//            if(usingFibers){
//                int fiber = subCells[i][j].getFiberNum();
//                Line A = singleDimensionFiberArr[fiber];
//                
//                if(A.getFiberType() == 1 && subCells[i][j].getDiskType() == 1){//Perc Disk on Perc Fiber
//                    subCells[i][j].setDiskType(2);
//                }
//            }
        }
    
    for(int i = 0; i < N; i++)
        for(int j = 0; j < subCells[i].size(); j++){
            if(usingFibers){
                int fiber = subCells[i][j].getFiberNum();
                Line A = singleDimensionFiberArr[fiber];

                if(A.getFiberType() == 1 && subCells[i][j].getDiskType() == 1)//Perc Disk on Perc Fiber
                    subCells[i][j].setDiskType(2);
            }
        }
}
void findClustersFibers(int rCell, int rId){
    //find all roots
    for(int i = 0; i < N_FIBERS; i++)
        for(int j = 0; j < fiberSubCells[i].size(); j++){
            findRoot(fiberSubCells[i][j], i, j);
            
            if(fiberSubCells[i][j].getRCellId() == rCell && fiberSubCells[i][j].getRootId() == rId)
                fiberSubCells[i][j].setFiberType(1);//Percolating cluster
            else
                fiberSubCells[i][j].setFiberType(0);//All other clusters
        }
}

void printAllDisks(){
    printf("rCellId\taId\tx\t\ty\t\tsize\tdx\tdy\tdiskType\tfiberNum\n");
    for(int i = 0; i < N; i++){
        for(int j = 0; j < subCells[i].size(); j++){
            subCells[i][j].print();
        }
    }
}
void printAllFibers(){
    printf("rCellId\taId\tx\t\ty\t\tsize\tdx\tdy\tpercClust\n");
    for(int i = 0; i < N_FIBERS; i++)
        for(int j = 0; j < fiberSubCells[i].size(); j++)
            fiberSubCells[i][j].print();
}

int percolateCircles(){
    
    bool wrap = 0;
    int count = 0;
    float x, y;
    
    //if DEBUGGING
//    std::vector<Point> diskChronoVec;

    
    if(importData)
        importPoints();
    
    //Generate PRNG and distribution
    std::random_device randDev;
    std::seed_seq seed{ randDev() };
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> dis(0, L); // [0,L)
    
    do
    {
        if(DEBUGGING){
            if(count > 10 * ETA * N * 4 /PI){ // If making far too many disks...
                //print list of disk and fibers
//                diskChronoDataForImport(diskChronoVec);
                printf("Lines: %i, Printing a line just to do set a breakpoint on this line.\n", num_fibers);
                exit(0);
            }
        }
        count++;
        int ptr[9];
        x = dis(gen);
        y = dis(gen);
        
        if(importData){
            if(2*count <= coorVec.size()){
                x = coorVec[2*(count - 1)];
                y = coorVec[2*count - 1];
            }else if(2*count <= coorVec.size() + 1){
                x = 0;
                y = 0;
            }else{
                printf("Trying to import points but not working.\n");
                exit(0);
            }
        }
        
        //        std::vector<std::vector<int>> vec(L);
        //        vec[0].push_back(8);
        //        printf("vec[0]: %i\n", vec[0][0]);
        
//        printf("count: %i ---------------------\n", count);
        /*
         switch (count)
         {
         case 1:
                 x = 1;
                 y = 1;
                 break;
         case 2:
                 x = 1.75;
                 y = 1;
                 break;
         case 3:
                 x = 2.5;
                 y = 1;
                 break;
         case 4:
                 x = 3.25;
                 y = 1;
                 break;
         case 5:
                x = 4;
                y = 1;
                break;
         case 6:
                x = 4.75;
                y = 1;
                break;
         case 7:
                 x = 5.5;
                 y = 1;
                 exit(0);
                 break;
         }
         //*/
        
        //        printf("count: %i\n", count);

        int fibNum = -1;
        if(usingFibers)
            fibNum = getPointOnFiber(x, y, count);//Use x and y = dis(gen) as seed for tubule disk coordinate
        
        int subC = (int)(x) + (int)(y) * L; //find the subcell which store x, y coordinate
        int aId = (int)subCells[subC].size();
        
        //Initialize Disk A
        Point A(x, y, subC, aId);
        
        if(usingFibers)
            A.setFiberNum(fibNum);
        
//        printf("%f,%f\n", x, y);
//        printf("x: %f, y: %f, count: %i\n", x, y, count);
//        printf("subC: %i, aId: %i\n", subC, aId);
        
        nineCells(subC, ptr);//find the 9 explicit subcells that are around A

        //Place A's contents into array-vector
        subCells[subC].push_back(A);//store Point A at the end of the subCells[s]
//        diskChronoVec.push_back(A);

        
        int j = 0;
        
        while(wrap == 0 && j < 9){//go through each cell of those 9 subcells
            while(subCells[ptr[j]].size() == 0 && j < 8)//Still executes while loop at j == 8, can't figure out
                j++;                                    // why program bugs out when tackling this inefficiency
            
            unsigned int i = 0;
            int k = (int)subCells[ptr[j]].size();
            
            //IS THIS BULLSHIT?
            if(j == 4)
                k--;//Saves us from comparing Point A to itself below
            
            //check all points in each of the subcell
            while (wrap == 0 && i < k) {
                
                int disX, disY;//get bond vector between point A and the neighbor subCell
                bondvector(j, disX, disY);
                
                Point B = subCells[ptr[j]][i];
                
                //If two points are amalagamable
                if (distance(A, B) <= d) {
                    
                    Point rB = findRoot(B,ptr[j],i);

                    Point rA = findRoot(subCells[subC][aId], subC, (int)aId);

                    int subCellrB = rB.getRCellId(), IdrB = rB.getRootId();
                    
                    if ((rA.getRCellId() != rB.getRCellId()) || (rA.getRootId() != rB.getRootId())) {//if A and B have dif roots
                        if (rA.getS() >= rB.getS()) {//compare negative size==> neighbor's cluster is larger or equal A's cluster
                            int newSize = rB.getS() + rA.getS();
                            int cellIDrootA = rA.getRCellId();
                            int rIDrootA = rA.getRootId();
                            rB.setS(newSize);
                            rA.setS(newSize);
                            
                            int newRCellId = rB.getRCellId();
                            int newRootId = rB.getRootId();//?
                            //update Point A since now A root is R
                            A.setId(newRCellId, newRootId);	//A's new root
                            rA.setId(newRCellId, newRootId);//change the old root of A to the new root
                            
                            int oldDisAx = A.getDisX();
                            int oldDisAy = A.getDisY();
                            int displacementAx = B.getDisX() - disX;
                            int displacementAy = B.getDisY() - disY;
                            A.setDis(displacementAx,displacementAy);//displacement to the newRoot
                            rA.setDis(displacementAx - oldDisAx, displacementAy - oldDisAy);//update displacement of the A's old root
                            
                            //update these point rA, A, B, rB back to the subCells array
                            subCells[cellIDrootA][rIDrootA] = rA;
                            subCells[B.getRCellId()][B.getRootId()] = rB;
                            subCells[subC][aId] = A;
                            subCells[ptr[j]][i] = B;
                            subCells[newRCellId][newRootId].setDis(0, 0);//since it is a new root, dx=dy=0
                            subCells[newRCellId][newRootId].setS(newSize);
                        }
                        else {//A is a larger cluster
                            //Point RootA = findRoot(A,subC,aId);
                            int newSize = rB.getS() + A.getS();
                            rB.setS(newSize);
                            rA.setS(newSize);
                            
                            //update rootID and displacement for point B and its root
                            int newRCellId = A.getRCellId();
                            int newRootId = A.getRootId();
                            B.setId(newRCellId, newRootId);
                            rB.setId(newRCellId, newRootId);
                            int oldBdx = B.getDisX();
                            int oldBdy = B.getDisY();
                            int newBdx = A.getDisX() + disX;
                            int newBdy = A.getDisY() + disY;
                            B.setDis(newBdx, newBdy);
                            rB.setDis(newBdx - oldBdx, newBdy - oldBdy);
                            
                            //Update point A, B, R back to the subCells array
                            subCells[subCellrB][IdrB]=rB;
                            subCells[subC][aId] = A;
                            subCells[ptr[j]][i] = B;
                        }
                    }
                    else {
                        wrap = detectWrapping(j,A, B);
                        if(wrap && printClusters)
                            findClustersDisks(rA.getRCellId(), rA.getRootId());
                    }
                }
                i++;
            }
            j++;
        }
    } while (!wrap);
    
    return count;
}
int percolateFibers(){
    
    bool wrap = 0;
    int count = 0;
    float x, y, theta;
    
    if(importData)
        importLines();
    
//    std::vector<Line> lineChronoVec;
    
    //Generate PRNG and distribution
    std::random_device randDev;
    std::seed_seq seed{ randDev() };
    std::mt19937 gen(seed);
    std::uniform_real_distribution<float> dis(0, L); // [0,L)
    
    do
    {
        //Increase count by 1 since we are as of now adding 1 fiber
        count++;
        
        // 16 stickSubCells to check instead of 9 because L_FIBERS % fiber_length is not always equal
        //      to a whole number so the edge cells are less than fiber_length in width/height.
        int ptr[16];
        
        // Generate x and y coordinates and theta angle for the fiber
        x = dis(gen);
        y = dis(gen);
        theta = 2 * (dis(gen) / L - 0.5) * dtheta;
        
        /*
        if(importData){
            int num = (int)coorVec.size();//Perhaps I need an Xcode update because this line and the next are a must
         
            if(3*count <= num){
                x = coorVec[3*count - 3];
                y = coorVec[3*count - 2];
                theta = coorVec[3*count - 1];
            }else if(3*count <= num + 3){
                x = 0;
                y = 0;
                theta = 0;
                printf("extra linedata being imported?\n");
            }else{
                printf("Trying to import lines but not working.\n");
                exit(0);
            }
        }
        //*/
        
        /*
         switch (count)
         {
         case 1:
             //left
             x = 2.9476254;
             y = 3.3958137;
             theta = 0.19234024;
             break;
         case 2:
             //bottom
             x = 3.7387991;
             y = 3.4841223;
             theta = -1.0173935;
             break;
         }
         //*/
//        printf("count: %i\n", count);
        
        int subC = (int)(x / fiber_length) + (int)(y / fiber_length) * L_FIBERS;
        int aId = (int)fiberSubCells[subC].size();
        
        Line A(x, y, theta, fiber_length);
        
        A.setId(subC, aId);
        
        fiberSubCells[subC].push_back(A);
        nineFiberCells(subC, ptr);//find the 9 explicit subcells that are around A
        
//        lineChronoVec.push_back(A);
        
        if(DEBUGGING){
            //Just didn't want to print these the last time I debugged
            //            printf("\ncount: %i, %f, %f, %f \n", count, x, y, theta);
            //            printf("subC: %i, aId: %i\n", subC, aId);
            
//            fiberChronoDataForPythonCode(lineChronoVec);
//            lineChronoDataForImport(lineChronoVec);
            
            if(count > int(5 * 5.6372858 * N / fiber_length / fiber_length)){
                printf("Too many sticks! count: %i\n", count);
                
                if(!importData){
//                    fiberChronoDataForPythonCode(lineChronoVec);
//                    lineChronoDataForImport(lineChronoVec);
                }
                exit(0);
            }
        }
        
        int cellsToCheck = 9; //16=corner, 12=edge, 9=neither
        int j = 0;
        
        int bondVecArrDisX[7];
        int bondVecArrDisY[7];
        
        //If last edges of L_FIBERS X L_FIBERS grid are not 1x1 squares
        if(L/fiber_length != L_FIBERS){
            
            //corners
            //top right corner
            if(subC == N_FIBERS - L_FIBERS - 2){//L_FIBERS * (L_FIBERS - 1) - 2
                cellsToCheck = 16;
                
                //above
                ptr[9] = L_FIBERS - 3;
                ptr[10] = L_FIBERS - 2;
                ptr[11] = L_FIBERS - 1;
                
                //right
                ptr[12] = subC - L_FIBERS - L_FIBERS + 2;
                ptr[13] = subC - L_FIBERS + 2;
                ptr[14] = subC + 2;
                
                //Corner
                ptr[15] = 0;
                
                //First 3 above, 3 to right, corner last
                int arrOfDisX[] = {1,0,-1,-2,-2,-2,-2};
                int arrOfDisY[] = {-2,-2,-2,1,0,-1,-2};
                
                for(int i = 0; i < sizeof(arrOfDisX)/sizeof(arrOfDisX[0]); i++){
                    bondVecArrDisX[i] = arrOfDisX[i];
                    bondVecArrDisY[i] = arrOfDisY[i];
                    //                    printf("bondvec: %i\n", bondVecArrDisX[i+9]);
                }
            }
            //top left corner
            else if(subC == N_FIBERS - L_FIBERS - L_FIBERS){//L_FIBERS * (L_FIBERS - 2)
                cellsToCheck = 16;
                
                //above
                ptr[9] = L_FIBERS - 1;
                ptr[10] = 0;
                ptr[11] = 1;
                //left
                ptr[12] = subC + L_FIBERS + L_FIBERS - 2;
                ptr[13] = subC + L_FIBERS - 2;
                ptr[14] = subC - 2;
                
                ptr[15] = L_FIBERS - 2;
                
                //First 3 above, 3 left, corner last
                int arrOfDisX[] = {1,0,-1,2,2,2,2};
                int arrOfDisY[] = {-2,-2,-2,1,0,-1,-2};
                
                for(int i = 0; i < sizeof(arrOfDisX)/sizeof(arrOfDisX[0]); i++){
                    bondVecArrDisX[i] = arrOfDisX[i];
                    bondVecArrDisY[i] = arrOfDisY[i];
                }
            }
            //bottom left corner
            else if(subC == 0){
                //below
                ptr[9] = N_FIBERS - L_FIBERS - 1;//L_FIBERS * (L_FIBERS - 1) - 1
                ptr[10] = N_FIBERS - L_FIBERS - L_FIBERS;//L_FIBERS * (L_FIBERS - 2)
                ptr[11] = N_FIBERS - L_FIBERS - L_FIBERS + 1;//L_FIBERS * (L_FIBERS - 2) + 1
                //left
                ptr[12] = L_FIBERS + L_FIBERS - 2;
                ptr[13] = L_FIBERS - 2;
                ptr[14] = N_FIBERS - 2;
                
                ptr[15] = N_FIBERS - L_FIBERS - 2;//L_FIBERS * (L_FIBERS - 1) - 2
                cellsToCheck = 16;
                
                //First 3 below, 3 left, corner last
                int arrOfDisX[] = {1,0,-1,2,2,2,2};
                int arrOfDisY[] = {2,2,2,-1,0,1,2};
                
                for(int i = 0; i < sizeof(arrOfDisX)/sizeof(arrOfDisX[0]); i++){
                    bondVecArrDisX[i] = arrOfDisX[i];
                    bondVecArrDisY[i] = arrOfDisY[i];
                }
            }
            //bottom right corner
            else if(subC == L_FIBERS - 2){
                //below
                ptr[9] = N_FIBERS - L_FIBERS - 3;//L_FIBERS * (L_FIBERS - 1) - 3
                ptr[10] = N_FIBERS - L_FIBERS - 2;//L_FIBERS * (L_FIBERS - 1) - 2
                ptr[11] = N_FIBERS - L_FIBERS - 1;//L_FIBERS * (L_FIBERS - 1) - 1
                //right
                ptr[12] = L_FIBERS;
                ptr[13] = 0;
                ptr[14] = N_FIBERS - L_FIBERS;//L_FIBERS * (L_FIBERS - 1)
                
                ptr[15] = N_FIBERS - L_FIBERS - L_FIBERS;//L_FIBERS * (L_FIBERS - 2)
                cellsToCheck = 16;
                
                //First 3 below, 3 to right, corner last of corresponding dx and dy values
                int arrOfDisX[] = {1,0,-1,-2,-2,-2,-2};
                int arrOfDisY[] = {2,2,2,1,0,-1,2};
                
                for(int i = 0; i < sizeof(arrOfDisX)/sizeof(arrOfDisX[0]); i++){
                    bondVecArrDisX[i] = arrOfDisX[i];
                    bondVecArrDisY[i] = arrOfDisY[i];
                }
            }
            
            //edges
            else if((subC + 2) % L_FIBERS == 0 && subC < N_FIBERS - L_FIBERS){
                //column to right
                ptr[9] = subC - L_FIBERS - L_FIBERS + 2;
                ptr[10] = subC - L_FIBERS + 2;
                ptr[11] = subC + 2;
                cellsToCheck = 12;
                
                //going bottom to top
                int arrOfDisX[] = {-2,-2,-2};
                int arrOfDisY[] = {1,0,-1};
                
                for(int i = 0; i < sizeof(arrOfDisX)/sizeof(arrOfDisX[0]); i++){
                    bondVecArrDisX[i] = arrOfDisX[i];
                    bondVecArrDisY[i] = arrOfDisY[i];
                }
            }
            else if(subC % L_FIBERS == 0 && subC < N_FIBERS - L_FIBERS){
                //column to left
                ptr[9] = subC + L_FIBERS + L_FIBERS - 2;
                ptr[10] = subC + L_FIBERS - 2;
                ptr[11] = subC - 2;
                cellsToCheck = 12;
                
                //going top to bottom
                int arrOfDisX[] = {2,2,2};
                int arrOfDisY[] = {-1,0,1};
                
                for(int i = 0; i < sizeof(arrOfDisX)/sizeof(arrOfDisX[0]); i++){
                    bondVecArrDisX[i] = arrOfDisX[i];
                    bondVecArrDisY[i] = arrOfDisY[i];
                }
            }
            else if(subC >= N_FIBERS - L_FIBERS - L_FIBERS && subC < N_FIBERS - L_FIBERS - 1){
                //row above
                ptr[9] = subC - (N_FIBERS - L_FIBERS - L_FIBERS) - 1;//L_FIBERS * (L_FIBERS - 2) - 1
                ptr[10] = subC - (N_FIBERS - L_FIBERS - L_FIBERS);
                ptr[11] = subC - (N_FIBERS - L_FIBERS - L_FIBERS) + 1;
                cellsToCheck = 12;
                
                //going left to right
                int arrOfDisX[] = {1,0,-1};
                int arrOfDisY[] = {-2,-2,-2};
                
                for(int i = 0; i < sizeof(arrOfDisX)/sizeof(arrOfDisX[0]); i++){
                    bondVecArrDisX[i] = arrOfDisX[i];
                    bondVecArrDisY[i] = arrOfDisY[i];
                }
            }
            else if(subC < L_FIBERS - 2){
                //row below
                ptr[9] = subC + N_FIBERS - L_FIBERS - L_FIBERS - 1;
                ptr[10] = subC + N_FIBERS - L_FIBERS - L_FIBERS;
                ptr[11] = subC + N_FIBERS - L_FIBERS - L_FIBERS + 1;
                cellsToCheck = 12;
                
                //going left to right
                int arrOfDisX[] = {1,0,-1};
                int arrOfDisY[] = {2,2,2};
                
                for(int i = 0; i < sizeof(arrOfDisX)/sizeof(arrOfDisX[0]); i++){
                    bondVecArrDisX[i] = arrOfDisX[i];
                    bondVecArrDisY[i] = arrOfDisY[i];
                }
            }
        }
        
        
        while(wrap==0 && j < cellsToCheck){//go through each cell of those 9 subcells
            
            //Next line has us skip variable assignments if no neighbor exists
            while(fiberSubCells[ptr[j]].size() == 0 && j < 8)//Still executes while loop at j == 8, can't figure out
                j++;                                    // why program bugs out when tackling this inefficiency
            
            unsigned int i = 0;//pos vector of subCell list
            int k = (int)fiberSubCells[ptr[j]].size();
            if(j == 4)
                k--;//Saves us from comparing Point A to itself below

            //check all points in each of the subcell
            while (wrap == 0 && i < k) {
                
                int disX, disY;//get bond vector between point A and the neighbor subCell
                // Here things get a little confusing. Since we assign values of j to certain neighboring cells,
                // these values of j are not consistent with an extra need to check any sides or corners.
                // These must be done by manually in each case as defined above.
                
                //                NEED TO ADD THE ABOVE
                
                
                bondvector(j, disX, disY);
                
                if(j < 9)
                    bondvector(j, disX, disY);
                else{
                    disX = bondVecArrDisX[j-9];
                    disY = bondVecArrDisY[j-9];
                }
                
                Line B = fiberSubCells[ptr[j]][i];
                
                //If two sticks are amalagamable)
                if(intersect(A, B)){
                    
                    Line rB = findRoot(B,ptr[j],i);
                    Line rA = findRoot(A,subC, aId);
                    int subCellrB = rB.getRCellId(), IdrB = rB.getRootId();
                    
                    if ((A.getRCellId() != rB.getRCellId()) || (A.getRootId() != rB.getRootId())) {
                        //if A and B have dif roots
                        if (rA.getS() >= rB.getS()) {//compare negative size==> neighbor's cluster is larger or equal A's cluster
                            int newSize = rB.getS() + rA.getS();
                            int cellIDrootA = rA.getRCellId();
                            int rIDrootA = rA.getRootId();
                            rB.setS(newSize);
                            rA.setS(newSize);
                            
                            int newRCellId = rB.getRCellId();
                            int newRootId = rB.getRootId();//?
                            //update Point A since now A root is R
                            A.setId(newRCellId, newRootId);	//A's new root
                            rA.setId(newRCellId, newRootId);//change the old root of A to the new root
                            
                            int oldDisAx = A.getDisX();
                            int oldDisAy = A.getDisY();
                            int displacementAx = B.getDisX() - disX;
                            int displacementAy = B.getDisY() - disY;
                            A.setDis(displacementAx,displacementAy);//displacement to the newRoot
                            rA.setDis(displacementAx - oldDisAx, displacementAy - oldDisAy);//update displacement of the A's old root
                            
                            //update these point rA, A, B, rB back to the subCells array
                            fiberSubCells[cellIDrootA][rIDrootA] = rA;
                            fiberSubCells[B.getRCellId()][B.getRootId()] = rB;
                            fiberSubCells[subC][aId] = A;
                            fiberSubCells[ptr[j]][i] = B;
                            fiberSubCells[newRCellId][newRootId].setDis(0, 0);//since it is a new root, dx=dy=0
                            fiberSubCells[newRCellId][newRootId].setS(newSize);
                        }
                        else {//A is a larger cluster
                            //Point RootA = findRoot(A,subC,aId);
                            int newSize = rB.getS() + A.getS();
                            rB.setS(newSize);
                            A.setS(newSize);
                            
                            //update rootID and displacement for point B and its root
                            int newRCellId = A.getRCellId();
                            int newRootId = A.getRootId();
                            B.setId(newRCellId, newRootId);
                            rB.setId(newRCellId, newRootId);
                            int oldBdx = B.getDisX();
                            int oldBdy = B.getDisY();
                            int newBdx = A.getDisX() + disX;
                            int newBdy = A.getDisY() + disY;
                            B.setDis(newBdx, newBdy);
                            rB.setDis(newBdx - oldBdx, newBdy - oldBdy);
                            
                            //Update point A, B, R back to the subCells array
                            fiberSubCells[subCellrB][IdrB]=rB;
                            fiberSubCells[subC][aId] = A;
                            fiberSubCells[ptr[j]][i] = B;
                        }
                    }
                    else {
                        wrap = detectWrapping(j,A, B, disX, disY);
                        if(wrap && printClusters)
                            findClustersFibers(rA.getRCellId(), rA.getRootId());
                        
                    }
                }
                i++;
            }
            j++;
            
        }
    } while (wrap == 0);
    
//    lineChronoDataForImport(lineChronoVec);
    
    return count;
}

int main()
{
    //Choose a directory that is easy to access all output files.
    chdir("/Users/AndrewFroning/Desktop/PercOutput");

    /*
     //generate L and fiber lengths to run
     float len = 2;
     
     //Run i fiber lengths
     for(int i = 1; i < 20; i++)
     for(int j = 0; j < 4; j++)
     printf("%i, %f\n", (int)pow(2, 5+j), i/len);
     //*/
    char buffer[200];
    getcwd(buffer, sizeof(buffer));
    printf("Current directory: %s\n", buffer);
    
    if(usingFibers)

        printf("L_fibers: %i, fiber_length: %f\n", L_FIBERS, fiber_length);
        if(percFibers && fiber_length < 1){
            printf("Cannot use fiber perc for < 1, too much computation.\n");
            exit(0);
        }
    
    clock_t startTime = clock();
    
    long int populationData[NUM_RUN] ,fiberPopulationData[NUM_RUN], fibersWithDisks[NUM_RUN] = {0};
    
//***********************************************************************
//                     BEGIN PERCOLATION SIMULATIONS
//***********************************************************************
    printf("L = %i, NUM_RUN = %i\n\nStarting percolation...\n", L, NUM_RUN);

    for (int i = 0; i < NUM_RUN; i++)
    {
        //The next three IF statements indicate how far along the program is.
        if(i == NUM_RUN * 1 / 4)
            printf("1/4 of trials completed at %.2fs\n", (float)(clock() - startTime)/CLOCKS_PER_SEC);
        else if(i == NUM_RUN / 2)
            printf("1/2 of trials completed at %.2fs\n", (float)(clock() - startTime)/CLOCKS_PER_SEC);
        else if(i == NUM_RUN * 3 / 4)
            printf("3/4 of trials completed at %.2fs\n", (float)(clock() - startTime)/CLOCKS_PER_SEC);
        
        //Determine what types of percolation will be run.
        if(!usingDisks && !usingFibers){
            printf("Not using either disks or fibers. Program will end.\n");
            exit(0);
        }else{
            if(usingDisks && usingFibers){
                if(percFibers){
                    //Both Disk and Fiber Percolation
                    
                    //First percolate fibers
                    fiberPopulationData[i] = num_fibers = percolateFibers();
                    
                    // if DEBUGGING
//                    plotFiberCoor();
                    
                    long int fiberIterator = 0;
                    for(int i = 0; i < N_FIBERS; i++)
                        for(int j = 0; j < fiberSubCells[i].size(); j++){
                            singleDimensionFiberArr[fiberIterator] = fiberSubCells[i][j];
                            fiberIterator++;
                        }

                    //Then percolate disks on those fibers
                    populationData[i] = percolateCircles();
                    printf("fiberData: %li, popData: %li\n", fiberPopulationData[i], populationData[i]);
                }
                else{
                    //Disk Percolation on top of Pre-determined fibers
                    
                    setFibers();
                    
                    populationData[i] = percolateCircles();
                    fibersWithDisks[i] = populatedFibers();
                    printf("popData: %li, numPopFibersPerTrial = %li\n", populationData[i], fibersWithDisks[i]);
                }
            }
            else if(!usingDisks && usingFibers){
                //Only Fiber Percolation
                
                fiberPopulationData[i] = percolateFibers();
                populationData[i] = fiberPopulationData[i]; //Function "printPercThresholds" takes "populationData" as input
                printf("fiberPopulationData = %li\n", fiberPopulationData[i]);
            }
            else if(usingDisks && !usingFibers){
                //Only Disk Percolation
                
                populationData[i] = percolateCircles();
                printf("populationData = %li\n", populationData[i]);
            }
        }
        ///////////////////////////////////
        //Print Clusters for Python Plotting
        if(printClusters/* && NUM_RUN == 1*/){
            
//            plotDiskCoor();
//            plotFiberCoor();
            
//            printAllDisks();
//            printAllFibers();
            
            clusterPercentages(populationData[i], fiberPopulationData[i], i);
//            clusterPercentageData();
        }
//        else if(printClusters)
//            printf("Change NUM_RUN = 1 to printClusters of a single run.\n");
        

        
        ////////////////////
        
        //Clear used space
        if(usingFibers)
            for (long int j = 0; j < fiberPopulationData[i]; j++)
                fiberSubCells[j].clear();
        if(usingDisks)
            for (int j = 0; j < N; j++)
                subCells[j].clear();	// This line "reset" the subCells array back to being empty
    }
    float time =(float)(clock() - startTime)/CLOCKS_PER_SEC;
    printf("\nRun time: %.3fs\n", time);
    
    if(printClusters)
        clusterPercentageData();
    
//    *******************************************************************
//                        DETERMINE RAW ETA VALUES
//    *******************************************************************
    float etaDisk = 0.0, etaFiber = 0.0, etaDiskSD = 0.0, etaFiberSD = 0.0;

    if(usingDisks && usingFibers){
        if(percFibers){
            //Both Disk and Fiber Percolation
            
            etaDisk = calcEtaDisk(populationData);
            etaDiskSD = calcStdDeviation(populationData) / N;
            etaFiber = calcEtaFiber(fiberPopulationData);
            etaFiberSD = calcStdDeviation(fiberPopulationData) / N_FIBERS;
            printf("etaDiskSD: %f +/- %f\n", etaDisk, etaDiskSD);
            printf("etaFiberSD: %f +/- %f\n", etaFiber, etaFiberSD);
        }
        else{
            //Disk Percolation on top of Pre-determined fibers
            
            etaDisk = calcEtaDisk(populationData);
            etaDiskSD = calcStdDeviation(populationData) / N;
            etaFiber = calcEtaFiber(fiberPopulationData);
            etaFiberSD = calcStdDeviation(fiberPopulationData) / N_FIBERS;
            printf("etaDiskSD: %f +/- %f\n", etaDisk, etaDiskSD);
            printf("etaFiberSD: %f +/- %f\n", etaFiber, etaFiberSD);
        }
    }
    else if(!usingDisks && usingFibers){
        //Only Fiber Percolation
        etaDisk = etaDiskSD = 0;
        etaFiber = calcEtaFiber(fiberPopulationData);
        etaFiberSD = calcStdDeviation(fiberPopulationData) / N_FIBERS;
        printf("etaFiberSD: %f +/- %f\n", etaFiber, etaFiberSD);
    }
    else if(usingDisks && !usingFibers){
        //Only Disk Percolation
        
        etaDisk = calcEtaDisk(populationData);
        etaDiskSD = calcStdDeviation(populationData) / N;
        etaFiber = etaFiberSD = 0;
        printf("etaDiskSD: %f +/- %f\n", etaDisk, etaDiskSD);
    }
    
//    *******************************************************************
//                              PRINTING FILES
//    *******************************************************************
    
    //Print the Perclation Thresholds to a file with a Unix time stamp
    printPercThresholds(populationData, etaDiskSD);

    
    //Write the density Standard Deviation to a file for plotting
//    plotStdDeviation(etaDiskSD);
    
    //Print to a Master Log of all conditions run to see program execution history
    printToMasterLog(time, etaDisk, etaFiber);
    
    return 0;
}


/*
 FILE OUTPUT:
  
 Run Num Key Site
 --- --- --------
 0   X   0
 1   X   0
 2   X   0
 3   X   0
 Prints all runs
 ///////////////////End of File Input
 
 Press any key to continue...
//////////////////////////////////////////////////////////////////////////////////////////////////
 
 4/9 mtg
 1. Fixed DistSquar Function using std::, abs() was not using the library std
 
 4/11
 1. Noticed a point in rawdata.txt held same value for x and y, but in sorteddata.txt x and y were different values. The algorithm was carried out with the different values.
  
 4/26 && 4/27
 I have modified Percolation code to import a list of points, specifically "Points.dat"
 need to
5/9 Notes:
 -20 runs, L=128, Mersenne Twister:     41.42, 39.57, 40.46,    39.71,  39.16,  39.34
 -                              etas:                       ,   1.116,  1.108,  1.113
 -20 runs, L=128, Subtract with Carry:  9.41,  8.93,  9.12,     9.18,   8.93,   9.01
 -                              etas:   x   ,  1.122, 1.1095,   1.1223,  1.111, 1.114
 -100 runs, L=256, Subtract with Carry: 184.43
 -                              etas:   1.122
 -           using tubules...    803.65s, eta: 1.1184,   diskEta: 1.768
 & 2 randevs in seed_seq:  487s, 529.25s, eta: 1.118056, diskEta: 1.766643
 & 1 randevs in seed_seq:        520.52s, eta: 1.116978, diskEta: 1.764066
  
 *************************************************
 Generating seed only once...using Linear Congruential, random points
 std::random_device randDev;
 std::seed_seq seed{ randDev() };
 
 std::minstd_rand gen(seed);//~5.5 times faster than Mersenne Twister(mt19937), Program runs ~4 times faster
 std::uniform_real_distribution<double> dis(0, L); // [0,L)
 
 // *******random generator ran once per trial
 -1000 runs, L=256, Linear Congruential: 194.99s,   199.02s,
 -                                etas:  1.121412,  1.120429
 -1000 runs, L=512, Linear Congruential: 1078.15s,  1089.75s,
 -                                etas:  1.123389,  1.124648,
 -1000 runs, L=1024, Linear Congruential:   4044.30s, 4354.99s, 3845.21s, 3779.00s, 3790.47s, 4174.49s,
 -                                etas:     1.130051, 1.129973, 1.130064, 1.129943, 1.129941, 1.129735,
 
 
 // *********Mersenne Twister once per trial, RANDOM POINTS************
 -10000 runs, L=256, Mersenne Twister:  1936.79s,
 -                               etas:  1.119938,
 -1000 runs, L=1024, Mersenne Twister:  4042.80s,   4139.15s,
 -                               etas:  1.125206,   1.125415,
 -200 runs, L=2048, Mersenne Twister:  3543.76s,   3667.49s,
 -                               etas:  1.126657,  1.126307,
 -1000 runs, L=2048, Mersenne Twister:  20293.33s,
 -                               etas:  0.322167,   <-----5/16 bug
 
 // *********Mersenne Twister once per trial, TUBULES************
 -1000 runs, L=256, Mersenne Twister:  1830.80s,  NEW TIMEs     524.70s
 -                          line etas:  1.118883,               1.120441
 -                          disk etas:  1.769219,               1.773270
 -200 runs, L=512, Mersenne Twister:   1545.31s,
 -                          line etas:  1.121229,
 -                          disk etas:  1.775181,
 -1000 runs, L=512, Mersenne Twister:   2545.14s,
 -                          line etas:  1.123473,
 -                          disk etas:  1.781332,
 Clearly Can't use anything but Mersenne Twister, as it provides the necessary random seeds
 
 
 -50 runs, L=1024, Linear Congruential: 200.95s,    201.27s,
 -                               etas:  1.130244,   1.131514
 -50 runs, L=2048, Linear Congruential: 902.69s,
 -                               etas:  1.139594
 -50 runs, L=2048, Mersenne Twister:    921.48s,
 -                               etas:  1.138019
  
 5/20 update:
 -Optimized program runtime (about 70-80 times faster) using by creating PRNG once per trial,
 and while using sticks, once for all stick characteristics (x,y,theta)
  
 TASKS:
 -increase sig figs to 7 or 8 for eta values. DONE
 
 */
