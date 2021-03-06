// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse
// Edited by: Camillo ARGUELLO
// Date:     2013/10/08 -> 2020/10

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}

/**
 * This function computes the given matches points into the Fundamental matrix
 * */
FMatrix<float,3,3> computeFundamentalMatrix(vector<Match>& matches){

    // Finding fundamental matrix
    FMatrix<float,9,9> A(-1.);
    float x1, y1, x2, y2;
    for (int i = 0; i < 8; i++) {
        Match match = matches[rand()%matches.size()];

        // Normalization of 8 given points 
        x1 = match.x1 * 10e-3; 
        y1 = match.y1 * 10e-3;
        x2 = match.x2 * 10e-3; 
        y2 = match.y2 * 10e-3;

        // Building the linear system equation
        A(i,0) = x1 * x2;
        A(i,1) = x1 * y2;
        A(i,2) = x1;
        
        A(i,3) = y1 * x2;
        A(i,4) = y1 * y2;
        A(i,5) = y1;
        
        A(i,6) = x2;
        A(i,7) = y2;
        A(i,8) = 1;
        
        // Fill of zeros
        A(8,i) = 0;
    }
    // Extra row
    A(8,8) = 0;

    // Solve linear system using svd
    FVector<float,9> S;
    FMatrix<float,9,9> U, V;

    // Compute SVD of A
    svd(A, U, S, V);

    // Extracting F out of V
    FMatrix<float,3,3> F;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            F(i,j) = V(8,3 * i + j);
        }
    }

    // Enforcing det(F) = 0 (Constraint of Rank)
    FVector<float,3> Sf;
    FMatrix<float,3,3> Uf, Vf;
    svd(F, Uf, Sf, Vf);
    Sf[2] = 0;
    F = Uf * Diagonal(Sf) * Vf;

    // Normalization for F
    FMatrix<float,3,3> N(0.0);
    N(0,0) = 10e-3; 
    N(1,1) = 10e-3; 
    N(2,2) = 1;
    F = N * F * N; 

    return F;
}

/**
 * This function calculates the epipolar Distance given a match point and the Fundamental Matrix
 * */
float epipolarDistance(Match m, FMatrix<float,3,3>& F) {
    float x1, x2, y1, y2;

    x1 = m.x1; 
    y1 = m.y1;
    x2 = m.x2; 
    y2 = m.y2;
    FloatPoint3 point, point1;
    point[0]=x1;
    point[1]=y1;
    point[2]=1;

    point1[0]=x2;
    point1[1]=y2;
    point1[2]=1;

    point = transpose(F) * point;
    float dist = abs(point1 * point);
    dist = dist / sqrt(point[0] * point[0] + point[1] * point[1]);
        
    return dist;
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {

    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    vector<int> bestInliers;

    int counter = 0;
    vector<int> inliers;
    int n = matches.size();
    cout << "Computing Fundamental Matrix..." << endl;
    while (counter < Niter){

        FMatrix<float,3,3> F = computeFundamentalMatrix(matches);

        // Get Epipolar Distance
        for (int i = 0; i < matches.size(); i++){
            if(epipolarDistance(matches[i],F) <= distMax) {
                cout << "the distance is " << epipolarDistance(matches[i],F) << endl;
                inliers.push_back(i);
            }
        }

        if (inliers.size() > bestInliers.size()){
            bestF = F;
            bestInliers = inliers;
            Niter = (int)(log(BETA) / log(1-(10e-3)-pow(float( inliers.size() ) / float(n), 8)));
        }

        counter++;
    }

    cout << "Iterations: " << counter << ", Inliers: " << inliers.size() << endl;
    
    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);

    cout << "The Fundamental Matrix was successfully calculated." << endl;
    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    cout << "Click to select a point." << endl;

    int image1 = I1.width();
    int image2 = I2.width();

    while(true) {
        int x,y;
        
        if(getMouse(x,y) == 3)
            break;

        if(getMouse(x,y) == 1) {
            IntPoint2 p0(x,y);
            drawCircle(p0, 4, GREEN, 2);
            
            DoublePoint3 p1, p2;
            FVector<float, 3> v;
            v[1] = y;
            v[2] = 1;
            IntPoint2 leftPoint, rightPoint;
            
            // LEFT IMAGE SCREEN
            if (x < image1) {
                v[0] = x;
                v = transpose(F) * v;
                v /= v[2];
                leftPoint[0] = image1;
                leftPoint[1] = (int) (-v[2] / v[1]);
                rightPoint[0] = image1 + image2;
                rightPoint[1] = (int) (-(v[2] + v[0] * image1) / v[1]);
            } else {
            // RIGHT IMAGE SCREEN
                v[0] = x - image1;
                v = F * v;
                v /= v[2];
                leftPoint[0] = 0;
                leftPoint[1] = -v[2] / v[1];
                rightPoint[0] = image1;
                rightPoint[1] = (int) (-(v[2] + v[0] * image1) / v[1]);

            }
            cout << "vector v=" << v << endl;
            drawLine(leftPoint, rightPoint, RED, 2);
        }
    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    cout << " matches: " << matches.size() << endl;
    click();
    
    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);        
    }
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
