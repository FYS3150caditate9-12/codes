#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>

#define N_MAX 10000.
#define T_MAX 100
#define EXPER 1

#define LEG 1.;

using namespace std;



double analitical_solution(double, double);
double rand(double, double);


int main(){
    //Monte Carlo in one dimension with Uniform distibution
    int n = 11;
    double histogram_uniform[n], temp[n], average[n][T_MAX];

    for(int i=0; i<n; i++)
    for(int t=0; t<T_MAX; t++)
        average[i][t]=0;

    histogram_uniform[0] = N_MAX;
    temp[0] = 0;

    int e=0;            //start N experiment on which make average
    while(e<EXPER){
          int t=0;

          for(int i=1; i<n; i++){      //Set variable equal to zero
              histogram_uniform[i] = 0;
              temp[i] =0;
          }

    int rando;

    print(histogram_uniform, n);
    while(t<T_MAX){

        for(int i=0; i<n-1; i++){                    //For each bar of the histogram
        for(int k=1; k<histogram_uniform[i]; k++ ){  //For each particle on one bar
            rando = rand() % 100;

            if(rando>= 50){
                if(i<n-2) temp[i+1]++;
                else      temp[i-1]++;
            }
        }

        for(int i=0; i<n; i++){
            histogram_uniform[i]= temp[i];           //Update the histogram
            temp[i]=0;                               //Set temporaney histogram to zero
            average[i][t] += histogram_uniform[i];   //Update average
        }

        histogram_uniform[0]= N_MAX;                 //boundary conditions

        print(histogram_uniform, n);
        t++;
        }
    e++;
    }

    cout << endl;
    for(int t=0; t<T_MAX; t+=1){
    for(int i=0; i<n; i++){
        cout <<  average[i][t]/EXPER << "  ";       //Print out hte normalizzed average of more experiments
    }
        cout << endl;
    }

    // Methods MonteCarlo one-two dimensions with normal distribution

    int dimension = 1;   //set 0 for one-dimension and 1 for two-dimension

    int N_particles=500;
    int N_delta_x=9;
    int T_total= 300 ;

    double legh=9;                        //Lehg of the interval X
    double delta_x = legh/N_delta_x;


    srand (time(NULL));

    clock_t start, finish;                  //calculate the time for each algorithm.
    start = clock();


    double istogramma[N_delta_x+2][(N_delta_x+2)][T_total];
    double istogramma_2[N_delta_x+2][(N_delta_x+2)][T_total];

    double normalizzation;
    if(dimension==0){
        normalizzation = N_particles;
    }else{
        normalizzation = pow(N_particles, 2.);
    }

    //Set the istogram to zero.
    for(int t=0; t<T_total; t++){
    for(int x=0; x<=N_delta_x+1; x++){
    for(int y=0; y<=N_delta_x+1; y++){
          istogramma[x][y][t] = 0;
        istogramma_2[x][y][t] = 0;
    }
    }
    }


    double sigma = 1/sqrt(2), density= 1./double(N_particles);

    for(int j=0; j<=(legh*N_particles)*dimension; j++) {

       for(int i=0; i< 3*N_particles; i++){

        double positio_x= -i*density; //Uniform distribution before zero;
        double positio_y=  j*density;

        if(dimension == 0)
            positio_y = 0; //to generate the istogram when we are in one dimension

        for (int t=0; t<T_total; t++) {

            if (dimension == 0){
                positio_x += sigma*randn(0, 1.);
            }else{
                double leght = sigma*randn(0,1.);
                double theta = (rand() % 1000)/1000.*2.*M_PI;
                positio_x += leght*cos(theta);
                positio_y += leght*sin(theta);
            }

            //set boudary condition
            if(positio_x <0 || positio_x >= legh){
                t = T_total+1;
            }
            if(positio_y < 0){
                 //t = T_total+1;              //boudary condition 0.
                positio_y = abs(positio_y);    //wall
            }

            if(positio_y > legh){
                 //t = T_total+1;
                positio_y = legh - (positio_y-legh);
            }


            if(t!=T_total+1){
               //syncronize istogram every Delta_t
               for(int k=0; k<=N_delta_x+1; k++){
               for(int z=0; z<=(N_delta_x+1)*dimension; z++){
                  if(positio_x>delta_x*(k) && positio_x<(k+1)*delta_x){
                  if(positio_y>=delta_x*(z) && positio_y<= (z+1)*delta_x){
                           istogramma[k][z][t] +=1; //map of the particles


                  }}
               }}
            }

    }
    }
    }

    cout << "time life of N particles" << endl;

    for(int t=0; t<T_total; t++){
        if(t%1){
        }else{
        for(int y=0; y<=(N_delta_x+1)*dimension ; y++){
        for(int x=0; x<N_delta_x+1 ; x++){
            cout  << istogramma[x][y][t] << "   ";
            }
        cout << endl;
           }
        cout << endl << endl;
        }
    }


    //Every Delta_t N particles enter in the system.
    for(int t=0; t<T_total; t++){
        for(int t_2=0; t_2<T_total-t; t_2++){

        for(int k=0; k<=N_delta_x+1; k++){
        for(int y=0; y<=(N_delta_x+1); y++){
            istogramma_2[k][y][t + t_2] += istogramma[k][y][t_2];
        }}
        }
    }

    cout << endl << "Diffusion problem:" << endl;

    for(int t=0; t<T_total; t++){
        if(t%5){

        }else{
        cout << "t= " << t << endl;
        for(int y=0; y<=(N_delta_x+1)*dimension ; y++){
            cout << "1  ";

        for(int k=0; k<=N_delta_x ; k++){

            cout << double(istogramma_2[k][y][t])/normalizzation << "  "  ;
            }

        cout << endl;

         }
        cout << endl;
        }

    }

    finish = clock();
    cout << double(finish-start)/CLOCKS_PER_SEC << " seconds to solve this algorithm" << endl;


   //STANDAR ALGORITHM TO SOLVE DIFFUSION
   int isto[N_delta_x+2][N_delta_x+2];

   vector <double>position_x;
   vector <double>position_y;

   int t=0;

   start = clock();

   while(t<T_total){
   t++;

   for(int j=0; j<=dimension*N_delta_x*N_particles;j++){
   for(int i=0; i<3*N_particles; i++){
       position_x.push_back(-i*density);             //we add particles with uniform
       position_y.push_back(j*density);              //distribution before x<0.
   }
   }

   for(int j=0; j<=(N_delta_x+1);j++){
   for(int k=0; k<=N_delta_x+1; k++){                //loop to set the istogram to zero
        isto[k][j]=0;
   }
   }

   for(int i=0; i<position_x.size(); i++){           //loop to move all the N particles

        if(dimension == 0){
            position_x[i] += randn(0,(1/sqrt(sqrt(2))));
        }else{
            double leght, theta;                     //in two dimensio we consider
            leght = randn(0,(1/sqrt(sqrt(2))));      //also the angle of the mouving
            theta = (rand() % 1000)/1000.*2.*M_PI;   //with unifrom distribution
                                                     //between 0 and 2pi
            position_x[i] += leght*cos(theta);
            position_y[i] += leght*sin(theta);

        }
   }


   for(int i=0; i<position_x.size(); i++){
      if(position_x[i] < 0 || position_x[i]>=legh || position_y[i] < 0 || position_y[i]>=legh ){

         position_x[i] = position_x[position_x.size()];   //If the particle goes outside the
         position_y[i] = position_y[position_y.size()];   //system becomes deleted.
         position_x.pop_back();                           //we invert the position of the last
         position_y.pop_back();                           //particle with the one that we want
         i--;                                             //delete and we use the function
      }                                                   //"pop.back()" to delete the last
   }                                                      //element of the vector

   for(int i=0; i<position_x.size(); i++)
     if(t==250)
        cout << position_x[i]<< "  " << position_y[i] << endl;

    //syncronize istogram
    for(int i=0; i<position_x.size(); i++){

       for(int k=0; k<=N_delta_x+1; k++){
       for(int z=0; z<=(N_delta_x+1); z++){
           if(position_x[i]>delta_x*(k-1) && position_x[i]<(k)*delta_x){
           if(position_y[i]>=delta_x*(z-1) && position_y[i]< z*delta_x){

                isto[k][z] +=1;                           //we look for each particles where
         }}                                               //is in side the istogram.
     }}
    }

    //loop to create the istogram
    for(int j=1; j<=(N_delta_x-1)*dimension +1; j++){
        cout << "1  ";
        for(int k=1; k<=N_delta_x+1; k++){
        cout << isto[k][j]/normalizzation<< "  ";       //print out the istogram
        }
        cout << endl;
    }
    cout << endl << endl;
 }

 finish = clock();
 cout << double(finish-start)/CLOCKS_PER_SEC << " seconds to solve the problem"<< endl;

 return 0;
}


double analitical_solution(double x, double t){

         double sum2=0;

         for(int n=1; n<100; n++){
                sum2 += (-2/(M_PI*n))*sin(M_PI*n*x)*exp(-M_PI*n*M_PI*n*t);
         }

return  1 - x + sum2;
}



double randn (double mu, double sigma)
{
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;

    if (call == 1)
    {
        call = !call;
        return (mu + sigma * (double) X2);
    }

    do
    {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);

    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;

    call = !call;

    return (mu + sigma * (double) X1);
}
