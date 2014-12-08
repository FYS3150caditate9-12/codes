#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <fstream>
#include <string>
#include <vector>

#include <armadillo>
#define LENGHT 11
using namespace arma;
using namespace std;

//two dimensional diffusion problem
double analitycal(double, double, double);
void baundary_conditions(double **, int);
void print(double **, int, double);
void explicit_solver(double **, double, int);

int main(){

    double **lattice = new double*[LENGHT];


    for(int i =0; i<LENGHT; i++){
        lattice[i]= new double[LENGHT];
    }

    //Initial conditons
    for(int x =0; x<LENGHT; x++)
    for(int y =0; y<LENGHT; y++)
        lattice[x][y] = 0;

    baundary_conditions(lattice, LENGHT);
    print(lattice, LENGHT, 0);

    double alpha =0.1;
    int time =0;

    while(time < 200){
          time++;
          explicit_solver(lattice, alpha, LENGHT);


          for(int x=1; x<LENGHT-1; x++){
          //  temp[x][0] = temp[x][1]; //boudary condition as a wall
          //  temp[x][LENGHT-1] = temp[x][LENGHT-2];
          }



          print(lattice, LENGHT, time*alpha*pow(1/(LENGHT+1),2.));
    }


    //Jacobi for laplace solution.
    int n=LENGHT-1;
    double h = 1/(n+1);

    mat u(n+1, n+1);
    mat u_guess(n+1, n+1);
    mat u_new(n+1,n+1);

    double a=2;

    u.zeros();
    u_new.zeros();
    u_guess.zeros();
    //boundary conditions
    for(int j =0; j<=n; j++){
        u(0,j) = 2 - j/double(n);
        u(n,j) = 1 - j/double(n);
        u(j,0) = 2 - j/double(n);
        u(j,n) = 1 - j/double(n);
    }

    u_new = u;
    u_guess = u;

    time = 0;

    while (time < 1000){
          time++;

          u = u_new;

          int iterations = 0, max_iter=1000000; double diff = 10.;

          while( (iterations <= max_iter) && ( diff > 0.00001) ){
                 diff=0;

                 for (int j = 1; j< n; j++){
                 for (int l = 1; l< n; l++){
                      u_new(j,l) = 1/(1+4*a)*(a*(u_guess(j+1,l)+u_guess(j-1,l)+ u_guess(j,l+1)+u_guess(j,l-1)) + u(j,l));

                      diff += fabs(u_new(j,l)-u_guess(j,l));
                  }}

                  iterations++;
                  diff /= pow((n),2.0);

                  u_guess = u_new;

                  //boudary condition with wall
                  for(int x=0; x<n;x++){
                     //u_guess(x,0) = u_guess(x,1);
                     //u_guess(x,n-1) = u_guess(x,n-2);
                  }



    }
    cout << "Number of iterations: "<< iterations << ", error: "<< diff << endl;
    cout << u_new << endl;


}
    return 0;
}
double analitycal(double x, double y, double t){
    double sum=0;
    for(int n=1; n<100 ; n++)
    for(int m=1; m<100 ; m++)
            sum += 1/(n*M_PI)*sin(n*M_PI*x)*1/(m*M_PI)*sin(m*M_PI*y)*exp(-M_PI*M_PI*(n*n+m*m)*t);

    return (1-x)+(1-y)-2*sum;
}

void baundary_conditions(double **lattice, int n){
    //Bourders conditions to compare with the analytical solution
    for(int j =0; j<n; j++){
        lattice[0][j] = 2 - j/double(n-1);
        lattice[n-1][j] = 1 - j/double(n-1);
        lattice[j][0] = 2 - j/double(n-1);
        lattice[j][n-1] = 1 - j/double(n-1);
    }
}


void print(double **u, int n, double t ){
    for(int i =0; i<n; i++){
    for(int j =0; j<n; j++){
         for(int x=0; x<=0; x++)
             cout << u[j][i] << "   ";
             //cout << u[j][i] - analitycal(double(j)/LENGHT,double(i)/LENGHT, t)<< "  " << endl;
         }
        cout << endl;
    }
cout << endl << endl;
}

void explicit_solver(double **lattice, double alpha, int n){
    double **temp = new double*[n];

    for(int i =0; i<n; i++){
        temp[i]= new double[n];
    }

    for(int x=1; x<n-1; x++){
    for(int y=1; y<n-1; y++){
        temp[x][y] = lattice[x][y];
        temp[x][y] += alpha*( lattice[x+1][y] + lattice[x-1][y] + lattice[x][y+1] + lattice[x][y-1] - 4*lattice[x][y]);

    }}

    for(int x=1; x<n-1; x++)
    for(int y=1; y<n-1; y++)
        lattice[x][y] = temp[x][y];
}
