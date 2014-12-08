#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <fstream>

using namespace std;

#define LEGH 1
#define TMAX 400


void tridiagonal_solver(double, double, double, double *, double *, int);
void print_3d(ofstream &, double * , int, int *, double, double);
void print_1d(ofstream &,  double * , int, double, double);
void explicit_solver(double *, int, double);
void implicit_solver(double *, int, double);
void cranic_solver(double *, double *, int, double);
double analitical_solution(double, double);
double condition_t_zero_dirichlet(double);
double dirichlet(double, double);



int main()
{
    int n = 10, print_way=1;
    double l=LEGH, delta_x;
    delta_x = l/double(n);

    double alpha = 0.4;  //minore di 0.5
    double delta_t = delta_x*delta_x*alpha;

    n++; //perchè l'array ha bis

    double y_es[n], y_im[n];
    double y_cn_1[n],y_cn[n];

    //Set boundary condition for t=0
    for(int i=0; i<n; i++){
        y_es[i] = y_im[i] =   y_cn_1[i]=  condition_t_zero_dirichlet(delta_x*double(i));
    }


    int t=0;

    //_____Printing_files_____
    char *filename = new char[1000];
    sprintf(filename, "diffusion_ex.dat");
    ofstream expli(filename);
    sprintf(filename, "diffusion_im.dat");
    ofstream impli(filename);
    sprintf(filename, "diffusion_Crank-Nicolson.dat");
    ofstream cra_nic(filename);
    sprintf(filename, "diffusion_Crank-Nicolson_3d.dat");
    ofstream cra_nic_3d(filename);
    // _______________________


    while( t < TMAX){
        t++;

        //explicit solution
        cout << endl << "explicit"<< endl;
        print_1d(expli, y_es, n, delta_x, delta_t*t);
        explicit_solver(y_es,n,alpha);

        //implicit solution
        cout << endl << "implicit"<< endl;
        print_1d(impli, y_im, n, delta_x, delta_t*t);
        implicit_solver(y_im, n, alpha);

        //Crank-Nicolson solution
        cout << endl << "cranic"<< endl;
        print_1d(cra_nic, y_cn_1, n, delta_x, delta_t*t);
        print_3d(cra_nic_3d, y_cn_1, n, &print_way, delta_x, t*delta_t);
        cranic_solver(y_cn_1, y_cn, n, alpha);
        for(int i=0; i<n; i++){
            y_cn_1[i] = y_cn[i] + dirichlet(i*delta_x, t*delta_t);
        }
        //for have variable boundary conditions
        y_cn_1[0] = dirichlet(0,(t+1)*delta_t);
        y_cn_1[n-1] = dirichlet(LEGH,(t+1)*delta_t);

        for(int i=0; i<n; i++){
            y_cn_1[i] -= dirichlet(i*delta_x,(t+1)*delta_t);
        }
        //______________________
    }

    expli.close();
    impli.close();
    cra_nic.close(); 


cout << "See files for data." << endl;

return 0;
}


double dirichlet(double x, double t){
    double a =1;
    //boundary condition that change in the time
    a= 0.5 + 0.5*cos(10*t + M_PI);
    double b = 0;
    double l = LEGH;

    return a + x*(b-a)/l;

}

void tridiagonal_solver(double a, double b, double c, double *u, double *solution, int n){

    double *coef = new double[n];
    double *u_primo = new double[n];

    coef[0] = c/b;

    u_primo[0] = u[0]/b;

    for(int i=1; i<n; i++){
        coef[i] = c/(b - a*coef[i-1]);
        u_primo[i] = (u[i] - a*u_primo[i-1])/(b - a*coef[i-1]);
    }

    solution[n-1] = u_primo[n-1];

    for(int i=n-2; i>=1; i--){
         solution[i] = u_primo[i] - coef[i]*solution[i+1];
    }

    solution[0] = u_primo[0] - (c/b)*solution[1];

}

void print_3d(ofstream &output, double *y , int n, int *print_way, double delta_x, double t){

    int i;

    if(*print_way==0){
       for(i=0; i<n ; i++){
          *print_way = 1;
           output <<std::scientific << i << "    " << t << "    " << y[i] + dirichlet((i)*delta_x, t) << endl;
         }
    }else{
      *print_way=0;

       for(i=n-1; i>=0; i--){
          output <<std::scientific << i << "    " << t << "    " << y[i] + dirichlet((i)*delta_x,t) << endl;
       }
    }

    /*Standard print
    for(int i=0; i<n ; i++){
         output <<std::scientific << y[i] + dirichlet((i)*delta_x, t) << "  ";
    }
    output << endl;
    */
}

void print_1d(ofstream &output,  double *y , int n, double delta_x, double t){

     output.precision(5);

     for(int i=0; i<n ; i++){
         output <<std::scientific  << y[i] + dirichlet(i*delta_x, t) << "    ";
         cout <<std::scientific  << y[i] + dirichlet(i*delta_x, t) << "    " << endl;
     }

     output << endl;
     cout << endl;
}

double condition_t_zero_dirichlet(double x){ //we set the diricet bandouri condition.

    double value_function = 0;

    if(x==0)  value_function = dirichlet(0,0); //this is the initial condition. (in our case we set a concentration non equale to zero only in x = 0)

    value_function -= dirichlet(x,0) ;  //diriclet boudary condition
    return value_function;

}

void explicit_solver(double *y, int n, double alpha){
    double y_mainos_one = 0;

    for(int i=1; i<n-1; i++){

        double temp =  alpha*(y_mainos_one +y[i+1]) + (1-2*alpha)*y[i];
        y_mainos_one = y[i];
        y[i] = temp;

    }
}
void implicit_solver(double *u, int n, double alpha){

     double b, temp1, temp2, *coef = new double[n];

     b=1+2*alpha;
     u[0] = u[n-1]=0;

     //set coefficients
     coef[1] = -alpha/b;
     u[1] = u[1]/b;

     for(int i=2; i<n-1; i++){
         coef[i] = -alpha/(b + alpha*coef[i-1]);
         u[i] = (u[i] + alpha*u[i-1])/(b +  alpha*coef[i-1]);
     }

     temp1 = 0;

     for(int i=n-2; i>=1; i--){
         temp2 = u[i] - coef[i]*temp1;
         temp1 = u[i];
         u[i] = temp2;
     }
}
void cranic_solver(double *v_1, double *v, int n, double alpha ){

    double temp, temp2, a, b;
    a = 2 - 2*alpha;
    b = alpha;

    temp= v_1[1];
    v_1[1] = a*v_1[1] + b*v_1[2];

    for(int i=2; i<n-2; i++){
        temp2 = v_1[i];
        v_1[i] = b*(temp + v_1[i+1]) + a*v_1[i];
        temp = temp2;
    }
    v_1[n-2] = b*temp2 + a*v_1[n-2];

    a= 2+2*alpha;
    b= -alpha;
    tridiagonal_solver(b,a,b,&v_1[1],&v[1], n-2);
}



double analitical_solution(double x, double t){

         double sum2=0;

         for(int n=1; n<100; n++){
                sum2 += (-2/(M_PI*n))*sin(M_PI*n*x)*exp(-M_PI*n*M_PI*n*t);
         }

return  1 - x + sum2;
}












