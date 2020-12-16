#ifndef SEIF_H_INCLUDED
#define SEIF_H_INCLUDED
#include <iostream>                                             //for cout
#include <fstream>                                              //for ofstream
#include <stdlib.h>                                             //srand, rand
#include <time.h>                                               //time
using namespace std;                                            //for Standard program
#define size_max 100                                            //Number of max neurons
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_____________________________________________________                   _________________________________________________________________
//_____________________________________________________     Variables     _________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
int Number_of_nodes=0;
float Possibility_of_connection=0;//Possibility of connection
float Matrix[size_max][size_max];//Because matrix start from 100
float External_electric_current=0;
float Length_of_steps=0;
int Final_time=0;
int Storage_start_time=0;
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//____________________________________________________                       ______________________________________________________________
//____________________________________________________   Adjacency matrix    ______________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//##############################################################
//####                                                      ####
//####            Create class Adjacency                    ####
//####                                                      ####
//##############################################################
class Adjacency
{
	private:                                                    //Fixed and variable values that we do not have access to outside the class
        int iSecret;                                            //Random probability between zero and one
public:                                                     //Fixed and variable values that we have access to from outside the class
	    ofstream temp1;                                         //create file for save data
	    ofstream temp2;                                         //create file for save data
        void Build_matrix(void);                                        //create Adjacency matrix
        void Show_Binary(void);                                        //show Adjacency matrix
        void Print_Network(void);                                       //print Adjacency matrix
        void Print_Binary(void);                                 //print Adjacency matrix
        int numb;                                            //Counting links
};
//________________________________ create Adjacency matrix __________________________________//
void Adjacency::Build_matrix(){
    srand (time(NULL));                                         //initialize random seed
    for (int i=0 ; i<Number_of_nodes;i++){
       for (int j=0 ; j<Number_of_nodes;j++){Matrix[i][j]=0;}}
    numb=0;
    for (int i=0 ; i<Number_of_nodes;i++){
        for (int j=i ; j<Number_of_nodes;j++){
            iSecret = rand() % 10 + 1;
            int pp=Possibility_of_connection*10;
            if (iSecret<=pp && i!=j ){
                Matrix[i][j]=1;
                numb=numb+1;
            }
            else{
               Matrix[i][j]=0;
            }
            Matrix[j][i]=Matrix[i][j];
        }
    }
}
//________________________________ show Adjacency matrix ____________________________________//
void Adjacency::Show_Binary(){
    for (int i=1 ; i<=Number_of_nodes;i++){cout<<'\t'<<i;}
    cout<<endl<<endl;
    for (int i=0 ; i<Number_of_nodes;i++){
        cout<<i+1<<'\t';
        for (int j=0 ; j<Number_of_nodes;j++){cout<<Matrix[i][j]<<'\t';}
        cout<<endl<<endl;
    }
    cout << "\nShow Binary Finished" << endl;
}
//_________________________________ print Adjacency matrix __________________________________//
 void Adjacency::Print_Network(){
    ofstream temp1("F:/Programs/C++/datas.txt");
    temp1<<"{\"nodes\":[{\"name\":\"0\",\"group\":1}";
    for (int i=1 ; i<Number_of_nodes;i++){temp1<<",{\"name\":\""<<i<<"\",\"group\":"<<i<<"}";}
    temp1<<"],\"links\":[";
    int number=0;
    for (int i=0 ; i<Number_of_nodes;i++){
        for (int j=i ; j<Number_of_nodes;j++){
            if (Matrix[i][j]==1){
                number=number+1;
                temp1<<"{\"source\":"<<i<<",\"target\":"<<j<<",\"value\":1}";
                if(number<numb){temp1<<",";}}}}
    temp1<<"]}";
    temp1.close();
    cout << "\nPrint Network Finished" << endl;}
//_________________________________ print binery Adjacency matrix __________________________________//
void Adjacency::Print_Binary(){
    ofstream temp2("F:/Programs/C++/bineryprint.txt");


    for (int i=0 ; i<Number_of_nodes;i++){
        for (int j=0 ; j<Number_of_nodes;j++){temp2<<Matrix[i][j]<<'\t';}
        temp2<<endl;
    }
    temp2.close();
    cout << "\nPrint Binary Finished" << endl;}
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_______________________________________________________                          ________________________________________________________
//_______________________________________________________      Hodgkin_Huxley      ________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//##############################################################
//####                                                      ####
//####          Create class Hodgkin_Huxley                 ####
//####                                                      ####
//##############################################################
#include <math.h>                                               //for pow()
class Hodgkin_Huxley
{
	private:                                                    //Fixed and variable values that we do not have access to outside the class
		double Gna = 120.0;                                     //sodiom_constant
		double Gk = 36.0;                                       //Potasiom_constant
		double Gl = 0.3;                                        //Leak_constant
		double Gsyn = 0.2;                                        //Leak_constant
        double Ena = 55.17;                                     //Votage_sodiom
		double Ek = -72.14;                                     //Votage_Potasiom
		double El = -49.42;                                     //Votage_Leak
		double Esyn = 20.0;                                     //Votage_Leak
        double Cm = 1.0;                                        //Capacitor_capacity
        double  k1, k2, k3, k4;                                 //define 4 point for calculate Runge-Kutta
	public:                                                     //Fixed and variable values that we have access to from outside the class
	    ofstream temp;                                          //create file for save data
	    //void Adjacency(void);
	    void Start(void);                                       //call start and Initial values
        void onedt(void);                                       //Run for one step dt
        void Print_Data(void);                                    //run for all time and print in file
        void Show_Data(void);                                     //run for all time and show in exe
		double alpha_n(double);                                 //calculate alpha n
		double beta_n(double);                                  //calculate beta n
	  	double alpha_m(double);                                 //calculate alpha m
	  	double beta_m(double);                                  //calculate beta m
	  	double alpha_h(double);                                 //calculate alpha h
	  	double beta_h(double);                                  //calculate beta h
	  	double alpha_r(double);                                 //calculate alpha r
	  	double beta_r(double);                                  //calculate beta r
	  	double n_inf(double);                                   //calculate infinitude n
	  	double m_inf(double);                                   //calculate infinitude m
	  	double h_inf(double);                                   //calculate infinitude h
        double r_inf(double);                                   //calculate infinitude r
	  	double INa1(double,double,double);                      //sodiom Current
	  	double IK1(double,double);                              //potasiom Current
	  	double Il1(double);                                     //leak Current
	  	double Isyn(double,int);                                     //leak Current
	  	double dvdt(double,double,double,double,double,int);        //Differential equation for voltage
	  	double dndt(double,double,double);                      //Differential equation for n
        double dmdt(double,double,double);                      //Differential equation for m
        double dhdt(double,double,double);                      //Differential equation for h
        double drdt(double,double,double);                      //Differential equation for r
	  	double rk4thOrder_v(double,double,double,double,double,double,int);//Runge-Kutta for voltage
	  	double rk4thOrder_n(double,double,double,double);       //Runge-Kutta for n
        double rk4thOrder_m(double,double,double,double);       //Runge-Kutta for m
	  	double rk4thOrder_h(double,double,double,double);       //Runge-Kutta for h
        double rk4thOrder_r(double,double,double,double);       //Runge-Kutta for r
        double dt = Length_of_steps;                                //length steps
        double V[size_max][2];                                 //define matrix for pre and post voltage of any neuron
        double N[size_max][2];                                 //define matrix for pre and post n of any neuron
        double H[size_max][2];                                 //define matrix for pre and post h of any neuron
        double M[size_max][2];                                 //define matrix for pre and post m of any neuron
        double R[size_max][2];                                 //define matrix for pre and post r of any neuron
        //double A[size_max][size_max];                              //Because matrix start from 0
	  	double t0;                                              //Pedometer
	  	int conter;
  		double v,n,m,h,r;};                                       //variable values for each step

//_________________________________Calculate alpha and betas_________________________________//

double Hodgkin_Huxley::alpha_n(double v){return   0.01*(v+50)/(1-exp(-(v+50)/10));}
double Hodgkin_Huxley::beta_n(double v){return   0.125*exp(-(v+60)/80);}
double Hodgkin_Huxley::alpha_m(double v){return   0.1*(v+35)/(1-exp(-(v+35)/10));}
double Hodgkin_Huxley::beta_m(double v){return   4.0*exp(-0.0556*(v+60));}
double Hodgkin_Huxley::alpha_h(double v){return   0.07*exp(-0.05*(v+60));}
double Hodgkin_Huxley::beta_h(double v){return   1/(1+exp(-(0.1)*(v+30)));}
double Hodgkin_Huxley::alpha_r(double v){return   1.875/(1+exp(-(v+20)));}
double Hodgkin_Huxley::beta_r(double v){return   0.125;}

//__________________________Calculate infinite activation variables__________________________//

double Hodgkin_Huxley::n_inf(double v){return alpha_n(v)/(alpha_n(v)+beta_n(v));}
double Hodgkin_Huxley::h_inf(double v){return alpha_h(v)/(alpha_h(v)+beta_h(v));}
double Hodgkin_Huxley::m_inf(double v){return alpha_m(v)/(alpha_m(v)+beta_m(v));}
double Hodgkin_Huxley::r_inf(double v){return alpha_r(v)/(alpha_r(v)+beta_r(v));}

//__________________________________Calculation of currents__________________________________//

double Hodgkin_Huxley::INa1(double v,double h,double m) {return Gna*h*pow(m,3)*(v-Ena);}
double Hodgkin_Huxley::IK1(double v,double n) {return Gk*pow(n,4)*(v-Ek);}
double Hodgkin_Huxley::Il1(double v) {return Gl*(v-El);}

double Hodgkin_Huxley::Isyn(double v,int conter) {

    double sigma=0;
    for(int i=0; i<Number_of_nodes ; i++){

            sigma=sigma+Matrix[conter][i]*R[i][0]*(Esyn-V[i][0]);
            //cout<<"Isyn"<<i<<'\t'<<A[conter][i]<<endl;
    }
    return Gsyn*sigma;
    }

//___________________________________Differential Equations__________________________________//

double Hodgkin_Huxley::dvdt(double t, double v,double n,double h,double m,int conter){return  (1/Cm)*(External_electric_current +Isyn(v,conter)-(INa1(v,h,m)+IK1(v,n)+Il1(v)));}
double Hodgkin_Huxley::dndt(double t,double n, double v){return  ((alpha_n(v)*(1-n))-beta_n(v)*n);}
double Hodgkin_Huxley::dhdt(double t, double h, double v){return   ((alpha_h(v)*(1-h))-beta_h(v)*h);}
double Hodgkin_Huxley::dmdt(double t, double m, double v){return   ((alpha_m(v)*(1-m))-beta_m(v)*m);}
double Hodgkin_Huxley::drdt(double t, double r, double v){return   ((alpha_r(v)*(1-r))-beta_r(v)*r);}

//__________________________________Runge-Kutta calculations_________________________________//

double Hodgkin_Huxley::rk4thOrder_v(double t0, double v, double dt,double n,double h,double m,int conter) {
            k1=     dt*dvdt(t0, v,n,h,m,conter);
            k2=     dt*dvdt((t0+dt/2), (v+k1/2),n,h,m,conter);
            k3=     dt*dvdt((t0+dt/2), (v+k2/2),n,h,m,conter);
            k4=     dt*dvdt((t0+dt), (v+k3),n,h,m,conter);
            v=      v+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
   return   v;}
double Hodgkin_Huxley::rk4thOrder_n(double t0, double v, double dt, double n) {
            k1=     dt*dndt(t0, n,v);
            k2=     dt*dndt((t0+dt/2), (n+k1/2),v);
            k3=     dt*dndt((t0+dt/2), (n+k2/2),v);
            k4=     dt*dndt((t0+dt), (n+k3),v);
            n=      n+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
   return   n;}
double Hodgkin_Huxley::rk4thOrder_h(double t0, double v, double dt,double h) {
            k1=     dt*dhdt(t0, h,v);
            k2=     dt*dhdt((t0+dt/2), (h+k1/2),v);
            k3=     dt*dhdt((t0+dt/2), (h+k2/2),v);
            k4=     dt*dhdt((t0+dt), (h+k3),v);
            h=      h+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
    return   h;}
double Hodgkin_Huxley::rk4thOrder_m(double t0, double v, double dt,double m) {
            k1=     dt*dmdt(t0, m,v);
            k2=     dt*dmdt((t0+dt/2), (m+k1/2),v);
            k3=     dt*dmdt((t0+dt/2), (m+k2/2),v);
            k4=     dt*dmdt((t0+dt), (m+k3),v);
            m=      m+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
    return   m;}

double Hodgkin_Huxley::rk4thOrder_r(double t0, double v, double dt,double r) {
            k1=     dt*dmdt(t0, r,v);
            k2=     dt*dmdt((t0+dt/2), (r+k1/2),v);
            k3=     dt*dmdt((t0+dt/2), (r+k2/2),v);
            k4=     dt*dmdt((t0+dt), (r+k3),v);
            r=      r+double((1.0/6.0)*(k1+2*k2+2*k3+k4));
    return   r;}
//_________________________________One step run calculations_________________________________//

void Hodgkin_Huxley::Start(){

    for (int i=0 ; i<Number_of_nodes;i++){
       for (int j=0 ; j<=1;j++){
            V[i][j]=0;
            N[i][j]=0;
            H[i][j]=0;
            M[i][j]=0;
            R[i][j]=0;}}
    v=-20.0;
    for(int i=0;i<Number_of_nodes;i++){
        V[i][1]=v;

        N[i][1]=n_inf(v);
        H[i][1]=h_inf(v);
        M[i][1]=m_inf(v);
        R[i][1]=r_inf(v);
        //cout<<i+1<<M[i][1]<<endl;
        }
}

void Hodgkin_Huxley::onedt(){
    v=rk4thOrder_v(t0, v, dt,n,h,m,conter);
    n=rk4thOrder_n(t0,v, dt ,n);
    h=rk4thOrder_h(t0, v, dt ,h);
    m=rk4thOrder_m(t0, v, dt ,m);
    r=rk4thOrder_r(t0, v, dt ,r);
    }

//__________________________run for all steps and print in file _____________________________//
void Hodgkin_Huxley::Print_Data(){
    //ofstream temp;
    ofstream temp("F:/Programs/C++/temp.txt");
    for (t0=dt ; t0<=Final_time ;t0=t0 + dt){
    conter=-1;
    //cout<<t0<<endl;
        for ( int i=0 ; i<Number_of_nodes ; i++){
            conter=conter+1;
            v=V[i][0];
            n=N[i][0];
            h=H[i][0];
            m=M[i][0];
            r=R[i][0];
            onedt();

            V[i][1]=v;
            N[i][1]=n;
            H[i][1]=h;
            M[i][1]=m;
            R[i][1]=r;

            V[i][0]=V[i][1];
            N[i][0]=N[i][1];
            H[i][0]=H[i][1];
            M[i][0]=M[i][1];
            R[i][0]=R[i][1];
        }
        if(t0>=2000.0){
        double Rsum=0;
        temp<<t0;
        for(int k=0;k<Number_of_nodes;k++){
            Rsum= Rsum+ R[k][1];
            temp<<'\t'<<V[k][1];
        }
        temp<<'\t'<<Rsum<<endl;
        }
    }
    temp.close();
    cout<<"\nPrint Data Finished"<<endl;
}
//___________________________ run for all steps and show in exe _____________________________//
void Hodgkin_Huxley::Show_Data(){
    for (t0=dt ; t0<=Final_time ;t0=t0 + dt){
    conter=-1;
    //cout<<t0<<endl;
        for ( int i=0 ; i<Number_of_nodes ; i++){
            conter=conter+1;
            v=V[i][0];
            n=N[i][0];
            h=H[i][0];
            m=M[i][0];
            r=R[i][0];
            onedt();

            V[i][1]=v;
            N[i][1]=n;
            H[i][1]=h;
            M[i][1]=m;
            R[i][1]=r;

            V[i][0]=V[i][1];
            N[i][0]=N[i][1];
            H[i][0]=H[i][1];
            M[i][0]=M[i][1];
            R[i][0]=R[i][1];
        }
        if(t0>=2000.0){
        double Rsum=0;
        cout<<t0;
        for(int k=0;k<Number_of_nodes;k++){
            Rsum= Rsum+ R[k][1];
            cout<<'\t'<<V[k][1];
        }
        cout<<'\t'<<Rsum<<endl;
        }
    }
    cout<<"\nShow Data Finished"<<endl;
}
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________________________




#endif // SEIF_H_INCLUDED
