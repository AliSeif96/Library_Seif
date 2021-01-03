#ifndef NEURONAL_MODELS_H_INCLUDED
#define NEURONAL_MODELS_H_INCLUDED

//{includes
#include <iostream>                                             //for cout
#include <fstream>                                              //for ofstream
#include <stdlib.h>                                             //srand, rand
#include <time.h>                                               //time
#include <sstream>                                              //for stringstream
#include <fstream>                                              //for ofstream
#include <math.h>                                               //for pow()
using namespace std;                                            //for Standard program
//}

/* Help
input

Number_of_nodes=----------->number of neurons from 0 to 100(size_max)
Matrix[size_max][size_max]=>Possibility of connection from 0 to 1
External_electric_current_HH=------------------->address save data
Length_of_steps=------------------->variable change
Final_time
Coupling_strength_HH

output

Beginning




*/

//{Variables
#define size_max 200                                            //Number of neurons
int Number_of_nodes_N;
int Matrix_N[size_max][size_max];                                       //Because matrix start from 100
float External_electric_current_HH;
double Length_of_steps;
int Final_time;
float Coupling_strength_HH;
float Changer_N;
float phi_wb;

//}

//{class Hodgkin_Huxley
//##############################################################
//####                                                      ####
//####          Create class Hodgkin_Huxley                 ####
//####                                                      ####
//##############################################################
class Hodgkin_Huxley
{
	private:                                                    //Fixed and variable values that we do not have access to outside the class
		float Gna = 120.0;                                     //sodiom_constant
		float Gk = 36.0;                                       //Potasiom_constant
		float Gl = 0.3;                                        //Leak_constant
		float Gsyn = Coupling_strength_HH;                                        //Leak_constant
        float Ena = 55.17;                                     //Votage_sodiom
		float Ek = -72.14;                                     //Votage_Potasiom
		float El = -49.42;                                     //Votage_Leak
		float Esyn = 20.0;                                     //Votage_Leak
        float Cm = 1.0;                                        //Capacitor_capacity
        float  k1, k2, k3, k4;                                 //define 4 point for calculate Runge-Kutta
	public:                                                     //Fixed and variable values that we have access to from outside the class
	    ofstream temp;                                          //create file for save data
	    void Beginning(void);                                       //call start and Initial values
        void onedt(void);                                       //Run for one step dt
        void Print_Data(void);                                    //run for all time and print in file
        void Show_Data(void);                                     //run for all time and show in exe
		float alpha_n(float);                                 //calculate alpha n
		float beta_n(float);                                  //calculate beta n
	  	float alpha_m(float);                                 //calculate alpha m
	  	float beta_m(float);                                  //calculate beta m
	  	float alpha_h(float);                                 //calculate alpha h
	  	float beta_h(float);                                  //calculate beta h
	  	float alpha_r(float);                                 //calculate alpha r
	  	float beta_r(float);                                  //calculate beta r
        float n_inf(float);                                   //calculate infinitude n
	  	float m_inf(float);                                   //calculate infinitude m
	  	float h_inf(float);                                   //calculate infinitude h
        float r_inf(float);                                   //calculate infinitude r
	  	float INa1(float,float,float);                      //sodiom Current
	  	float IK1(float,float);                              //potasiom Current
	  	float Il1(float);                                     //leak Current
	  	float Isyn(float,int);                                     //leak Current
	  	float dvdt(float,float,float,float,float,int);        //Differential equation for voltage
	  	float dndt(float,float,float);                      //Differential equation for n
        float dmdt(float,float,float);                      //Differential equation for m
        float dhdt(float,float,float);                      //Differential equation for h
        float drdt(float,float,float);                      //Differential equation for r
	  	float rk4thOrder_v(float,float,float,float,float,float,int);//Runge-Kutta for voltage
	  	float rk4thOrder_n(float,float,float,float);       //Runge-Kutta for n
        float rk4thOrder_m(float,float,float,float);       //Runge-Kutta for m
	  	float rk4thOrder_h(float,float,float,float);       //Runge-Kutta for h
        float rk4thOrder_r(float,float,float,float);       //Runge-Kutta for r
        float dt = Length_of_steps;                                //length steps
        float V[size_max][2];                                 //define matrix for pre and post voltage of any neuron
        float N[size_max][2];                                 //define matrix for pre and post n of any neuron
        float H[size_max][2];                                 //define matrix for pre and post h of any neuron
        float M[size_max][2];                                 //define matrix for pre and post m of any neuron
        float R[size_max][2];                                 //define matrix for pre and post r of any neuron
        //float A[size_max][size_max];                              //Because matrix start from 0
	  	double t0;                                              //Pedometer
	  	int conter;
  		float v,n,m,h,r;                                       //variable values for each step
};

//_________________________________Calculate alpha and betas_________________________________//

float Hodgkin_Huxley::alpha_n(float v){return   0.01*(v+50)/(1-exp(-(v+50)/10));}
float Hodgkin_Huxley::beta_n(float v){return   0.125*exp(-(v+60)/80);}
float Hodgkin_Huxley::alpha_m(float v){return   0.1*(v+35)/(1-exp(-(v+35)/10));}
float Hodgkin_Huxley::beta_m(float v){return   4.0*exp(-0.0556*(v+60));}
float Hodgkin_Huxley::alpha_h(float v){return   0.07*exp(-0.05*(v+60));}
float Hodgkin_Huxley::beta_h(float v){return   1/(1+exp(-(0.1)*(v+30)));}
float Hodgkin_Huxley::alpha_r(float v){return   1.875/(1+exp(-(v+20)));}
float Hodgkin_Huxley::beta_r(float v){return   0.125;}

//__________________________Calculate infinite activation variables__________________________//

float Hodgkin_Huxley::n_inf(float v){return alpha_n(v)/(alpha_n(v)+beta_n(v));}
float Hodgkin_Huxley::h_inf(float v){return alpha_h(v)/(alpha_h(v)+beta_h(v));}
float Hodgkin_Huxley::m_inf(float v){return alpha_m(v)/(alpha_m(v)+beta_m(v));}
float Hodgkin_Huxley::r_inf(float v){return alpha_r(v)/(alpha_r(v)+beta_r(v));}

//__________________________________Calculation of currents__________________________________//

float Hodgkin_Huxley::INa1(float v,float h,float m) {return Gna*h*pow(m,3)*(v-Ena);}
float Hodgkin_Huxley::IK1(float v,float n) {return Gk*pow(n,4)*(v-Ek);}
float Hodgkin_Huxley::Il1(float v) {return Gl*(v-El);}

float Hodgkin_Huxley::Isyn(float v,int conter) {

    float sigma=0;
    for(int i=0; i<Number_of_nodes_N ; i++){

            sigma=sigma+Matrix_N[conter][i]*R[i][0]*(Esyn-V[i][0]);
            //cout<<"Isyn"<<i<<'\t'<<A[conter][i]<<endl;
    }
    return Gsyn*sigma;
    }

//___________________________________Differential Equations__________________________________//

float Hodgkin_Huxley::dvdt(float t, float v,float n,float h,float m,int conter){
    return  (1/Cm)*(External_electric_current_HH +Isyn(v,conter)-(INa1(v,h,m)+IK1(v,n)+Il1(v)));}
float Hodgkin_Huxley::dndt(float t,float n, float v){return  ((alpha_n(v)*(1-n))-beta_n(v)*n);}
float Hodgkin_Huxley::dhdt(float t, float h, float v){return   ((alpha_h(v)*(1-h))-beta_h(v)*h);}
float Hodgkin_Huxley::dmdt(float t, float m, float v){return   ((alpha_m(v)*(1-m))-beta_m(v)*m);}
float Hodgkin_Huxley::drdt(float t, float r, float v){return   ((alpha_r(v)*(1-r))-beta_r(v)*r);}

//__________________________________Runge-Kutta calculations_________________________________//

float Hodgkin_Huxley::rk4thOrder_v(float t0, float v, float dt,float n,float h,float m,int conter) {
            k1=     dt*dvdt(t0, v,n,h,m,conter);
            k2=     dt*dvdt((t0+dt/2), (v+k1/2),n,h,m,conter);
            k3=     dt*dvdt((t0+dt/2), (v+k2/2),n,h,m,conter);
            k4=     dt*dvdt((t0+dt), (v+k3),n,h,m,conter);
            v=      v+float((1.0/6.0)*(k1+2*k2+2*k3+k4));
   return   v;}
float Hodgkin_Huxley::rk4thOrder_n(float t0, float v, float dt, float n) {
            k1=     dt*dndt(t0, n,v);
            k2=     dt*dndt((t0+dt/2), (n+k1/2),v);
            k3=     dt*dndt((t0+dt/2), (n+k2/2),v);
            k4=     dt*dndt((t0+dt), (n+k3),v);
            n=      n+float((1.0/6.0)*(k1+2*k2+2*k3+k4));
   return   n;}
float Hodgkin_Huxley::rk4thOrder_h(float t0, float v, float dt,float h) {
            k1=     dt*dhdt(t0, h,v);
            k2=     dt*dhdt((t0+dt/2), (h+k1/2),v);
            k3=     dt*dhdt((t0+dt/2), (h+k2/2),v);
            k4=     dt*dhdt((t0+dt), (h+k3),v);
            h=      h+float((1.0/6.0)*(k1+2*k2+2*k3+k4));
    return   h;}
float Hodgkin_Huxley::rk4thOrder_m(float t0, float v, float dt,float m) {
            k1=     dt*dmdt(t0, m,v);
            k2=     dt*dmdt((t0+dt/2), (m+k1/2),v);
            k3=     dt*dmdt((t0+dt/2), (m+k2/2),v);
            k4=     dt*dmdt((t0+dt), (m+k3),v);
            m=      m+float((1.0/6.0)*(k1+2*k2+2*k3+k4));
    return   m;}

float Hodgkin_Huxley::rk4thOrder_r(float t0, float v, float dt,float r) {
            k1=     dt*dmdt(t0, r,v);
            k2=     dt*dmdt((t0+dt/2), (r+k1/2),v);
            k3=     dt*dmdt((t0+dt/2), (r+k2/2),v);
            k4=     dt*dmdt((t0+dt), (r+k3),v);
            r=      r+float((1.0/6.0)*(k1+2*k2+2*k3+k4));
    return   r;}
//_________________________________One step run calculations_________________________________//

void Hodgkin_Huxley::Beginning(){


    for (int i=0 ; i<Number_of_nodes_N;i++){
       for (int j=0 ; j<=1;j++){
            V[i][j]=0;
            N[i][j]=0;
            H[i][j]=0;
            M[i][j]=0;
            R[i][j]=0;}}
    v=-20.0;
    for(int i=0;i<Number_of_nodes_N;i++){
        V[i][1]=v;

        N[i][1]=n_inf(v);
        H[i][1]=h_inf(v);
        M[i][1]=m_inf(v);
        R[i][1]=r_inf(v);
        cout<<i+1<<'\t'<<V[i][1]<<endl;
        }
    cout << "\nBeginning HH Finish" << endl;

}

//__________________________run for all steps and print in file _____________________________//
void Hodgkin_Huxley::onedt(){
    v=rk4thOrder_v(t0, v, dt,n,h,m,conter);
    n=rk4thOrder_n(t0,v, dt ,n);
    h=rk4thOrder_h(t0, v, dt ,h);
    m=rk4thOrder_m(t0, v, dt ,m);
    r=rk4thOrder_r(t0, v, dt ,r);
    }

//__________________________run for all steps and print in file _____________________________//


void Hodgkin_Huxley::Print_Data(){

    std::string s;
	std::stringstream ss;
	ss << Changer_N;
	ss >> s;
    std::string scheme ("E:/programs/sync/read_data_aida/data/Hodgkin_Huxley_");
    std::string hostname;
    std::string url;
    hostname = s+".txt" ;
    url = scheme+ hostname;
    char* char_arr;
    string str_obj(url);
    char_arr = &str_obj[0];
    char * filename = char_arr;
    ofstream temp(filename);
    cout<<filename;
    for(t0=0;t0<Final_time;t0+=0.1){
        //cout<<t0<<"\t"<<v<<endl;
        conter=-1;
        for ( int i=0 ; i<Number_of_nodes_N ; i++){
        conter=conter+1;
            v=V[i][0];
            n=N[i][0];
            h=H[i][0];
            m=M[i][0];
            r=R[i][0];
            //cout<<v;
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
        if(t0>=5.0){
        float Rsum=0;
        temp << t0;
        for(int k=0;k<Number_of_nodes_N;k++){
            Rsum= Rsum+ R[k][1];
            temp<<'\t'<<V[k][1];
        }
        temp<<'\t'<<Rsum<<endl;
        }
    }
    temp.close();
    cout << "\nPrint Data HH Finish" << endl;
}
//}

//{class Wang_Buzsaki
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//_________________________________________________________________________________________________________________________________________     //$$$
//_________________________________________________________________________________________________________________________________________     //$$$
//_________________________________________________________________________________________________________________________________________     //$$$
//_________________________________________________________________________________________________________________________________________     //$$$
//__________________________________________________                   ____________________________________________________________________     //$$$
//__________________________________________________    Wang_Buzsaki   ____________________________________________________________________     //$$$
//_________________________________________________________________________________________________________________________________________     //$$$
//_________________________________________________________________________________________________________________________________________     //$$$
//_________________________________________________________________________________________________________________________________________     //$$$
//_________________________________________________________________________________________________________________________________________     //$$$
//##############################################################                                                                                //$$$
//####                                                      ####                                                                                //$$$
//####            Create class Wang_Buzsaki                 ####                                                                                //$$$
//####                                                      ####                                                                                //$$$
//##############################################################                                                                                //$$$
class Wang_Buzsaki                                                                                                                              //$$$
{                                                                                                                                               //$$$
	private:                                                    //Fixed and variable values that we do not have access to outside the class     //$$$
		int Gna = 35;                                           //sodiom_constant                                                               //$$$
		int Gk = 9;                                             //Potasiom_constant                                                             //$$$
		float Gl = 0.1 ;                                        //Leak_constant                                                                 //$$$
		float Gsyn = 0.1;                                       //synaps_constant                                                               //$$$
        int Ena = 55;                                           //Votage_sodiom                                                                 //$$$
		int Ek = -90;                                           //Votage_Potasiom                                                               //$$$
		int El = -65;                                           //Votage_Leak                                                                   //$$$
		float Esyn = -75;                                       //Votage_synaps                                                                 //$$$
        float Cm = 1.0;                                         //Capacitor_capacity                                                            //$$$
        float phi=phi_wb;                                       //phi                                                                           //$$$
        float K=Coupling_strength_HH;                           //Coupling_strength                                                             //$$$
        float  k1, k2, k3, k4;                                  //define 4 point for calculate Runge-Kutta                                      //$$$
	public:                                                     //Fixed and variable values that we have access to from outside the class       //$$$
        ofstream temp2;                                         //create file for save data                                                     //$$$
	    void Beginning(void);                                   //call start and Initial values                                                 //$$$
        void onedt(void);                                       //Run for one step dt                                                           //$$$
        void Print_Data(void);                                  //run for all time and print in file                                            //$$$
		float alpha_n(float);                                   //calculate alpha n                                                             //$$$
		float beta_n(float);                                    //calculate beta n                                                              //$$$
	  	float alpha_m(float);                                   //calculate alpha m                                                             //$$$
	  	float beta_m(float);                                    //calculate beta m                                                              //$$$
	  	float alpha_h(float);                                   //calculate alpha h                                                             //$$$
	  	float beta_h(float);                                    //calculate beta h                                                              //$$$
        float F(float);                                         //calculate F                                                                   //$$$
        float n_inf(float);                                     //calculate infinitude n                                                        //$$$
	  	float m_inf(float);                                     //calculate infinitude m                                                        //$$$
	  	float h_inf(float);                                     //calculate infinitude h                                                        //$$$
        float s_inf(float);                                     //calculate infinitude s                                                        //$$$
	  	float INa1(float,float);                                //sodiom Current                                                                //$$$
	  	float IK1(float,float);                                 //potasiom Current                                                              //$$$
	  	float Il1(float);                                       //leak Current                                                                  //$$$
	  	float Isyn(float,int);                                  //post synaptic Current                                                         //$$$
	  	float dvdt(float,float,float,float,int);                //Differential equation for voltage                                             //$$$
	  	float dndt(float,float,float);                          //Differential equation for n                                                   //$$$
        float dhdt(float,float,float);                          //Differential equation for h                                                   //$$$
        float dsdt(float,float,float);                          //Differential equation for s                                                   //$$$
	  	float rk4thOrder_v(float,float,float,float,float,int);  //Runge-Kutta for voltage                                                       //$$$
	  	float rk4thOrder_n(float,float,float,float);            //Runge-Kutta for n                                                             //$$$
	  	float rk4thOrder_h(float,float,float,float);            //Runge-Kutta for h                                                             //$$$
        float rk4thOrder_s(float,float,float,float);            //Runge-Kutta for s                                                             //$$$
        float dt = Length_of_steps;                             //length steps                                                                  //$$$
        float V[size_max][2];                                   //define matrix for pre and post voltage of any neuron                          //$$$
        float N[size_max][2];                                   //define matrix for pre and post n of any neuron                                //$$$
        float H[size_max][2];                                   //define matrix for pre and post h of any neuron                                //$$$
        float S[size_max][2];                                   //define matrix for pre and post s of any neuron                                //$$$
        //float A[size_max][size_max];                          //Because matrix start from 0                                                   //$$$
	  	double t0;                                              //Pedometer                                                                     //$$$
	  	int conter;                                             //Counter                                                                       //$$$
  		float v,n,h,s;                                          //variable values for each step                                                 //$$$
};                                                              //variable values for each step                                                 //$$$
//_________________________________________________________________________________________________________Calculate alpha and betas____________//$$$
float Wang_Buzsaki::alpha_n(float v){return   -0.01*(v+34.0)/(exp(-0.1*(v+34.0))-1.0);}                                                         //$$$
float Wang_Buzsaki::beta_n(float v){return   0.125*exp(-(v+44.0)/80.0);}                                                                        //$$$
float Wang_Buzsaki::alpha_m(float v){return   -0.1*(v+35.0)/(exp(-0.1*(v+35.0))-1.0);}                                                          //$$$
float Wang_Buzsaki::beta_m(float v){return   4.0*exp(-(v+60.0)/18.0);}                                                                          //$$$
float Wang_Buzsaki::alpha_h(float v){return   0.07*exp(-(v+58.0)/20.0);}                                                                        //$$$
float Wang_Buzsaki::beta_h(float v){return   1.0/(exp(-0.1*(v+28.0))+1.0);}                                                                     //$$$
float Wang_Buzsaki::F(float v){return   1.0/(1.0+exp(-0.5*(v-0.0)));}                                                                           //$$$
//_________________________________________________________________________________________________Calculate infinite activation variables______//$$$
float Wang_Buzsaki::n_inf(float v){return alpha_n(v)/(alpha_n(v)+beta_n(v));}                                                                   //$$$
float Wang_Buzsaki::h_inf(float v){return alpha_h(v)/(alpha_h(v)+beta_h(v));}                                                                   //$$$
float Wang_Buzsaki::m_inf(float v){return alpha_m(v)/(alpha_m(v)+beta_m(v));}                                                                   //$$$
float Wang_Buzsaki::s_inf(float v){return 12.0/(12.0+0.1);}                                                                                     //$$$
//________________________________________________________________________________________________________________Calculation of currents_______//$$$
float Wang_Buzsaki::INa1(float v,float h) {return Gna*h*pow(m_inf(v),3)*(v-Ena);}                                                               //$$$
float Wang_Buzsaki::IK1(float v,float n) {return Gk*pow(n,4)*(v-Ek);}                                                                           //$$$
float Wang_Buzsaki::Il1(float v) {return Gl*(v-El);}                                                                                            //$$$
float Wang_Buzsaki::Isyn(float v,int conter) {                                                                                                  //$$$
    float sigma=0;                                                                                                                              //$$$
    for(int i=0; i<Number_of_nodes_N ; i++){                                                                                                    //$$$
            sigma=sigma+Matrix_N[conter][i]*Gsyn*S[i][0]*(V[i][0]-Esyn);                                                                        //$$$
    }                                                                                                                                           //$$$
    return (K/Number_of_nodes_N)*sigma;                                                                                                         //$$$
    }                                                                                                                                           //$$$
//____________________________________________________________________________________________________________________Differential Equations____//$$$
float Wang_Buzsaki::dvdt(float t, float v,float n,float h,int conter){                                                                          //$$$
    return  (1/Cm)*(External_electric_current_HH -Isyn(v,conter)-(INa1(v,h)+IK1(v,n)+Il1(v)));}                                                 //$$$
float Wang_Buzsaki::dndt(float t,float n, float v){return  phi*((alpha_n(v)*(1-n))-beta_n(v)*n);}                                               //$$$
float Wang_Buzsaki::dhdt(float t, float h, float v){return   phi*((alpha_h(v)*(1-h))-beta_h(v)*h);}                                             //$$$
float Wang_Buzsaki::dsdt(float t, float s, float v){return   ((12.0*F(v)*(1-s))-0.1*s);}                                                        //$$$
//_________________________________________________________________________________________________________________Runge-Kutta calculations_____//$$$
float Wang_Buzsaki::rk4thOrder_v(float t0, float v, float dt,float n,float h,int conter) {                                                      //$$$
            k1=     dt*dvdt(t0, v,n,h,conter);                                                                                                  //$$$
            k2=     dt*dvdt((t0+dt/2), (v+k1/2),n,h,conter);                                                                                    //$$$
            k3=     dt*dvdt((t0+dt/2), (v+k2/2),n,h,conter);                                                                                    //$$$
            k4=     dt*dvdt((t0+dt), (v+k3),n,h,conter);                                                                                        //$$$
            v=      v+float((1.0/6.0)*(k1+2*k2+2*k3+k4));                                                                                       //$$$
   return   v;}                                                                                                                                 //$$$
float Wang_Buzsaki::rk4thOrder_n(float t0, float v, float dt, float n) {                                                                        //$$$
            k1=     dt*dndt(t0, n,v);                                                                                                           //$$$
            k2=     dt*dndt((t0+dt/2), (n+k1/2),v);                                                                                             //$$$
            k3=     dt*dndt((t0+dt/2), (n+k2/2),v);                                                                                             //$$$
            k4=     dt*dndt((t0+dt), (n+k3),v);                                                                                                 //$$$
            n=      n+float((1.0/6.0)*(k1+2*k2+2*k3+k4));                                                                                       //$$$
   return   n;}                                                                                                                                 //$$$
float Wang_Buzsaki::rk4thOrder_h(float t0, float v, float dt,float h) {                                                                         //$$$
            k1=     dt*dhdt(t0, h,v);                                                                                                           //$$$
            k2=     dt*dhdt((t0+dt/2), (h+k1/2),v);                                                                                             //$$$
            k3=     dt*dhdt((t0+dt/2), (h+k2/2),v);                                                                                             //$$$
            k4=     dt*dhdt((t0+dt), (h+k3),v);                                                                                                 //$$$
            h=      h+float((1.0/6.0)*(k1+2*k2+2*k3+k4));                                                                                       //$$$
    return   h;}                                                                                                                                //$$$
float Wang_Buzsaki::rk4thOrder_s(float t0, float v, float dt,float s) {                                                                         //$$$
            k1=     dt*dsdt(t0, s,v);                                                                                                           //$$$
            k2=     dt*dsdt((t0+dt/2), (s+k1/2),v);                                                                                             //$$$
            k3=     dt*dsdt((t0+dt/2), (s+k2/2),v);                                                                                             //$$$
            k4=     dt*dsdt((t0+dt), (s+k3),v);                                                                                                 //$$$
            s=      s+float((1.0/6.0)*(k1+2*k2+2*k3+k4));                                                                                       //$$$
    return   s;}                                                                                                                                //$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//##############################################################                                                                                //$$$
//####                                                      ####                                                                                //$$$
//####                     Beginning                        ####                                                                                //$$$
//####             One step run calculations                ####                                                                                //$$$
//##############################################################                                                                                //$$$
void Wang_Buzsaki::Beginning(){                                                                                                                 //$$$
    cout<<endl;                                                                                                                                 //$$$
    for (int i=0 ; i<Number_of_nodes_N;i++){                                                                                                    //$$$
       for (int j=0 ; j<=1;j++){                                                                                                                //$$$
            V[i][j]=0;                                                                                                                          //$$$
            N[i][j]=0;                                                                                                                          //$$$
            H[i][j]=0;                                                                                                                          //$$$
            S[i][j]=0;}}                                                                                                                        //$$$
    v=-63.0;                                                                                                                                    //$$$
    for(int i=0;i<Number_of_nodes_N;i++){                                                                                                       //$$$
        V[i][1]=v;                                                                                                                              //$$$
        N[i][1]=n_inf(v);                                                                                                                       //$$$
        H[i][1]=h_inf(v);                                                                                                                       //$$$
        S[i][1]=s_inf(v);                                                                                                                       //$$$
        cout<<i+1<<'\t'<<V[i][1]<<endl;                                                                                                         //$$$
        }                                                                                                                                       //$$$
    cout << "\nBeginning WB Finish" << endl;                                                                                                    //$$$
}                                                                                                                                               //$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//##############################################################                                                                                //$$$
//####                                                      ####                                                                                //$$$
//####                       onedt                          ####                                                                                //$$$
//####                                                      ####                                                                                //$$$
//##############################################################                                                                                //$$$
void Wang_Buzsaki::onedt(){                                                                                                                     //$$$
    v=rk4thOrder_v(t0, v, dt,n,h,conter);                                                                                                       //$$$
    n=rk4thOrder_n(t0,v, dt ,n);                                                                                                                //$$$
    h=rk4thOrder_h(t0, v, dt ,h);                                                                                                               //$$$
    s=rk4thOrder_s(t0, v, dt ,s);                                                                                                               //$$$
    }                                                                                                                                           //$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//##############################################################                                                                                //$$$
//####                                                      ####                                                                                //$$$
//####         run for all steps and print in file          ####                                                                                //$$$
//####                                                      ####                                                                                //$$$
//##############################################################                                                                                //$$$
void Wang_Buzsaki::Print_Data(){                                                                                                                //$$$
    std::string ssss;                                                                                                                           //$$$
	std::stringstream ss;                                                                                                                       //$$$
	ss << Changer_N;                                                                                                                            //$$$
	ss >> ssss;                                                                                                                                 //$$$
    std::string scheme ("E:/programs/sync/read_data_aida/data/Wang_Buzsaki_");                                                                  //$$$
    std::string hostname;                                                                                                                       //$$$
    std::string url;                                                                                                                            //$$$
    hostname = ssss+".txt" ;                                                                                                                    //$$$
    url = scheme+ hostname;                                                                                                                     //$$$
    char* char_arr;                                                                                                                             //$$$
    string str_obj(url);                                                                                                                        //$$$
    char_arr = &str_obj[0];                                                                                                                     //$$$
    char * filename1 = char_arr;                                                                                                                //$$$
    ofstream temp2(filename1);                                                                                                                  //$$$
    cout<<filename1;                                                                                                                            //$$$
    cout<<"asdf"<<endl;                                                                                                                         //$$$
    for (t0=dt ; t0<=Final_time ;t0=t0 + dt){                                                                                                   //$$$
    conter=-1;                                                                                                                                  //$$$
    //cout<<t0<<endl;                                                                                                                           //$$$
        for ( int i=0 ; i<Number_of_nodes_N ; i++){                                                                                             //$$$
            conter=conter+1;                                                                                                                    //$$$
            v=V[i][0];                                                                                                                          //$$$
            n=N[i][0];                                                                                                                          //$$$
            h=H[i][0];                                                                                                                          //$$$
            s=S[i][0];                                                                                                                          //$$$
            onedt();                                                                                                                            //$$$
            V[i][1]=v;                                                                                                                          //$$$
            N[i][1]=n;                                                                                                                          //$$$
            H[i][1]=h;                                                                                                                          //$$$
            S[i][1]=s;                                                                                                                          //$$$
            V[i][0]=V[i][1];                                                                                                                    //$$$
            N[i][0]=N[i][1];                                                                                                                    //$$$
            H[i][0]=H[i][1];                                                                                                                    //$$$
            S[i][0]=S[i][1];                                                                                                                    //$$$
        }                                                                                                                                       //$$$
        float Ssum=0;                                                                                                                           //$$$
        temp2<<t0;                                                                                                                              //$$$
        for(int k=0;k<Number_of_nodes_N;k++){                                                                                                   //$$$
            Ssum= Ssum+ S[k][1];                                                                                                                //$$$
            temp2<<'\t'<<V[k][1];                                                                                                               //$$$
        }                                                                                                                                       //$$$
        temp2<<'\t'<<Ssum<<endl;                                                                                                                //$$$
    }                                                                                                                                           //$$$
    temp2.close();                                                                                                                              //$$$
    cout<<"\nPrint Data Wang_Buzsaki Finished"<<endl;                                                                                           //$$$
}                                                                                                                                               //$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//}

#endif // NEURONAL_MODELS_H_INCLUDED
