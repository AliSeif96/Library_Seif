#ifndef ADJACENCY_MATRIX_H_INCLUDED
#define ADJACENCY_MATRIX_H_INCLUDED

#include <iostream>                                             //for cout
#include <fstream>                                              //for ofstream
#include <stdlib.h>                                             //srand, rand
#include <time.h>                                               //time
#include <sstream>                                              //for stringstream
#include <fstream>                                              //for ofstream
using namespace std;                                            //for Standard program

/*
input

Number_of_nodes=----------->number of neurons from 0 to 100(size_max)
Possibility_of_connection=->Possibility of connection from 0 to 1
Address=------------------->address save data
Changer=------------------->variable change


output

Build_matrix_zero=--------->create matrix that all nodes is zero
Erdos_Renyi=--------------->create Erdos_Renyi network connection matrix that nodes is 0 or 1
Show_Binary=--------------->show binary connection in monitor
Print_Network=------------->3d network print for python
Print_Binary=-------------->binery print for python
Build_array_rand=---------->create array that first 100 rand numbers

*/

#define size_max 100                                            //Number of neurons
float Matrix[size_max][size_max];                                       //Because matrix start from 100
int Number_of_nodes;
float Possibility_of_connection;                              //Possibility of connection
float Changer;
string Address;
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
        void Build_matrix_zero(void);
        void Erdos_Renyi(void);                                        //create Adjacency matrix
        void Show_Binary(void);                                        //show Adjacency matrix
        void Print_Network(void);                                       //print Adjacency matrix
        void Print_Binary(void);                                 //print Adjacency matrix
        void Build_array_rand(void);
        void Build_onelink_step(void);

        int arra[size_max*size_max];
        int numb;
};                                             //Counting links


//________________________________ Build_matrix_zero __________________________________//

void Adjacency::Build_matrix_zero(){
    for (int i=0 ; i<Number_of_nodes;i++){
       for (int j=0 ; j<Number_of_nodes;j++){
            Matrix[i][j]=0;
       }
    }
    cout << "\nBuild matrix zero Finish" << endl;
}
//________________________________ create Adjacency matrix __________________________________//

void Adjacency::Erdos_Renyi(){
    srand (time(NULL));                                         //initialize random seed
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
            Matrix[j][i]=Matrix[i][j];}}
    cout << "\nErdos-Renyi Finish" << endl;
}
//________________________________ show Adjacency matrix ____________________________________//

void Adjacency::Show_Binary(){
    for (int i=1 ; i<=Number_of_nodes;i++){cout<<'\t'<<i;}
    cout<<endl<<endl;
    for (int i=0 ; i<Number_of_nodes;i++){
        cout<<i+1<<'\t';
        for (int j=0 ; j<Number_of_nodes;j++){

            cout<<Matrix[i][j]<<'\t';
        }
        cout<<endl<<endl;
    }
    //cout << "\nShow Binary Finished" << endl;
}
//_________________________________ print Adjacency matrix __________________________________//

void Adjacency::Print_Network(){

    std::string s;
	std::stringstream ss;
	ss << Changer;
	ss >> s;
    std::string scheme ("Network/network_HH_");
    std::string hostname;
    std::string url;
    hostname = s+".txt" ;
    url = Address+scheme+ hostname;
    char* char_arr;
    string str_obj(url);
    char_arr = &str_obj[0];
    char * filename = char_arr;

    ofstream temp1(filename);
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
    cout << "\nPrint Network Finish" << endl;}


//_________________________________ print binery Adjacency matrix __________________________________//
void Adjacency::Print_Binary(){
    std::string s;
	std::stringstream ss;
	ss << Changer;
	ss >> s;
    std::string scheme ("Binery/bineryprint_HH_");
    std::string hostname;
    std::string url;
    hostname = s+".txt" ;
    url = Address+scheme+ hostname;
    char* char_arr;
    string str_obj(url);
    char_arr = &str_obj[0];
    char * filename = char_arr;
    ofstream temp2(filename);
    for (int i=0 ; i<Number_of_nodes;i++){
        for (int j=0 ; j<Number_of_nodes;j++){temp2<<Matrix[i][j]<<'\t';}
        temp2<<endl;}
    temp2.close();
    cout << "\nPrint Binary Finish" << endl;}
//________________________________ Build_array_rand __________________________________//
void Adjacency::Build_array_rand(){
    srand (time(NULL));
    int j;

    int i = 0;
    while( i < (Number_of_nodes*Number_of_nodes) ){
        arra[i] = (rand() % (Number_of_nodes*Number_of_nodes));
        for (j = 0; j < i; j++)
            if (arra[j] == arra[i])
                break;
            if ( ! (j < i) )
            i++;
    }
        for(j = 0; j <(Number_of_nodes*Number_of_nodes); j++){
       cout<<arra[j]<<endl;
    }
        cout<<"__________________"<<endl;
}
//________________________________ Build_onelink_step __________________________________//
void Adjacency::Build_onelink_step(){
        for(int i=Changer*100 ; i<((Changer+1)*100) ; i++){
        int x1=arra[i]%10;
        int x2=((arra[i]-x1)/10)%10;
        int a=((10*x2)+x1);
        int b=(arra[i]-((10*x2)+x1))/100;
        if(b==a){
        cout <<  Matrix[b][a]<<'\t'<<b<<'\t'<<a<<"||||||||||||||||||||||||||||||||||||||||||"<< endl;
            Matrix[b][a]=0;
        }else{

            Matrix[b][a]=1;

        cout <<  Matrix[b][a]<<'\t'<<b<<'\t'<<a<< endl;
        }
    }

}
#endif // ADJACENCY_MATRIX_H_INCLUDED
