#include <iostream>
#include<Adjacency_matrix.h>
#include<Neuronal_models.h>
using namespace std;

int main()
{
    Number_of_nodes=10;
    Possibility_of_connection=0.5;
    Address="E:/programs/add corent/first/";
    Changer=5;
    class Adjacency adj;
    adj.Build_matrix_zero();
    adj.Erdos_Renyi();
    adj.Show_Binary();
    adj.Print_Network();
    adj.Print_Binary();
    //adj.Build_array_rand();
    //adj.Build_onelink_step();

    for (int i=0 ; i<Number_of_nodes;i++){
       for (int j=0 ; j<Number_of_nodes;j++){
            Matrix_N[i][j]=Matrix[i][j];
       }
    }

    Number_of_nodes_N=10;
    Changer_N=5;
    Address="E:/programs/add corent/first/";
    Length_of_steps=0.1;                                            //(0,1] ms
    Final_time=2500;
    External_electric_current_HH=9.741;//9.741                                  //[9.7,]
    Coupling_strength_HH=0.001;//                                         //[0.0003,]
    class Hodgkin_Huxley hh;
    hh.Beginning();
    //hh.openfile();
    hh.Print_Data();

    cout << "\nrun done" << endl;
    return 0;
}
