#include <iostream>
#include "Seif.h"

using namespace std;

int main()
{


    Number_of_nodes=4;                      //[1,100]
    Possibility_of_connection=0.9;          //[0,1]
    External_electric_current=6.5;          //[0,7]
    Storage_start_time=2000;                //The execution time of the program is from zero
    Length_of_steps=0.1;                    //(0,1] ms
    Final_time=2500;                        //ms



    class Adjacency adj;
    adj.Build_matrix();
    //adj.Show_Binary();
    adj.Print_Binary();
    adj.Print_Network();


    class Hodgkin_Huxley HH;

    HH.Start();
    HH.Print_Data();
    //HH.Show_Data();



    return 0;
}
