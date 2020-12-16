#include <iostream>
#include "Seif.h"

using namespace std;

int main()
{

    //##############################################################
    //####                                                      ####
    //####               class Adjacency matrix                 ####
    //####                                                      ####
    //##############################################################
    Number_of_nodes=10;                                              //[1,100]
    Possibility_of_connection=0.9;                                  //[0,1]
    class Adjacency adj;
    adj.Build_matrix();
    //adj.Show_Binary();
    adj.Print_Binary();
    adj.Print_Network();
    //##############################################################
    //####                                                      ####
    //####            class Hodgkin_Huxley Excitatory           ####
    //####                                                      ####
    //##############################################################
    External_electric_current_HH=9.7;                                  //[9.7,]
    Storage_start_time=2000;                                        //The execution time of the program is from zero
    Length_of_steps=0.1;                                            //(0,1] ms
    Final_time=2500;                                                //ms
    Coupling_strength_HH=0.05;                                         //[0.0003,]
    class Hodgkin_Huxley HH;
    HH.Start();
    HH.Print_Data();
    //HH.Show_Data();
    //##############################################################
    //####                                                      ####
    //####             class Wang-Buzsaki Inhibitory            ####
    //####                                                      ####
    //##############################################################
    External_electric_current_WB=5;                                  //[0.16,]
    Coupling_strength_WB=5;                                         //[0,0.3]
    class Wang_Buzsaki WB;
    WB.Start();
    WB.Print_Data();
    //WB.Show_Data();

    return 0;
}
