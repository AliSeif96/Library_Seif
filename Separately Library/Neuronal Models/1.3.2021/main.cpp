//{about code
/************************************************************************************************/
/*** Topic: Read a 2D matrix from a .txt file and save it in a 2D array in C++                ***/
/***	                                                                             Ali-Seif ***/
/*** Version Release 17.12 rev 11256                                                          ***/
/*** Date: 1/2/2021                                                                          ***/
/*** Code implemented in CodeBlocks C++ compiler (v. 17.12),                                  ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
//}

//{include
#include <iostream>                                                             //for cout
#include <sstream>                                                              //for stringstream
#include <fstream>                                                              //for ofstream
#include <string>                                                               //for std::string
#include <cstdlib>                                                              //for exit(1)
#include <vector>                                                               //for vector <double> column;
#include <Neuronal_models.h>
using namespace std;                                                            //for Standard program
//}

//{defines
//##############################################################
//####                                                      ####
//####                 size of matrix                       ####
//####                                                      ####
//##############################################################
#define COLUMNS 200
//}

int main()
{

//{Read a 2D matrix from a .txt file and save it in a 2D array in C++


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//##############################################################                                                                //$$$
//####                                                      ####                                                                //$$$
//####        read file and create file for print           ####                                                                //$$$
//####                                                      ####                                                                //$$$
//##############################################################                                                                //$$$
    //ofstream temp("E:/programs/sync/read_data_aida/cout.txt");                  //address for print                           //$$$
    float Changer=20;                                                                                                           //$$$
    std::string s;                                                                                                              //$$$
	std::stringstream ss;                                                                                                       //$$$
	ss << Changer;                                                                                                              //$$$
	ss >> s;                                                                                                                    //$$$
    std::string scheme ("E:/programs/sync/read_data_aida/matrix-SF-DAG1src-"); //read                                           //$$$
    std::string hostname;                                                                                                       //$$$
    std::string url;                                                                                                            //$$$
    hostname = s+".txt" ;                                                                                                       //$$$
    url =scheme+ hostname;                                                                                                      //$$$
    char* char_arr;                                                                                                             //$$$
    string str_obj(url);                                                                                                        //$$$
    char_arr = &str_obj[0];                                                                                                     //$$$
    char * filename1 = char_arr;                                                                                                //$$$
    vector< vector <double> > data;                                                                                             //$$$
    string filename = filename1;                                                                                                //$$$
    ifstream ifile(filename.c_str());                                                                                           //$$$
//##############################################################                                                                //$$$
//####                                                      ####                                                                //$$$
//####                   read file .txt                     ####                                                                //$$$
//####                                                      ####                                                                //$$$
//##############################################################                                                                //$$$
    if (ifile.is_open()) {                                                                                                      //$$$
        double num;                                                                                                             //$$$
        vector <double> numbers_in_line;                                                                                        //$$$
        while (ifile >> num) {                                                                                                  //$$$
            numbers_in_line.push_back(num);                                                                                     //$$$
            if (numbers_in_line.size() == COLUMNS) {                                                                            //$$$
                data.push_back(numbers_in_line);                                                                                //$$$
                numbers_in_line.clear();                                                                                        //$$$
            }                                                                                                                   //$$$
        }                                                                                                                       //$$$
    }                                                                                                                           //$$$
    else {                                                                                                                      //$$$
        cerr << "There was an error opening the input file!\n";                                                                 //$$$
        exit(1);                                                                                                                //$$$
    }                                                                                                                           //$$$
//##############################################################                                                                //$$$
//####                                                      ####                                                                //$$$
//####            print or cout one column                  ####                                                                //$$$
//####                                                      ####                                                                //$$$
//##############################################################                                                                //$$$
    for(int j=1;j<=COLUMNS;j++){                                                                                                //$$$
        vector <double> column;                                                                                                 //$$$
        int col = j;            //[1 200];//example: the 2nd column                                                             //$$$
        for (int i = 0; i < data.size(); ++i) {                                                                                 //$$$
            column.push_back(data[i][col - 1]);                                                                                 //$$$
            Matrix_N[col - 1][i]=column[i];                                             //Matrix_N[-][|]                        //$$$
            //temp<<j<<'\t'<<i<<'\t'<< column[i] << endl;                                                                       //$$$
        }                                                                                                                       //$$$
    }                                                                                                                           //$$$
    ifile.close();                                                                                                              //$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//}

//{example of matrix that read


//===================================================================================================================================
//##############################################################                                                                //===
//####                                                      ####                                                                //===
//####            example of matrix that read               ####                                                                //===
//####                                                      ####                                                                //===
//##############################################################                                                                //===
     for(int j=0;j<10;j++){                                                                                                     //===
        for (int i = 0; i < 10; ++i) {                                                                                          //===
            cout<<Matrix_N[i][j];               //Matrix_N[-][|]                                                                //===
        }                                                                                                                       //===
        cout<<endl;                                                                                                             //===
    }                                                                                                                           //===
//===================================================================================================================================
//}

//{conditions


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Length_of_steps=0.1;                                            //(0,1] ms                                                  //+++
    Final_time=2500;                                                                                                            //+++
    External_electric_current_HH=0.3;                                                                                           //+++
    phi_wb=7;                                                                                                                   //+++
    Coupling_strength_HH=0.15;//                                         //[0.0003,]                                            //+++
    Number_of_nodes_N=200;                                                                                                      //+++
    Changer_N=1;                                                                                                                //+++
    //class Hodgkin_Huxley hh;                                                                                                  //+++
    //hh.Beginning();                                                                                                           //+++
    //hh.Print_Data();                                                                                                          //+++
    class Wang_Buzsaki wb;                                                                                                      //+++
    wb.Beginning();                                                                                                             //+++
    wb.Print_Data();                                                                                                            //+++
    cout << "\n2assssssssssssssssssss56a" << endl;                                                                              //+++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//}
    return 0;
}
