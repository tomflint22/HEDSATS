#include <vector>
#include <algorithm>
#include<cmath>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<cctype>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include "double_exponential_transform_header.h"
#include "Faddeeva.cpp"
#include "Analytical_Functions.cpp"

#include <omp.h>



using namespace ThermalA;
using namespace std;

/////////////////////////////////////////////FUNCTION DEFENITIONS//////////////////////////////////////////////////
int number_of_lines_in_data_file(string file_name);
int number_of_columns_in_data_file(string file_name);
void assign_array_from_TC_file(string file_name, double **Experimental_array,int number_of_thermocouples);
void assign_TC_coordinates(string file_name, double **coordinate_array,int number_of_thermocouples);
int peak_array_index(double **array_to_index, int size_of_said_array,int array_column_identifier);
void populate_analytical_array_for_a_TC(double **array_to_populate,int array_TC_index_to_populate,int number_of_elements_in_array_to_populate, double t_gap,double **TC_coords_array,double a_DE_in,double b_DE_in,double cr_DE_in,double cf_DE_in,double a_EB_in,double b_EB_in,double cr_EB_in,double cf_EB_in,
double thermal_conductivity,double density,double specific_heat,double v_arc, double T_0, double current, double voltage,
double eff, double B, double D, double L,double b_g, double d_g,double y_i,double portion_of_heat_in_DE,
double t_start_ramp,double t_finish_ramp,double ramp_rate_start,double ramp_rate_end);
double Error_indication(double **experimental_array,double **calculated_array, int size_of_arrays,int num_TCs);
void Print_Final_Results(string file_path,double **Analytically_computed_array,double **Experimental_data_array,double **Thermocouple_coordinates_array,int num_TCs,
                int number_of_elements_in_data_arrays);
void evaluate (double **Analytically_computed_array,double **Experimental_data_array,double **Thermocouple_coordinates_array,int num_TCs,
                int number_of_elements_in_data_arrays,double t_gap,double aDE,double bDE,double c_rDE,double c_fDE,double aEB,double bEB,double c_rEB,double c_fEB,double d_g,double PortionOfHeatInSurfaceDE,
                double thermal_conductivity,double density,double specific_heat,double v_arc, double T_0, double current, double voltage,double eff, double B, double D, double L,double b_g,
                double y_i,double t_start_ramp,double t_finish_ramp,double ramp_rate_start,double ramp_rate_end);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{

int nthreads, tid;


/* Fork a team of threads giving them their own copies of variables */
#pragma omp parallel private(nthreads, tid)
  {

  /* Obtain thread number */
  tid = omp_get_thread_num();
  /* Only master thread does this */
  if (tid == 0)
    {
    nthreads = omp_get_num_threads();
    cout<<"HEDSATS is Running on "<<nthreads<<" threads"<<endl;
    }

  }  /* All threads join master thread and disband */





  ofstream variableEvolution,efficiency,cool_out;
string path;

char char_path_choice;


cout<<"Enter path to Working Directory (' ./ ' for current folder): "<<endl;
cin>>path;


do{
cout<<"Working Directory is set to: \t"<<path<<"\t \n Is this correct? (y/n)"<<endl;
cin>>char_path_choice;
if(char_path_choice!='y'&&char_path_choice!='Y'){
cout<<"Enter New Path :"<<endl;
cin>>path;
}
}
while(char_path_choice!='y'&&char_path_choice!='Y');
cout<<"Working Directory is set to: \t"<<path<<"\n"<<endl;

string Thermocouple_data=path+"DATA.dat";
string Thermocouple_coordinates=path+"Thermocouple_coordinates.dat";
//string variables_path=path+"varEvo.dat";
//string efficiency_path=path+"efficiency.dat";
string coolout_path=path+"cooling.dat";

//variableEvolution.open(variables_path);
//efficiency.open(efficiency_path);
cool_out.open(coolout_path);

int N_THERMOCOUPLES=number_of_columns_in_data_file(Thermocouple_coordinates);
cout<<"NUMBER OF THERMOCOUPLES IS : "<<number_of_columns_in_data_file(Thermocouple_coordinates)<<"\n"<<endl;

int num_lines_exp_1_i=number_of_lines_in_data_file(Thermocouple_data);//find number of temperature points
//int initial_exp_1_array_size=num_lines_exp_1_i-1;

cout<<"NUMBER OF DATA LINES IS "<<num_lines_exp_1_i<<endl;



double **Thermocouple_DATA;
Thermocouple_DATA=new double *[N_THERMOCOUPLES+1];
for(int i=0;i<N_THERMOCOUPLES+1;i++){
Thermocouple_DATA[i]=new double[num_lines_exp_1_i-1];
}
//  DECLARE EXPERIMENTAL DATA ARRAYS

assign_array_from_TC_file(Thermocouple_data,Thermocouple_DATA,N_THERMOCOUPLES);
////////Thermocouple_DATA[0 is time column 1-N are thermocouple data columns][rows of data]
double t_gap=Thermocouple_DATA[0][10]-Thermocouple_DATA[0][9];

double T_0=Thermocouple_DATA[1][0];  //  COULD HAVE DIFFERENT T_0 FOR EACH TC HERE
cout<<"INITIAL TEMPERATURE: "<<T_0<<endl;
///////////////////////////      ANALYSIS DATA - CONSTANTS ///////////////////////////////////////////

double B=170;
double D=30.0;
double L=300;
char plate_dims_choice;
do{
cout<<"Plate Dimensions Are set to: \n B = "<<B<<" mm \nD = "<<D<<" mm \nL = "<<L<<" mm \n Is this correct? (y/n)"<<endl;
cin>>plate_dims_choice;
if(plate_dims_choice!='y'&&plate_dims_choice!='Y'){
cout<<"Enter New Plate Dimensions :"<<endl;
cout<<"B :";
cin>>B;
cout<<"D :";
cin>>D;
cout<<"L :";
cin>>L;
}
}
while(plate_dims_choice!='y'&&plate_dims_choice!='Y');
cout<<"Plate Dimensions Are set to: \n B = "<<B<<" mm \nD = "<<D<<" mm \nL = "<<L<<" mm \n"<<endl;


double voltage=150;
double current=90;                                                                                 //
double v_arc=3.333333333;

char weld_params_choice;
do{
cout<<"Welding Parameters Are set to: \n voltage = "<<voltage<<" V \ncurrent = "<<current<<" A \nwelding travel velocity = "<<v_arc<<" mm/s \n Is this correct? (y/n)"<<endl;
cin>>weld_params_choice;
if(weld_params_choice!='y'&&weld_params_choice!='Y'){
cout<<"Enter New Welding Parameters :"<<endl;
cout<<"voltage :";
cin>>voltage;
cout<<"current :";
cin>>current;
cout<<"welding travel velocity :";
cin>>v_arc;
}
}
while(weld_params_choice!='y'&&weld_params_choice!='Y');
cout<<"Welding Parameters Are set to: \n voltage = "<<voltage<<" V \ncurrent = "<<current<<" A \nwelding travel velocity = "<<v_arc<<" mm/s \n"<<endl;

double time_ramp_start=0.0;
double time_ramp_end=L/v_arc;
double ramp_up_rate=8.0;
double ramp_down_rate=8.0;

//int number_of_elements_in_analytical_array=static_cast<int>(L/v_arc);
int number_of_elements_in_analytical_array=250;


double **Computed_Results;
Computed_Results=new double *[N_THERMOCOUPLES+1];
for(int i=0;i<N_THERMOCOUPLES+1;i++){
Computed_Results[i]=new double[number_of_elements_in_analytical_array];
}





double y_i=0;


double thermal_conductivity=0.0334;//0.0396;                                                   //
double specific_heat=605;//605;//480;                                                          //
double Density=7690*(1E-9);//7690*(1E-9);

char thermal_props_choice;
do{
cout<<"Material thermal properties Are set to: \n k = "<<thermal_conductivity<<" W/mm K \ncp = "<<specific_heat<<" J/kg K \nDensity = "<<Density<<" kg/mm^3 \n Is this correct? (y/n)"<<endl;
cin>>thermal_props_choice;
if(thermal_props_choice!='y'&&thermal_props_choice!='Y'){
cout<<"Enter New Welding Parameters :"<<endl;
cout<<"k :";
cin>>thermal_conductivity;
cout<<"specific_heat :";
cin>>specific_heat;
cout<<"Density :";
cin>>Density;
}
}
while(thermal_props_choice!='y'&&thermal_props_choice!='Y');
cout<<"Material thermal properties Are set to: \n k = "<<thermal_conductivity<<" W/mm K \ncp = "<<specific_heat<<" J/kg K \nDensity = "<<Density<<" kg/mm^3 \n"<<endl;



double b_g=84.9;//B/2;//85;


char b_g_choice;
do{
cout<<"b_g is set to: \n b_g = "<<b_g<<" mm \n Is this correct? (y/n)"<<endl;
cin>>b_g_choice;
if(b_g_choice!='y'&&b_g_choice!='Y'){
cout<<"Enter New Welding Parameters :"<<endl;
cout<<"b_g :";
cin>>b_g;
}
}
while(b_g_choice!='y'&&b_g_choice!='Y');
cout<<"b_g is set to: \n k = "<<b_g<<" mm \n"<<endl;


double **TC_coords;
TC_coords=new double *[N_THERMOCOUPLES];
for(int i=0;i<N_THERMOCOUPLES;i++){
TC_coords[i]=new double[3];
}


assign_TC_coordinates(Thermocouple_coordinates,TC_coords,N_THERMOCOUPLES);



int closest_TC;
cout<<"Enter the index of the closest thermocouple, (For temporal normalisation) ";
cin>>closest_TC;
int index_of_closest_TC=closest_TC;
//int index_of_closest_TC=2;

populate_analytical_array_for_a_TC(Computed_Results,index_of_closest_TC,number_of_elements_in_analytical_array,t_gap,TC_coords,0.5,0.5,
0.5,0.5,0.5,0.25,0.5,0.5,thermal_conductivity,Density,specific_heat,
v_arc,T_0,current,voltage,0.5,B,D,L,b_g,D*0.9,y_i,0.0,time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate);

int analytical_peak_index_closest=peak_array_index(Computed_Results,number_of_elements_in_analytical_array,index_of_closest_TC);
int experimental_peak_index_closest=peak_array_index(Thermocouple_DATA, num_lines_exp_1_i-1,index_of_closest_TC);
int peak_shift=experimental_peak_index_closest-analytical_peak_index_closest;

cout<<"experimental_peak: "<<experimental_peak_index_closest<<endl;
cout<<"analytical_peak: "<<analytical_peak_index_closest<<endl;
//cout<<"peak shift 1: "<<peak_shift_1<<endl;
cout<<"peak shift: "<<peak_shift<<endl;

cout<<"peak in analytical data is at array index "<<analytical_peak_index_closest<<endl;//find peak in experimental data
cout<<"corresponding to a time of "<<analytical_peak_index_closest*t_gap<<" seconds"<<endl;//find peak in experimental data


double **Thermocouple_DATA_shifted;
Thermocouple_DATA_shifted=new double *[N_THERMOCOUPLES+1];
for(int i=0;i<N_THERMOCOUPLES+1;i++){
Thermocouple_DATA_shifted[i]=new double[number_of_elements_in_analytical_array];
}


for(int k=0;k<number_of_elements_in_analytical_array-1;k+=1){
Thermocouple_DATA_shifted[0][k]=Thermocouple_DATA[0][k];
}
for(int i=1;i<N_THERMOCOUPLES+1;i++){//shift arrays by same time offset!!!! due to closest thermocouple on face where beam is incident gives delta T

for(int j=0;j<number_of_elements_in_analytical_array-1;j+=1){
Thermocouple_DATA_shifted[i][j]=Thermocouple_DATA[i][j+peak_shift];//-Thermocouple_DATA[i][peak_shift];
}
    }


double t_end=number_of_elements_in_analytical_array*t_gap;
cout<<"end time for analysis"<<t_end<<endl;
//double T_0=0.0;//EB_exp_data_2_shifted[0];    //  COULD HAVE DIFFERENT T_0 FOR EACH TC HERE
cout<<"peak in shifted experimental data is at array index "<<peak_array_index(Thermocouple_DATA_shifted, number_of_elements_in_analytical_array,index_of_closest_TC)<<endl;//find peak in experimental data

//cout<<"TEST"<<Thermocouple_DATA_shifted[5][0];
int farthest_TC;
cout<<"Enter the index of the farthest thermocouple from the heat source, (For efficiency calculation) ";
cin>>farthest_TC;
int index_of_farthest_TC=farthest_TC;

//int index_of_farthest_TC=5;

populate_analytical_array_for_a_TC(Computed_Results,index_of_farthest_TC,number_of_elements_in_analytical_array,t_gap,TC_coords,0.5,0.5,
0.5,0.5,0.5,0.25,0.5,0.5,thermal_conductivity,Density,specific_heat,
v_arc,T_0,current,voltage,1.0,B,D,L,b_g,D,y_i,0.0,time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate);


//for(int y=0;y<number_of_elements_in_analytical_array-1;y++){
//cout<<Computed_Results[1][y]<<"\t"<<Thermocouple_DATA_shifted[1][y]<<endl;
//}


double Temperature_rise_at_farthest_thermocouple_analytical=Computed_Results[index_of_farthest_TC][peak_array_index(Computed_Results,number_of_elements_in_analytical_array,index_of_farthest_TC)]-T_0;

double Temperature_rise_at_farthest_thermocouple_experimental=Thermocouple_DATA_shifted[index_of_farthest_TC][peak_array_index(Thermocouple_DATA_shifted,number_of_elements_in_analytical_array,index_of_farthest_TC)]-T_0;
cout<<"deltaT far TC= "<<Temperature_rise_at_farthest_thermocouple_experimental<<endl;

cout<<"Temperature Rise at farthest thermocouple is:  "<<Temperature_rise_at_farthest_thermocouple_experimental<<endl;
cout<<"Temperature Rise at farthest thermocouple Analytical  is:  "<<Temperature_rise_at_farthest_thermocouple_analytical<<endl;

//  FIND FARTHEST THERMOCOUPLE AND THEN FIND EFFICIENCY
double thermal_efficiency=Temperature_rise_at_farthest_thermocouple_experimental/Temperature_rise_at_farthest_thermocouple_analytical;

cout<<"THERMAL EFFICIENCY THEREFORE IS "<<thermal_efficiency<<endl;

//efficiency<<"thermal_efficiency: "<<thermal_efficiency<<endl;




double aDE = 3.0;
double bDE = 1.0;
double c_rDE = 2.0;
double c_fDE = 2.0;
double aEB = 1.2;
double bEB = 2.0;
double c_rEB = 3.0;
double c_fEB = 1.0;
double d_g = 28.0;
double PortionOfHeatInSurfaceDE = 0.005;



char Heat_Source_Geometry_choice;
do{
cout<<"Heat Source Parameters Are set to:"<<endl;
cout<<"aDE = "<<aDE<<" mm"<<endl;
cout<<"bDE = "<<bDE<<" mm"<<endl;
cout<<"c_rDE = "<<c_rDE<<" mm"<<endl;
cout<<"c_fDE = "<<c_fDE<<" mm"<<endl;
cout<<"aEB = "<<aEB<<" mm"<<endl;
cout<<"bEB = "<<bEB<<" mm"<<endl;
cout<<"c_rEB = "<<c_rEB<<" mm"<<endl;
cout<<"c_fEB = "<<c_fEB<<" mm"<<endl;
cout<<"d_g = "<<d_g<<" mm"<<endl;
cout<<"f [portion of heat in surface double ellipsoidal component] = "<<PortionOfHeatInSurfaceDE<<endl;




  evaluate (Computed_Results,Thermocouple_DATA_shifted,TC_coords,N_THERMOCOUPLES,number_of_elements_in_analytical_array,t_gap,
  aDE,bDE,c_rDE,c_fDE,aEB,bEB,c_rEB,c_fEB,d_g,PortionOfHeatInSurfaceDE,
  thermal_conductivity,Density,specific_heat,v_arc,T_0,current,voltage,thermal_efficiency,B,D,L,b_g,y_i,time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate);


Print_Final_Results(path,Computed_Results,Thermocouple_DATA_shifted,TC_coords,N_THERMOCOUPLES,
                number_of_elements_in_analytical_array);


cout<<"Are You Happy With The fit? (y/n)"<<endl;
cin>>Heat_Source_Geometry_choice;
if(Heat_Source_Geometry_choice!='y'&&Heat_Source_Geometry_choice!='Y'){
cout<<"Enter New Parameters :"<<endl;
cout<<"aDE :";
cin>>aDE;
cout<<"bDE :";
cin>>bDE;
cout<<"c_rDE :";
cin>>c_rDE;
cout<<"c_fDE :";
cin>>c_fDE;
cout<<"aEB :";
cin>>aEB;
cout<<"bEB :";
cin>>bEB;
cout<<"c_rEB :";
cin>>c_rEB;
cout<<"c_fEB :";
cin>>c_fEB;
cout<<"d_g :";
cin>>d_g;
cout<<"f [portion of heat in surface double ellipsoidal component] :";
cin>>PortionOfHeatInSurfaceDE;

}
}
while(Heat_Source_Geometry_choice!='y'&&Heat_Source_Geometry_choice!='Y');



cout<<"Writing Cooling Profile"<<endl;

cool_out<<"time";
for(double x_out=b_g;x_out<=B;x_out+=1.0){
cool_out<<"\t"<<"x="<<x_out;
}
cool_out<<endl;



for(double time=0;time<=60;time+=1.0){
cool_out<<time;
for(double x_out=b_g;x_out<=B;x_out+=1.0){

double TempTemp=Temperature_from_combined_DECbeam_DE_heat_source(x_out,D/2.0,L/2.0,aDE,bDE,c_rDE,c_fDE,
aEB,bEB,c_rEB,c_fEB,time,thermal_conductivity,Density,specific_heat,
v_arc,T_0,current,voltage,thermal_efficiency,B,D,L,b_g,d_g,y_i,PortionOfHeatInSurfaceDE,time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate,21,22,22);

cool_out<<"\t"<<TempTemp;


}
cool_out<<"\n";

}










delete[] Thermocouple_DATA;
delete[] Computed_Results;
delete[] Thermocouple_DATA_shifted;
delete[] TC_coords;

cool_out.close();
//efficiency.close();
variableEvolution.close();


return 0;
}


//FUNCTIONS

int number_of_lines_in_data_file(string file_name){
string line;
  fstream file_data;
  cout<<"\t \t ******************** \n  counting number of lines in "<<file_name<<"\n \t \t ********************"<<endl;
  file_data.open(file_name, std::fstream::in);
int line_counter=0;
  while ( getline (file_data,line) )
    {
    line_counter+=1;
    }
    file_data.close();

return line_counter;
}

int number_of_columns_in_data_file(string file_name){
string line,temp;
  fstream file_data;
  cout<<"\t \t ******************** \n  counting number of columns in "<<file_name<<"\n \t \t ********************"<<endl;
  file_data.open(file_name, std::fstream::in);

stringstream ss;
int ncols=0;
getline(file_data, line);//discard first line

      getline(file_data, line);
    ss.clear();
    ss << line;

     while (ss >> temp)
        {
        ncols++;
         //cout << temp << " ";
        }

    file_data.close();

return ncols;
}

void assign_array_from_TC_file(string file_name, double **Experimental_array,int number_of_thermocouples){
cout<<"\t \t ********************"<<endl;
cout<<"Assigning thermocouple data to array construct from data file "<<file_name<<endl;
cout<<"\t \t ********************"<<endl;
string line,dummy_line;
string::size_type sz;
fstream file_data;
file_data.open(file_name,std::fstream::in);
getline(file_data,dummy_line);

int array_counter=0;
    while ( getline (file_data,line) )
    {
istringstream ss(line);


      //  double time=stod(line,&sz);
      //  double value=stod(line.substr(sz));
      for(int T=0;T<number_of_thermocouples+1;T++){
  ss>>Experimental_array[T][array_counter];
}

array_counter+=1;
      }

//for(int y=0;y<array_counter;y++){
//cout<<Experimental_array[9][y]<<endl;
//}


file_data.close();

}
void assign_TC_coordinates(string file_name, double **coordinate_array,int number_of_thermocouples){
cout<<"\t \t ********************"<<endl;
cout<<"Reading Thermocouple Co-Ordinates from file: "<<file_name<<endl;
cout<<"\t \t ********************"<<endl;
string line,dummy_line;
string::size_type sz;
fstream file_data;
file_data.open(file_name,std::fstream::in);
getline(file_data,dummy_line);

int array_counter=0;
    while ( getline (file_data,line) )
    {
istringstream ss(line);


      //  double time=stod(line,&sz);
      //  double value=stod(line.substr(sz));
      for(int T=0;T<number_of_thermocouples;T++){
  ss>>coordinate_array[T][array_counter];
}

array_counter+=1;
      }

//for(int y=0;y<array_counter;y++){
//cout<<Experimental_array[9][y]<<endl;
//}


file_data.close();

}
int peak_array_index(double **array_to_index, int size_of_said_array,int array_column_identifier){//returns array index of peak value
int index=0;
double max_value=0;
//    double current_max_value=array_to_index[0];
for(int i=0;i<size_of_said_array-1;i+=1){

double current_value=array_to_index[array_column_identifier][i];
double temp_value2=array_to_index[array_column_identifier][i+1];
    if(current_value>=max_value){
    max_value=current_value;
    index=i;
    }
//cout<<max_value<<endl;

//cout<<index<<endl;
}
return index;
}
//  Populate array from t=0
void populate_analytical_array_for_a_TC(double **array_to_populate,int array_TC_index_to_populate,int number_of_elements_in_array_to_populate, double t_gap,double **TC_coords_array,double a_DE_in,double b_DE_in,double cr_DE_in,double cf_DE_in,double a_EB_in,double b_EB_in,double cr_EB_in,double cf_EB_in,
double thermal_conductivity,double density,double specific_heat,double v_arc, double T_0, double current, double voltage,
double eff, double B, double D, double L,double b_g, double d_g,double y_i,double portion_of_heat_in_DE,
double t_start_ramp,double t_finish_ramp,double ramp_rate_start,double ramp_rate_end){

int k=0;
cout<<"\n Populating Analytical Array for Coordinates: ("<<TC_coords_array[array_TC_index_to_populate-1][0]<<","<<TC_coords_array[array_TC_index_to_populate-1][1]<<","<<
TC_coords_array[array_TC_index_to_populate-1][2]<<") \n"<<endl;
//////      assign analytical values

for(double time=0.0;time<number_of_elements_in_array_to_populate*t_gap;time+=t_gap){

//EB_analytical_time[k]=time;

array_to_populate[array_TC_index_to_populate][k]=Temperature_from_combined_DECbeam_DE_heat_source(TC_coords_array[array_TC_index_to_populate-1][0],TC_coords_array[array_TC_index_to_populate-1][1],
TC_coords_array[array_TC_index_to_populate-1][2],a_DE_in,b_DE_in,cr_DE_in,cf_DE_in,a_EB_in,b_EB_in,cr_EB_in,cf_EB_in,time,thermal_conductivity,density,specific_heat,v_arc,T_0,current,
voltage,eff,B,D,L,b_g,d_g,y_i,portion_of_heat_in_DE,t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end,21,22,22);

k++;
}

}
double Error_indication(double **experimental_array,double **calculated_array, int number_of_elements_in_data_arrays,int num_TCs){

double running_total=0.0;
for(int i=0;i<number_of_elements_in_data_arrays;i++){

for(int track=1;track<num_TCs;track++){
running_total+=abs(calculated_array[track][i]-experimental_array[track][i]);
}

}

return running_total;
}
void Print_Final_Results(string file_path,double **Analytically_computed_array,double **Experimental_data_array,double **Thermocouple_coordinates_array,int num_TCs,
                int number_of_elements_in_data_arrays){
ofstream final_output;
string output_path=file_path+"Output.dat";
final_output.open(output_path);
cout<<"OUTPUTTING RESULTS TO FILE"<<endl;

final_output<<"time \t";

for(int title=0;title<num_TCs;title++){
final_output<<"Experimental("<<Thermocouple_coordinates_array[title][0]<<","<<Thermocouple_coordinates_array[title][1]<<","<<Thermocouple_coordinates_array[title][2]<<")\t"<<
        "Computed("<<Thermocouple_coordinates_array[title][0]<<","<<Thermocouple_coordinates_array[title][1]<<","<<Thermocouple_coordinates_array[title][2]<<")\t";

}
final_output<<endl;
    int k=0;
//////      assign analytical values
for(int time=0;time<number_of_elements_in_data_arrays-1;time++){
final_output<<Experimental_data_array[0][time];
for(int track=1;track<=num_TCs;track++){

final_output<<"\t"<<Experimental_data_array[track][k]<<"\t"<<Analytically_computed_array[track][k];


}
final_output<<endl;


k++;
}






final_output.close();
}
void evaluate (double **Analytically_computed_array,double **Experimental_data_array,double **Thermocouple_coordinates_array,int num_TCs,
                int number_of_elements_in_data_arrays,double t_gap, double aDE,double bDE,double c_rDE,double c_fDE,double aEB,double bEB,double c_rEB,double c_fEB,double d_g,double PortionOfHeatInSurfaceDE,
                double thermal_conductivity,double density,double specific_heat,double v_arc, double T_0, double current, double voltage,double eff, double B, double D, double L,double b_g,
                double y_i,double t_start_ramp,double t_finish_ramp,double ramp_rate_start,double ramp_rate_end)

//****************************************************************************
{

cout<<"HEDSATS is Evaluating Semi-Analytical Solution"<<endl;
//  int member;
//  int i;
//  double x[NVARS+1];
//  for ( member = 0; member < POPSIZE; member++ )
//  {
//    for ( i = 0; i < NVARS; i++ )
//    {
//      x[i+1] = population[member].gene[i];
//    }

//cout<<"x1: "<<x[1]<<"\nx2: "<<x[2]<<"\nx3: "<<x[3]<<"\nx4: "<<x[4]<<"\nx5: "<<x[5]<<"\nx6: "<<x[6]<<"\nx7: "<<x[7]<<"\nx8: "<<x[8]<<"\nx9: "<<x[9]<<"\nx10: "<<x[10]<<endl;
    int k=0;
//////      assign analytical values
for(double time=0;time<number_of_elements_in_data_arrays*t_gap;time+=t_gap){
#pragma omp parallel for
for(int track=1;track<=num_TCs;track++){

Analytically_computed_array[track][k]=Temperature_from_combined_DECbeam_DE_heat_source(Thermocouple_coordinates_array[track-1][0],Thermocouple_coordinates_array[track-1][1],Thermocouple_coordinates_array[track-1][2],
aDE,bDE,c_rDE,c_fDE,aEB,bEB,c_rEB,c_fEB,time,thermal_conductivity,density,specific_heat,v_arc,T_0,current,voltage,eff,B,D,L,b_g,d_g,y_i,PortionOfHeatInSurfaceDE,t_start_ramp,
t_finish_ramp,ramp_rate_start,ramp_rate_end,21,22,22);

//double aDE,double bDE,double c_rDE,double c_fDE,double aEB,double bEB,double c_rEB,double c_fEB,double d_g,double PortionOfHeatInSurfaceDE
//cout<<"\n  TEST \n"<<Analytically_computed_array[track][k]<<"\n"<<endl;

}
//
////array_to_populate1[k]=Temperature_from_combined_DECbeam_DE_heat_source(TC1x_in,TC1y_in,TC1z_in,x[5],x[6],x[7],
////x[8],x[1],x[2],x[3],x[4],time,thermal_conductivity,density,specific_heat,v_arc,T_0,current,
////voltage,eff,B,D,L,b_g,x[10],y_i,x[9],t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end);
//


k++;
}

//double running_total=0.0;
//for(int i=0;i<number_of_elements_in_data_arrays;i++){
//
//for(int track=1;track<num_TCs;track++){
//running_total+=abs(Analytically_computed_array[track][i]-Experimental_data_array[track][i]);
//}
////running_total+=(pow(array_to_populate2[i]-experimental_array2[i],2)/array_to_populate2[i])+(pow(array_to_populate3[i]-experimental_array3[i],2)/array_to_populate3[i])+
////                (pow(array_to_populate6[i]-experimental_array6[i],2)/array_to_populate6[i])+(pow(array_to_populate7[i]-experimental_array7[i],2)/array_to_populate7[i]);
//
//
////cout<<"\n  TEST \n"<<running_total<<"\n"<<endl;
//}
cout<<"Least Squares Sum = "<<Error_indication(Experimental_data_array,Analytically_computed_array,number_of_elements_in_data_arrays,num_TCs)<<endl;;
//cout<<" Something Like Least Squares Data Fit: "<<running_total<<endl;

//    population[member].fitness =1/running_total;
//  }
  return;
}

