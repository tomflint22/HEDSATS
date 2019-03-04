
#include<cmath>
#include<iostream>
#include<fstream>

#include "double_exponential_transform_header.h"
#include "Faddeeva.cpp"
#include "Analytical_Functions.cpp"




using namespace std;

int main()
{

ThermalA::Welcome();






double T_0=273.15;

double B=50.0;
double D=20.0;
double L=100.0;

double thermal_conductivity=0.0334;//0.0396;                                                   //
double specific_heat=605;//605;//480;                                                          //
double Density=7690*(1E-9);//7690*(1E-9);



double voltage=10.5;
double current=260.0;                                                                                 //
double v_arc=2.0;

double time_ramp_start=0.0;
double time_ramp_end=L/v_arc;
double ramp_up_rate=10.0;
double ramp_down_rate=10.0;

double y_i=0;

double aDE = 3.0;
double bDE = 2.5;
double c_rDE = 4.0;
double c_fDE = 1.0;
double d_g = 0.0;


double b_g=(B/2.0)-5.0;//B/2;//85;
double thermal_efficiency=0.6;


//Location to Probe Temperature

double probe_x = b_g-5.0;
double probe_y = 0.0;
double probe_z = 2.0*(L/8.0);




ofstream data;
data.open("data.dat");


cout<<
"For a domain size of "<<B<<" mm in the x-direction, "
<<D<<" mm in the y-direction"<<
"and "<<L<<" mm in the z-direction \n"<<
"HEDSATS will now compare the thermal response due to a DEC heat source with different boundary condition combinations. \n \n"<<endl;

cout<<"The output will be written to data.dat"<<endl;

data<<"time \t X22Y22Z22 \t X33Y33Z33 \t X11Y22Z11 \t X21Y22Z33 \t X21Y22Z11"<<endl;

for(double time=0.0;time<400.0;time+=1.0){

cout<<time<<" s"<<endl;

data<<time<<"\t"<<

ThermalA::Temp_Finite_DE(probe_x,probe_y,probe_z,aDE,bDE,c_rDE,c_fDE,time,thermal_conductivity,Density,specific_heat,v_arc,T_0,current,voltage,thermal_efficiency,B,D,L,b_g,d_g,
time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate,22,22,22)
<<"\t"<<
ThermalA::Temp_Finite_DE(probe_x,probe_y,probe_z,aDE,bDE,c_rDE,c_fDE,time,thermal_conductivity,Density,specific_heat,v_arc,T_0,current,voltage,thermal_efficiency,B,D,L,b_g,d_g,
time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate,33,33,33)
<<"\t"<<
ThermalA::Temp_Finite_DE(probe_x,probe_y,probe_z,aDE,bDE,c_rDE,c_fDE,time,thermal_conductivity,Density,specific_heat,v_arc,T_0,current,voltage,thermal_efficiency,B,D,L,b_g,d_g,
time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate,11,22,11)
<<"\t"<<
ThermalA::Temp_Finite_DE(probe_x,probe_y,probe_z,aDE,bDE,c_rDE,c_fDE,time,thermal_conductivity,Density,specific_heat,v_arc,T_0,current,voltage,thermal_efficiency,B,D,L,b_g,d_g,
time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate,21,22,33)
<<"\t"<<
ThermalA::Temp_Finite_DE(probe_x,probe_y,probe_z,aDE,bDE,c_rDE,c_fDE,time,thermal_conductivity,Density,specific_heat,v_arc,T_0,current,voltage,thermal_efficiency,B,D,L,b_g,d_g,
time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate,21,22,11)
//<<"\t"<<
//ThermalA::Temperature_convection_on_all_surfaces(probe_x,probe_y,probe_z,aDE,bDE,c_rDE,c_fDE,time,thermal_conductivity,Density,specific_heat,v_arc,T_0,current,voltage,
//thermal_efficiency,B,D,L,time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate,b_g,0.0)
<<endl;

}

//std::cout<<Temperature<<std::endl;


data.close();

return 0;
}



//#include "../../HEDSATSlib/Faddeeva.cpp"
//#include "../../HEDSATSlib/Analytical_Functions.cpp"
