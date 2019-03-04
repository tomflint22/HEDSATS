
#include<cmath>
#include<iostream>

//#include <omp.h>

#include "double_exponential_transform_header.h"
#include "Faddeeva.cpp"
#include "Analytical_Functions.cpp"



using namespace std;

int main()
{
ThermalA::Welcome();






double T_0=105.0;

double B=170;
double D=30.0;
double L=300;

double thermal_conductivity=0.0334;//0.0396;                                                   //
double specific_heat=605;//605;//480;                                                          //
double Density=7690*(1E-9);//7690*(1E-9);



double voltage=150;
double current=90;                                                                                 //
double v_arc=3.333333333;

double time_ramp_start=0.0;
double time_ramp_end=L/v_arc;
double ramp_up_rate=8.0;
double ramp_down_rate=8.0;

double y_i=0;

double aDE = 3.0;
double bDE = 1.0;
double c_rDE = 2.0;
double c_fDE = 2.0;
double aEB = 1.2;
double bEB = 3.5;
double c_rEB = 4.0;
double c_fEB = 1.0;
double d_g = 28.0;
double PortionOfHeatInSurfaceDE = 0.02;

double b_g=84.5;//B/2;//85;
double thermal_efficiency=0.58;


//double probe_x = b_g;
//double probe_y = 0.0;
//double probe_z = 0.0;
//double time=2000.0;


double Temperature=0.0;//


int xBoundary=22;
int yBoundary=22;
int zBoundary=22;

double probe_x,probe_y,probe_z,time;

cout<<"In this example the HEDSATS librray will be utilised to probe the temperature in a domain experiencing EB welding. \n"<<endl;
cout<<"Domain geometry \n ******************"<<endl;
cout<<"B= "<<B<<" mm"<<endl;
cout<<"D= "<<D<<" mm"<<endl;
cout<<"L= "<<L<<" mm"<<endl;
cout<<"******************"<<endl;

cout<<"enter x-coordinate"<<endl;
cin>>probe_x;
cout<<"enter y-coordinate"<<endl;
cin>>probe_y;
cout<<"enter z-coordinate"<<endl;
cin>>probe_z;
cout<<"time"<<endl;
cin>>time;
//Return Temperature due to combined DEC - DE model, where DE is additionally applied at surface to represent surface effects

Temperature=ThermalA::Temperature_from_combined_DECbeam_DE_heat_source(probe_x,probe_y,probe_z,aDE,bDE,c_rDE,c_fDE,
aEB,bEB,c_rEB,c_fEB,time,thermal_conductivity,Density,specific_heat,
v_arc,T_0,current,voltage,thermal_efficiency,B,D,L,b_g,d_g,y_i,PortionOfHeatInSurfaceDE,time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate,
xBoundary,yBoundary,zBoundary);



//Return Temperature due to DEC model, where a_i==a, cr_i==cr, cfi==cf etc for heat source orthogonality
//ThermalA::Temp_Finite_EB(probe_x,probe_y,probe_z,aDE,bDE,c_rDE,c_fDE,time,thermal_conductivity,Density,specific_heat,v_arc,T_0,current,voltage,thermal_efficiency,B,D,L,b_g,d_g,0.0,
//time_ramp_start,time_ramp_end,ramp_up_rate,ramp_down_rate,xBoundary,yBoundary,zBoundary);


std::cout<<"Temperature at ("<<probe_x<<","<<probe_y<<","<<probe_z<<") and "<<time<<" seconds ="<<Temperature<<std::endl;




return 0;
}
