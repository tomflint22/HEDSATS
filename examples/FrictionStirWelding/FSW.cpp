
#include<cmath>
#include<iostream>

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
double L=200;

double thermal_conductivity=0.0334;//0.0396;                                                   //
double specific_heat=605;//605;//480;                                                          //
double Density=7690*(1E-9);//7690*(1E-9);



double voltage=150;
double current=90;                                                                                 //
double vel=5.0;


double y_i=0.0;

double a = 3.0;
double b = 1.0;
double rr = 7.0;
double rf = 5.0;

double d_g = 18.0;
double l_g = 10.0;

double b_g=B/2.0;//B/2;//85;
double thermal_efficiency=0.5;


double time_ramp_start=2.0;
double time_ramp_end=L/vel;
double ramp_up_rate=8.0;
double ramp_down_rate=8.0;






double Temperature=0.0;//
double Temperature_dwell=0.0;//
double Temperature_motion=0.0;//


cout<<"In this example HEDSATS is used to determine the temperature due to a FSW process where there is an initial 2s dwell period. Showing how to use multiple heat source models"<<endl;
cout<<"After the 2s dwell, the 'dwell' heat source ramps out at the same rate as the FSW heat source in motion"<<endl;


cout<<"Domain geometry \n ******************"<<endl;
cout<<"B= "<<B<<" mm"<<endl;
cout<<"D= "<<D<<" mm"<<endl;
cout<<"L= "<<L<<" mm"<<endl;
cout<<"******************"<<endl;

cout<<"Tool Input Power = "<<voltage*current*thermal_efficiency<<endl;


double probe_x,probe_y,probe_z,time;

cout<<"In this example the HEDSATS librray will be utilised to probe the temperature in a domain experiencing EB welding."<<endl;
cout<<"enter x-coordinate"<<endl;
cin>>probe_x;
cout<<"enter y-coordinate"<<endl;
cin>>probe_y;
cout<<"enter z-coordinate"<<endl;
cin>>probe_z;
cout<<"time"<<endl;
cin>>time;


Temperature_dwell=ThermalA::Temperature_Finite_4Q_FSW_DWELL(probe_x,probe_y,probe_z,a,b,rr,rf,time,thermal_conductivity,Density,specific_heat,
0.0,current,voltage,thermal_efficiency,B,D,L,b_g,d_g,l_g,y_i,0.0,2.0,ramp_up_rate,ramp_down_rate);


Temperature_motion=ThermalA::Temperature_Finite_4Q_FSW(probe_x,probe_y,probe_z,a,b,rr,rf,time,thermal_conductivity,Density,specific_heat,
vel,0.0,current,voltage,thermal_efficiency,B,D,L,b_g,d_g,y_i,2.0,time_ramp_end,ramp_up_rate,ramp_down_rate);

Temperature=T_0+Temperature_dwell+Temperature_motion;





std::cout<<Temperature<<std::endl;




return 0;
}

