

#include "Analytical_Functions.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <complex>

using namespace std;


#define pi 3.141592654
#define eulers 2.718281828459045
#define kb 1.380648E-23
#define Na 6.022E23
#define stef_B_const 5.670373E-8


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


                                 //               Class  Function  Defenitions                  //

namespace ThermalA {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FINITE_plate_to_be_integrated_DE_Heat_source::set_values (double x_in,double y_in,double z_in,double a_in,double b_in,double cr_in,double cf_in,double t_in,double k_conductivity_in,
double density_in,double cp_steel_in,double v_arc_in,double B_in, double D_in, double L_in,double b_g_in, double d_g_in,double t_start_ramp_in,double t_finish_ramp_in,
double ramp_rate_start_in,double ramp_rate_end_in,int xBC_in,int yBC_in,int zBC_in) {
  x = x_in;
  y = y_in;
  z = z_in;
  a = a_in;
  b = b_in;
  cr = cr_in;
  cf = cf_in;
  t = t_in;
  k_conductivity = k_conductivity_in;
  density = density_in;
  cp_steel = cp_steel_in;
  v_arc = v_arc_in;
  B=B_in;
	  D=D_in;
	  L=L_in;
    b_g=b_g_in;
    d_g=d_g_in;
    t_start_ramp=t_start_ramp_in;
    t_finish_ramp=t_finish_ramp_in;
    ramp_rate_start=ramp_rate_start_in;
    ramp_rate_end=ramp_rate_end_in;
    xBC=xBC_in;
    yBC=yBC_in;
    zBC=zBC_in;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void FINITE_plate_to_be_integrated_EB::set_values (double x_in,double y_in,double z_in,double a_in,double b_in,double cr_in,double cf_in,double t_in,double k_conductivity_in,
 double density_in,double cp_steel_in,double v_arc_in,double B_in, double D_in, double L_in,double b_g_in,double d_g_in,double y_i_in,double t_start_ramp_in,double t_finish_ramp_in,
 double ramp_rate_start_in,double ramp_rate_end_in,int xBC_in,int yBC_in,int zBC_in) {
  x = x_in;
  y = y_in;
  z = z_in;
  a = a_in;
  b = b_in;
  cr = cr_in;
  cf = cf_in;
  t = t_in;
  k_conductivity = k_conductivity_in;
  density = density_in;
  cp_steel = cp_steel_in;
  v_arc = v_arc_in;
  B=B_in;
	  D=D_in;
	  L=L_in;
	b_g=b_g_in;
	d_g=d_g_in;
	y_i=y_i_in;
	t_start_ramp=t_start_ramp_in;
    t_finish_ramp=t_finish_ramp_in;
    ramp_rate_start=ramp_rate_start_in;
    ramp_rate_end=ramp_rate_end_in;
        xBC=xBC_in;
    yBC=yBC_in;
    zBC=zBC_in;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FINITE_plate_to_be_integrated_FSW::set_values (double x_in,double y_in,double z_in,double r_in,double b_in,double t_in,double k_conductivity_in,double density_in,double cp_steel_in,double v_arc_in,double B_in, double D_in, double L_in,double b_g_in,double d_g_in,double y_i_in,double t_start_point_RAMP_in,double t_end_point_RAMP_in,double uprate_RAMP_in, double downrate_RAMP_in){
  x = x_in;
  y = y_in;
  z = z_in;
  r = r_in;
  b = b_in;
  t = t_in;
  k_conductivity = k_conductivity_in;
  density = density_in;
  cp_steel = cp_steel_in;
  v_arc = v_arc_in;
  B=B_in;
    D=D_in;
	  L=L_in;
	b_g=b_g_in;
	d_g=d_g_in;
	y_i=y_i_in;
	t_start_point_RAMP=t_start_point_RAMP_in;
	t_end_point_RAMP=t_end_point_RAMP_in;
	uprate_RAMP=uprate_RAMP_in;
	downrate_RAMP=downrate_RAMP_in;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FINITE_plate_to_be_integrated_FSW_4Q::set_values(double x_in,double y_in,double z_in,double a_in,double b_in,double r_r_in,
    double r_f_in,double t_in,double k_conductivity_in,double density_in,double cp_steel_in,double travel_velocity_in,double B_in,
    double D_in,double L_in,double b_g_in, double d_g_in,double y_i_in,double t_start_point_RAMP_in,double t_end_point_RAMP_in,
    double uprate_RAMP_in,double downrate_RAMP_in){
//,double velocity_ramp_uprate_in,double velocity_ramp_downdate_in,double time_2_in,double time_3_in

  x = x_in;
  y = y_in;
  z = z_in;
  a = a_in;
  b = b_in;
  r_r=r_r_in;
  r_f=r_f_in;
  t = t_in;
  k_conductivity = k_conductivity_in;
  density = density_in;
  cp_steel = cp_steel_in;
  travel_velocity=travel_velocity_in;
  B=B_in;
  D=D_in;
  L=L_in;
  b_g=b_g_in;
  d_g=d_g_in;
  y_i=y_i_in;
  t_start_point_RAMP=t_start_point_RAMP_in;
  t_end_point_RAMP=t_end_point_RAMP_in;
  uprate_RAMP=uprate_RAMP_in;
  downrate_RAMP=downrate_RAMP_in;
 // velocity_ramp_uprate=velocity_ramp_uprate_in;
  //velocity_ramp_downdate=velocity_ramp_downdate_in;
  //time_2=time_2_in;
  //time_3=time_3_in;

}

void FINITE_plate_to_be_integrated_FSW_4Q_DWELL::set_values(double x_in,double y_in,double z_in,double a_in,double b_in,double r_r_in,
    double r_f_in,double t_in,double k_conductivity_in,double density_in,double cp_steel_in,double B_in,
    double D_in,double L_in,double b_g_in, double d_g_in,double l_g_in,double y_i_in,double t_start_point_RAMP_in,
    double t_end_point_RAMP_in,double uprate_RAMP_in,double downrate_RAMP_in){

  x = x_in;
  y = y_in;
  z = z_in;
  a = a_in;
  b = b_in;
  r_r=r_r_in;
  r_f=r_f_in;
  t = t_in;
  k_conductivity = k_conductivity_in;
  density = density_in;
  cp_steel = cp_steel_in;
  B=B_in;
  D=D_in;
  L=L_in;
  b_g=b_g_in;
  d_g=d_g_in;
  l_g=l_g_in;
  y_i=y_i_in;
  t_start_point_RAMP=t_start_point_RAMP_in;
  t_end_point_RAMP=t_end_point_RAMP_in;
  uprate_RAMP=uprate_RAMP_in;
  downrate_RAMP=downrate_RAMP_in;
}






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FINITE_plate_to_be_integrated_with_ramping_CONVECTION_ON_ALL_SURFACES::set_values (double x_in,double y_in,double z_in,double a_in,double b_in,double cr_in,double cf_in,double t_in,double k_conductivity_in,double density_in,double cp_steel_in,double v_arc_in,double B_in, double D_in, double L_in,double ts_in,double tf_in,double r_start_in,double r_fall_in,double b_g_in,double d_g_in) {
  x = x_in;
  y = y_in;
  z = z_in;
  a = a_in;
  b = b_in;
  cr = cr_in;
  cf = cf_in;
  t = t_in;
  k_conductivity = k_conductivity_in;
  density = density_in;
  cp_steel = cp_steel_in;
  v_arc = v_arc_in;
  B=B_in;
	  D=D_in;
	  L=L_in;
	ts=ts_in;
	tf=tf_in;
	r_start=r_start_in;
	r_fall=r_fall_in;
	b_g=b_g_in;
	d_g=d_g_in;
}

//////////////////////////////////////////////////////FUNCTION DEFENITIONS////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Welcome(){
std::cout<<"******************************************"<<std::endl;
std::cout<<"******************************************"<<std::endl;
std::cout<<"Welcome to HEDSATS - High Energy Density Semi-Analytical Thermal Solution"<<std::endl;
std::cout<<"******************************************"<<std::endl;
std::cout<<"******************************************"<<std::endl;
std::cout<<"This software is used to solve the heat equation subject to the application of heat source models representative of Arc, Friction Stir Welding and Electron Beam Welding."<<std::endl;
std::cout<<"The spatial computation is an exact solution with the only numerical operation being the temporal integral of the heat kernels."<<std::endl;
std::cout<<"A veriety of boundary conditions are provided. The library of solutions is growing"<<std::endl;
std::cout<<"******************************************"<<std::endl;
std::cout<<"******************************************"<<std::endl;
}

double erf(double x){
	double magxerf=sqrt(ThermalA::intpow(x,2));
 double y = 1.0 / ( 1.0 + 0.3275911 * magxerf);
 double erf= 1 - (((((
        + 1.061405429  * y
        - 1.453152027) * y
        + 1.421413741) * y
        - 0.284496736) * y
        + 0.254829592) * y)
        * exp (-magxerf * magxerf);
 if(0<x){
return erf;
 }
 else{
return -erf;
 }
}


double inverf(double x){

	double MAXDUB=4;//5.1969857; /* made this up , is there a macro? */ incidentally controlls peak temp for finite plate model
	double sign;

	if( x >= 0.0 ){
sign=1.00;
	}
	else{
sign=-1.00;
	}


  double ax, t, ret;

  ax = fabs(x);

/* This approximation, taken from Table 10 of Blair et al., is valid
   for |x|<=0.75 and has a maximum relative error of 4.47 x 10^-8. */

  if( ax <= 0.75 ) {

    double p[] = {-13.0959967422,26.785225760,-9.289057635};
    double q[] = {-12.0749426297,30.960614529,-17.149977991,1.00000000};

    t = x*x-0.75*0.75;
    ret = x*(p[0]+t*(p[1]+t*p[2]))/(q[0]+t*(q[1]+t*(q[2]+t*q[3])));

  } else if( ax >= 0.75 && ax <= 0.9375 ) {
    double p[] = {-.12402565221,1.0688059574,-1.9594556078,.4230581357};
    double q[] = {-.08827697997,.8900743359,-2.1757031196,1.0000000000};

/* This approximation, taken from Table 29 of Blair et al., is valid
   for .75<=|x|<=.9375 and has a maximum relative error of 4.17 x 10^-8. */

    t = x*x - 0.9375*0.9375;
    ret = x*(p[0]+t*(p[1]+t*(p[2]+t*p[3])))/
         (q[0]+t*(q[1]+t*(q[2]+t*q[3])));

  } else if( ax >= 0.9375 && ax <= (1.0-1.0e-100) ){
    double p[] = {.1550470003116,1.382719649631,.690969348887,
       -1.128081391617, .680544246825,-.16444156791};
    double q[]={.155024849822,1.385228141995,1.000000000000};

/* This approximation, taken from Table 50 of Blair et al., is valid
   for .9375<=|x|<=1-10^-100 and has a maximum relative error of 2.45 x 10^-8. */

    t=1.0/sqrt(-log(1.0-ax));
    ret = sign*(p[0]/t+p[1]+t*(p[2]+t*(p[3]+t*(p[4]+t*p[5]))))/
          (q[0]+t*(q[1]+t*(q[2])));
    } else  {
      ret = sign*MAXDUB;
    }

  if(fabs(ret)>MAXDUB){
ret=sign*MAXDUB;
  }

   return(ret);
}

double intpow( double base, int exponent )
{
int i;
double out = base;
for( i=1 ; i < exponent ; i++ )
{
out *= base;
}
if(exponent==0){
out=1;
}
return out;
}


double power_ramping_function(double t_prime,double ts, double tf, double r_start, double r_fall){

double ramping=((1/(1+(exp(-r_start*(t_prime-(ts+(log(99.0)/r_start)))))))*((1/(1+(exp(r_fall*(t_prime-(tf-(log(99.0)/r_fall)))))))));
return ramping;

}



//double part_past_integral_infinite(double x,double y,double z,double a,double b,double cr,double cf,double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc){
//
//	double kappa=(k_conductivity/(cp_steel*density));
//
//	double rf=((2*cf)/(cr+cf));
//	double rr=((2*cr)/(cr+cf));
//
//	double xx=ThermalA::intpow(x,2);
//	double yy=ThermalA::intpow(y,2);
//	double aa=ThermalA::intpow(a,2);
//	double bb=ThermalA::intpow(b,2);
//	double crcr=ThermalA::intpow(cr,2);
//	double cfcf=ThermalA::intpow(cf,2);
//	double t_minus_t_prime=(sqrt(ThermalA::intpow((t-t_prime),2)));
//double kappa12ttp=(12*kappa*(t_minus_t_prime));
//	//double kappa12ttp=(12*kappa*(t-t_prime));
//
//
//	double A_r=(1/sqrt((kappa12ttp+crcr)))*(exp(-3*((ThermalA::intpow((z-(v_arc*t_prime)),2))/(kappa12ttp+crcr))));
//	double A_f=(1/sqrt((kappa12ttp+cfcf)))*(exp(-3*((ThermalA::intpow((z-(v_arc*t_prime)),2))/(kappa12ttp+cfcf))));
//
//	double B_r= ThermalA::erf((cr/2)*((z-(v_arc*t_prime))/(sqrt(kappa12ttp/12)*sqrt(kappa12ttp+crcr))));
//	double B_f= ThermalA::erf((cf/2)*((z-(v_arc*t_prime))/(sqrt(kappa12ttp/12)*sqrt(kappa12ttp+cfcf))));
////cout<<B_f<<endl;
//double part_past_integral_infinite=((1/(sqrt(kappa12ttp+aa)*sqrt(kappa12ttp+bb)))*(exp(-3*((xx/(kappa12ttp+aa))+(yy/(kappa12ttp+bb)))))*((rr*A_r*(1-B_r))+(rf*A_f*(1+B_f))));
//return part_past_integral_infinite;
//}



//void matlab_write(int county, int countx, double B, double D, double L, double v_arc, double z, double tstart, double tgap, double tend, double xstart, double xend, double yend, double max_scale, string loc1, string loc2, string loc3, string folder){
//
//
//
//
//	ofstream matlab_file;
//string file_name_matlab=loc1+loc2+loc3;
//	matlab_file.open(file_name_matlab.c_str());
//
//cout<<"WRITING COMPILATION FILE"<<"\t"<<file_name_matlab<<endl;
//
//matlab_file<<"clc; \n clear all;"<<endl;
//
//matlab_file<<"range= "<<"[0 0 "<<county-1+0<<" "<<countx-1+0<<"];"<<endl;
//matlab_file<<"aviobjinfinite = VideoWriter(strcat('"<<loc1<<"Videos\\"<<"_B="<<B<<"D="<<D<<"L="<<L<<"arc_velocity="<<v_arc<<", "<<z<<"mm along welding axis','.avi'));"<<endl;
//matlab_file<<"open(aviobjinfinite);"<<endl;
//matlab_file<<"for genericcounter="<<tstart<<":"<<tgap<<":"<<tend<<endl;
//matlab_file<<"addressstartinfinite='"<<loc1+folder<<"';"<<endl;
//
//
//matlab_file<<"addressendinfinite='s';"<<endl;
//matlab_file<<"addressstringinfinite=strcat(addressstartinfinite,num2str(genericcounter),addressendinfinite,'"<<"','.txt');"<<endl;
//matlab_file<<"cinf=dlmread(addressstringinfinite,'\\t',range);"<<endl;
//matlab_file<<"infinlabel=strcat(num2str(genericcounter),' s');"<<endl;
//matlab_file<<"INFINITEfig=figure('Visible','off');"<<endl;
//matlab_file<<"annotation('textbox',[0 0 1 1],'string',infinlabel)"<<endl;
//matlab_file<<"imagesc("<<xstart<<":"<<xend<<",0:"<<fabs(yend)<<",cinf,[0 "<<max_scale<<"])"<<endl;
//matlab_file<<"colorbar"<<endl;
//matlab_file<<"writelocationINFINITEstart='E:\\Green\\Images\\"<<"';"<<endl;
//matlab_file<<"finalwritelocationINFINITE=strcat(writelocationINFINITEstart,num2str(genericcounter),'s','.png');"<<endl;
//matlab_file<<"saveas(INFINITEfig,finalwritelocationINFINITE)"<<endl;
//matlab_file<<"infinframe=getframe(INFINITEfig);"<<endl;
//matlab_file<<"close(INFINITEfig)"<<endl;
//matlab_file<<"writeVideo(aviobjinfinite,infinframe);"<<endl;
//matlab_file<<"end"<<endl;
////matlab_file<<"end"<<endl;
//matlab_file<<"close(aviobjinfinite);"<<endl;
//
//matlab_file.close();
//}
//void data_write(double **cross_section, int array_sizeoverx, int array_sizeovery, string folder_path,string file_identifier){
////cout<<"WRITING FILE - "<<loc3<<" - "<<Temp_or_flux<<endl;
//ofstream datafile;
//string file_name=folder_path+file_identifier;
////cout<<file_name<<endl;
//datafile.open(file_name.c_str());
//
//	for(int cy=0;cy<=array_sizeovery-1;cy++){
//	for(int cx=0;cx<=array_sizeoverx-1;cx++){
////            if(cross_section[cx][cy]==0){
////
////              datafile<<"\t"<<"\t";
////
////            }
////            else{
//datafile<<cross_section[cx][cy]<<"\t";
//  //          }
//	}
//datafile<<endl;
//	}
//
//	datafile.close();
//
//}

//void data_write_3D(double ***Volume_Array, int array_sizeoverx, int array_sizeovery,int array_sizeoverz, string folder_path,string file_identifier){
//
//ofstream datafile;
//string file_name=folder_path+file_identifier;
//datafile.open(file_name.c_str());
//
//    for(int cx=0;cx<=array_sizeoverx-1;cx++){
//    for(int cy=0;cy<=array_sizeovery-1;cy++){
//    for(int cz=0;cz<=array_sizeoverz-1;cz++){
//
//
//
//datafile<<cx<<"\t"<<cy<<"\t"<<cz<<"\t"<<Volume_Array[cx][cy][cz]<<endl;
//
//	}
////	datafile<<endl;
//	}
//	datafile<<endl;
//    }
//
//
//
//
//
//
//
//
//
//	datafile.close();
//}

//double y_integrated_between_dg_and_D(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g){//2/D inside
//
//double Beta_m[48]={//for B=0.05
//  0.22176,3.15743,6.29113,9.43008,
//  12.5703,15.7111,18.8522,21.9934,
//  25.1347,28.2761,31.4175,34.559,
//  37.7004,40.8419,43.9834,47.125,
//  50.2665,53.408,56.5496,59.6911,
//  62.8326,65.9742,69.1158,72.2573,
//  75.3989,78.5405,81.682,84.8236,
//  87.9652,91.1067,94.2483,97.3899,
//  100.531,103.673,106.815,109.956,
//  113.098,116.239,119.381,122.523,
//  125.664,128.806,131.947,135.089,
//  138.23,141.372,144.514,147.655};
//
//double B_Nuss=0.05;
//	double alpha=(k_conductivity/(cp_steel*density));
//	//cout<<"alpha"<<"\t"<<alpha<<endl;
//double Y33=0;
//for(int m=0;m<48;m++){//48
//
//double argument_of_trig_fns1=Beta_m[m]*(-1+(d_g/D));
//double argument_of_trig_fns2=Beta_m[m]-((Beta_m[m]*d_g)/(D));
//complex<double> A_y_arg;
//complex<double> B_y_arg;
//complex<double> C_y_arg;
//complex<double> D_y_arg;
//complex<double> E_y_arg;
//complex<double> F_y_arg;
//
//double first_part=((ThermalA::intpow(b,2)*Beta_m[m])/(D))/(2*sqrt(3.00)*b);//for
//double second_part=(6*(D-d_g))/(2*sqrt(3.00)*b);//for
//double third_part=(6*(d_g-D))/(2*sqrt(3.00)*b);//for
//double fourth_part=(b*Beta_m[m])/(2*sqrt(3.00)*D);//for ThermalA::erf
//
//
//  A_y_arg = complex<double>(first_part,second_part);
//  B_y_arg = complex<double>(first_part,third_part);
//  C_y_arg = complex<double>(fourth_part,0);
//  D_y_arg = complex<double>(first_part,second_part);
//  E_y_arg = complex<double>(fourth_part,0);
//  F_y_arg = complex<double>(first_part,third_part);
//
//double A_y=cos(argument_of_trig_fns1)*imag(Faddeeva::erfi(A_y_arg));
//double B_y=cos(argument_of_trig_fns2)*imag(Faddeeva::erfi(B_y_arg));
//double C_y=sin(argument_of_trig_fns1)*real(Faddeeva::erfi(C_y_arg));
//double D_y=sin(argument_of_trig_fns1)*real(Faddeeva::erfi(D_y_arg));
//double E_y=sin(argument_of_trig_fns2)*real(Faddeeva::erfi(E_y_arg));
//double F_y=sin(argument_of_trig_fns2)*real(Faddeeva::erfi(F_y_arg));
//
//
//double Omegla=(A_y)-(B_y)-(C_y)+(D_y)+(E_y)-(F_y);
//double denominator=(4*D*(B_Nuss+ThermalA::intpow(B_Nuss,2)+ThermalA::intpow(Beta_m[m],2)));
////double Omegla=(2*cos((Beta_m[m]*(ThermalA::intpow(D,2)-(D*d_g)))/(ThermalA::intpow(D,2)))*ThermalA::erf((6.00*(ThermalA::intpow(D,2)-(D*d_g)))/(2*sqrt(3.00)*b*D)))+(2*Faddeeva::erfi((b*Beta_m[m])/(2*sqrt(3.00)*D))*sin((Beta_m[m]*(ThermalA::intpow(D,2)-(D*d_g)))/(ThermalA::intpow(D,2))));
////double numerator=2*b*(exp(-1*((pow(Beta_m[m],2.00)*(pow(b,2.00)+(12*alpha*(t-t_prime))))/(12*pow(D,2.00)))))*(sqrt(pi/3.00))*(ThermalA::intpow(B_Nuss,2)+ThermalA::intpow(Beta_m[m],2))*(cos(Beta_m[m]-((y*Beta_m[m])/D)));
//double numerator=2*b*exp(-((ThermalA::intpow(Beta_m[m],2)*(ThermalA::intpow(b,2)+(12*alpha*(t-t_prime))))/(12*ThermalA::intpow(D,2))))*sqrt(pi/3.00)*(ThermalA::intpow(B_Nuss,2)+ThermalA::intpow(Beta_m[m],2))*cos(Beta_m[m]-((y*Beta_m[m])/(D)))*Omegla;
//
//
//Y33=Y33+((numerator)/denominator);
////cout<<Y33<<endl;
//}
//return Y33;
//}
//double y_integrated_between_0_and_D(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g){
/////////////////DO THIS CHANGE!!!!
//double Beta_m[48]={//for B=0.05
//  0.22176,3.15743,6.29113,9.43008,
//  12.5703,15.7111,18.8522,21.9934,
//  25.1347,28.2761,31.4175,34.559,
//  37.7004,40.8419,43.9834,47.125,
//  50.2665,53.408,56.5496,59.6911,
//  62.8326,65.9742,69.1158,72.2573,
//  75.3989,78.5405,81.682,84.8236,
//  87.9652,91.1067,94.2483,97.3899,
//  100.531,103.673,106.815,109.956,
//  113.098,116.239,119.381,122.523,
//  125.664,128.806,131.947,135.089,
//  138.23,141.372,144.514,147.655};
//double B_Nuss=0.05;
//	double alpha=(k_conductivity/(cp_steel*density));
//	//cout<<"alpha"<<"\t"<<alpha<<endl;
//double Y33=0;
//for(int m=0;m<48;m++){//48
//
//double argument_of_trig_fns=((Beta_m[m]*(ThermalA::intpow(D,2)-(D*d_g)))/(ThermalA::intpow(D,2)));
//complex<double> A_y_arg;
//complex<double> B_y_arg;
//complex<double> C_y_arg;
//complex<double> D_y_arg;
//complex<double> E_y_arg;
//complex<double> F_y_arg;
//complex<double> G_y_arg;
//complex<double> H_y_arg;
//
//double first_part=((ThermalA::intpow(b,2)*Beta_m[m])/(2*sqrt(3.00)*b*D));//for erfi
//double second_part=((6*D*d_g)/(2*sqrt(3.00)*b*D));//for erfi
//double third_part=((6*(ThermalA::intpow(D,2)-(D*d_g)))/(2*sqrt(3.00)*b*D));//for ThermalA::erf
////double fourth_part=(()/(2*sqrt(3.00)*b*D));//for ThermalA::erf
//
//
//  A_y_arg = complex<double>(first_part,-second_part);
//  B_y_arg = complex<double>(first_part,second_part);
//  C_y_arg = complex<double>(third_part,-first_part);
//  D_y_arg = complex<double>(third_part,first_part);
//  E_y_arg = complex<double>(third_part,-first_part);
//  F_y_arg = complex<double>(first_part,-second_part);
//  G_y_arg = complex<double>(third_part,first_part);
//
//  H_y_arg = complex<double>(first_part,second_part);
//
//
//double A_y=(cos(argument_of_trig_fns)*imag(Faddeeva::erfi(A_y_arg)));
//double B_y=(cos(argument_of_trig_fns)*imag(Faddeeva::erfi(B_y_arg)));
//double C_y=(cos(argument_of_trig_fns)*real(Faddeeva::erf(C_y_arg)));
//double D_y=(cos(argument_of_trig_fns)*real(Faddeeva::erf(D_y_arg)));
//double E_y=(sin(argument_of_trig_fns)*imag(Faddeeva::erf(E_y_arg)));
//double F_y=(sin(argument_of_trig_fns)*real(Faddeeva::erfi(F_y_arg)));
//double G_y=(sin(argument_of_trig_fns)*imag(Faddeeva::erf(G_y_arg)));
//double H_y=(sin(argument_of_trig_fns)*real(Faddeeva::erfi(H_y_arg)));
//
//double Omegla=(-1*(A_y))+(B_y)+(C_y)+(D_y)+(E_y)+(F_y)-(G_y)+(H_y);
//
//double denominator=(D*4*(B_Nuss+ThermalA::intpow(B_Nuss,2)+ThermalA::intpow(Beta_m[m],2)));
//double numerator=2*b*exp(-((ThermalA::intpow(Beta_m[m],2)*((ThermalA::intpow(b,2))+(12*alpha*(t-t_prime))))/(12*ThermalA::intpow(D,2))))*sqrt(pi/3.00)*(ThermalA::intpow(B_Nuss,2)+ThermalA::intpow(Beta_m[m],2))*cos(Beta_m[m]-((y*Beta_m[m])/(D)))*Omegla;
//Y33=Y33+((numerator)/denominator);
////cout<<"to test"<<"\t"<<Y33<<endl;
//
//}
//return Y33;
//}

//double y_insulating_integrated_between_d_g_and_D(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g){
//double alpha=(k_conductivity/(cp_steel*density));
////cout<<"alpha"<<alpha<<endl;
//
//////////////////    1'  GF between d_g and D, insulating at d_g and D, heat source at d_g
//
//
//double Y22=0;
//for(int n=-30;n<=30;n++){
//
//double erf_denom=2*b*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(b,2)+(12*alpha*(t-t_prime)));
//double numerator=b*exp(-((3*ThermalA::intpow(((2*D*n)+y-(d_g*(1+(2*n)))),2))/(ThermalA::intpow(b,2)+(12*alpha*(t-t_prime)))))*sqrt(alpha)*(ThermalA::erf(((ThermalA::intpow(b,2)*(D+(2*D*n)+y))+(12*D*t*alpha)-(2*d_g*((ThermalA::intpow(b,2)*(1+n))+(6*t*alpha)))+(12*alpha*t_prime*(-D+d_g)))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(b,2)*(D-(2*D*n)-y))+(12*D*t*alpha)-(12*D*t_prime*alpha)+(2*d_g*((ThermalA::intpow(b,2)*n)-(6*t*alpha)+(6*t_prime*alpha))))/(erf_denom)))*sqrt(t-t_prime);
//double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(b,2)+(12*alpha*(t-t_prime)));
//
//
//Y22+=numerator/denominator;
////Y22+=PART_2;
////cout<<X22<<endl;
//}
//return Y22;
//}


////////////////CONVECTION AT BOTH SURFACES////////////////////

double y_component_convection_at_D_convection_at_0_hs_at_dg(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g){

//  SHOULD THE 2/L BE INSIDE OR OUTSIDE THE SUM     THINK ABOUT IT

//double Beta_m[64]={//for B=0.05
// 0,0.38412,3.18862,6.30697,9.44067,12.5783,15.7175,18.8575,21.998,25.1387,
//   28.2796,31.4207,34.5619,37.7031,40.8444,43.9857,47.1271,50.2685,53.4099,
//   56.5513,59.6928,62.8342,65.9757,69.1172,72.2587,75.4002,78.5417,81.6832,
//   84.8248,87.9663,91.1078,94.2494,97.3909,100.532,103.674,106.816,109.957,
//   113.099,116.24,119.382,122.523,125.665,128.806,131.948,135.09,138.231,
//   141.373,144.514,147.656,150.797,153.939,157.081,160.222,163.364,166.505,
//   169.647,172.788,175.93,179.072,182.213,185.355,188.496,191.638,194.78};
//
//   //,197.921
//
//double B_Nuss_0=0.1;
//double B_Nuss_D=0.05;
double Beta_m[637]={//for B=0.1

 0.0,0.443521,3.20399,6.31485,9.44595,12.5823,15.7207,18.8602,22.0002,25.1407,28.2814,31.4223,
 34.5633,37.7044,40.8456,43.9868,47.1281,50.2695,53.4108,56.5522,59.6936,62.835,65.9765,69.1179,
 72.2594,75.4009,78.5424,81.6839,84.8254,87.9669,91.1084,94.2499,97.3914,100.533,103.674,106.816,
 109.958,113.099,116.241,119.382,122.524,125.665,128.807,131.948,135.09,138.232,141.373,144.515,
 147.656,150.798,153.939,157.081,160.222,163.364,166.506,169.647,172.789,175.93,179.072,182.213,
 185.355,188.497,191.638,194.78,197.921,201.063,204.205,207.346,210.488,213.629,216.771,219.912,
 223.054,226.196,229.337,232.479,235.62,238.762,241.903,245.045,248.187,251.328,254.47,257.611,
 260.753,263.895,267.036,270.178,273.319,276.461,279.602,282.744,285.886,289.027,292.169,295.31,
 298.452,301.594,304.735,307.877,311.018,314.16,317.301,320.443,323.585,326.726,329.868,333.009,
 336.151,339.293,342.434,345.576,348.717,351.859,355.001,358.142,361.284,364.425,367.567,370.708,
 373.85,376.992,380.133,383.275,386.416,389.558,392.7,395.841,398.983,402.124,405.266,408.408,
 411.549,414.691,417.832,420.974,424.115,427.257,430.399,433.54,436.682,439.823,442.965,446.107,
 449.248,452.39,455.531,458.673,461.815,464.956,468.098,471.239,474.381,477.523,480.664,483.806,
 486.947,490.089,493.23,496.372,499.514,502.655,505.797,508.938,512.08,515.222,518.363,521.505,
 524.646,527.788,530.93,534.071,537.213,540.354,543.496,546.637,549.779,552.921,556.062,559.204,
 562.345,565.487,568.629,571.77,574.912,578.053,581.195,584.337,587.478,590.62,593.761,596.903,
 600.045,603.186,606.328,609.469,612.611,615.752,618.894,622.036,625.177,628.319,631.46,634.602,
 637.744,640.885,644.027,647.168,650.31,653.452,656.593,659.735,662.876,666.018,669.16,672.301,
 675.443,678.584,681.726,684.867,688.009,691.151,694.292,697.434,700.575,703.717,706.859,710.0,
 713.142,716.283,719.425,722.567,725.708,728.85,731.991,735.133,738.275,741.416,744.558,747.699,
 750.841,753.983,757.124,760.266,763.407,766.549,769.69,772.832,775.974,779.115,782.257,785.398,
 788.54,791.682,794.823,797.965,801.106,804.248,807.39,810.531,813.673,816.814,819.956,823.098,
 826.239,829.381,832.522,835.664,838.805,841.947,845.089,848.23,851.372,854.513,857.655,860.797,
 863.938,867.08,870.221,873.363,876.505,879.646,882.788,885.929,889.071,892.213,895.354,898.496,
 901.637,904.779,907.92,911.062,914.204,917.345,920.487,923.628,926.77,929.912,933.053,936.195,
 939.336,942.478,945.62,948.761,951.903,955.044,958.186,961.328,964.469,967.611,970.752,973.894,
 977.036,980.177,983.319,986.46,989.602,992.743,995.885,999.027,1002.17,1005.31,1008.45,1011.59,
 1014.73,1017.88,1021.02,1024.16,1027.3,1030.44,1033.58,1036.73,1039.87,1043.01,1046.15,1049.29,
 1052.43,1055.58,1058.72,1061.86,1065.,1068.14,1071.28,1074.42,1077.57,1080.71,1083.85,1086.99,
 1090.13,1093.27,1096.42,1099.56,1102.7,1105.84,1108.98,1112.12,1115.27,1118.41,1121.55,1124.69,
 1127.83,1130.97,1134.12,1137.26,1140.4,1143.54,1146.68,1149.82,1152.96,1156.11,1159.25,1162.39,
 1165.53,1168.67,1171.81,1174.96,1178.1,1181.24,1184.38,1187.52,1190.66,1193.81,1196.95,1200.09,
 1203.23,1206.37,1209.51,1212.65,1215.8,1218.94,1222.08,1225.22,1228.36,1231.5,1234.65,1237.79,
 1240.93,1244.07,1247.21,1250.35,1253.5,1256.64,1259.78,1262.92,1266.06,1269.2,1272.35,1275.49,
 1278.63,1281.77,1284.91,1288.05,1291.19,1294.34,1297.48,1300.62,1303.76,1306.9,1310.04,1313.19,
 1316.33,1319.47,1322.61,1325.75,1328.89,1332.04,1335.18,1338.32,1341.46,1344.6,1347.74,1350.88,
 1354.03,1357.17,1360.31,1363.45,1366.59,1369.73,1372.88,1376.02,1379.16,1382.3,1385.44,1388.58,
 1391.73,1394.87,1398.01,1401.15,1404.29,1407.43,1410.58,1413.72,1416.86,1420.,1423.14,1426.28,
 1429.42,1432.57,1435.71,1438.85,1441.99,1445.13,1448.27,1451.42,1454.56,1457.7,1460.84,1463.98,
 1467.12,1470.27,1473.41,1476.55,1479.69,1482.83,1485.97,1489.12,1492.26,1495.4,1498.54,1501.68,
 1504.82,1507.96,1511.11,1514.25,1517.39,1520.53,1523.67,1526.81,1529.96,1533.1,1536.24,1539.38,
 1542.52,1545.66,1548.81,1551.95,1555.09,1558.23,1561.37,1564.51,1567.65,1570.8,1577.08,1580.22,
 1583.36,1586.5,1589.65,1592.79,1595.93,1599.07,1602.21,1605.35,1608.5,1611.64,1614.78,1617.92,
 1621.06,1624.2,1627.35,1630.49,1633.63,1636.77,1639.91,1643.05,1646.19,1649.34,1652.48,1655.62,
 1658.76,1661.9,1665.04,1668.19,1671.33,1674.47,1677.61,1680.75,1683.89,1687.04,1690.18,1693.32,
 1696.46,1699.6,1702.74,1705.88,1709.03,1712.17,1715.31,1718.45,1721.59,1724.73,1727.88,1731.02,
 1734.16,1737.3,1740.44,1743.58,1746.73,1749.87,1753.01,1756.15,1759.29,1762.43,1765.58,1768.72,
 1771.86,1775.,1778.14,1781.28,1784.42,1787.57,1790.71,1793.85,1796.99,1800.13,1803.27,1806.42,
 1809.56,1812.7,1815.84,1818.98,1822.12,1825.27,1828.41,1831.55,1834.69,1837.83,1840.97,1844.11,
 1847.26,1850.4,1853.54,1856.68,1859.82,1862.96,1866.11,1869.25,1872.39,1875.53,1878.67,1881.81,
 1884.96,1888.1,1891.24,1894.38,1897.52,1900.66,1903.81,1906.95,1910.09,1913.23,1916.37,1919.51,
 1922.65,1925.8,1928.94,1932.08,1935.22,1938.36,1941.5,1944.65,1947.79,1950.93,1954.07,1957.21,
 1960.35,1963.5,1966.64,1969.78,1972.92,1976.06,1979.2,1982.35,1985.49,1988.63,1991.77,1994.91,1998.05

 };

double B_Nuss_0=0.1;
double B_Nuss_D=0.1;

	double alpha=(k_conductivity/(cp_steel*density));
	//cout<<"alpha"<<"\t"<<alpha<<endl;
	//cout<<"y"<<endl;
double Y33=0;
for(int m=0;m<228;m++){//64

double argument_of_trig_fns=Beta_m[m]*(d_g/D);
//double argument_of_trig_fns2=Beta_m[m]-((Beta_m[m]*d_g)/(D));

complex<double> A_y_arg;
complex<double> B_y_arg;
complex<double> C_y_arg;
complex<double> D_y_arg;


double first_part=(d_g*sqrt(3.00))/(b);//for
double second_part=(b*Beta_m[m])/(2*sqrt(3.00)*D);//for
double third_part=(6*(ThermalA::intpow(D,2)-(D*d_g)))/(2*sqrt(3.00)*b*D);//for

  A_y_arg = complex<double>(first_part,-second_part);
  B_y_arg = complex<double>(first_part,second_part);
  C_y_arg = complex<double>(third_part,second_part);
  D_y_arg = complex<double>(third_part,-second_part);


double A_y=cos(argument_of_trig_fns)*imag(Faddeeva::erf(A_y_arg))*B_Nuss_0;
double B_y=cos(argument_of_trig_fns)*imag(Faddeeva::erf(B_y_arg))*B_Nuss_0;
double C_y=sin(argument_of_trig_fns)*real(Faddeeva::erf(C_y_arg))*B_Nuss_0;
double D_y=sin(argument_of_trig_fns)*real(Faddeeva::erf(A_y_arg))*B_Nuss_0;
double E_y=sin(argument_of_trig_fns)*real(Faddeeva::erf(D_y_arg))*B_Nuss_0;
double F_y=sin(argument_of_trig_fns)*real(Faddeeva::erf(B_y_arg))*B_Nuss_0;
double G_y=cos(argument_of_trig_fns)*real(Faddeeva::erf(D_y_arg))*Beta_m[m];
double H_y=cos(argument_of_trig_fns)*real(Faddeeva::erf(C_y_arg))*Beta_m[m];
double I_y=cos(argument_of_trig_fns)*real(Faddeeva::erf(A_y_arg))*Beta_m[m];
double J_y=cos(argument_of_trig_fns)*real(Faddeeva::erf(B_y_arg))*Beta_m[m];
double K_y=sin(argument_of_trig_fns)*imag(Faddeeva::erf(A_y_arg))*Beta_m[m];
double L_y=sin(argument_of_trig_fns)*imag(Faddeeva::erf(B_y_arg))*Beta_m[m];


double Omegla=(-A_y)+(B_y)+(C_y)+(D_y)+(E_y)+(F_y)+(G_y)+(H_y)+(I_y)+(J_y)+(K_y)-(L_y)+(imag(Faddeeva::erf(C_y_arg))*((-1*cos(argument_of_trig_fns)*B_Nuss_0)+(sin(argument_of_trig_fns)*Beta_m[m])))+(imag(Faddeeva::erf(D_y_arg))*((cos(argument_of_trig_fns)*B_Nuss_0)-(sin(argument_of_trig_fns)*Beta_m[m])));
double denominator=(2*D*((B_Nuss_0*(ThermalA::intpow(B_Nuss_D,2)+ThermalA::intpow(Beta_m[m],2)))+(ThermalA::intpow(B_Nuss_0,2)*(B_Nuss_D+ThermalA::intpow(B_Nuss_D,2)+ThermalA::intpow(Beta_m[m],2)))+(ThermalA::intpow(Beta_m[m],2)*(B_Nuss_D+ThermalA::intpow(B_Nuss_D,2)+ThermalA::intpow(Beta_m[m],2)))));
//double Omegla=(2*cos((Beta_m[m]*(ThermalA::intpow(D,2)-(D*d_g)))/(ThermalA::intpow(D,2)))*ThermalA::erf((6.00*(ThermalA::intpow(D,2)-(D*d_g)))/(2*sqrt(3.00)*b*D)))+(2*Faddeeva::erfi((b*Beta_m[m])/(2*sqrt(3.00)*D))*sin((Beta_m[m]*(ThermalA::intpow(D,2)-(D*d_g)))/(ThermalA::intpow(D,2))));
//double numerator=2*b*(exp(-1*((pow(Beta_m[m],2.00)*(pow(b,2.00)+(12*alpha*(t-t_prime))))/(12*pow(D,2.00)))))*(sqrt(pi/3.00))*(ThermalA::intpow(B_Nuss,2)+ThermalA::intpow(Beta_m[m],2))*(cos(Beta_m[m]-((y*Beta_m[m])/D)));
double numerator=b*exp(-((ThermalA::intpow(Beta_m[m],2)*(ThermalA::intpow(b,2)+(12*alpha*(t-t_prime))))/(12*ThermalA::intpow(D,2))))*sqrt(pi/3.00)*((sin((y*Beta_m[m])/(D))*B_Nuss_0)+(cos((y*Beta_m[m])/(D))*Beta_m[m]))*(ThermalA::intpow(B_Nuss_D,2)+ThermalA::intpow(Beta_m[m],2))*Omegla;
double inc=numerator/denominator;

Y33+=((numerator)/denominator);
//cout<<Y33<<endl;
}
//if(isnan(Y33)==1){
//Y33=0.0;
//}


return Y33;
}


double x_component_convection_at_B_convection_at_0_hs_at_bg(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g){

//  SHOULD THE 2/L BE INSIDE OR OUTSIDE THE SUM     THINK ABOUT IT

double Beta_m[637]={//for B=0.05

 0.0,0.443521,3.20399,6.31485,9.44595,12.5823,15.7207,18.8602,22.0002,25.1407,28.2814,31.4223,
 34.5633,37.7044,40.8456,43.9868,47.1281,50.2695,53.4108,56.5522,59.6936,62.835,65.9765,69.1179,
 72.2594,75.4009,78.5424,81.6839,84.8254,87.9669,91.1084,94.2499,97.3914,100.533,103.674,106.816,
 109.958,113.099,116.241,119.382,122.524,125.665,128.807,131.948,135.09,138.232,141.373,144.515,
 147.656,150.798,153.939,157.081,160.222,163.364,166.506,169.647,172.789,175.93,179.072,182.213,
 185.355,188.497,191.638,194.78,197.921,201.063,204.205,207.346,210.488,213.629,216.771,219.912,
 223.054,226.196,229.337,232.479,235.62,238.762,241.903,245.045,248.187,251.328,254.47,257.611,
 260.753,263.895,267.036,270.178,273.319,276.461,279.602,282.744,285.886,289.027,292.169,295.31,
 298.452,301.594,304.735,307.877,311.018,314.16,317.301,320.443,323.585,326.726,329.868,333.009,
 336.151,339.293,342.434,345.576,348.717,351.859,355.001,358.142,361.284,364.425,367.567,370.708,
 373.85,376.992,380.133,383.275,386.416,389.558,392.7,395.841,398.983,402.124,405.266,408.408,
 411.549,414.691,417.832,420.974,424.115,427.257,430.399,433.54,436.682,439.823,442.965,446.107,
 449.248,452.39,455.531,458.673,461.815,464.956,468.098,471.239,474.381,477.523,480.664,483.806,
 486.947,490.089,493.23,496.372,499.514,502.655,505.797,508.938,512.08,515.222,518.363,521.505,
 524.646,527.788,530.93,534.071,537.213,540.354,543.496,546.637,549.779,552.921,556.062,559.204,
 562.345,565.487,568.629,571.77,574.912,578.053,581.195,584.337,587.478,590.62,593.761,596.903,
 600.045,603.186,606.328,609.469,612.611,615.752,618.894,622.036,625.177,628.319,631.46,634.602,
 637.744,640.885,644.027,647.168,650.31,653.452,656.593,659.735,662.876,666.018,669.16,672.301,
 675.443,678.584,681.726,684.867,688.009,691.151,694.292,697.434,700.575,703.717,706.859,710.0,
 713.142,716.283,719.425,722.567,725.708,728.85,731.991,735.133,738.275,741.416,744.558,747.699,
 750.841,753.983,757.124,760.266,763.407,766.549,769.69,772.832,775.974,779.115,782.257,785.398,
 788.54,791.682,794.823,797.965,801.106,804.248,807.39,810.531,813.673,816.814,819.956,823.098,
 826.239,829.381,832.522,835.664,838.805,841.947,845.089,848.23,851.372,854.513,857.655,860.797,
 863.938,867.08,870.221,873.363,876.505,879.646,882.788,885.929,889.071,892.213,895.354,898.496,
 901.637,904.779,907.92,911.062,914.204,917.345,920.487,923.628,926.77,929.912,933.053,936.195,
 939.336,942.478,945.62,948.761,951.903,955.044,958.186,961.328,964.469,967.611,970.752,973.894,
 977.036,980.177,983.319,986.46,989.602,992.743,995.885,999.027,1002.17,1005.31,1008.45,1011.59,
 1014.73,1017.88,1021.02,1024.16,1027.3,1030.44,1033.58,1036.73,1039.87,1043.01,1046.15,1049.29,
 1052.43,1055.58,1058.72,1061.86,1065.,1068.14,1071.28,1074.42,1077.57,1080.71,1083.85,1086.99,
 1090.13,1093.27,1096.42,1099.56,1102.7,1105.84,1108.98,1112.12,1115.27,1118.41,1121.55,1124.69,
 1127.83,1130.97,1134.12,1137.26,1140.4,1143.54,1146.68,1149.82,1152.96,1156.11,1159.25,1162.39,
 1165.53,1168.67,1171.81,1174.96,1178.1,1181.24,1184.38,1187.52,1190.66,1193.81,1196.95,1200.09,
 1203.23,1206.37,1209.51,1212.65,1215.8,1218.94,1222.08,1225.22,1228.36,1231.5,1234.65,1237.79,
 1240.93,1244.07,1247.21,1250.35,1253.5,1256.64,1259.78,1262.92,1266.06,1269.2,1272.35,1275.49,
 1278.63,1281.77,1284.91,1288.05,1291.19,1294.34,1297.48,1300.62,1303.76,1306.9,1310.04,1313.19,
 1316.33,1319.47,1322.61,1325.75,1328.89,1332.04,1335.18,1338.32,1341.46,1344.6,1347.74,1350.88,
 1354.03,1357.17,1360.31,1363.45,1366.59,1369.73,1372.88,1376.02,1379.16,1382.3,1385.44,1388.58,
 1391.73,1394.87,1398.01,1401.15,1404.29,1407.43,1410.58,1413.72,1416.86,1420.,1423.14,1426.28,
 1429.42,1432.57,1435.71,1438.85,1441.99,1445.13,1448.27,1451.42,1454.56,1457.7,1460.84,1463.98,
 1467.12,1470.27,1473.41,1476.55,1479.69,1482.83,1485.97,1489.12,1492.26,1495.4,1498.54,1501.68,
 1504.82,1507.96,1511.11,1514.25,1517.39,1520.53,1523.67,1526.81,1529.96,1533.1,1536.24,1539.38,
 1542.52,1545.66,1548.81,1551.95,1555.09,1558.23,1561.37,1564.51,1567.65,1570.8,1577.08,1580.22,
 1583.36,1586.5,1589.65,1592.79,1595.93,1599.07,1602.21,1605.35,1608.5,1611.64,1614.78,1617.92,
 1621.06,1624.2,1627.35,1630.49,1633.63,1636.77,1639.91,1643.05,1646.19,1649.34,1652.48,1655.62,
 1658.76,1661.9,1665.04,1668.19,1671.33,1674.47,1677.61,1680.75,1683.89,1687.04,1690.18,1693.32,
 1696.46,1699.6,1702.74,1705.88,1709.03,1712.17,1715.31,1718.45,1721.59,1724.73,1727.88,1731.02,
 1734.16,1737.3,1740.44,1743.58,1746.73,1749.87,1753.01,1756.15,1759.29,1762.43,1765.58,1768.72,
 1771.86,1775.,1778.14,1781.28,1784.42,1787.57,1790.71,1793.85,1796.99,1800.13,1803.27,1806.42,
 1809.56,1812.7,1815.84,1818.98,1822.12,1825.27,1828.41,1831.55,1834.69,1837.83,1840.97,1844.11,
 1847.26,1850.4,1853.54,1856.68,1859.82,1862.96,1866.11,1869.25,1872.39,1875.53,1878.67,1881.81,
 1884.96,1888.1,1891.24,1894.38,1897.52,1900.66,1903.81,1906.95,1910.09,1913.23,1916.37,1919.51,
 1922.65,1925.8,1928.94,1932.08,1935.22,1938.36,1941.5,1944.65,1947.79,1950.93,1954.07,1957.21,
 1960.35,1963.5,1966.64,1969.78,1972.92,1976.06,1979.2,1982.35,1985.49,1988.63,1991.77,1994.91,1998.05

 };

double B_Nuss_0=0.1;
double B_Nuss_B=0.1;
	double alpha=(k_conductivity/(cp_steel*density));
	//cout<<"alpha"<<"\t"<<alpha<<endl;
	//cout<<"x"<<endl;
double X33=0;
for(int m=0;m<228;m++){//64

double argument_of_trig_fns=Beta_m[m]*(b_g/B);
//double argument_of_trig_fns2=Beta_m[m]-((Beta_m[m]*d_g)/(D));

complex<double> A_y_arg;
complex<double> B_y_arg;
complex<double> C_y_arg;
complex<double> D_y_arg;


double first_part=(b_g*sqrt(3.00))/(a);//for
double second_part=(a*Beta_m[m])/(2*sqrt(3.00)*B);//for
double third_part=(6*(ThermalA::intpow(B,2)-(B*b_g)))/(2*sqrt(3.00)*a*B);//for

  A_y_arg = complex<double>(first_part,-second_part);
  B_y_arg = complex<double>(first_part,second_part);
  C_y_arg = complex<double>(third_part,second_part);
  D_y_arg = complex<double>(third_part,-second_part);


double A_y=cos(argument_of_trig_fns)*imag(Faddeeva::erf(A_y_arg))*B_Nuss_0;
double B_y=cos(argument_of_trig_fns)*imag(Faddeeva::erf(B_y_arg))*B_Nuss_0;
double C_y=sin(argument_of_trig_fns)*real(Faddeeva::erf(C_y_arg))*B_Nuss_0;
double D_y=sin(argument_of_trig_fns)*real(Faddeeva::erf(A_y_arg))*B_Nuss_0;
double E_y=sin(argument_of_trig_fns)*real(Faddeeva::erf(D_y_arg))*B_Nuss_0;
double F_y=sin(argument_of_trig_fns)*real(Faddeeva::erf(B_y_arg))*B_Nuss_0;
double G_y=cos(argument_of_trig_fns)*real(Faddeeva::erf(D_y_arg))*Beta_m[m];
double H_y=cos(argument_of_trig_fns)*real(Faddeeva::erf(C_y_arg))*Beta_m[m];
double I_y=cos(argument_of_trig_fns)*real(Faddeeva::erf(A_y_arg))*Beta_m[m];
double J_y=cos(argument_of_trig_fns)*real(Faddeeva::erf(B_y_arg))*Beta_m[m];
double K_y=sin(argument_of_trig_fns)*imag(Faddeeva::erf(A_y_arg))*Beta_m[m];
double L_y=sin(argument_of_trig_fns)*imag(Faddeeva::erf(B_y_arg))*Beta_m[m];


double Omegla=(-A_y)+(B_y)+(C_y)+(D_y)+(E_y)+(F_y)+(G_y)+(H_y)+(I_y)+(J_y)+(K_y)-(L_y)+(imag(Faddeeva::erf(C_y_arg))*((-1*cos(argument_of_trig_fns)*B_Nuss_0)+(sin(argument_of_trig_fns)*Beta_m[m])))+(imag(Faddeeva::erf(D_y_arg))*((cos(argument_of_trig_fns)*B_Nuss_0)-(sin(argument_of_trig_fns)*Beta_m[m])));
double denominator=(2*B*((B_Nuss_0*(ThermalA::intpow(B_Nuss_B,2)+ThermalA::intpow(Beta_m[m],2)))+(ThermalA::intpow(B_Nuss_0,2)*(B_Nuss_B+ThermalA::intpow(B_Nuss_B,2)+ThermalA::intpow(Beta_m[m],2)))+(ThermalA::intpow(Beta_m[m],2)*(B_Nuss_B+ThermalA::intpow(B_Nuss_B,2)+ThermalA::intpow(Beta_m[m],2)))));
//double Omegla=(2*cos((Beta_m[m]*(ThermalA::intpow(D,2)-(D*d_g)))/(ThermalA::intpow(D,2)))*ThermalA::erf((6.00*(ThermalA::intpow(D,2)-(D*d_g)))/(2*sqrt(3.00)*b*D)))+(2*Faddeeva::erfi((b*Beta_m[m])/(2*sqrt(3.00)*D))*sin((Beta_m[m]*(ThermalA::intpow(D,2)-(D*d_g)))/(ThermalA::intpow(D,2))));
//double numerator=2*b*(exp(-1*((pow(Beta_m[m],2.00)*(pow(b,2.00)+(12*alpha*(t-t_prime))))/(12*pow(D,2.00)))))*(sqrt(pi/3.00))*(ThermalA::intpow(B_Nuss,2)+ThermalA::intpow(Beta_m[m],2))*(cos(Beta_m[m]-((y*Beta_m[m])/D)));
double numerator=a*exp(-((ThermalA::intpow(Beta_m[m],2)*(ThermalA::intpow(a,2)+(12*alpha*(t-t_prime))))/(12*ThermalA::intpow(B,2))))*sqrt(pi/3.00)*((sin((x*Beta_m[m])/(B))*B_Nuss_0)+(cos((x*Beta_m[m])/(B))*Beta_m[m]))*(ThermalA::intpow(B_Nuss_B,2)+ThermalA::intpow(Beta_m[m],2))*Omegla;


X33=X33+((numerator)/denominator);
//cout<<Y33<<endl;
}

//if(isnan(X33)==1){
//X33=0.0;
//}

return X33;
}

double z_rear_component_convection_at_L_convection_at_0_hs_at_vt(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double cr,double v_arc){

//  SHOULD THE 2/L BE INSIDE OR OUTSIDE THE SUM     THINK ABOUT IT

double Beta_m[637]={//for B=0.05

 0.0,0.443521,3.20399,6.31485,9.44595,12.5823,15.7207,18.8602,22.0002,25.1407,28.2814,31.4223,
 34.5633,37.7044,40.8456,43.9868,47.1281,50.2695,53.4108,56.5522,59.6936,62.835,65.9765,69.1179,
 72.2594,75.4009,78.5424,81.6839,84.8254,87.9669,91.1084,94.2499,97.3914,100.533,103.674,106.816,
 109.958,113.099,116.241,119.382,122.524,125.665,128.807,131.948,135.09,138.232,141.373,144.515,
 147.656,150.798,153.939,157.081,160.222,163.364,166.506,169.647,172.789,175.93,179.072,182.213,
 185.355,188.497,191.638,194.78,197.921,201.063,204.205,207.346,210.488,213.629,216.771,219.912,
 223.054,226.196,229.337,232.479,235.62,238.762,241.903,245.045,248.187,251.328,254.47,257.611,
 260.753,263.895,267.036,270.178,273.319,276.461,279.602,282.744,285.886,289.027,292.169,295.31,
 298.452,301.594,304.735,307.877,311.018,314.16,317.301,320.443,323.585,326.726,329.868,333.009,
 336.151,339.293,342.434,345.576,348.717,351.859,355.001,358.142,361.284,364.425,367.567,370.708,
 373.85,376.992,380.133,383.275,386.416,389.558,392.7,395.841,398.983,402.124,405.266,408.408,
 411.549,414.691,417.832,420.974,424.115,427.257,430.399,433.54,436.682,439.823,442.965,446.107,
 449.248,452.39,455.531,458.673,461.815,464.956,468.098,471.239,474.381,477.523,480.664,483.806,
 486.947,490.089,493.23,496.372,499.514,502.655,505.797,508.938,512.08,515.222,518.363,521.505,
 524.646,527.788,530.93,534.071,537.213,540.354,543.496,546.637,549.779,552.921,556.062,559.204,
 562.345,565.487,568.629,571.77,574.912,578.053,581.195,584.337,587.478,590.62,593.761,596.903,
 600.045,603.186,606.328,609.469,612.611,615.752,618.894,622.036,625.177,628.319,631.46,634.602,
 637.744,640.885,644.027,647.168,650.31,653.452,656.593,659.735,662.876,666.018,669.16,672.301,
 675.443,678.584,681.726,684.867,688.009,691.151,694.292,697.434,700.575,703.717,706.859,710.0,
 713.142,716.283,719.425,722.567,725.708,728.85,731.991,735.133,738.275,741.416,744.558,747.699,
 750.841,753.983,757.124,760.266,763.407,766.549,769.69,772.832,775.974,779.115,782.257,785.398,
 788.54,791.682,794.823,797.965,801.106,804.248,807.39,810.531,813.673,816.814,819.956,823.098,
 826.239,829.381,832.522,835.664,838.805,841.947,845.089,848.23,851.372,854.513,857.655,860.797,
 863.938,867.08,870.221,873.363,876.505,879.646,882.788,885.929,889.071,892.213,895.354,898.496,
 901.637,904.779,907.92,911.062,914.204,917.345,920.487,923.628,926.77,929.912,933.053,936.195,
 939.336,942.478,945.62,948.761,951.903,955.044,958.186,961.328,964.469,967.611,970.752,973.894,
 977.036,980.177,983.319,986.46,989.602,992.743,995.885,999.027,1002.17,1005.31,1008.45,1011.59,
 1014.73,1017.88,1021.02,1024.16,1027.3,1030.44,1033.58,1036.73,1039.87,1043.01,1046.15,1049.29,
 1052.43,1055.58,1058.72,1061.86,1065.,1068.14,1071.28,1074.42,1077.57,1080.71,1083.85,1086.99,
 1090.13,1093.27,1096.42,1099.56,1102.7,1105.84,1108.98,1112.12,1115.27,1118.41,1121.55,1124.69,
 1127.83,1130.97,1134.12,1137.26,1140.4,1143.54,1146.68,1149.82,1152.96,1156.11,1159.25,1162.39,
 1165.53,1168.67,1171.81,1174.96,1178.1,1181.24,1184.38,1187.52,1190.66,1193.81,1196.95,1200.09,
 1203.23,1206.37,1209.51,1212.65,1215.8,1218.94,1222.08,1225.22,1228.36,1231.5,1234.65,1237.79,
 1240.93,1244.07,1247.21,1250.35,1253.5,1256.64,1259.78,1262.92,1266.06,1269.2,1272.35,1275.49,
 1278.63,1281.77,1284.91,1288.05,1291.19,1294.34,1297.48,1300.62,1303.76,1306.9,1310.04,1313.19,
 1316.33,1319.47,1322.61,1325.75,1328.89,1332.04,1335.18,1338.32,1341.46,1344.6,1347.74,1350.88,
 1354.03,1357.17,1360.31,1363.45,1366.59,1369.73,1372.88,1376.02,1379.16,1382.3,1385.44,1388.58,
 1391.73,1394.87,1398.01,1401.15,1404.29,1407.43,1410.58,1413.72,1416.86,1420.,1423.14,1426.28,
 1429.42,1432.57,1435.71,1438.85,1441.99,1445.13,1448.27,1451.42,1454.56,1457.7,1460.84,1463.98,
 1467.12,1470.27,1473.41,1476.55,1479.69,1482.83,1485.97,1489.12,1492.26,1495.4,1498.54,1501.68,
 1504.82,1507.96,1511.11,1514.25,1517.39,1520.53,1523.67,1526.81,1529.96,1533.1,1536.24,1539.38,
 1542.52,1545.66,1548.81,1551.95,1555.09,1558.23,1561.37,1564.51,1567.65,1570.8,1577.08,1580.22,
 1583.36,1586.5,1589.65,1592.79,1595.93,1599.07,1602.21,1605.35,1608.5,1611.64,1614.78,1617.92,
 1621.06,1624.2,1627.35,1630.49,1633.63,1636.77,1639.91,1643.05,1646.19,1649.34,1652.48,1655.62,
 1658.76,1661.9,1665.04,1668.19,1671.33,1674.47,1677.61,1680.75,1683.89,1687.04,1690.18,1693.32,
 1696.46,1699.6,1702.74,1705.88,1709.03,1712.17,1715.31,1718.45,1721.59,1724.73,1727.88,1731.02,
 1734.16,1737.3,1740.44,1743.58,1746.73,1749.87,1753.01,1756.15,1759.29,1762.43,1765.58,1768.72,
 1771.86,1775.,1778.14,1781.28,1784.42,1787.57,1790.71,1793.85,1796.99,1800.13,1803.27,1806.42,
 1809.56,1812.7,1815.84,1818.98,1822.12,1825.27,1828.41,1831.55,1834.69,1837.83,1840.97,1844.11,
 1847.26,1850.4,1853.54,1856.68,1859.82,1862.96,1866.11,1869.25,1872.39,1875.53,1878.67,1881.81,
 1884.96,1888.1,1891.24,1894.38,1897.52,1900.66,1903.81,1906.95,1910.09,1913.23,1916.37,1919.51,
 1922.65,1925.8,1928.94,1932.08,1935.22,1938.36,1941.5,1944.65,1947.79,1950.93,1954.07,1957.21,
 1960.35,1963.5,1966.64,1969.78,1972.92,1976.06,1979.2,1982.35,1985.49,1988.63,1991.77,1994.91,1998.05

 };

double B_Nuss_0=0.1;
double B_Nuss_L=0.1;
	double alpha=(k_conductivity/(cp_steel*density));
	//cout<<"alpha"<<"\t"<<alpha<<endl;
	//cout<<"zrear"<<endl;
double Z33=0;
for(int m=0;m<228;m++){//64
//cout<<"only one term in sum"<<endl;
double argument_of_trig_fn=Beta_m[m]*((v_arc*t_prime)/L);
//double argument_of_trig_fns2=Beta_m[m]-((Beta_m[m]*d_g)/(D));

complex<double> A_z_arg;
complex<double> B_z_arg;
complex<double> C_z_arg;



double first_part=(6*t_prime*v_arc)/(2.0*sqrt(3.00)*cr);//for
double second_part=(cr*Beta_m[m])/(2.0*sqrt(3.00)*L);//for


  A_z_arg = complex<double>(second_part,-first_part);
  B_z_arg = complex<double>(second_part,0);
  C_z_arg = complex<double>(first_part,-second_part);


double A_z=cos(argument_of_trig_fn)*real(Faddeeva::erfi(A_z_arg))*B_Nuss_0;
double B_z=sin(argument_of_trig_fn)*imag(Faddeeva::erfi(B_z_arg))*B_Nuss_0;
double C_z=sin(argument_of_trig_fn)*real(Faddeeva::erf(C_z_arg))*B_Nuss_0;
double D_z=sin(argument_of_trig_fn)*imag(Faddeeva::erfi(B_z_arg))*B_Nuss_0;
double E_z=sin(argument_of_trig_fn)*imag(Faddeeva::erfi(A_z_arg))*B_Nuss_0;
double F_z=cos(argument_of_trig_fn)*imag(Faddeeva::erfi(B_z_arg))*Beta_m[m];
double G_z=cos(argument_of_trig_fn)*imag(Faddeeva::erfi(B_z_arg))*Beta_m[m];
double H_z=cos(argument_of_trig_fn)*imag(Faddeeva::erfi(A_z_arg))*Beta_m[m];
double I_z=cos(argument_of_trig_fn)*real(Faddeeva::erf(C_z_arg))*Beta_m[m];
double J_z=sin(argument_of_trig_fn)*real(Faddeeva::erfi(A_z_arg))*Beta_m[m];



double Omegla=(A_z-B_z+C_z+D_z-E_z-F_z+G_z-H_z+I_z-J_z)+(imag(Faddeeva::erf(C_z_arg))*((sin(argument_of_trig_fn)*Beta_m[m])-(cos(argument_of_trig_fn)*B_Nuss_0)))+((real(Faddeeva::erfi(B_z_arg)))*((2*sin(argument_of_trig_fn)*Beta_m[m])-(2*cos(argument_of_trig_fn)*B_Nuss_0)));
double denominator=(2*L*((B_Nuss_0*(ThermalA::intpow(B_Nuss_L,2)+ThermalA::intpow(Beta_m[m],2)))+(ThermalA::intpow(B_Nuss_0,2)*(B_Nuss_L+ThermalA::intpow(B_Nuss_L,2)+ThermalA::intpow(Beta_m[m],2)))+(ThermalA::intpow(Beta_m[m],2)*(B_Nuss_L+ThermalA::intpow(B_Nuss_L,2)+ThermalA::intpow(Beta_m[m],2)))));
double numerator=(exp(-((ThermalA::intpow(Beta_m[m],2)*((ThermalA::intpow(cr,2))+(12*alpha*(t-t_prime))))/(12*ThermalA::intpow(L,2))))*sqrt(pi/3.00)*cr*((sin((z*Beta_m[m])/(L))*B_Nuss_0)+(cos((z*Beta_m[m])/(L))*Beta_m[m]))*(ThermalA::intpow(B_Nuss_L,2)+ThermalA::intpow(Beta_m[m],2)))*Omegla;



Z33=Z33+((numerator)/denominator);
//Z33=exp(-((ThermalA::intpow(Beta_m[m],2)*((ThermalA::intpow(cr,2))+(12*alpha*(t-t_prime))))/(12*ThermalA::intpow(L,2))));
//cout<<Y33<<endl;
}

//if(isnan(Z33)==1){
//Z33=0.0;
//}

return Z33;
}


double z_front_component_convection_at_L_convection_at_0_hs_at_vt(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double cf,double v_arc){

//  SHOULD THE 2/L BE INSIDE OR OUTSIDE THE SUM     THINK ABOUT IT

double Beta_m[637]={//for B=0.05

 0.0,0.443521,3.20399,6.31485,9.44595,12.5823,15.7207,18.8602,22.0002,25.1407,28.2814,31.4223,
 34.5633,37.7044,40.8456,43.9868,47.1281,50.2695,53.4108,56.5522,59.6936,62.835,65.9765,69.1179,
 72.2594,75.4009,78.5424,81.6839,84.8254,87.9669,91.1084,94.2499,97.3914,100.533,103.674,106.816,
 109.958,113.099,116.241,119.382,122.524,125.665,128.807,131.948,135.09,138.232,141.373,144.515,
 147.656,150.798,153.939,157.081,160.222,163.364,166.506,169.647,172.789,175.93,179.072,182.213,
 185.355,188.497,191.638,194.78,197.921,201.063,204.205,207.346,210.488,213.629,216.771,219.912,
 223.054,226.196,229.337,232.479,235.62,238.762,241.903,245.045,248.187,251.328,254.47,257.611,
 260.753,263.895,267.036,270.178,273.319,276.461,279.602,282.744,285.886,289.027,292.169,295.31,
 298.452,301.594,304.735,307.877,311.018,314.16,317.301,320.443,323.585,326.726,329.868,333.009,
 336.151,339.293,342.434,345.576,348.717,351.859,355.001,358.142,361.284,364.425,367.567,370.708,
 373.85,376.992,380.133,383.275,386.416,389.558,392.7,395.841,398.983,402.124,405.266,408.408,
 411.549,414.691,417.832,420.974,424.115,427.257,430.399,433.54,436.682,439.823,442.965,446.107,
 449.248,452.39,455.531,458.673,461.815,464.956,468.098,471.239,474.381,477.523,480.664,483.806,
 486.947,490.089,493.23,496.372,499.514,502.655,505.797,508.938,512.08,515.222,518.363,521.505,
 524.646,527.788,530.93,534.071,537.213,540.354,543.496,546.637,549.779,552.921,556.062,559.204,
 562.345,565.487,568.629,571.77,574.912,578.053,581.195,584.337,587.478,590.62,593.761,596.903,
 600.045,603.186,606.328,609.469,612.611,615.752,618.894,622.036,625.177,628.319,631.46,634.602,
 637.744,640.885,644.027,647.168,650.31,653.452,656.593,659.735,662.876,666.018,669.16,672.301,
 675.443,678.584,681.726,684.867,688.009,691.151,694.292,697.434,700.575,703.717,706.859,710.0,
 713.142,716.283,719.425,722.567,725.708,728.85,731.991,735.133,738.275,741.416,744.558,747.699,
 750.841,753.983,757.124,760.266,763.407,766.549,769.69,772.832,775.974,779.115,782.257,785.398,
 788.54,791.682,794.823,797.965,801.106,804.248,807.39,810.531,813.673,816.814,819.956,823.098,
 826.239,829.381,832.522,835.664,838.805,841.947,845.089,848.23,851.372,854.513,857.655,860.797,
 863.938,867.08,870.221,873.363,876.505,879.646,882.788,885.929,889.071,892.213,895.354,898.496,
 901.637,904.779,907.92,911.062,914.204,917.345,920.487,923.628,926.77,929.912,933.053,936.195,
 939.336,942.478,945.62,948.761,951.903,955.044,958.186,961.328,964.469,967.611,970.752,973.894,
 977.036,980.177,983.319,986.46,989.602,992.743,995.885,999.027,1002.17,1005.31,1008.45,1011.59,
 1014.73,1017.88,1021.02,1024.16,1027.3,1030.44,1033.58,1036.73,1039.87,1043.01,1046.15,1049.29,
 1052.43,1055.58,1058.72,1061.86,1065.,1068.14,1071.28,1074.42,1077.57,1080.71,1083.85,1086.99,
 1090.13,1093.27,1096.42,1099.56,1102.7,1105.84,1108.98,1112.12,1115.27,1118.41,1121.55,1124.69,
 1127.83,1130.97,1134.12,1137.26,1140.4,1143.54,1146.68,1149.82,1152.96,1156.11,1159.25,1162.39,
 1165.53,1168.67,1171.81,1174.96,1178.1,1181.24,1184.38,1187.52,1190.66,1193.81,1196.95,1200.09,
 1203.23,1206.37,1209.51,1212.65,1215.8,1218.94,1222.08,1225.22,1228.36,1231.5,1234.65,1237.79,
 1240.93,1244.07,1247.21,1250.35,1253.5,1256.64,1259.78,1262.92,1266.06,1269.2,1272.35,1275.49,
 1278.63,1281.77,1284.91,1288.05,1291.19,1294.34,1297.48,1300.62,1303.76,1306.9,1310.04,1313.19,
 1316.33,1319.47,1322.61,1325.75,1328.89,1332.04,1335.18,1338.32,1341.46,1344.6,1347.74,1350.88,
 1354.03,1357.17,1360.31,1363.45,1366.59,1369.73,1372.88,1376.02,1379.16,1382.3,1385.44,1388.58,
 1391.73,1394.87,1398.01,1401.15,1404.29,1407.43,1410.58,1413.72,1416.86,1420.,1423.14,1426.28,
 1429.42,1432.57,1435.71,1438.85,1441.99,1445.13,1448.27,1451.42,1454.56,1457.7,1460.84,1463.98,
 1467.12,1470.27,1473.41,1476.55,1479.69,1482.83,1485.97,1489.12,1492.26,1495.4,1498.54,1501.68,
 1504.82,1507.96,1511.11,1514.25,1517.39,1520.53,1523.67,1526.81,1529.96,1533.1,1536.24,1539.38,
 1542.52,1545.66,1548.81,1551.95,1555.09,1558.23,1561.37,1564.51,1567.65,1570.8,1577.08,1580.22,
 1583.36,1586.5,1589.65,1592.79,1595.93,1599.07,1602.21,1605.35,1608.5,1611.64,1614.78,1617.92,
 1621.06,1624.2,1627.35,1630.49,1633.63,1636.77,1639.91,1643.05,1646.19,1649.34,1652.48,1655.62,
 1658.76,1661.9,1665.04,1668.19,1671.33,1674.47,1677.61,1680.75,1683.89,1687.04,1690.18,1693.32,
 1696.46,1699.6,1702.74,1705.88,1709.03,1712.17,1715.31,1718.45,1721.59,1724.73,1727.88,1731.02,
 1734.16,1737.3,1740.44,1743.58,1746.73,1749.87,1753.01,1756.15,1759.29,1762.43,1765.58,1768.72,
 1771.86,1775.,1778.14,1781.28,1784.42,1787.57,1790.71,1793.85,1796.99,1800.13,1803.27,1806.42,
 1809.56,1812.7,1815.84,1818.98,1822.12,1825.27,1828.41,1831.55,1834.69,1837.83,1840.97,1844.11,
 1847.26,1850.4,1853.54,1856.68,1859.82,1862.96,1866.11,1869.25,1872.39,1875.53,1878.67,1881.81,
 1884.96,1888.1,1891.24,1894.38,1897.52,1900.66,1903.81,1906.95,1910.09,1913.23,1916.37,1919.51,
 1922.65,1925.8,1928.94,1932.08,1935.22,1938.36,1941.5,1944.65,1947.79,1950.93,1954.07,1957.21,
 1960.35,1963.5,1966.64,1969.78,1972.92,1976.06,1979.2,1982.35,1985.49,1988.63,1991.77,1994.91,1998.05

 };

double B_Nuss_0=0.1;
double B_Nuss_L=0.1;
	double alpha=(k_conductivity/(cp_steel*density));
	//cout<<"alpha"<<"\t"<<alpha<<endl;
	//cout<<"zfront"<<endl;
double Z33=0;
for(int m=0;m<228;m++){//64
//cout<<"only one term in sum"<<endl;
double argument_of_trig_fn=Beta_m[m]*((v_arc*t_prime)/L);
//double argument_of_trig_fns2=Beta_m[m]-((Beta_m[m]*d_g)/(D));

complex<double> A_z_arg;
complex<double> B_z_arg;
complex<double> C_z_arg;



double first_part=(6*(L-(v_arc*t_prime)))/(2.0*sqrt(3.00)*cf);//for
double second_part=(cf*Beta_m[m])/(2.0*sqrt(3.00)*L);//for


  A_z_arg = complex<double>(second_part,0);
  B_z_arg = complex<double>(first_part,second_part);
  C_z_arg = complex<double>(first_part,-second_part);


double A_z=cos(argument_of_trig_fn)*real(Faddeeva::erfi(A_z_arg))*B_Nuss_0;
double B_z=cos(argument_of_trig_fn)*real(Faddeeva::erfi(A_z_arg))*B_Nuss_0;
double C_z=sin(argument_of_trig_fn)*imag(Faddeeva::erfi(A_z_arg))*B_Nuss_0;
double D_z=sin(argument_of_trig_fn)*real(Faddeeva::erf(B_z_arg))*B_Nuss_0;

double E_z=sin(argument_of_trig_fn)*imag(Faddeeva::erfi(A_z_arg))*B_Nuss_0;
double F_z=sin(argument_of_trig_fn)*real(Faddeeva::erf(C_z_arg))*B_Nuss_0;
double G_z=cos(argument_of_trig_fn)*imag(Faddeeva::erfi(A_z_arg))*Beta_m[m];
double H_z=cos(argument_of_trig_fn)*imag(Faddeeva::erfi(A_z_arg))*Beta_m[m];

double I_z=cos(argument_of_trig_fn)*real(Faddeeva::erf(C_z_arg))*Beta_m[m];
double J_z=cos(argument_of_trig_fn)*real(Faddeeva::erf(B_z_arg))*Beta_m[m];
double K_z=sin(argument_of_trig_fn)*real(Faddeeva::erfi(A_z_arg))*Beta_m[m];

double L_z=sin(argument_of_trig_fn)*real(Faddeeva::erfi(A_z_arg))*Beta_m[m];


double Omegla=(A_z+B_z+C_z+D_z-E_z+F_z+G_z-H_z+I_z+J_z-K_z-L_z)+(imag(Faddeeva::erf(B_z_arg))*((sin(argument_of_trig_fn)*Beta_m[m])-(cos(argument_of_trig_fn)*B_Nuss_0)))+((imag(Faddeeva::erf(C_z_arg)))*((cos(argument_of_trig_fn)*B_Nuss_0)-(sin(argument_of_trig_fn)*Beta_m[m])));
double denominator=(2*L*((B_Nuss_0*(ThermalA::intpow(B_Nuss_L,2)+ThermalA::intpow(Beta_m[m],2)))+(ThermalA::intpow(B_Nuss_0,2)*(B_Nuss_L+ThermalA::intpow(B_Nuss_L,2)+ThermalA::intpow(Beta_m[m],2)))+(ThermalA::intpow(Beta_m[m],2)*(B_Nuss_L+ThermalA::intpow(B_Nuss_L,2)+ThermalA::intpow(Beta_m[m],2)))));
double numerator=(exp(-((ThermalA::intpow(Beta_m[m],2)*((ThermalA::intpow(cf,2))+(12*alpha*(t-t_prime))))/(12*ThermalA::intpow(L,2))))*sqrt(pi/3.00)*cf*((sin((z*Beta_m[m])/(L))*B_Nuss_0)+(cos((z*Beta_m[m])/(L))*Beta_m[m]))*(ThermalA::intpow(B_Nuss_L,2)+ThermalA::intpow(Beta_m[m],2)))*Omegla;



Z33=Z33+((numerator)/denominator);
//Z33=exp(-((ThermalA::intpow(Beta_m[m],2)*((12*alpha*(t-t_prime))+(ThermalA::intpow(cf,2))))/(12*ThermalA::intpow(L,2))));
//cout<<Y33<<endl;
}
//if(isnan(Z33)==1){
//Z33=0.0;
//}
return Z33;
}
double part_past_time_integral_convection_at_all_faces_offset_in_x_and_y(double x,double y,double z,double a,double b,double cr,double cf,double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc, double B, double D, double L, double ts, double tf, double r_start, double r_fall, double b_g, double d_g){



	double rf=((2*cf)/(cr+cf));
	double rr=((2*cr)/(cr+cf));


double x_comp=x_component_convection_at_B_convection_at_0_hs_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
double y_comp=y_component_convection_at_D_convection_at_0_hs_at_dg(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);
double z_rear_comp=z_rear_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cr,v_arc);
double z_front_comp=z_front_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cf,v_arc);
////////////////////////////////////
////////////////////////////////////    Do We need to do this, does this make sense??????

double z_total_comp=(((rr/cr)*z_rear_comp)+((rf/cf)*z_front_comp));
//cout<<z_front_comp<<endl;
//Z components are going negative
double ramping=power_ramping_function(t_prime,ts,tf,r_start,r_fall);

double part_past_integral_all_convec_and_ramping=ramping*x_comp*y_comp*z_total_comp;

return part_past_integral_all_convec_and_ramping;

}


double Temperature_convection_on_all_surfaces(double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,
double cp_steel,double v_arc, double T_0, double current, double voltage, double eff, double B, double D, double L,double ts,double tf,double r_start,double r_fall,double b_g, double d_g){

if(t==0){
return T_0;
}


double integral=0;
//	int n_ROOT=0;
FINITE_plate_to_be_integrated_with_ramping_CONVECTION_ON_ALL_SURFACES W_piece1;//
W_piece1.set_values(x,y,z,a,b,cr,cf,t,k_conductivity,density,cp_steel,v_arc,B,D,L,ts,tf,r_start,r_fall,b_g,d_g);
integral=DEIntegrator<FINITE_plate_to_be_integrated_with_ramping_CONVECTION_ON_ALL_SURFACES>::Integrate(W_piece1, 0, t, 1e-20);

double Temperature_convec_all_surfaces=(T_0+(integral*((voltage*current*eff*6*sqrt(3.00))/(a*b*pi*sqrt(pi)*density*cp_steel))));

return Temperature_convec_all_surfaces;
}




















double x_DE_part_HS_at_bg(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g){
double alpha=((k_conductivity)/(cp_steel*density));
double X22=0;

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha));
double erf_denom=2*a*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha));
for(int n=-5;n<=5;n+=1){



double portion_1=exp(-((3*ThermalA::intpow((-b_g+(2*B*n)+x),2))/(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha))))*(-ThermalA::erf(((ThermalA::intpow(a,2)*((B*(-1+(2*n)))+x))+(12*(-B+b_g)*t*alpha)+(12*(B-b_g)*alpha*t_prime))/erf_denom)+ThermalA::erf(((ThermalA::intpow(a,2)*((2*B*n)+x))+(12*b_g*alpha*(t-t_prime)))/erf_denom));
double portion_2=exp(-(3*ThermalA::intpow((b_g+(2*B*n)+x),2))/(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha)))*(ThermalA::erf(((ThermalA::intpow(a,2)*(B+(2+B+n)+x))+(12*(B-b_g)*t*alpha)-(12*(B-b_g)*alpha*t_prime))/erf_denom)-ThermalA::erf(((ThermalA::intpow(a,2)*((2*B*n)+x))-(12*b_g*alpha*(t-t_prime)))/erf_denom));
double numerator=a*sqrt(alpha)*sqrt(t-t_prime)*(portion_1+portion_2);


double int_X22=numerator/denominator;

X22=X22+int_X22;
//cout<<Y22<<endl;

}
return X22;
}

double x_DE_part_HS_at_bg_insulating_at_0_dirichlet_at_B(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g){
double alpha=((k_conductivity)/(cp_steel*density));
double X22=0;

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha));
double erf_denom=2*a*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha));
for(int n=-5;n<=5;n+=1){



double portion_1=exp(-((3*ThermalA::intpow((-b_g+(2*B*n)+x),2))/(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha))))*(-ThermalA::erf(((ThermalA::intpow(a,2)*((B*(-1+(2*n)))+x))+(12*(-B+b_g)*t*alpha)+(12*(B-b_g)*alpha*t_prime))/erf_denom)+ThermalA::erf(((ThermalA::intpow(a,2)*((2*B*n)+x))+(12*b_g*alpha*(t-t_prime)))/erf_denom));
double portion_2=exp(-(3*ThermalA::intpow((b_g+(2*B*n)+x),2))/(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha)))*(ThermalA::erf(((ThermalA::intpow(a,2)*(B+(2+B+n)+x))+(12*(B-b_g)*t*alpha)-(12*(B-b_g)*alpha*t_prime))/erf_denom)-ThermalA::erf(((ThermalA::intpow(a,2)*((2*B*n)+x))-(12*b_g*alpha*(t-t_prime)))/erf_denom));
double numerator=(pow(-1.0,n))*a*sqrt(alpha)*sqrt(t-t_prime)*(portion_1+portion_2);


double int_X22=numerator/denominator;

X22=X22+int_X22;
//cout<<Y22<<endl;

}
return X22;
}
double x_DE_part_HS_at_bg_dirichlet_at_0_dirichlet_at_B(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g){
double alpha=((k_conductivity)/(cp_steel*density));
double X22=0;

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha));
double erf_denom=2*a*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha));
for(int n=-5;n<=5;n+=1){



double portion_1=exp(-((3*ThermalA::intpow((-b_g+(2*B*n)+x),2))/(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha))))*(-ThermalA::erf(((ThermalA::intpow(a,2)*((B*(-1+(2*n)))+x))+(12*(-B+b_g)*t*alpha)+(12*(B-b_g)*alpha*t_prime))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(a,2)*((2*B*n)+x))+(12*b_g*t*alpha)-(12*b_g*t_prime*alpha))/(erf_denom)));
double portion_2=exp(-(3*ThermalA::intpow((b_g+(2*B*n)+x),2))/(ThermalA::intpow(a,2)+(12*t*alpha)-(12*t_prime*alpha)))*(-ThermalA::erf(((ThermalA::intpow(a,2)*(B+(2*B*n)+x))+(12*(B-b_g)*t*alpha)-(12*(B-b_g)*alpha*t_prime))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(a,2)*((2*B*n)+x))-(12*b_g*t*alpha)+(12*b_g*t_prime*alpha))/(erf_denom)));
double numerator=a*sqrt(alpha)*sqrt(t-t_prime)*(portion_1+portion_2);


double int_X22=numerator/denominator;

X22=X22+int_X22;
//cout<<Y22<<endl;

}
return X22;
}
double y_DE_part_HS_at_dg_integrated_dg_D(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g){
double alpha=((k_conductivity)/(cp_steel*density));
double Y22=0;

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha));

double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha));
for(int n=-5;n<=5;n+=1){



double portion_1=exp(-((3*ThermalA::intpow((-d_g+(2*D*n)+y),2))/(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha))))*(ThermalA::erf((b*(d_g-(2*D*n)-y))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(b,2)*((D*(-1+(2*n)))+y))+(12*(d_g-D)*alpha*(t-t_prime)))/(erf_denom*b)));
double portion_2=exp(-(3*ThermalA::intpow((d_g+(2*D*n)+y),2))/(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha)))*(-ThermalA::erf((b*(d_g+(2*D*n)+y))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(b,2)*(D+(2*D*n)+y))+(12*(D-d_g)*t*alpha)+(12*(d_g-D)*alpha*t_prime))/(erf_denom*b)));
double numerator=b*sqrt(alpha)*sqrt(t-t_prime)*(-portion_1+portion_2);

double int_Y22=numerator/denominator;

Y22=Y22+int_Y22;

}
return Y22;
}
double z_DE_part_rear_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double v_arc, double cr){
double alpha=((k_conductivity)/(cp_steel*density));
double Z22=0;

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(cr,2)+(12*t*alpha)-(12*t_prime*alpha));

double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(cr,2)+(12*t*alpha)-(12*t_prime*alpha));
for(int n=-5;n<=5;n+=1){



double portion_1=exp(-((3*ThermalA::intpow(((2*L*n)+z-(v_arc*t_prime)),2))/(ThermalA::intpow(cr,2)+(12*t*alpha)-(12*t_prime*alpha))))*(-ThermalA::erf((cr*((2*L*n)+z-(v_arc*t_prime)))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(cr,2)*((2*L*n)+z))+(12*v_arc*alpha*(t-t_prime)*t_prime))/(erf_denom*cr)));
double portion_2=exp(-(3*ThermalA::intpow(((2*L*n)+z+(v_arc*t_prime)),2))/(ThermalA::intpow(cr,2)+(12*t*alpha)-(12*t_prime*alpha)))*(ThermalA::erf((cr*((2*L*n)+z+(v_arc*t_prime)))/(erf_denom))-ThermalA::erf(((ThermalA::intpow(cr,2)*((2*L*n)+z))-(12*t*v_arc*alpha*t_prime)+(12*v_arc*alpha*ThermalA::intpow(t_prime,2)))/(erf_denom*cr)));
double numerator=cr*sqrt(alpha)*sqrt(t-t_prime)*(portion_1+portion_2);


double int_Z22=numerator/denominator;

Z22=Z22+int_Z22;

}
return Z22;
}
double z_DE_part_front_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double v_arc, double cf){
double alpha=((k_conductivity)/(cp_steel*density));
double Z22=0;

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(cf,2)+(12*t*alpha)-(12*t_prime*alpha));

double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(cf,2)+(12*t*alpha)-(12*t_prime*alpha));
for(int n=-5;n<=5;n+=1){



double portion_1=exp(-((3*ThermalA::intpow(((2*L*n)+z-(v_arc*t_prime)),2))/(ThermalA::intpow(cf,2)+(12*t*alpha)-(12*t_prime*alpha))))*(ThermalA::erf((cf*((2*L*n)+z-(v_arc*t_prime)))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(cf,2)*(L-(2*L*n)-z))+(12*L*t*alpha)-(12*(L+(t*v_arc))*alpha*t_prime)+(12*v_arc*alpha*ThermalA::intpow(t_prime,2)))/(erf_denom*cf)));
double portion_2=exp(-(3*ThermalA::intpow(((2*L*n)+z+(v_arc*t_prime)),2))/(ThermalA::intpow(cf,2)+(12*t*alpha)-(12*t_prime*alpha)))*(-ThermalA::erf((cf*((2*L*n)+z+(v_arc*t_prime)))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(cf,2)*(L+(2*L*n)+z))+(12*L*t*alpha)-(12*alpha*t_prime*(L+(t*v_arc)-(v_arc*t_prime))))/(erf_denom*cf)));
double numerator=cf*sqrt(alpha)*sqrt(t-t_prime)*(portion_1+portion_2);


double int_Z22=numerator/denominator;

Z22=Z22+int_Z22;

}
return Z22;
}

double z_DE_part_rear_Z11(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double v_arc, double cr){
double alpha=((k_conductivity)/(cp_steel*density));
double Z11=0;

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(cr,2)+(12*t*alpha)-(12*t_prime*alpha));

double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(cr,2)+(12*t*alpha)-(12*t_prime*alpha));
for(int n=-5;n<=5;n+=1){



double portion_1=-exp(-((3*ThermalA::intpow(((2*L*n)+z-(v_arc*t_prime)),2))/(ThermalA::intpow(cr,2)+(12*t*alpha)-(12*t_prime*alpha))))*(ThermalA::erf((cr*((2*L*n)+z-(v_arc*t_prime)))/(erf_denom))
        +ThermalA::erf(((-1*ThermalA::intpow(cr,2)*((2*L*n)+z))+(12*v_arc*alpha*(-t+t_prime)*t_prime))/(erf_denom*cr)));

double portion_2=exp(-(3*ThermalA::intpow(((2*L*n)+z+(v_arc*t_prime)),2))/(ThermalA::intpow(cr,2)+(12*t*alpha)-(12*t_prime*alpha)))*(-ThermalA::erf((cr*((2*L*n)+z+(v_arc*t_prime)))/(erf_denom))
        +ThermalA::erf(((ThermalA::intpow(cr,2)*((2*L*n)+z))+(12*t_prime*v_arc*alpha*(-t+t_prime)))/(erf_denom*cr)));
double numerator=cr*sqrt(alpha)*sqrt(t-t_prime)*(portion_1+portion_2);


double int_Z11=numerator/denominator;

Z11=Z11+int_Z11;

}

return Z11;
}
double z_DE_part_front_Z11(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double v_arc, double cf){
double alpha=((k_conductivity)/(cp_steel*density));
double Z11=0;

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(cf,2)+(12*t*alpha)-(12*t_prime*alpha));

double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(cf,2)+(12*t*alpha)-(12*t_prime*alpha));
for(int n=-5;n<=5;n+=1){



double portion_1=exp(-((3*ThermalA::intpow(((2*L*n)+z-(v_arc*t_prime)),2))/(ThermalA::intpow(cf,2)+(12*t*alpha)-(12*t_prime*alpha))))*(ThermalA::erf((cf*((2*L*n)+z-(v_arc*t_prime)))/(erf_denom))
        +ThermalA::erf(((ThermalA::intpow(cf,2)*(L-(2*L*n)-z))+(12*L*t*alpha)-(12*alpha*t_prime*(L+(t*v_arc)-(v_arc*t_prime))))/(erf_denom*cf)));
double portion_2=exp(-(3*ThermalA::intpow(((2*L*n)+z+(v_arc*t_prime)),2))/(ThermalA::intpow(cf,2)+(12*t*alpha)-(12*t_prime*alpha)))*(ThermalA::erf((cf*((2*L*n)+z+(v_arc*t_prime)))/(erf_denom))
        -ThermalA::erf(((ThermalA::intpow(cf,2)*(L+(2*L*n)+z))+(12*L*t*alpha)-(12*alpha*t_prime*(L+(t*v_arc)-(v_arc*t_prime))))/(erf_denom*cf)));
double numerator=cf*sqrt(alpha)*sqrt(t-t_prime)*(portion_1+portion_2);


double int_Z11=numerator/denominator;

Z11=Z11+int_Z11;

}
return Z11;
}

double y_EB_conical_part(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y,double d_g,double y_i){
double alpha=((k_conductivity)/(cp_steel*density));
double Y22=0;
double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime);
for(int n=-5;n<=5;n+=1){


double int_Y22=(sqrt(alpha)*(ThermalA::erf((d_g-(2*D*n)-y)/erf_denom)+ThermalA::erf((d_g+(2*D*n)+y)/erf_denom)+ThermalA::erf(((2*D*n)+y-y_i)/erf_denom)-ThermalA::erf(((2*D*n)+y+y_i)/erf_denom))*sqrt(t-t_prime))/(2*sqrt(alpha*(t-t_prime)));
Y22=Y22+int_Y22;
//cout<<Y22<<endl;

}


return Y22;
}

double partpastintegral_DEC_EB_WELD_HEAT_SOURCE(double x,double y,double z,double a,double b,double cr,double cf,
double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc, double B,double D,
double L, double b_g, double d_g, double y_i,double t_start_ramp,double t_finish_ramp,double ramp_rate_start,
double ramp_rate_end,int xBC,int yBC,int zBC){


double r_r=(6*b*cr*ThermalA::intpow(eulers,3))/((cr+cf)*((3*b*ThermalA::intpow(eulers,3))+(2*(ThermalA::intpow(eulers,3)-1)*sqrt(3*pi)*(d_g-y_i))));

double r_f=(6*b*cf*ThermalA::intpow(eulers,3))/((cr+cf)*((3*b*ThermalA::intpow(eulers,3))+(2*(ThermalA::intpow(eulers,3)-1)*sqrt(3*pi)*(d_g-y_i))));

double rc_r=(4*cr*(ThermalA::intpow(eulers,3)-1)*(d_g-y_i))/((cf+cr)*((-2*d_g)+(ThermalA::intpow(eulers,3)*((2*d_g)+(b*sqrt(3.0/pi))-(2*y_i)))+(2*y_i)));

double rc_f=(4*cf*(ThermalA::intpow(eulers,3)-1)*(d_g-y_i))/((cf+cr)*((-2*d_g)+(ThermalA::intpow(eulers,3)*((2*d_g)+(b*sqrt(3.0/pi))-(2*y_i)))+(2*y_i)));


//x_component_convection_at_B_convection_at_0_hs_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
//y_component_convection_at_D_convection_at_0_hs_at_dg(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);//y_i=0 in this case!
//z_rear_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cr,v_arc);
//z_front_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cf,v_arc);

double x_DE,y_DE,z_r_DE,z_f_DE,x_C,y_C,z_r_C,z_f_C;

if(xBC==11){
x_DE=x_DE_part_HS_at_bg_dirichlet_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
x_C=x_DE_part_HS_at_bg_dirichlet_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
}
else if(xBC==21){
x_DE=x_DE_part_HS_at_bg_insulating_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
x_C=x_DE_part_HS_at_bg_insulating_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
}
else if(xBC==22){
x_DE=x_DE_part_HS_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
x_C=x_DE_part_HS_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
}
else if(xBC==33){
x_DE=x_component_convection_at_B_convection_at_0_hs_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
x_C=x_component_convection_at_B_convection_at_0_hs_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
}
else{
cerr<<"The Boundary Conditions "<<xBC<<" are not yet included in HEDSATS. Please contact the author if you would like them"<<endl;
}


if(yBC==22){
y_DE=y_DE_part_HS_at_dg_integrated_dg_D(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);
y_C=y_EB_conical_part(k_conductivity,cp_steel,density,t,t_prime,D,y,d_g,y_i);
}
else if(yBC==33){
y_DE=y_component_convection_at_D_convection_at_0_hs_at_dg(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);//y_i=0 in this case!
y_C=y_component_convection_at_D_convection_at_0_hs_at_dg(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);//y_i=0 in this case!
}
else{
cerr<<"The Boundary Conditions "<<yBC<<" are not yet included in HEDSATS. Please contact the author if you would like them"<<endl;
}


if(zBC==11){
z_r_DE=z_DE_part_rear_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
z_f_DE=z_DE_part_front_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);
z_r_C=z_DE_part_rear_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
z_f_C=z_DE_part_front_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);
}
else if(zBC==22){
z_r_DE=z_DE_part_rear_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
z_f_DE=z_DE_part_front_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);
z_r_C=z_DE_part_rear_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
z_f_C=z_DE_part_front_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);
}
else if(zBC==33){
z_r_DE=z_rear_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cr,v_arc);
z_f_DE=z_front_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cf,v_arc);
z_r_C=z_rear_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cr,v_arc);
z_f_C=z_front_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cf,v_arc);
}
else{
cerr<<"The Boundary Conditions "<<zBC<<" are not yet included in HEDSATS. Please contact the author"<<endl;
}

//double x_DE=x_DE_part_HS_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
//double x_DE=x_DE_part_HS_at_bg_insulating_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
//double x_DE=x_DE_part_HS_at_bg_dirichlet_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
//double y_DE=y_DE_part_HS_at_dg_integrated_dg_D(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);
//double y_DE=y_component_convection_at_D_convection_at_0_hs_at_dg(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);//y_i=0 in this case!


//double z_r_DE=z_DE_part_rear_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
//double z_f_DE=z_DE_part_front_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);
//double z_r_DE=z_DE_part_rear_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
//double z_f_DE=z_DE_part_front_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);


//double x_C=x_DE_part_HS_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
//double x_C=x_DE_part_HS_at_bg_insulating_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
//double x_C=x_DE_part_HS_at_bg_dirichlet_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
//double y_C=y_EB_conical_part(k_conductivity,cp_steel,density,t,t_prime,D,y,d_g,y_i);
//
//double z_r_C=z_DE_part_rear_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
//double z_f_C=z_DE_part_front_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);
//double z_r_C=z_DE_part_rear_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
//double z_f_C=z_DE_part_front_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);


//double z_DE_part_rear_Z11()
//double z_DE_part_front_Z11()

double rear_conical_part=((54.0*ThermalA::intpow(eulers,3)*rc_r)/(ThermalA::intpow(pi,2)*(ThermalA::intpow(eulers,3)-1)*(2*cr*3*a)*(d_g-y_i)))*(x_C*y_C*z_r_C);
double front_conical_part=((54.0*ThermalA::intpow(eulers,3)*rc_f)/(ThermalA::intpow(pi,2)*(ThermalA::intpow(eulers,3)-1)*(2*cf*3*a)*(d_g-y_i)))*(x_C*y_C*z_f_C);
double rear_double_ellipsoidal_part=((6.0*sqrt(3.0)*r_r)/(a*cr*b*pi*sqrt(pi)))*(x_DE*y_DE*z_r_DE);
double front_double_ellipsoidal_part=((6.0*sqrt(3.0)*r_f)/(a*cf*b*pi*sqrt(pi)))*(x_DE*y_DE*z_f_DE);


//cout<<"rear conical part= "<<rear_conical_part<<"\t"<<"front conical part= "<<front_conical_part<<"\t"<<"rear DE part= "<<rear_double_ellipsoidal_part<<"\t"<<"front DE part= "<<front_double_ellipsoidal_part<<endl;

double ramping=power_ramping_function(t_prime,t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end);
////////                    have different ramping rates here
//if(t_prime>(L/v_arc)){
//ramping=0;
//}
//else{
//ramping=1;
//}

double combine_portions;

combine_portions=(1/(density*cp_steel))*(rear_conical_part+front_conical_part+rear_double_ellipsoidal_part+front_double_ellipsoidal_part);

return (combine_portions*ramping);

}
double part_past_integral_DE_HEAT_SOURCE(double x,double y,double z,double a,double b,double cr,double cf,
double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc, double B, double D,
double L, double b_g, double d_g,double t_start_ramp,double t_finish_ramp,
double ramp_rate_start,double ramp_rate_end,int xBC,int yBC,int zBC){

	double rf=((2*cf)/(cr+cf));
	double rr=((2*cr)/(cr+cf));

//x_component_convection_at_B_convection_at_0_hs_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
//y_component_convection_at_D_convection_at_0_hs_at_dg(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);//y_i=0 in this case!
//z_rear_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cr,v_arc);
//z_front_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cf,v_arc);
double x_DE,y_DE,z_r_DE,z_f_DE,x_C,y_C,z_r_C,z_f_C;

if(xBC==11){
x_DE=x_DE_part_HS_at_bg_dirichlet_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
}
else if(xBC==21){
x_DE=x_DE_part_HS_at_bg_insulating_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
}
else if(xBC==22){
x_DE=x_DE_part_HS_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
}
else if(xBC==33){
x_DE=x_component_convection_at_B_convection_at_0_hs_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
}
else{
cerr<<"The Boundary Conditions "<<xBC<<" are not yet included in HEDSATS. Please contact the author if you would like them"<<endl;
}


if(yBC==22){
y_DE=y_DE_part_HS_at_dg_integrated_dg_D(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);
}
else if(yBC==33){
//y_DE=x_component_convection_at_B_convection_at_0_hs_at_bg(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);//
y_DE=y_component_convection_at_D_convection_at_0_hs_at_dg(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);//y_i=0 in this case!
}
else{
cerr<<"The Boundary Conditions "<<yBC<<" are not yet included in HEDSATS. Please contact the author if you would like them"<<endl;
}


if(zBC==11){
z_r_DE=z_DE_part_rear_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
z_f_DE=z_DE_part_front_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);
}
else if(zBC==22){
z_r_DE=z_DE_part_rear_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
z_f_DE=z_DE_part_front_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);
}
else if(zBC==33){
z_r_DE=z_rear_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cr,v_arc);
z_f_DE=z_front_component_convection_at_L_convection_at_0_hs_at_vt(k_conductivity,cp_steel,density,t,t_prime,L,z,cf,v_arc);
}
else{
cerr<<"The Boundary Conditions "<<zBC<<" are not yet included in HEDSATS. Please contact the author"<<endl;
}


//double x_DE=x_DE_part_HS_at_bg(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
////double x_DE=x_DE_part_HS_at_bg_insulating_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
////double x_DE=x_DE_part_HS_at_bg_dirichlet_at_0_dirichlet_at_B(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
////double y_DE=y_DE_part_HS_at_dg_integrated_dg_D(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);
//double y_DE=y_component_convection_at_D_convection_at_0_hs_at_dg(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);//y_i=0 in this case!
//
//double z_r_DE=z_DE_part_rear_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
//double z_f_DE=z_DE_part_front_component(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);
//double z_r_DE=z_DE_part_rear_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cr);
//double z_f_DE=z_DE_part_front_Z11(k_conductivity,cp_steel,density,t,t_prime,L,z,v_arc,cf);


double rear_part_DE=((6.0*sqrt(3.0)*rr)/(a*b*pi*sqrt(pi)*cr))*(x_DE*y_DE*z_r_DE);
double front_part_DE=((6.0*sqrt(3.0)*rf)/(a*b*pi*sqrt(pi)*cf))*(x_DE*y_DE*z_f_DE);

double ramping=power_ramping_function(t_prime,t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end);
////////                    have different ramping rates here
//if(t_prime>(L/v_arc)){
//ramping=0;
//}
//else{
//ramping=1;
//}


double part_past_integral_finite=(1/(density*cp_steel))*(rear_part_DE+front_part_DE);


return (part_past_integral_finite*ramping);

}
//double Temperature_simpsons_integrated(double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,double cp_steel,double v_arc, double T_0, double B, double D, double L, double b_g, double d_g,double voltage, double current, double eff, double y_i){
//
//
//double htrap=5E-5;
//double k_it;
////double lowlimitt=0;////////////
////double highlimitt=t;
//double tot=0;
////cout<<tot<<"\t"<<t<<endl;
//for(k_it=0;k_it<t;k_it+=htrap){
////cout<<k_it<<endl;
//	double a_trap=k_it;
//double b_trap=a_trap+htrap;
//
////double slice=((b_trap-a_trap)/6)*((part_past_integral_taking_arrays(x,y,z,a,b,cr,cf,t,a_trap,k_of_t_array,density,cp_of_t_array,v_arc,B,D,L,ts,tf,r_start,r_fall,tstart,tgap,array_sizeovert))+(4*part_past_integral_taking_arrays(x,y,z,a,b,cr,cf,t,((a_trap+b_trap)/2),k_of_t_array,density,cp_of_t_array,v_arc,B,D,L,ts,tf,r_start,r_fall,tstart,tgap,array_sizeovert))+(part_past_integral_taking_arrays(x,y,z,a,b,cr,cf,t,b_trap,k_of_t_array,density,cp_of_t_array,v_arc,B,D,L,ts,tf,r_start,r_fall,tstart,tgap,array_sizeovert)));
//
//
//
//double slice=((b_trap-a_trap)/6)*((partpastintegral_DEC_EB_WELD_HEAT_SOURCE(x,y,z,a,b,cr,cf,t,a_trap,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i))+(4*partpastintegral_DEC_EB_WELD_HEAT_SOURCE(x,y,z,a,b,cr,cf,t,((a_trap+b_trap)/2),k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i))+(partpastintegral_DEC_EB_WELD_HEAT_SOURCE(x,y,z,a,b,cr,cf,t,b_trap,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i)));
//if(isnan(slice)==1){
//slice=0;
//}
//
////partpastintegral_DEC_EB_WELD_HEAT_SOURCE(x,y,z,a,b,cr,cf,t,a_trap,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i)
////partpastintegral_DEC_EB_WELD_HEAT_SOURCE(x,y,z,a,b,cr,cf,t,((a_trap+b_trap)/2),k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i)
////partpastintegral_DEC_EB_WELD_HEAT_SOURCE(x,y,z,a,b,cr,cf,t,b_trap,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i)
//
////cout<<slice<<"\t"<<tot<<endl;
//
//tot=tot+slice;
//}
//cout<<tot<<endl;
//return (T_0+(tot*voltage*current*eff));
//}
double Temp_Finite_EB(double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,double cp_steel,
double v_arc, double T_0, double current, double voltage, double eff, double B, double D, double L,double b_g, double d_g,double y_i,double t_start_ramp,
double t_finish_ramp,double ramp_rate_start,double ramp_rate_end,int xBC,int yBC,int zBC){

if(t==0){
return T_0;
}
/*
if(t>0&&t<=600.0){
double integral=0;
FINITE_plate_to_be_integrated_EB W_piece1;//
W_piece1.set_values(x,y,z,a,b,cr,cf,t,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i);

integral=DEIntegrator<FINITE_plate_to_be_integrated_EB>::Integrate(W_piece1, 0, t, 1e-20);
double Temp_Finite_EB=(T_0+((integral)*voltage*current*eff));
return Temp_Finite_EB;
}
*/
//if(t>600.0){
double integral1=0;
double integral2=0;
double integral3=0;
double integral4=0;
double integral5=0;

FINITE_plate_to_be_integrated_EB W_piece1;//
W_piece1.set_values(x,y,z,a,b,cr,cf,t,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i,t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end,xBC,yBC,zBC);

integral1=DEIntegrator<FINITE_plate_to_be_integrated_EB>::Integrate(W_piece1, 0, t/16, 1e-20);
integral2=DEIntegrator<FINITE_plate_to_be_integrated_EB>::Integrate(W_piece1, t/16, (2*t)/16, 1e-20);
integral3=DEIntegrator<FINITE_plate_to_be_integrated_EB>::Integrate(W_piece1, (2*t)/16, (3*t)/16, 1e-20);
integral4=DEIntegrator<FINITE_plate_to_be_integrated_EB>::Integrate(W_piece1, (3*t)/16,t/2, 1e-20);
integral5=DEIntegrator<FINITE_plate_to_be_integrated_EB>::Integrate(W_piece1, t/2, t, 1e-20);

double Temp_Finite_EB=(T_0+((integral1+integral2+integral3+integral4+integral5)*voltage*current*eff));
return Temp_Finite_EB;
//}


}

double Temp_Finite_DE(double x,double y,double z,double a,double b,double cr,double cf,double t,
double k_conductivity,double density,double cp_steel,double v_arc, double T_0, double current, double voltage, double eff,
 double B, double D, double L,double b_g, double d_g,double t_start_ramp,double t_finish_ramp,double ramp_rate_start,double ramp_rate_end,int xBC,int yBC,int zBC){

if(t==0){
return T_0;
}
/*
if(t>0&&t<=600.0){
double integral=0;
FINITE_plate_to_be_integrated_DE_Heat_source W_piece1;//
W_piece1.set_values(x,y,z,a,b,cr,cf,t,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g);
integral=DEIntegrator<FINITE_plate_to_be_integrated_DE_Heat_source>::Integrate(W_piece1, 0, t, 1e-20);
double Temp_Finite_DE=(T_0+((integral)*voltage*current*eff));
return Temp_Finite_DE;
}
*/
//if(t>600.0){
double integral1=0;
double integral2=0;
double integral3=0;
double integral4=0;
double integral5=0;

FINITE_plate_to_be_integrated_DE_Heat_source W_piece1;//
W_piece1.set_values(x,y,z,a,b,cr,cf,t,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end,xBC,yBC,zBC);

integral1=DEIntegrator<FINITE_plate_to_be_integrated_DE_Heat_source>::Integrate(W_piece1, 0, t/16, 1e-20);
integral2=DEIntegrator<FINITE_plate_to_be_integrated_DE_Heat_source>::Integrate(W_piece1, t/16, t/8, 1e-20);
integral3=DEIntegrator<FINITE_plate_to_be_integrated_DE_Heat_source>::Integrate(W_piece1, t/8, (3*t)/16, 1e-20);
integral4=DEIntegrator<FINITE_plate_to_be_integrated_DE_Heat_source>::Integrate(W_piece1, (3*t)/16, t/2, 1e-20);
integral5=DEIntegrator<FINITE_plate_to_be_integrated_DE_Heat_source>::Integrate(W_piece1, t/2, t, 1e-20);

double Temp_Finite_DE=(T_0+((integral1+integral2+integral3+integral4+integral5)*voltage*current*eff));
return Temp_Finite_DE;
//}

}

double Temperature_from_combined_DECbeam_DE_heat_source(double x,double y,double z,double a_DE_in,double b_DE_in,double cr_DE_in,double cf_DE_in,double a_EB_in,double b_EB_in,
double cr_EB_in,double cf_EB_in,double t,double k_conductivity,double density,double cp_steel,double v_arc, double T_0, double current, double voltage, double eff, double B, double D,
 double L,double b_g, double d_g,double y_i,double portion_of_heat_in_DE,double t_start_ramp,double t_finish_ramp,double ramp_rate_start,double ramp_rate_end,int xBC,int yBC,int zBC){

double portion_of_heat_in_EB=1.0-portion_of_heat_in_DE;

double DE_component=portion_of_heat_in_DE*ThermalA::Temp_Finite_DE(x,y,z,a_DE_in,b_DE_in,cr_DE_in,cf_DE_in,t,k_conductivity,density,cp_steel,v_arc,0,current,voltage,eff,B,D,L,b_g,0.0,t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end,xBC,yBC,zBC);
double EB_component=portion_of_heat_in_EB*ThermalA::Temp_Finite_EB(x,y,z,a_EB_in,b_EB_in,cr_EB_in,cf_EB_in,t,k_conductivity,density,cp_steel,v_arc,0,current,voltage,eff,B,D,L,b_g,d_g,y_i,t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end,xBC,yBC,zBC);
double total_component=T_0+(DE_component+EB_component);
return total_component;
}









double FSW_Heat_Source_Component_x_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double r,double b_g){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double X_22=0;


double erf_denom=2*pi*r*sqrt(alpha)*sqrt(t-t_prime)*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha));
for(int n=-10;n<=10;n+=1){

double numerator=pi*r*sqrt(alpha)*sqrt(t-t_prime)*((exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow((-b_g+(2*B*n)+x),2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(r,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha))))*(-ThermalA::erf(((ThermalA::intpow(pi*r,2)*((B*(-1+(2*n)))+x))+(108*(-B+b_g)*ThermalA::intpow(k_p,2)*t*alpha)+(108*(B-b_g)*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)
+ThermalA::erf(((ThermalA::intpow(pi*r,2)*((2*B*n)+x))+(108*b_g*ThermalA::intpow(k_p,2)*t*alpha)-(108*b_g*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)))+(exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow((b_g+(2*B*n)+x),2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(r,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha))))*
(-ThermalA::erf(((ThermalA::intpow(pi*r,2)*((2*B*n)+x))-(108*b_g*ThermalA::intpow(k_p,2)*t*alpha)+(108*b_g*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)+ThermalA::erf(((ThermalA::intpow(pi*r,2)*(B+(2*B*n)+x))+(108*(B-b_g)*ThermalA::intpow(k_p,2)*t*alpha)+(108*(-B+b_g)*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom))));

double int_X_22=numerator/denominator;

X_22=X_22+int_X_22;

}
return X_22;

}

double FSW_Heat_Source_Component_x_NOT_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double r,double b_g){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double X_22=0;


double erf_denom=2*r*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(r,2)+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt((ThermalA::intpow(r,2))+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha));
for(int n=-10;n<=10;n+=1){

double numerator=r*sqrt(alpha)*sqrt(t-t_prime)*((exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow((-b_g+(2*B*n)+x),2))/((ThermalA::intpow(r,2))+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha))))*(-ThermalA::erf(((ThermalA::intpow(r,2)*((B*(-1+(2*n)))+x))+(12*(-B+b_g)*ThermalA::intpow(k_p,2)*t*alpha)+(12*(B-b_g)*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)
+ThermalA::erf(((ThermalA::intpow(r,2)*((2*B*n)+x))+(12*b_g*ThermalA::intpow(k_p,2)*t*alpha)-(12*b_g*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)))+(exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow((b_g+(2*B*n)+x),2))/((ThermalA::intpow(r,2))+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha))))*
(-ThermalA::erf(((ThermalA::intpow(r,2)*((2*B*n)+x))-(12*b_g*ThermalA::intpow(k_p,2)*t*alpha)+(12*b_g*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)+ThermalA::erf(((ThermalA::intpow(r,2)*(B+(2*B*n)+x))+(12*(B-b_g)*ThermalA::intpow(k_p,2)*t*alpha)+(12*(-B+b_g)*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom))));

double int_X_22=numerator/denominator;

X_22=X_22+int_X_22;

}
return X_22;

}

double FSW_Heat_Source_Component_z_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r,double v_arc){
//
//
//FOR THE CASE WHERE r_r==r_f
//
//
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double Z_22=0;


double erf_denom=2*pi*r*sqrt(alpha)*sqrt(t-t_prime)*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha));
for(int n=-10;n<=10;n+=1){

double first_exponential=exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z-(v_arc*t_prime),2))/(((ThermalA::intpow(pi,2)*ThermalA::intpow(r,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha)))));
double second_exponential=exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z+(v_arc*t_prime),2))/(((ThermalA::intpow(pi,2)*ThermalA::intpow(r,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha)))));

double numerator=pi*r*sqrt(alpha)*sqrt(t-t_prime)*((first_exponential*(ThermalA::erf(((ThermalA::intpow(pi*r,2)*((2*L*n)+z))+(108*ThermalA::intpow(k_p,2)*v_arc*(t-t_prime)*t_prime*alpha))/erf_denom)
+ThermalA::erf(((ThermalA::intpow(pi*r,2)*(L-(2*L*n)-z))+(108*L*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha*(L+(v_arc*t)-(v_arc*t_prime))))/erf_denom)))+(second_exponential*
(-ThermalA::erf(((ThermalA::intpow(pi*r,2)*((2*L*n)+z))+(108*ThermalA::intpow(k_p,2)*v_arc*t_prime*(-t+t_prime)*alpha))/erf_denom)+ThermalA::erf(((ThermalA::intpow(pi*r,2)*(L+(2*L*n)+z))+(108*ThermalA::intpow(k_p,2)*L*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha*(L+(v_arc*t)-(v_arc*t_prime))))/erf_denom))));

double int_Z_22=numerator/denominator;


Z_22=Z_22+int_Z_22;

}
return Z_22;

}


double FSW_Heat_Source_Component_z_NOT_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r,double v_arc){
//
//
//FOR THE CASE WHERE r_r==r_f
//
//
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double Z_22=0;


double erf_denom=2*r*sqrt(alpha)*sqrt(t-t_prime)*sqrt((ThermalA::intpow(r,2))+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha));

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt((ThermalA::intpow(r,2))+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha));
for(int n=-10;n<=10;n+=1){

double first_exponential=exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z-(v_arc*t_prime),2))/(((ThermalA::intpow(r,2))+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha)))));
double second_exponential=exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z+(v_arc*t_prime),2))/(((ThermalA::intpow(r,2))+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha)))));

double numerator=r*sqrt(alpha)*sqrt(t-t_prime)*((first_exponential*(ThermalA::erf(((ThermalA::intpow(r,2)*((2*L*n)+z))+(12*ThermalA::intpow(k_p,2)*v_arc*(t-t_prime)*t_prime*alpha))/erf_denom)
+ThermalA::erf(((ThermalA::intpow(r,2)*(L-(2*L*n)-z))+(12*L*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha*(L+(v_arc*t)-(v_arc*t_prime))))/erf_denom)))+(second_exponential*
(-ThermalA::erf(((ThermalA::intpow(r,2)*((2*L*n)+z))+(12*ThermalA::intpow(k_p,2)*v_arc*t_prime*(-t+t_prime)*alpha))/erf_denom)+ThermalA::erf(((ThermalA::intpow(r,2)*(L+(2*L*n)+z))+(12*ThermalA::intpow(k_p,2)*L*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha*(L+(v_arc*t)-(v_arc*t_prime))))/erf_denom))));

double int_Z_22=numerator/denominator;


Z_22=Z_22+int_Z_22;

}
return Z_22;

}

double FSW_Heat_Source_y_column_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y,double d_g,double y_i){

double alpha=((k_conductivity)/(cp_steel*density));
double Y22=0;
double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime);
for(int n=-10;n<=10;n+=1){


double int_Y22=(sqrt(alpha)*(ThermalA::erf((d_g-(2*D*n)-y)/erf_denom)+ThermalA::erf((d_g+(2*D*n)+y)/erf_denom)+ThermalA::erf(((2*D*n)+y-y_i)/erf_denom)-ThermalA::erf(((2*D*n)+y+y_i)/erf_denom))*sqrt(t-t_prime))/(2*sqrt(alpha*(t-t_prime)));
Y22=Y22+int_Y22;
//cout<<Y22<<endl;

}


return Y22;
}

double FSW_Heat_Source_y_doughnut_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g){
//
//
//FOR THE CASE WHERE r_r==r_f
//
//
double alpha=((k_conductivity)/(cp_steel*density));
double Y22=0;

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha));

double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha));
for(int n=-10;n<=10;n+=1){



double portion_1=exp(-((3*ThermalA::intpow((-d_g+(2*D*n)+y),2))/(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha))))*(ThermalA::erf((b*(d_g-(2*D*n)-y))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(b,2)*((D*(-1+(2*n)))+y))+(12*(d_g-D)*alpha*(t-t_prime)))/(erf_denom*b)));
double portion_2=exp(-(3*ThermalA::intpow((d_g+(2*D*n)+y),2))/(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha)))*(-ThermalA::erf((b*(d_g+(2*D*n)+y))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(b,2)*(D+(2*D*n)+y))+(12*(D-d_g)*t*alpha)+(12*(d_g-D)*alpha*t_prime))/(erf_denom*b)));
double numerator=b*sqrt(alpha)*sqrt(t-t_prime)*(-portion_1+portion_2);

double int_Y22=numerator/denominator;

Y22=Y22+int_Y22;

}
return Y22;
}







double FSW_4Q_Heat_Source_Component_x_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double X_22=0;


double erf_denom=2*pi*a*sqrt(alpha)*sqrt(t-t_prime)*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(a,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(a,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha));
for(int n=-10;n<=10;n+=1){

double numerator=pi*a*sqrt(alpha)*sqrt(t-t_prime)*((exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow((-b_g+(2*B*n)+x),2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(a,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha))))*(-ThermalA::erf(((ThermalA::intpow(pi*a,2)*((B*(-1+(2*n)))+x))+(108*(-B+b_g)*ThermalA::intpow(k_p,2)*t*alpha)+(108*(B-b_g)*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)
+ThermalA::erf(((ThermalA::intpow(pi*a,2)*((2*B*n)+x))+(108*b_g*ThermalA::intpow(k_p,2)*t*alpha)-(108*b_g*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)))+(exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow((b_g+(2*B*n)+x),2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(a,2))+(108*ThermalA::intpow(k_p,2)*t*alpha)-(108*ThermalA::intpow(k_p,2)*t_prime*alpha))))*
(-ThermalA::erf(((ThermalA::intpow(pi*a,2)*((2*B*n)+x))-(108*b_g*ThermalA::intpow(k_p,2)*t*alpha)+(108*b_g*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)+ThermalA::erf(((ThermalA::intpow(pi*a,2)*(B+(2*B*n)+x))+(108*(B-b_g)*ThermalA::intpow(k_p,2)*t*alpha)+(108*(-B+b_g)*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom))));

double int_X_22=numerator/denominator;

X_22=X_22+int_X_22;

}
return X_22;

}
double FSW_4Q_Heat_Source_Component_x_NOT_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double X_22=0;


double erf_denom=2*a*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(a,2)+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt((ThermalA::intpow(a,2))+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha));
for(int n=-10;n<=10;n+=1){

double numerator=a*sqrt(alpha)*sqrt(t-t_prime)*((exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow((-b_g+(2*B*n)+x),2))/((ThermalA::intpow(a,2))+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha))))*(-ThermalA::erf(((ThermalA::intpow(a,2)*((B*(-1+(2*n)))+x))+(12*(-B+b_g)*ThermalA::intpow(k_p,2)*t*alpha)+(12*(B-b_g)*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)
+ThermalA::erf(((ThermalA::intpow(a,2)*((2*B*n)+x))+(12*b_g*ThermalA::intpow(k_p,2)*t*alpha)-(12*b_g*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)))+(exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow((b_g+(2*B*n)+x),2))/((ThermalA::intpow(a,2))+(12*ThermalA::intpow(k_p,2)*t*alpha)-(12*ThermalA::intpow(k_p,2)*t_prime*alpha))))*
(-ThermalA::erf(((ThermalA::intpow(a,2)*((2*B*n)+x))-(12*b_g*ThermalA::intpow(k_p,2)*t*alpha)+(12*b_g*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom)+ThermalA::erf(((ThermalA::intpow(a,2)*(B+(2*B*n)+x))+(12*(B-b_g)*ThermalA::intpow(k_p,2)*t*alpha)+(12*(-B+b_g)*ThermalA::intpow(k_p,2)*t_prime*alpha))/erf_denom))));

double int_X_22=numerator/denominator;

X_22=X_22+int_X_22;

}
return X_22;
}
double FSW_4Q_Heat_Source_Component_z_P_Front_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_f,double v_arc){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double Z_22=0;


double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
for(int n=-10;n<=10;n+=1){

double first_exponential=exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z-(v_arc*t_prime),2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));
double second_exponential=exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z+(v_arc*t_prime),2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));

double numerator=pi*r_f*sqrt(alpha)*sqrt(t-t_prime)*((first_exponential*(ThermalA::erf((pi*r_f*((2*L*n)+z-(v_arc*t_prime)))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2)*(L-(2*L*n)-z))+(108*ThermalA::intpow(k_p,2)*L*t*alpha)-(108*ThermalA::intpow(k_p,2)*alpha*t_prime*(L+(v_arc*t)-(v_arc*t_prime))))/(pi*r_f*erf_denom))))+
(second_exponential*(-ThermalA::erf((pi*r_f*((2*L*n)+z+(v_arc*t_prime)))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2)*(L+(2*L*n)+z))+(108*ThermalA::intpow(k_p,2)*L*t*alpha)-(108*ThermalA::intpow(k_p,2)*alpha*t_prime*(L+(v_arc*t)-(v_arc*t_prime))))/(pi*r_f*erf_denom)))));


double int_Z_22=numerator/denominator;


Z_22=Z_22+int_Z_22;

}
return Z_22;
}
double FSW_4Q_Heat_Source_Component_z_NOT_P_Front_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_f,double v_arc){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double Z_22=0;


double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(r_f,2)+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(r_f,2)+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
for(int n=-10;n<=10;n+=1){

double first_exponential=exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z-(v_arc*t_prime),2))/((ThermalA::intpow(r_f,2))+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));
double second_exponential=exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z+(v_arc*t_prime),2))/((ThermalA::intpow(r_f,2))+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));

double numerator=r_f*sqrt(alpha)*sqrt(t-t_prime)*((first_exponential*(ThermalA::erf((r_f*((2*L*n)+z-(v_arc*t_prime)))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(r_f,2)*(L-(2*L*n)-z))+(12*ThermalA::intpow(k_p,2)*L*t*alpha)-(12*ThermalA::intpow(k_p,2)*alpha*t_prime*(L+(v_arc*t)-(v_arc*t_prime))))/(r_f*erf_denom))))+
(second_exponential*(-ThermalA::erf((r_f*((2*L*n)+z+(v_arc*t_prime)))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(r_f,2)*(L+(2*L*n)+z))+(12*ThermalA::intpow(k_p,2)*L*t*alpha)-(12*ThermalA::intpow(k_p,2)*alpha*t_prime*(L+(v_arc*t)-(v_arc*t_prime))))/(r_f*erf_denom)))));


double int_Z_22=numerator/denominator;


Z_22=Z_22+int_Z_22;

}
return Z_22;
}
double FSW_4Q_Heat_Source_Component_z_P_Rear_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_r,double v_arc){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double Z_22=0;


double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
for(int n=-10;n<=10;n+=1){

double first_exponential=exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z-(v_arc*t_prime),2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));
double second_exponential=exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z+(v_arc*t_prime),2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));

double numerator=pi*r_r*sqrt(alpha)*sqrt(t-t_prime)*((first_exponential*(-ThermalA::erf((pi*r_r*((2*L*n)+z-(v_arc*t_prime)))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2)*((2*L*n)+z))+(108*ThermalA::intpow(k_p,2)*v_arc*alpha*(t-t_prime)*t_prime))/(pi*r_r*erf_denom))))+
(second_exponential*(ThermalA::erf((pi*r_r*((2*L*n)+z+(v_arc*t_prime)))/(erf_denom))
-ThermalA::erf(((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2)*((2*L*n)+z))+(108*ThermalA::intpow(k_p,2)*v_arc*alpha*t_prime*(t_prime-t)))/(pi*r_r*erf_denom)))));


double int_Z_22=numerator/denominator;


Z_22=Z_22+int_Z_22;

}
return Z_22;
}
double FSW_4Q_Heat_Source_Component_z_NOT_P_Rear_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_r,double v_arc){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double Z_22=0;


double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(r_r,2)+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(r_r,2)+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
for(int n=-10;n<=10;n+=1){

double first_exponential=exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z-(v_arc*t_prime),2))/((ThermalA::intpow(r_r,2))+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));
double second_exponential=exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow((2*L*n)+z+(v_arc*t_prime),2))/((ThermalA::intpow(r_r,2))+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));

double numerator=r_r*sqrt(alpha)*sqrt(t-t_prime)*((first_exponential*(-ThermalA::erf((r_r*((2*L*n)+z-(v_arc*t_prime)))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(r_r,2)*((2*L*n)+z))+(12*ThermalA::intpow(k_p,2)*v_arc*alpha*(t-t_prime)*t_prime))/(r_r*erf_denom))))+
(second_exponential*(ThermalA::erf((r_r*((2*L*n)+z+(v_arc*t_prime)))/(erf_denom))
-ThermalA::erf(((ThermalA::intpow(r_r,2)*((2*L*n)+z))+(12*ThermalA::intpow(k_p,2)*v_arc*alpha*t_prime*(t_prime-t)))/(r_r*erf_denom)))));


double int_Z_22=numerator/denominator;


Z_22=Z_22+int_Z_22;

}
return Z_22;
}

double FSW_4Q_Heat_Source_y_column_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y,double d_g,double y_i){
double alpha=((k_conductivity)/(cp_steel*density));
double Y22=0;
double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime);
for(int n=-10;n<=10;n+=1){

double int_Y22=(sqrt(alpha)*(ThermalA::erf((d_g-(2*D*n)-y)/erf_denom)+ThermalA::erf((d_g+(2*D*n)+y)/erf_denom)+ThermalA::erf(((2*D*n)+y-y_i)/erf_denom)-ThermalA::erf(((2*D*n)+y+y_i)/erf_denom))*sqrt(t-t_prime))/(2*sqrt(alpha*(t-t_prime)));
Y22=Y22+int_Y22;
//cout<<Y22<<endl;

}
return Y22;
}
double FSW_4Q_Heat_Source_y_column_component_insulating_at_0_dirichlet_at_D(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y,double d_g,double y_i){
double alpha=((k_conductivity)/(cp_steel*density));
double Y22=0;
double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime);
for(int n=-10;n<=10;n+=1){

double int_Y22=(pow(-1.0,n))*(sqrt(alpha)*(ThermalA::erf((d_g-(2*D*n)-y)/erf_denom)+ThermalA::erf((d_g+(2*D*n)+y)/erf_denom)+ThermalA::erf(((2*D*n)+y-y_i)/erf_denom)-ThermalA::erf(((2*D*n)+y+y_i)/erf_denom))*sqrt(t-t_prime))/(2*sqrt(alpha*(t-t_prime)));
Y22=Y22+int_Y22;
//cout<<Y22<<endl;

}
return Y22;
}


////(pow(-1.0,n))*
double FSW_4Q_Heat_Source_y_doughnut_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g){
double alpha=((k_conductivity)/(cp_steel*density));
double Y22=0;

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha));

double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha));
for(int n=-10;n<=10;n+=1){

double portion_1=exp(-((3*ThermalA::intpow((-d_g+(2*D*n)+y),2))/(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha))))*(ThermalA::erf((b*(d_g-(2*D*n)-y))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(b,2)*((D*(-1+(2*n)))+y))+(12*(d_g-D)*alpha*(t-t_prime)))/(erf_denom*b)));
double portion_2=exp(-(3*ThermalA::intpow((d_g+(2*D*n)+y),2))/(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha)))*(-ThermalA::erf((b*(d_g+(2*D*n)+y))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(b,2)*(D+(2*D*n)+y))+(12*(D-d_g)*t*alpha)+(12*(d_g-D)*alpha*t_prime))/(erf_denom*b)));
double numerator=b*sqrt(alpha)*sqrt(t-t_prime)*(-portion_1+portion_2);
double int_Y22=numerator/denominator;
Y22=Y22+int_Y22;
}
return Y22;
}

double FSW_4Q_Heat_Source_y_doughnut_component_insulating_at_0_dirichlet_at_D(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g){
double alpha=((k_conductivity)/(cp_steel*density));
double Y22=0;

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha));

double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha));
for(int n=-10;n<=10;n+=1){

double portion_1=exp(-((3*ThermalA::intpow((-d_g+(2*D*n)+y),2))/(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha))))*(ThermalA::erf((b*(d_g-(2*D*n)-y))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(b,2)*((D*(-1+(2*n)))+y))+(12*(d_g-D)*alpha*(t-t_prime)))/(erf_denom*b)));
double portion_2=exp(-(3*ThermalA::intpow((d_g+(2*D*n)+y),2))/(ThermalA::intpow(b,2)+(12*t*alpha)-(12*t_prime*alpha)))*(-ThermalA::erf((b*(d_g+(2*D*n)+y))/(erf_denom))+ThermalA::erf(((ThermalA::intpow(b,2)*(D+(2*D*n)+y))+(12*(D-d_g)*t*alpha)+(12*(d_g-D)*alpha*t_prime))/(erf_denom*b)));
double numerator=(pow(-1.0,n))*b*sqrt(alpha)*sqrt(t-t_prime)*(-portion_1+portion_2);
double int_Y22=numerator/denominator;
Y22=Y22+int_Y22;
}
return Y22;
}




double part_past_integral_4Q_FSW_HEAT_SOURCE(double x,double y,double z,double a,double b,double r_r, double r_f,double t,double t_prime,
double k_conductivity,double density,double cp_steel,double travel_velocity, double B, double D, double L, double b_g, double d_g,
double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP){
//,double velocity_ramp_uprate,double velocity_ramp_downdate,double time_2, double time_3)

double k_p=0.5907132527;

///////////////////SINE
double power_ramping=ThermalA::power_ramping_function(t_prime,t_start_point_RAMP,t_end_point_RAMP,uprate_RAMP,downrate_RAMP);
//double power_ramping=0.5+(0.5*sin(2*pi*0.25*t_prime));
/*
if(t_prime>(L/v_arc)){
ramping=0;
}
else{
ramping=1;
}
*/

//          t_start_point_RAMP=pin down time

//      time_2 - time_1 is hold time        therefore ramping up from time_1 to time_2
double v_arc=travel_velocity;//velocity_ramping_function(maximum_travel_velocity,t_prime,t_start_point_RAMP,time_2,time_3,t_end_point_RAMP,velocity_ramp_uprate,velocity_ramp_downdate);


double X_P=ThermalA::FSW_4Q_Heat_Source_Component_x_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
double X_not_P=ThermalA::FSW_4Q_Heat_Source_Component_x_NOT_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
double Z_P_Front=ThermalA::FSW_4Q_Heat_Source_Component_z_P_Front_component_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r_f,v_arc);
double Z_not_P_Front=ThermalA::FSW_4Q_Heat_Source_Component_z_NOT_P_Front_component_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r_f,v_arc);
double Z_P_Rear=ThermalA::FSW_4Q_Heat_Source_Component_z_P_Rear_component_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r_r,v_arc);
double Z_not_P_Rear=ThermalA::FSW_4Q_Heat_Source_Component_z_NOT_P_Rear_component_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r_r,v_arc);

double Y_Doughnut=ThermalA::FSW_4Q_Heat_Source_y_doughnut_component(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);
double Y_Column=ThermalA::FSW_4Q_Heat_Source_y_column_component(k_conductivity,cp_steel,density,t,t_prime,D,y,d_g,y_i);



/*
double Y_Doughnut=((FSW_4Q_Heat_Source_y_doughnut_component(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g)+
FSW_4Q_Heat_Source_y_doughnut_component_insulating_at_0_dirichlet_at_D(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g))/2);
double Y_Column=((FSW_4Q_Heat_Source_y_column_component(k_conductivity,cp_steel,density,t,t_prime,D,y,d_g,y_i)+
FSW_4Q_Heat_Source_y_column_component_insulating_at_0_dirichlet_at_D(k_conductivity,cp_steel,density,t,t_prime,D,y,d_g,y_i))/2);
*/


double const_1=((b*r_f)/((r_f+r_r)*(b+(2*sqrt(3/pi)*(d_g-y_i)))));
double const_2=((b*r_r)/((r_f+r_r)*(b+(2*sqrt(3/pi)*(d_g-y_i)))));
double const_3=((6*r_f*(d_g-y_i))/((r_f+r_r)*((6*d_g)+(b*sqrt(3*pi))-(6*y_i))));
double const_4=((6*r_r*(d_g-y_i))/((r_f+r_r)*((6*d_g)+(b*sqrt(3*pi))-(6*y_i))));

double Quad_1_factor=const_1*((108*sqrt(3)*ThermalA::intpow(k_p,2))/(a*b*pi*sqrt(pi)*r_f*(ThermalA::intpow(pi,2)-9)));
double Quad_2_factor=const_2*((108*sqrt(3)*ThermalA::intpow(k_p,2))/(a*b*pi*sqrt(pi)*r_r*(ThermalA::intpow(pi,2)-9)));
double Quad_3_factor=const_3*((54*ThermalA::intpow(k_p,2))/(a*pi*r_f*(d_g-y_i)*(ThermalA::intpow(pi,2)-9)));
double Quad_4_factor=const_4*((54*ThermalA::intpow(k_p,2))/(a*pi*r_r*(d_g-y_i)*(ThermalA::intpow(pi,2)-9)));


double Quad_1_Greens=((X_P*Z_P_Front)-(X_not_P*Z_not_P_Front))*(Y_Doughnut);
double Quad_2_Greens=((X_P*Z_P_Rear)-(X_not_P*Z_not_P_Rear))*(Y_Doughnut);
double Quad_3_Greens=((X_P*Z_P_Front)-(X_not_P*Z_not_P_Front))*(Y_Column);;
double Quad_4_Greens=((X_P*Z_P_Rear)-(X_not_P*Z_not_P_Rear))*(Y_Column);


double part_past_integral_finite=(1/(density*cp_steel))*((Quad_1_factor*Quad_1_Greens)+(Quad_2_factor*Quad_2_Greens)+(Quad_3_factor*Quad_3_Greens)+(Quad_4_factor*Quad_4_Greens));

//((column_component_factor*column_component)+(doughnut_component_factor*doughnut_component));


return part_past_integral_finite*power_ramping;
}

double FSW_4Q_Heat_Source_Component_z_P_Front_component_DWELL_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_f,double l_g){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double Z_22=0;


double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));

double denominator=2*sqrt(alpha*(t-t_prime))*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));

for(int n=-10;n<=10;n+=1){

double first_exponential=-exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow(-l_g+(2*L*n)+z,2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));
double second_exponential=exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow(l_g+(2*L*n)+z,2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));

double numerator=pi*r_f*sqrt(alpha)*sqrt(t-t_prime)*((first_exponential*(ThermalA::erf((pi*r_f*(l_g-(2*L*n)-z))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2)*((L*(-1+(2*n)))+z))+(108*ThermalA::intpow(k_p,2)*(-L+l_g)*t*alpha)+(108*ThermalA::intpow(k_p,2)*(L-l_g)*alpha*t_prime))/(pi*r_f*erf_denom))))+
(second_exponential*(-ThermalA::erf((pi*r_f*(l_g+(2*L*n)+z))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(pi,2)*ThermalA::intpow(r_f,2)*(L+(2+L+n)+z))+(108*ThermalA::intpow(k_p,2)*(L-l_g)*t*alpha)+(108*ThermalA::intpow(k_p,2)*(-L+l_g)*t_prime*alpha))/(pi*r_f*erf_denom)))));


double int_Z_22=numerator/denominator;


Z_22=Z_22+int_Z_22;

}
return Z_22;
}
double FSW_4Q_Heat_Source_Component_z_NOT_P_Front_component_DWELL_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_f,double l_g){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double Z_22=0;


double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(r_f,2)+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(r_f,2)+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
for(int n=-10;n<=10;n+=1){

double first_exponential=-exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow(-l_g+(2*L*n)+z,2))/((ThermalA::intpow(r_f,2))+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));
double second_exponential=exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow(l_g+(2*L*n)+z,2))/((ThermalA::intpow(r_f,2))+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));

double numerator=r_f*sqrt(alpha)*sqrt(t-t_prime)*((first_exponential*(ThermalA::erf((r_f*(l_g-(2*L*n)-z))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(r_f,2)*((L*(-1+(2*n)))+z))+(12*ThermalA::intpow(k_p,2)*(-L+l_g)*t*alpha)+(12*ThermalA::intpow(k_p,2)*(L-l_g)*t_prime*alpha))/(r_f*erf_denom))))+
(second_exponential*(-ThermalA::erf((r_f*(l_g+(2*L*n)+z))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(r_f,2)*(L+(2*L*n)+z))+(12*ThermalA::intpow(k_p,2)*(L-l_g)*t*alpha)+(12*ThermalA::intpow(k_p,2)*(-L+l_g)*t_prime*alpha))/(r_f*erf_denom)))));


double int_Z_22=numerator/denominator;


Z_22=Z_22+int_Z_22;

}
return Z_22;
}
double FSW_4Q_Heat_Source_Component_z_P_Rear_component_DWELL_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_r,double l_g){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double Z_22=0;


double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
for(int n=-10;n<=10;n+=1){

double first_exponential=exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow(-l_g+(2*L*n)+z,2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));
double second_exponential=exp(-((27*ThermalA::intpow(k_p,2)*ThermalA::intpow(l_g+(2*L*n)+z,2))/((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2))+(108*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));

double numerator=pi*r_r*sqrt(alpha)*sqrt(t-t_prime)*((first_exponential*(ThermalA::erf((pi*r_r*(l_g-(2*L*n)-z))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2)*((2*L*n)+z))+(108*ThermalA::intpow(k_p,2)*l_g*alpha*(t-t_prime)))/(pi*r_r*erf_denom))))+
(second_exponential*(ThermalA::erf((pi*r_r*(l_g+(2*L*n)+z))/(erf_denom))
-ThermalA::erf(((ThermalA::intpow(pi,2)*ThermalA::intpow(r_r,2)*((2*L*n)+z))-(108*ThermalA::intpow(k_p,2)*l_g*alpha*t)+(108*ThermalA::intpow(k_p,2)*l_g*alpha*t_prime))/(pi*r_r*erf_denom)))));


double int_Z_22=numerator/denominator;


Z_22=Z_22+int_Z_22;

}
return Z_22;
}
double FSW_4Q_Heat_Source_Component_z_NOT_P_Rear_component_DWELL_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_r,double l_g){
double alpha=((k_conductivity)/(cp_steel*density));
double k_p=0.5907132527;

double Z_22=0;


double erf_denom=2*sqrt(alpha)*sqrt(t-t_prime)*sqrt(ThermalA::intpow(r_r,2)+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
double denominator=2*sqrt(alpha*(t-t_prime))*sqrt(ThermalA::intpow(r_r,2)+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)));
for(int n=-10;n<=10;n+=1){

double first_exponential=exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow(-l_g+(2*L*n)+z,2))/((ThermalA::intpow(r_r,2))+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));
double second_exponential=exp(-((3*ThermalA::intpow(k_p,2)*ThermalA::intpow(l_g+(2*L*n)+z,2))/((ThermalA::intpow(r_r,2))+(12*ThermalA::intpow(k_p,2)*alpha*(t-t_prime)))));

double numerator=r_r*sqrt(alpha)*sqrt(t-t_prime)*((first_exponential*(ThermalA::erf((r_r*(l_g-(2*L*n)-z))/(erf_denom))
+ThermalA::erf(((ThermalA::intpow(r_r,2)*2*L*n)+(ThermalA::intpow(r_r,2)*z)+(12*ThermalA::intpow(k_p,2)*l_g*alpha*(t-t_prime)))/(r_r*erf_denom))))+
(second_exponential*(ThermalA::erf((r_r*(l_g+(2*L*n)+z))/(erf_denom))
-ThermalA::erf(((ThermalA::intpow(r_r,2)*2*L*n)+(ThermalA::intpow(r_r,2)*z)-(12*ThermalA::intpow(k_p,2)*l_g*alpha*t)+(12*ThermalA::intpow(k_p,2)*l_g*alpha*t_prime))/(r_r*erf_denom)))));


double int_Z_22=numerator/denominator;


Z_22=Z_22+int_Z_22;

}
return Z_22;
}

double part_past_integral_4Q_FSW_HEAT_SOURCE_DWELL(double x,double y,double z,double a,double b,double r_r, double r_f,double t,double t_prime,
double k_conductivity,double density,double cp_steel, double B, double D, double L, double b_g, double d_g, double l_g,
double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP){


double k_p=0.5907132527;


double power_ramping=ThermalA::power_ramping_function(t_prime,t_start_point_RAMP,t_end_point_RAMP,uprate_RAMP,downrate_RAMP);

/*
if(t_prime>(L/v_arc)){
ramping=0;
}
else{
ramping=1;
}
*/

//          t_start_point_RAMP=pin down time

//      time_2 - time_1 is hold time        therefore ramping up from time_1 to time_2


double X_P=ThermalA::FSW_4Q_Heat_Source_Component_x_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
double X_not_P=ThermalA::FSW_4Q_Heat_Source_Component_x_NOT_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,B,x,a,b_g);
double Z_P_Front=ThermalA::FSW_4Q_Heat_Source_Component_z_P_Front_component_DWELL_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r_f,l_g);
double Z_not_P_Front=ThermalA::FSW_4Q_Heat_Source_Component_z_NOT_P_Front_component_DWELL_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r_f,l_g);
double Z_P_Rear=ThermalA::FSW_4Q_Heat_Source_Component_z_P_Rear_component_DWELL_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r_r,l_g);
double Z_not_P_Rear=ThermalA::FSW_4Q_Heat_Source_Component_z_NOT_P_Rear_component_DWELL_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r_r,l_g);
double Y_Doughnut=ThermalA::FSW_4Q_Heat_Source_y_doughnut_component(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g);
double Y_Column=ThermalA::FSW_4Q_Heat_Source_y_column_component(k_conductivity,cp_steel,density,t,t_prime,D,y,d_g,y_i);



double const_1=((b*r_f)/((r_f+r_r)*(b+(2*sqrt(3/pi)*(d_g-y_i)))));
double const_2=((b*r_r)/((r_f+r_r)*(b+(2*sqrt(3/pi)*(d_g-y_i)))));
double const_3=((6*r_f*(d_g-y_i))/((r_f+r_r)*((6*d_g)+(b*sqrt(3*pi))-(6*y_i))));
double const_4=((6*r_r*(d_g-y_i))/((r_f+r_r)*((6*d_g)+(b*sqrt(3*pi))-(6*y_i))));

double Quad_1_factor=const_1*((108*sqrt(3)*ThermalA::intpow(k_p,2))/(a*b*pi*sqrt(pi)*r_f*(ThermalA::intpow(pi,2)-9)));
double Quad_2_factor=const_2*((108*sqrt(3)*ThermalA::intpow(k_p,2))/(a*b*pi*sqrt(pi)*r_r*(ThermalA::intpow(pi,2)-9)));
double Quad_3_factor=const_3*((54*ThermalA::intpow(k_p,2))/(a*pi*r_f*(d_g-y_i)*(ThermalA::intpow(pi,2)-9)));
double Quad_4_factor=const_4*((54*ThermalA::intpow(k_p,2))/(a*pi*r_r*(d_g-y_i)*(ThermalA::intpow(pi,2)-9)));


double Quad_1_Greens=((X_P*Z_P_Front)-(X_not_P*Z_not_P_Front))*(Y_Doughnut);
double Quad_2_Greens=((X_P*Z_P_Rear)-(X_not_P*Z_not_P_Rear))*(Y_Doughnut);
double Quad_3_Greens=((X_P*Z_P_Front)-(X_not_P*Z_not_P_Front))*(Y_Column);;
double Quad_4_Greens=((X_P*Z_P_Rear)-(X_not_P*Z_not_P_Rear))*(Y_Column);


double part_past_integral_finite=(1/(density*cp_steel))*((Quad_1_factor*Quad_1_Greens)+(Quad_2_factor*Quad_2_Greens)+(Quad_3_factor*Quad_3_Greens)+(Quad_4_factor*Quad_4_Greens));

//((column_component_factor*column_component)+(doughnut_component_factor*doughnut_component));


return part_past_integral_finite*power_ramping;
}




double Temperature_Finite_4Q_FSW(double x,double y,double z,double a,double b,double r_r, double r_f,double t,double k_conductivity,
double density,double cp_steel,double travel_velocity, double T_0, double current, double voltage, double eff, double B, double D, double L,
double b_g, double d_g,double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP){
//,double velocity_ramp_uprate,double velocity_ramp_downdate,double time_2, double time_3){

if(t==0){
return T_0;
}

//if(t>600.0){
double integral1=0;
double integral2=0;
double integral3=0;
double integral4=0;
double integral5=0;
double integral6=0;
double integral7=0;
double integral8=0;
double integral9=0;
double integral10=0;
double integral11=0;
double integral12=0;
double integral13=0;
double integral14=0;
double integral15=0;
double integral33=0;

FINITE_plate_to_be_integrated_FSW_4Q W_piece1;//
W_piece1.set_values(x,y,z,a,b,r_r,r_f,t,k_conductivity,density,cp_steel,travel_velocity,B,D,L,b_g,d_g,y_i,t_start_point_RAMP,
                    t_end_point_RAMP,uprate_RAMP,downrate_RAMP); //,velocity_ramp_uprate,velocity_ramp_downdate,time_2,time_3

integral1=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, 0, t/256, 1e-50);
integral2=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, t/256, t/128, 1e-50);

integral3=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, t/128, (2*t)/128, 1e-50);
integral4=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (2*t)/128, (3*t)/128, 1e-50);
integral5=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (3*t)/128, (4*t)/128, 1e-50);

integral33=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (4*t)/128, (5*t)/128, 1e-50);

integral6=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (5*t)/128, (6*t)/128, 1e-50);
integral7=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (6*t)/128,(8*t)/128, 1e-50);
integral8=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (8*t)/128, (12*t)/128, 1e-50);
integral9=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (12*t)/128, (16*t)/128, 1e-50);
integral10=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (16*t)/128, (20*t)/128, 1e-50);
integral11=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (20*t)/128, (24*t)/128, 1e-50);
integral12=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (24*t)/128, (32*t)/128, 1e-50);
integral13=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (32*t)/128, (64*t)/128, 1e-50);
integral14=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (64*t)/128, (96*t)/128, 1e-50);
integral15=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q>::Integrate(W_piece1, (96*t)/128, t, 1e-50);

double Temperature_4Q_FSW=(T_0+((integral1+integral2+integral3+integral4+integral5+integral6+integral7
+integral8+integral9+integral10+integral11+integral12+integral13+integral14+integral15+integral33)*voltage*current*eff));

return Temperature_4Q_FSW;
//}


}

double Temperature_Finite_4Q_FSW_DWELL(double x,double y,double z,double a,double b,double r_r, double r_f,double t,double k_conductivity,
double density,double cp_steel, double T_0, double current, double voltage, double eff, double B, double D, double L,
double b_g, double d_g, double l_g,double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP){
//,double velocity_ramp_uprate,double velocity_ramp_downdate,double time_2, double time_3){
//
//
//FOR THE CASE WHERE FSW heat source velocity is zero but is l_g along z-axis
//
//
if(t==0){
return T_0;
}

//if(t>600.0){
double integral1=0;
double integral2=0;
double integral3=0;
double integral4=0;
double integral5=0;
double integral6=0;
double integral7=0;
double integral8=0;
double integral9=0;
double integral10=0;
double integral11=0;
double integral12=0;
double integral13=0;
double integral14=0;
double integral15=0;
double integral33=0;

FINITE_plate_to_be_integrated_FSW_4Q_DWELL W_piece1;//
W_piece1.set_values(x,y,z,a,b,r_r,r_f,t,k_conductivity,density,cp_steel,B,D,L,b_g,d_g,l_g,y_i,t_start_point_RAMP,
                    t_end_point_RAMP,uprate_RAMP,downrate_RAMP); //,velocity_ramp_uprate,velocity_ramp_downdate,time_2,time_3

integral1=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, 0, t/256, 1e-50);
integral2=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, t/256, t/128, 1e-50);

integral3=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, t/128, (2*t)/128, 1e-50);
integral4=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (2*t)/128, (3*t)/128, 1e-50);
integral5=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (3*t)/128, (4*t)/128, 1e-50);

integral33=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (4*t)/128, (5*t)/128, 1e-50);

integral6=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (5*t)/128, (6*t)/128, 1e-50);
integral7=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (6*t)/128,(8*t)/128, 1e-50);
integral8=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (8*t)/128, (12*t)/128, 1e-50);
integral9=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (12*t)/128, (16*t)/128, 1e-50);
integral10=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (16*t)/128, (20*t)/128, 1e-50);
integral11=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (20*t)/128, (24*t)/128, 1e-50);
integral12=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (24*t)/128, (32*t)/128, 1e-50);
integral13=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (32*t)/128, (64*t)/128, 1e-50);
integral14=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (64*t)/128, (96*t)/128, 1e-50);
integral15=DEIntegrator<FINITE_plate_to_be_integrated_FSW_4Q_DWELL>::Integrate(W_piece1, (96*t)/128, t, 1e-50);

double Temperature_4Q_FSW_DWELL=(T_0+((integral1+integral2+integral3+integral4+integral5+integral6+integral7
+integral8+integral9+integral10+integral11+integral12+integral13+integral14+integral15+integral33)*voltage*current*eff));

return Temperature_4Q_FSW_DWELL;
//}


}




double part_past_integral_FSW_HEAT_SOURCE(double x,double y,double z,double r,double b,double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc, double B, double D, double L, double b_g, double d_g,double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP){
//
//
//FOR THE CASE WHERE r_r==r_f
//
//

double k_p=0.5907132527;

double column_component=((ThermalA::FSW_Heat_Source_Component_x_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,B,x,r,b_g)*ThermalA::FSW_Heat_Source_Component_z_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r,v_arc))-
(ThermalA::FSW_Heat_Source_Component_x_NOT_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,B,x,r,b_g)*ThermalA::FSW_Heat_Source_Component_z_NOT_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r,v_arc)))*(ThermalA::FSW_Heat_Source_y_column_component(k_conductivity,cp_steel,density,t,t_prime,D,y,d_g,y_i));

double TC=1/(1+((27*b*pow(pi,3/2))/(54*pi*(d_g-y_i)*sqrt(3.0))));
double BC=1-TC;

double column_component_factor=((27*ThermalA::intpow(k_p,2))/(pi*(ThermalA::intpow(pi,2)-9)*ThermalA::intpow(r,2)*(d_g-y_i)))*(TC);

double doughnut_component=((ThermalA::FSW_Heat_Source_Component_x_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,B,x,r,b_g)*ThermalA::FSW_Heat_Source_Component_z_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r,v_arc))-
(ThermalA::FSW_Heat_Source_Component_x_NOT_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,B,x,r,b_g)*ThermalA::FSW_Heat_Source_Component_z_NOT_P_component_insulating(k_conductivity,cp_steel,density,t,t_prime,L,z,r,v_arc)))*(ThermalA::FSW_Heat_Source_y_doughnut_component(k_conductivity,cp_steel,density,t,t_prime,D,y,b,d_g));

double doughnut_component_factor=((54*sqrt(3.0)*ThermalA::intpow(k_p,2))/(b*pow(pi,3/2)*(ThermalA::intpow(pi,2)-9)*ThermalA::intpow(r,2)))*(BC);



double ramping=ThermalA::power_ramping_function(t_prime,t_start_point_RAMP,t_end_point_RAMP,uprate_RAMP,downrate_RAMP);

/*
if(t_prime>(L/v_arc)){
ramping=0;
}
else{
ramping=1;
}
*/

double part_past_integral_finite=(1/(density*cp_steel))*((column_component_factor*column_component)+(doughnut_component_factor*doughnut_component));


return part_past_integral_finite*ramping;

}


double Temp_Finite_FSW(double x,double y,double z,double r,double b,double t,double k_conductivity,double density,
double cp_steel,double v_arc, double T_0, double current, double voltage, double eff, double B, double D, double L,double b_g,
 double d_g,double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP){
//
//
//FOR THE CASE WHERE r_r==r_f
//
//
if(t==0){
return T_0;
}

//if(t>600.0){
double integral1=0;
double integral2=0;
double integral3=0;
double integral4=0;
double integral5=0;
double integral6=0;
double integral7=0;
double integral8=0;
double integral9=0;
double integral10=0;
double integral11=0;
double integral12=0;
double integral13=0;
double integral14=0;
double integral15=0;
double integral33=0;

FINITE_plate_to_be_integrated_FSW W_piece1;//
W_piece1.set_values(x,y,z,r,b,t,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i,t_start_point_RAMP,t_end_point_RAMP,uprate_RAMP,downrate_RAMP);

integral1=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, 0, t/256, 1e-50);
integral2=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, t/256, t/128, 1e-50);

integral3=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, t/128, (2*t)/128, 1e-50);
integral4=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (2*t)/128, (3*t)/128, 1e-50);
integral5=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (3*t)/128, (4*t)/128, 1e-50);

integral33=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (4*t)/128, (5*t)/128, 1e-50);

integral6=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (5*t)/128, (6*t)/128, 1e-50);
integral7=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (6*t)/128,(8*t)/128, 1e-50);
integral8=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (8*t)/128, (12*t)/128, 1e-50);
integral9=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (12*t)/128, (16*t)/128, 1e-50);
integral10=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (16*t)/128, (20*t)/128, 1e-50);
integral11=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (20*t)/128, (24*t)/128, 1e-50);
integral12=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (24*t)/128, (32*t)/128, 1e-50);
integral13=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (32*t)/128, (64*t)/128, 1e-50);
integral14=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (64*t)/128, (96*t)/128, 1e-50);
integral15=DEIntegrator<FINITE_plate_to_be_integrated_FSW>::Integrate(W_piece1, (96*t)/128, t, 1e-50);

double Temp_Finite_FSW=(T_0+((integral1+integral2+integral3+integral4+integral5+integral6+integral7
+integral8+integral9+integral10+integral11+integral12+integral13+integral14+integral15+integral33)*voltage*current*eff));

return Temp_Finite_FSW;
//}


}








}
