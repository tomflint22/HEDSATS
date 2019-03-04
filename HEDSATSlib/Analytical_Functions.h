#ifndef ANALYTICAL_FUNCTIONS_H_INCLUDED
#define ANALYTICAL_FUNCTIONS_H_INCLUDED



//using namespace std;
namespace ThermalA {

void Welcome();
double erf(double x);
double inverf(double x);
double intpow(double base, int exponent);

//double part_past_integral_finite_with_ramping_and_convection_at_surface(double x,double y,double z,double a,double b,double cr,double cf,double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc, double B, double D, double L, double ts, double tf, double r_start, double r_fall,double B_Nuss);
double power_ramping_function(double t_prime,double ts, double tf, double r_start, double r_fall);
//double Temperature_finite(double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,double cp_steel,double v_arc, double T_0, double current, double voltage, double eff, double B, double D, double L);
//double Temp_Finite_R_nogroove_convec_at_surf(double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,double cp_steel,double v_arc, double T_0, double current, double voltage, double eff, double B, double D, double L,double ts,double tf,double r_start,double r_fall,double B_Nuss);
//double Temp_w_array(double x,double y,double z,double a,double b,double cr,double cf,double t,double k_of_t_array[],double density,double cp_of_t_array[],double v_arc, double T_0, double current, double voltage, double eff, double B, double D, double L,double ts,double tf,double r_start,double r_fall,double tstart,double tgap,int array_sizeovert);
//double part_past_integral_infinite(double x,double y,double z,double a,double b,double cr,double cf,double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc);
double Temperature_infinite(double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,double cp_steel,double v_arc, double T_0, double current, double voltage, double eff);
double Flux_Double_Ellipsoid(double x, double y, double z, double t, double a, double b, double cr, double cf, double eff, double voltage, double current, double v_arc);


//double Temperature_new(double x,double y,double z,double a,double b,double cr,double cf,double t,double k_of_t_array[],double density,double cp_of_t_array[],double v_arc, double B, double D, double L, double ts, double tf, double r_start, double r_fall,double tstart,double tgap,int array_sizeovert,double T_0,double lowlimitt);

//MODIFIED DEC HEAT SOURCE FOR EB WELDING//

//GOLDAK Parts - Different Boundary Conditions
double x_DE_part_HS_at_bg(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g);
double x_DE_part_HS_at_bg_insulating_at_0_dirichlet_at_B(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g);
double x_DE_part_HS_at_bg_dirichlet_at_0_dirichlet_at_B(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g);
double y_DE_part_HS_at_dg_integrated_dg_D(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g);
double z_DE_part_rear_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double v_arc, double cr);
double z_DE_part_front_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double v_arc, double cf);
double z_DE_part_rear_Z11(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double v_arc, double cr);
double z_DE_part_front_Z11(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double v_arc, double cf);
//Conical Parts
double y_EB_conical_part(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y,double d_g,double y_i);

//Combine//

double partpastintegral_DEC_EB_WELD_HEAT_SOURCE(double x,double y,double z,double a,double b,double cr,double cf,
double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc, double B,double D,
double L, double b_g, double d_g, double y_i,double t_start_ramp,double t_finish_ramp,double ramp_rate_start,
double ramp_rate_end,int xBC,int yBC,int zBC);

double part_past_integral_DE_HEAT_SOURCE(double x,double y,double z,double a,double b,double cr,double cf,
double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc, double B, double D,
double L, double b_g, double d_g,double t_start_ramp,double t_finish_ramp,
double ramp_rate_start,double ramp_rate_end,int xBC,int yBC,int zBC);

//double part_past_integral_DE_HEAT_SOURCE(double x,double y,double z,double a,double b,double cr,double cf,double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc, double B, double D, double L, double b_g, double d_g);
//CALCULATE
double Temp_Finite_EB(double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,double cp_steel,double v_arc, double T_0,
double current, double voltage, double eff, double B, double D, double L,double b_g, double d_g,double y_i,double t_start_ramp,double t_finish_ramp,double ramp_rate_start,double ramp_rate_end,
int xBC,int yBC,int zBC);

double Temp_Finite_DE(double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,double cp_steel,double v_arc, double T_0,
double current, double voltage, double eff, double B, double D, double L,double b_g, double d_g,double t_start_ramp,double t_finish_ramp,double ramp_rate_start,double ramp_rate_end,
int xBC,int yBC,int zBC);

double Temperature_from_combined_DECbeam_DE_heat_source(double x,double y,double z,double a_DE_in,double b_DE_in,double cr_DE_in,double cf_DE_in,double a_EB_in,double b_EB_in,
double cr_EB_in,double cf_EB_in,double t,double k_conductivity,double density,double cp_steel,double v_arc, double T_0, double current, double voltage, double eff, double B, double D,
 double L,double b_g, double d_g,double y_i,double portion_of_heat_in_DE,double t_start_ramp,double t_finish_ramp,double ramp_rate_start,double ramp_rate_end,
 int xBC,int yBC,int zBC);



////////////////CONVECTION AT BOTH SURFACES////////////////////
double y_component_convection_at_D_convection_at_0_hs_at_dg(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g);
double x_component_convection_at_B_convection_at_0_hs_at_bg(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g);
double z_rear_component_convection_at_L_convection_at_0_hs_at_vt(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double cr,double v_arc);
double z_front_component_convection_at_L_convection_at_0_hs_at_vt(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double cf,double v_arc);
double part_past_time_integral_convection_at_all_faces_offset_in_x_and_y(double x,double y,double z,double a,double b,double cr,double cf,double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc, double B, double D, double L, double ts, double tf, double r_start, double r_fall, double b_g, double d_g);
double Temperature_convection_on_all_surfaces(double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,double cp_steel,double v_arc, double T_0, double current, double voltage, double eff, double B, double D, double L,double ts,double tf,double r_start,double r_fall,double b_g, double d_g);
//////////////////////////////////////////////////////////////////////////////////



double FSW_Heat_Source_Component_x_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double r,double b_g);
double FSW_Heat_Source_Component_x_NOT_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double r,double b_g);
double FSW_Heat_Source_Component_z_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r,double v_arc);
double FSW_Heat_Source_Component_z_NOT_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r,double v_arc);
double FSW_Heat_Source_y_column_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y,double d_g,double y_i);
double FSW_Heat_Source_y_doughnut_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g);



double FSW_4Q_Heat_Source_Component_x_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g);
double FSW_4Q_Heat_Source_Component_x_NOT_P_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double B, double x, double a,double b_g);

double FSW_4Q_Heat_Source_Component_z_P_Front_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_f,double v_arc);
double FSW_4Q_Heat_Source_Component_z_NOT_P_Front_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_f,double v_arc);
double FSW_4Q_Heat_Source_Component_z_P_Rear_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_r,double v_arc);
double FSW_4Q_Heat_Source_Component_z_NOT_P_Rear_component_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_r,double v_arc);

double FSW_4Q_Heat_Source_y_column_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y,double d_g,double y_i);
double FSW_4Q_Heat_Source_y_doughnut_component(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g);

double FSW_4Q_Heat_Source_y_column_component_insulating_at_0_dirichlet_at_D(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y,double d_g,double y_i);
double FSW_4Q_Heat_Source_y_doughnut_component_insulating_at_0_dirichlet_at_D(double k_conductivity, double cp_steel, double density, double t, double t_prime, double D, double y, double b,double d_g);

double part_past_integral_4Q_FSW_HEAT_SOURCE(double x,double y,double z,double a,double b,double r_r, double r_f,double t,double t_prime,
double k_conductivity,double density,double cp_steel,double travel_velocity, double B, double D, double L, double b_g, double d_g,
double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP);
//double velocity_ramp_uprate,double velocity_ramp_downdate,double time_2, double time_3);

double FSW_4Q_Heat_Source_Component_z_P_Front_component_DWELL_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_f,double l_g);
double FSW_4Q_Heat_Source_Component_z_NOT_P_Front_component_DWELL_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_f,double l_g);
double FSW_4Q_Heat_Source_Component_z_P_Rear_component_DWELL_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_r,double l_g);
double FSW_4Q_Heat_Source_Component_z_NOT_P_Rear_component_DWELL_insulating(double k_conductivity, double cp_steel, double density, double t, double t_prime, double L, double z, double r_r,double l_g);


double part_past_integral_4Q_FSW_HEAT_SOURCE_DWELL(double x,double y,double z,double a,double b,double r_r, double r_f,double t,double t_prime,
double k_conductivity,double density,double cp_steel, double B, double D, double L, double b_g, double d_g, double l_g,
double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP);

double Temperature_Finite_4Q_FSW(double x,double y,double z,double a,double b,double r_r, double r_f,double t,double k_conductivity,
double density,double cp_steel,double travel_velocity, double T_0, double current, double voltage, double eff, double B, double D, double L,
double b_g, double d_g,double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP);
//,double velocity_ramp_uprate,double velocity_ramp_downdate,double time_2, double time_3);

double Temperature_Finite_4Q_FSW_DWELL(double x,double y,double z,double a,double b,double r_r, double r_f,double t,double k_conductivity,
double density,double cp_steel, double T_0, double current, double voltage, double eff, double B, double D, double L,
double b_g, double d_g, double l_g,double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP);


double part_past_integral_FSW_HEAT_SOURCE(double x,double y,double z,double r,double b,double t,double t_prime,double k_conductivity,double density,double cp_steel,double v_arc, double B, double D, double L, double b_g, double d_g,double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP);
double Temp_Finite_FSW(double x,double y,double z,double r,double b,double t,double k_conductivity,double density,double cp_steel,double v_arc, double T_0, double current, double voltage, double eff, double B, double D, double L,double b_g, double d_g,double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP);

/////////////////////////////////////Temperature dependent mechanical properties/////////////////////////////////////////////



class FINITE_plate_to_be_integrated_DE_Heat_source
{
	double x,y,z,a,b,cr,cf,t,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end;
	int xBC,yBC,zBC;

public:
	void set_values (double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,double cp_steel,double v_arc,double B,double D, double L,
	double b_g, double d_g,double t_start_ramp,double t_finish_ramp,double ramp_rate_start,double ramp_rate_end,int xBC,int yBC,int zBC);
	double operator()(double t_prime) const
	{
return part_past_integral_DE_HEAT_SOURCE(x,y,z,a,b,cr,cf,t,t_prime,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end,xBC,yBC,zBC);

}
};


class FINITE_plate_to_be_integrated_EB
{
	double x,y,z,a,b,cr,cf,t,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i,t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end;
	int xBC,yBC,zBC;
	//int n_ROOT;
public:
	void set_values (double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,double cp_steel,double v_arc,double B,double D, double L,
	double b_g, double d_g,double y_i,double t_start_ramp,double t_finish_ramp,double ramp_rate_start,double ramp_rate_end,int xBC,int yBC,int zBC);
	double operator()(double t_prime) const
	{
return partpastintegral_DEC_EB_WELD_HEAT_SOURCE(x,y,z,a,b,cr,cf,t,t_prime,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i,t_start_ramp,t_finish_ramp,ramp_rate_start,ramp_rate_end,xBC,yBC,zBC);

}
};




class FINITE_plate_to_be_integrated_FSW
{
	double x,y,z,r,b,t,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i,t_start_point_RAMP,t_end_point_RAMP,uprate_RAMP,downrate_RAMP;
	//int n_ROOT;
public:
	void set_values (double x,double y,double z,double r,double b,double t,double k_conductivity,double density,double cp_steel,double v_arc,double B,double D, double L,double b_g, double d_g,double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP);
	double operator()(double t_prime) const
	{
return part_past_integral_FSW_HEAT_SOURCE(x,y,z,r,b,t,t_prime,k_conductivity,density,cp_steel,v_arc,B,D,L,b_g,d_g,y_i,t_start_point_RAMP,t_end_point_RAMP,uprate_RAMP,downrate_RAMP);

}
};



class FINITE_plate_to_be_integrated_FSW_4Q
{
	double x,y,z,a,b,r_r,r_f,t,k_conductivity,density,cp_steel,travel_velocity,B,D,L,b_g,d_g,y_i,t_start_point_RAMP,t_end_point_RAMP,
	uprate_RAMP,downrate_RAMP;
	//,velocity_ramp_uprate,velocity_ramp_downdate,time_2, time_3;
	//int n_ROOT;

public:
	void set_values (double x,double y,double z,double a,double b,double r_r, double r_f,double t,double k_conductivity,
        double density,double cp_steel,double travel_velocity, double B, double D, double L, double b_g, double d_g,
        double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP);
        //double velocity_ramp_uprate,double velocity_ramp_downdate,double time_2, double time_3);

	double operator()(double t_prime) const
	{
return part_past_integral_4Q_FSW_HEAT_SOURCE(x,y,z,a,b,r_r,r_f,t,t_prime,
        k_conductivity, density, cp_steel,travel_velocity,B,D,L,b_g,d_g,
        y_i,t_start_point_RAMP,t_end_point_RAMP,uprate_RAMP,downrate_RAMP);
        //velocity_ramp_uprate,velocity_ramp_downdate,time_2,time_3);

}
};

class FINITE_plate_to_be_integrated_FSW_4Q_DWELL
{
	double x,y,z,a,b,r_r,r_f,t,k_conductivity,density,cp_steel,B,D,L,b_g,d_g,l_g,y_i,t_start_point_RAMP,t_end_point_RAMP,
	uprate_RAMP,downrate_RAMP;
	//,velocity_ramp_uprate,velocity_ramp_downdate,time_2, time_3;
	//int n_ROOT;

public:
	void set_values (double x,double y,double z,double a,double b,double r_r, double r_f,double t,double k_conductivity,
        double density,double cp_steel, double B, double D, double L, double b_g, double d_g, double l_g,
        double y_i,double t_start_point_RAMP,double t_end_point_RAMP,double uprate_RAMP, double downrate_RAMP);
        //double velocity_ramp_uprate,double velocity_ramp_downdate,double time_2, double time_3);

	double operator()(double t_prime) const
	{
return part_past_integral_4Q_FSW_HEAT_SOURCE_DWELL(x,y,z,a,b,r_r,r_f,t,t_prime,
        k_conductivity, density, cp_steel,B,D,L,b_g,d_g,l_g,
        y_i,t_start_point_RAMP,t_end_point_RAMP,uprate_RAMP,downrate_RAMP);
        //velocity_ramp_uprate,velocity_ramp_downdate,time_2,time_3);

}
};



class FINITE_plate_to_be_integrated_with_ramping_CONVECTION_ON_ALL_SURFACES
{
	double x,y,z,a,b,cr,cf,t,k_conductivity,density,cp_steel,v_arc,B,D,L,ts,tf,r_start,r_fall,b_g,d_g;
	//int n_ROOT;
public:
	void set_values (double x,double y,double z,double a,double b,double cr,double cf,double t,double k_conductivity,double density,double cp_steel,double v_arc,double B,double D, double L,double ts,double tf,double r_start,double r_fall,double b_g, double d_g);
	double operator()(double t_prime) const
	{
return part_past_time_integral_convection_at_all_faces_offset_in_x_and_y(x, y, z, a, b, cr, cf, t, t_prime, k_conductivity, density, cp_steel, v_arc,B,D,L,ts,tf,r_start,r_fall,b_g,d_g);

}
};





}

#endif // ANALYTICAL_FUNCTIONS_H_INCLUDED
