
struct LongitudinalCommand
{

};

struct AutopilotParams
{
	double airspeed_pitch_kp;
	double airspeed_pitch_ki;
	double theta_max;
	double theta_min;
	double airspeed_throttle_ki;
	double airspeed_throttle_kp;
	double throttle_trim;
	double altitude_ki;
	double altitude_kp;
	double pitch_kp;
	double pitch_kd;
	double delta_e_max;
	double Ts;
};


int main(int argc, char** argv){

}

double pitch_hold(double theta_c, double theta, double q, AutopilotParams P)
{
	double error = theta_c - theta;
	double up = error * P.pitch_kp;
	double ud = P.pitch_kd * q;
	double delta_e_unsat = -up + ud;
	double delta_e = saturate(delta_e_unsat, P.delta_e_max, -P.delta_e_max);

	return delta_e;
}

double altitude_hold(double h_c, double h, AutopilotParams P)
{
	static double integrator = 0;
	static double error_d1 = 0;

	double error = h_c - h;

	double inc = (P.Ts / 2)*(error + error_d1);
	integrator += inc;

	error_d1 = error;

	double ui = P.altitude_ki * integrator;
	double up = P.altitude_kp * error;

	double theta_c_unsat = ui + up;
	double theta_c = saturate(theta_c_unsat, P.theta_max, P.theta_min);

	if (theta_c_unsat != theta_c)
		integrator -= inc;

	return theta_c;
}


// Returns Theta_c from Va_c and current Va
double airspeed_with_pitch_hold(double Va_c, double Va, AutopilotParams P)
{
	static double integrator = 0;
	static double error_d1 = 0;

	double error = Va_c - Va;

	double inc = (P.Ts / 2)*(error + error_d1);
	integrator += inc;

	//proportional term
	double up = P.airspeed_pitch_kp*error;

	//integral term
	double ui = P.airspeed_pitch_ki*integrator;

	double theta_c_unsat = up + ui;
	double theta_c = saturate(theta_c_unsat, P.theta_max, P.theta_min);

	error_d1 = error;

	//Integrator Anti-Windup
	if (theta_c_unsat != theta_c)
		integrator -= inc;

	return theta_c;
}

double airspeed_with_throttle_hold(double Va_c, double Va, AutopilotParams P)
{
	static double integrator = 0;
	static double error_d1 = 0;

	double error = Va_c - Va;

	double inc = (P.Ts / 2)*(error + error_d1);
	integrator = integrator + inc;
	error_d1 = error;

	double delta_t_star = P.throttle_trim; //Trimmed Throttle

	double ui = P.airspeed_throttle_ki * integrator;
	double up = P.airspeed_throttle_kp * error;

	double delta_t_unsat = delta_t_star + ui + up;
	double delta_t = saturate(delta_t_unsat, 1, 0);

	if (delta_t_unsat != delta_t)
		integrator -= inc;

	return delta_t;
}

double saturate(double value, double max, double min)
{
	if (value > max)
		return max;
	else if (value < min)
		return min;
	else
		return value;
}