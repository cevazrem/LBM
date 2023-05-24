#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

const double pi = 3.14159265;

//Global attrs
int m = 101;			
double alpha = 1;
double left_temp = 1;
double right_temp = 0;
double len = 1;
double total_time = 0.2;
double dt, dx;
int time_steps;

//Struct for the func needed in LBM
struct all_data {
	double x;					// Coordinate
	double f0, f1, f2;			// Distribution functions
	double q0, q1, q2;			// Heat source distribution functions
	double feq0, feq1, feq2;	// Equilibrium distribution functions
	double rho;					// Curent temperature
	double ex;					// Exact solution
};

class lb_method {
private:
	std::vector <all_data> data;	// Array for LBM struct
	double dt2;						// Dimensionless time step
	double cur_time;				// Current solution time
	double tau, omega;				// To account for relaxation time
public:
	lb_method() {
		data.resize(m);

		// De - dimensioning
		cur_time = 0;
		dt2 = 1;

		// Set coordinates grid
		data[0].x = 0;
		for (int i = 1; i < m; ++i) {
			data[i].x = data[i - 1].x + dx;
		}

		// Normalizing due to Chapman-Enskog
		tau = alpha * 3 * dt / dx / dx + 0.5;	// Relaxation time
		omega = dt2 / tau;						// Simplification for operations

		for (int i = 0; i < m; ++i) {
			// Start temperature distribution
			data[i].rho = -data[i].x + 1;//data[i].x * (1 - data[i].x) * cosh(data[i].x);

			// Start values of distribution functions
			data[i].f0 = 2.0 / 3.0 * data[i].rho;
			data[i].f1 = 1.0 / 6.0 * data[i].rho;
			data[i].f2 = 1.0 / 6.0 * data[i].rho;
		}
	}

	//q(x, t) - Function of heat sources
	double right_part(double x) {
		return 0;//2 * exp(-cur_time) * (cosh(x) * (x * x - x + 1) - sinh(x)*(1 - 2 * x));
	}

	// Collision step
	void collision() {
		for (int i = 0; i < m; ++i) {
			// Calculating heat sources distribution functions

			data[i].q0 = 2.0 / 3.0 * right_part(data[i].x) * (1.0 - 1.0 / (2.0 * tau));
			data[i].q1 = 1.0 / 6.0 * right_part(data[i].x) * (1.0 - 1.0 / (2.0 * tau));
			data[i].q2 = 1.0 / 6.0 * right_part(data[i].x) * (1.0 - 1.0 / (2.0 * tau));

			// Calculating current temperature at the point
			data[i].rho = data[i].f0 + data[i].f1 + data[i].f2 + dt / 2.0 * right_part(data[i].x);

			// Calculating equilibrium distribution functions
			data[i].feq0 = 2.0 / 3.0 * data[i].rho;
			data[i].feq1 = 1.0 / 6.0 * data[i].rho;
			data[i].feq2 = 1.0 / 6.0 * data[i].rho;

			// Calculating new distribution functions
			data[i].f0 = (1 - omega) * data[i].f0 + omega * data[i].feq0 + dt * data[i].q0;
			data[i].f1 = (1 - omega) * data[i].f1 + omega * data[i].feq1 + dt * data[i].q1;
			data[i].f2 = (1 - omega) * data[i].f2 + omega * data[i].feq2 + dt * data[i].q2;

		}
	}

	//Streaming step
	void streaming() {
		for (int i = 0; i < m - 1; ++i) {
			// Transfer to the right
			data[m - i - 1].f1 = data[m - i - 2].f1;

			// Transfer to the left
			data[i].f2 = data[i + 1].f2;
		}
	}

	//Boundary conditions
	void bound_cond() {
		// Left temperature
		data[0].f1 = left_temp - data[0].f2 - data[0].f0;

		// Right temperature
		data[m - 1].f2 = right_temp - data[m - 1].f1 - data[m - 1].f0;
	}

	// For result file generation
	void write_to_tecplot(int step, double time)
	{
		std::string filename = "out\\out_" + std::to_string(step) + ".dat";
		std::string filename_ex = "out\\out_" + std::to_string(step) + "_ex.dat";

		std::ofstream myStream;
		std::ofstream myStream_ex;
		int k1 = m;

		myStream.fixed;
		myStream.precision(10);
		myStream_ex.fixed;
		myStream_ex.precision(10);

		myStream.open(filename, 'w');
		myStream << "TITLE=\"OUT\"" << std::endl;
		myStream << "VARIABLES=\"X\",\"T\"" << std::endl;
		myStream << "ZONE T=\"D1Q3\", I = " << k1 << ", ZONETYPE=\"ORDERED\", DATAPACKING=\"BLOCK\"" << std::endl;
		myStream << "STRANDID=1" << std::endl;
		myStream << "SOLUTIONTIME=" << time << std::endl << std::endl;

		myStream << "# X:" << std::endl;
		for (int i = 0; i < m; ++i) {
			myStream << data[i].x << std::endl;
		}

		myStream << std::endl << "# T:" << std::endl;
		for (int i = 0; i < m; ++i) {
			myStream << data[i].rho << std::endl;
		}

		myStream.close();

		myStream_ex.open(filename_ex, 'w');
		myStream_ex << "TITLE=\"OUT_EX\"" << std::endl;
		myStream_ex << "VARIABLES=\"X\",\"T\"" << std::endl;
		myStream_ex << "ZONE T=\"Analytical\", I = " << k1 << ", ZONETYPE=\"ORDERED\", DATAPACKING=\"BLOCK\"" << std::endl;
		myStream_ex << "STRANDID=2" << std::endl;
		myStream_ex << "SOLUTIONTIME=" << time << std::endl << std::endl;

		myStream_ex << "# X:" << std::endl;
		for (int i = 0; i < m; ++i) {
			myStream_ex << data[i].x << std::endl;
		}

		myStream_ex << std::endl << "# T:" << std::endl;
		for (int i = 0; i < m; ++i) {
			myStream_ex << data[i].ex << std::endl;
		}

		myStream_ex.close();
	}

	// Error calculation
	double calc_c_norm() {
		double res = 0.0;
		for (int i = 0; i < m; ++i) {
			if (abs(data[i].rho - data[i].ex) > res)
				res = abs(data[i].rho - data[i].ex);
		}

		return res;
	}

	// Error calculation
	double calc_l2_norm() {
		double res = 0.0;
		for (int i = 1; i < m - 1; ++i)
			res += pow((data[i].rho - data[i].ex), 2);

		return sqrt(res * dx);
	}

	// Main loop
	void calculate() {
		for (int i = 0; i < time_steps; ++i) {
			cur_time += dt;
			collision();
			streaming();
			bound_cond();

			// Calculating exact solution
			for (int j = 0; j < m; ++j) {
				data[j].ex = -data[j].x + 1;//data[j].x * (1 - data[j].x) * cosh(data[j].x) * exp(-cur_time);
			}
		}

		// Result output 
		write_to_tecplot(cur_time/dt, cur_time);

		// Errors output
		std::cout << "C_norm = " << calc_c_norm() << std::endl;
		std::cout << "L2_norm = " << calc_l2_norm() << std::endl;
	}
};

int main(int argc, char** argv) {
	std::cout << std::fixed << std::setprecision(10); // For the accuracy of the output

	// Set coord/time step
	dx = len / (m - 1);
	dt = dx * dx / 4 / alpha;

	// Total main loop steps needed
	time_steps = total_time / dt;

	// Initialize LBM
	lb_method solution_lb = lb_method();

	// Colculate sulution
	solution_lb.calculate();

	return 0;
}