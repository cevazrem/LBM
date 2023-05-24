#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

const double pi = 3.14159265;

//Global attrs
int n = 101;
int m = 101;
int time_steps;
double alpha = 1;
double len = pi;
double total_time = 0.1;
double dt, dx, dy;

//Struct for the func needed in LBM
struct all_data {
	double x, y;				// Coordinates
	double rho;					// Curent temperature
	double ex;					// Exact solution
	std::vector <double> f;		// Distribution functions
	std::vector <double> feq;	// Equilibrium distribution functions
	std::vector <double> q;		// Heat source distribution functions
};

//q(x, y, t) - Function of heat sources
double right_part(double x, double y, double t) {
	return sin(x) * sin(y) * exp(-t);
}

class lb_method {
private:
	std::vector <std::vector <all_data>> data;	// Array for LBM struct
	std::vector <double> weights;				// Weights array
	double dt2;									// Dimensionless time step
	double cur_time;							// Current solution time
	double tau, omega;							// To account for relaxation time
public:
	lb_method() {
		// De - dimensioning
		cur_time = 0;
		dt2 = 1;

		// Normalizing due to Chapman-Enskog
		tau = alpha * 3 * dt / dx / dx + 0.5;	// Relaxation time
		// Simplification for operations
		omega = dt2 / tau;

		// Set weights array
		weights.resize(9);
		weights[0] = static_cast<double>(4) / 9;

		for (int i = 1; i < 5; ++i) {
			weights[i] = static_cast<double>(1) / 9;
		}
		for (int i = 5; i < 9; ++i) {
			weights[i] = static_cast<double>(1) / 36;
		}

		// Set all arrays sizes
		data.resize(n);
		for (int i = 0; i < n; ++i) {
			data[i].resize(m);
			for (int j = 0; j < m; ++j) {
				data[i][j].f.resize(9);
				data[i][j].feq.resize(9);
				data[i][j].q.resize(9);
			}
		}

		// Set coordinates grid
		for (int i = 0; i < n; ++i) {
			data[i][0].x = 0;
			for (int j = 1; j < m; ++j) {
				data[i][j].x = data[i][j - 1].x + dx;
			}
		}

		for (int j = 0; j < m; ++j) {
			data[0][j].y = 0;
			for (int i = 1; i < n; ++i) {
				data[i][j].y = data[i - 1][j].y + dy;
			}
		}

		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				// Start temperature distribution
				data[i][j].rho = sin(data[i][j].x) * sin(data[i][j].y);
				
				// Start values of distribution functions
				for (int k = 0; k < 9; ++k) {
					data[i][j].f[k] = data[i][j].rho * weights[k];
				}
			}
		}
	}

	// Collision step
	void collision() {
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				// Calculating heat sources distribution functions
				for (int k = 0; k < 9; ++k)
					data[i][j].q[k] = weights[k] * right_part(data[i][j].x, data[i][j].y, cur_time) * (1.0 - 1.0 / (2.0 * tau));

				// Calculating current temperature at the point
				data[i][j].rho = 0;
				for (int k = 0; k < 9; ++k)
					data[i][j].rho += data[i][j].f[k];

				data[i][j].rho += dt / 2.0 * right_part(data[i][j].x, data[i][j].y, cur_time);

				for (int k = 0; k < 9; ++k) {
					// Calculating equilibrium distribution functions
					data[i][j].feq[k] = weights[k] * data[i][j].rho;
					// Calculating new distribution functions
					data[i][j].f[k] = (1 - omega) * data[i][j].f[k] + omega * data[i][j].feq[k] + dt * data[i][j].q[k];;
				}
			}
		}
	}

	// Streaming step
	void streaming() {
		for (int j = m - 1; j > 0; --j) {
			for (int i = 0; i < n; ++i) {
				data[i][j].f[2] = data[i][j - 1].f[2];
			}
		}

		for (int j = m - 1; j > 0; --j) {
			for (int i = 0; i < n - 1; ++i) {
				data[i][j].f[6] = data[i + 1][j - 1].f[6];
			}
		}

		for (int j = m - 1; j >= 0; --j) {
			for (int i = n - 1; i > 0; --i) {
				data[i][j].f[1] = data[i - 1][j].f[1];
			}
		}

		for (int j = m - 1; j > 0; --j) {
			for (int i = n - 1; i > 0; --i) {
				data[i][j].f[5] = data[i - 1][j - 1].f[5];
			}
		}

		for (int j = 0; j < m - 1; ++j) {
			for (int i = n - 1; i >= 0; --i) {
				data[i][j].f[4] = data[i][j + 1].f[4];
			}
		}

		for (int j = 0; j < m - 1; ++j) {
			for (int i = n - 1; i > 0; --i) {
				data[i][j].f[8] = data[i - 1][j + 1].f[8];
			}
		}

		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n - 1; ++i) {
				data[i][j].f[3] = data[i + 1][j].f[3];
			}
		}

		for (int j = 0; j < m - 1; ++j) {
			for (int i = 0; i < n - 1; ++i) {
				data[i][j].f[7] = data[i + 1][j + 1].f[7];
			}
		}
	}

	// Boundary conditions
	void bound_cond() {
		for (int j = 0; j < m; ++j) {
			data[0][j].f[3] = data[0][j].f[1];
			data[0][j].f[6] = data[0][j].f[8];
			data[0][j].f[7] = data[0][j].f[5];

			data[n - 1][j].f[1] = data[n - 1][j].f[3];
			data[n - 1][j].f[5] = data[n - 1][j].f[7];
			data[n - 1][j].f[8] = data[n - 1][j].f[6];
		}

		for (int i = 0; i < n; ++i) {
			data[i][0].f[4] = data[i][0].f[2];
			data[i][0].f[7] = data[i][0].f[5];
			data[i][0].f[8] = data[i][0].f[6];

			data[i][m - 1].f[2] = data[i][m - 1].f[4];
			data[i][m - 1].f[5] = data[i][m - 1].f[7];
			data[i][m - 1].f[6] = data[i][m - 1].f[8];
		}
	}

	// For result file generation
	void write_to_tecplot(int step, double time)
	{
		std::string filename = "out\\out_" + std::to_string(step) + ".dat";
		std::string filename_ex = "out\\out_" + std::to_string(step) + "_ex.dat";

		std::ofstream myStream;
		std::ofstream myStream_ex;
		int k1 = n;
		int k2 = m;

		myStream.fixed;
		myStream.precision(10);
		myStream_ex.fixed;
		myStream_ex.precision(10);

		myStream.open(filename, 'w');
		myStream << "TITLE=\"OUT\"" << std::endl;
		myStream << "VARIABLES=\"X\",\"Y\",\"T\"" << std::endl;
		myStream << "ZONE T=\"D2Q9\", I = " << k1 << ", J = " << k2 << ", ZONETYPE=\"ORDERED\", DATAPACKING=\"BLOCK\"" << std::endl;

		myStream << "STRANDID=1" << std::endl;
		myStream << "SOLUTIONTIME=" << time << std::endl << std::endl;

		myStream << "# X:" << std::endl;
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				myStream << data[i][j].x << std::endl;
			}
		}

		myStream << "# Y:" << std::endl;
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				myStream << data[i][j].y << std::endl;
			}
		}

		myStream << std::endl << "# T:" << std::endl;
		myStream << "# X:" << std::endl;
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				myStream << data[i][j].rho << std::endl;
			}
		}

		myStream.close();

		myStream_ex.open(filename_ex, 'w');
		myStream_ex << "TITLE=\"OUT_EX\"" << std::endl;
		myStream_ex << "VARIABLES=\"X\",\"Y\",\"T\"" << std::endl;
		myStream_ex << "ZONE T=\"ANALYTICAL\", I = " << k1 << ", J = " << k2 << ", ZONETYPE = \"ORDERED\", DATAPACKING=\"BLOCK\"" << std::endl;

		myStream_ex << "STRANDID=2" << std::endl;
		myStream_ex << "SOLUTIONTIME=" << time << std::endl << std::endl;

		myStream_ex << "# X:" << std::endl;
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				myStream_ex << data[i][j].x << std::endl;
			}
		}

		myStream_ex << "# Y:" << std::endl;
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				myStream_ex << data[i][j].y << std::endl;
			}
		}

		myStream_ex << std::endl << "# T:" << std::endl;
		myStream_ex << "# X:" << std::endl;
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				myStream_ex << data[i][j].ex << std::endl;
			}
		}

		myStream_ex.close();
	}

	// Error calculation
	double calc_c_norm() {
		double res = 0.0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				if (abs(data[i][j].rho - data[i][j].ex) > res)
					res = abs(data[i][j].rho - data[i][j].ex);

		return res;
	}

	// Error calculation
	double calc_l2_norm() {
		double res = 0.0;
		for (int i = 1; i < n - 1; ++i)
			for (int j = 1; j < m - 1; ++j)
				res += pow((data[i][j].rho - data[i][j].ex), 2);

		return sqrt(res * dx * dy);
	}

	// Main loop
	void calculate() {
		for (int kk = 0; kk < time_steps + 1; ++kk) {
			cur_time += dt;

			collision();
			streaming();
			bound_cond();

			// Calculating exact solution
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < n; ++i) {
					data[i][j].ex = sin(data[i][j].x) * sin(data[i][j].y) * exp(-cur_time);
				}
			}
		}

		// Result output 
		write_to_tecplot(cur_time / dt, cur_time);

		// Errors output
		std::cout << "c_norm=" << calc_c_norm() << std::endl;
		std::cout << "l2_norm=" << calc_l2_norm() << std::endl;
	}
};

int main(int argc, char** argv) {
	// For the accuracy of the output
	std::cout << std::fixed << std::setprecision(10);

	// Set coord/time step
	dx = len / (n - 1);
	dy = len / (m - 1);
	dt = dx * dx / 4 / alpha;

	// Total main loop steps needed
	time_steps = total_time / dt;

	// Initialize LBM
	lb_method solution = lb_method();

	// Colculate sulution
	solution.calculate();

	return 0;
}