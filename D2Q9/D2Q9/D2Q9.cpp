#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

const double pi = 3.14159265;
const double eps = 0.000001;

//0.0073829229
//0.0042361723
//0.0023070146
//0.0012024834

//0.0022583871
//0.0005648835
//0.0001411882
//0.0000352982

int n = 101;
int m = 101;
double alpha = 1;
double len = pi;
double total_time = 0.1;
double dt, dx, dy;
int time_steps;

struct all_data {
	std::vector <double> f, feq, q;
	double rho, x, y, ex;
};


double right_part(double x, double y, double t) {
	return sin(x) * sin(y) * exp(-t);
}


class lb_method {
private:
	std::vector <std::vector <all_data>> data;
	std::vector <double> weights;
	double dt2, dx2, dy2, csq, omega, tau, alpha2, cur_time;
	double sigma_0;
public:
	lb_method() {
		//normalize
		alpha2 = 1;
		dt2 = 1;
		dx2 = 1;
		dy2 = 1;
		sigma_0 = 0.04;
		cur_time = 0;

		std::cout << "dt=" << dt << std::endl;
		std::cout << "dx=" << dx << std::endl;
		std::cout << "alpha=" << alpha << std::endl;
		tau = alpha * 3 * dt / dx / dx + 0.5;
		std::cout << "tau=" << tau << std::endl;
		omega = dt2 / tau;
		std::cout << "omega=" << omega << std::endl;

		//set weights
		weights.resize(9);
		weights[0] = static_cast<double>(4) / 9;

		for (int i = 1; i < 5; ++i) {
			weights[i] = static_cast<double>(1) / 9;
		}
		for (int i = 5; i < 9; ++i) {
			weights[i] = static_cast<double>(1) / 36;
		}

		//set size of arrs
		data.resize(n);
		for (int i = 0; i < n; ++i) {
			data[i].resize(m);
			for (int j = 0; j < m; ++j) {
				data[i][j].f.resize(9);
				data[i][j].feq.resize(9);
				data[i][j].q.resize(9);
			}
		}

		//Set arrs x/y
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

		//Start distr
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				data[i][j].rho = sin(data[i][j].x) * sin(data[i][j].y);
				//data[i][j].rho = data[i][j].x * (pi - data[i][j].x) * sin(data[i][j].y);
				//data[i][j].rho = exp(-((data[i][j].x - 0.5) * (data[i][j].x - 0.5) + (data[i][j].y - 0.5) * (data[i][j].y - 0.5)) / (2 * sigma_0 * sigma_0));
				for (int k = 0; k < 9; ++k) {
					data[i][j].f[k] = data[i][j].rho * weights[k];
				}
			}
		}
	}


	void collision() {
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				for (int k = 0; k < 9; ++k)
					data[i][j].q[k] = weights[k] * right_part(data[i][j].x, data[i][j].y, cur_time) * (1.0 - 1.0 / (2.0 * tau));

				data[i][j].rho = 0;
				for (int k = 0; k < 9; ++k)
					data[i][j].rho += data[i][j].f[k];
				data[i][j].rho += dt / 2.0 * right_part(data[i][j].x, data[i][j].y, cur_time);

				for (int k = 0; k < 9; ++k) {
					data[i][j].feq[k] = weights[k] * data[i][j].rho;
					data[i][j].f[k] = (1 - omega) * data[i][j].f[k] + omega * data[i][j].feq[k] + dt * data[i][j].q[k];
				}	
			}
		}
	}

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
		/* grad = 0
		for (int j = 0; j < m; ++j) {
			data[0][j].f[1] = data[1][j].f[1];
			data[0][j].f[5] = data[1][j].f[5];
			data[0][j].f[8] = data[1][j].f[8];

			data[n - 1][j].f[3] = data[n - 2][j].f[3];
			data[n - 1][j].f[6] = data[n - 2][j].f[6];
			data[n - 1][j].f[7] = data[n - 2][j].f[7];
		}

		for (int i = 0; i < n; ++i) {
			data[i][0].f[2] = data[i][1].f[2];
			data[i][0].f[5] = data[i][1].f[5];
			data[i][0].f[6] = data[i][1].f[6];

			data[i][m - 1].f[4] = data[i][m - 2].f[4];
			data[i][m - 1].f[7] = data[i][m - 2].f[7];
			data[i][m - 1].f[8] = data[i][m - 2].f[8];
		}
		*/
	}

	void Write_To_TecPlot_Cons(int step, double time)
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
		myStream << "ZONE T=\"D2Q9\", I = " << k1 <<  ", J = " << k2 << ", ZONETYPE=\"ORDERED\", DATAPACKING=\"BLOCK\"" << std::endl;

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

	double calc_l2_norm() {
		double res = 0.0;
		for (int i = 1; i < n - 1; ++i)
			for(int j = 1; j < m - 1; ++j)
				res += abs((data[i][j].rho - data[i][j].ex) * (data[i][j].rho - data[i][j].ex));

		return sqrt(res * dx * dy);
	}

	double calc_c_norm() {
		double res = 0.0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				if (abs(data[i][j].rho - data[i][j].ex) > res)
					res = abs(data[i][j].rho - data[i][j].ex);

		return res;
	}

	void calculate() {
		double sigma_D;
		
		for (int kk = 0; kk < time_steps + 1; ++kk) {
			cur_time += dt;

			collision();
			streaming();
			bound_cond();

			sigma_D = sqrt(2 * kk * dt);
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < n; ++i) {
					data[i][j].ex = sin(data[i][j].x) * sin(data[i][j].y) * exp(-cur_time);
					//data[i][j].ex = data[i][j].x * (pi - data[i][j].x) * sin(data[i][j].y) * exp(-cur_time);
					//data[i][j].ex = (sigma_0 * sigma_0) / (sigma_0 * sigma_0 + sigma_D * sigma_D) * exp(-((data[i][j].x - 0.5) * (data[i][j].x - 0.5) + (data[i][j].y - 0.5) * (data[i][j].y - 0.5)) / (2 * (sigma_0 * sigma_0 + sigma_D * sigma_D)));
				}
			}

			//if ((kk % 10) == 0)
				//std::cout << kk << std::endl;
				//Write_To_TecPlot_Cons(kk, kk * dt);
		}

		//Write_To_TecPlot_Cons(1, total_time);
		std::cout <<  "c_norm=" << calc_c_norm() << std::endl;
		std::cout << "l2_norm=" << calc_l2_norm() << std::endl;
	}
};

class lb_method_2DQ5 {
private:
	std::vector <std::vector <all_data>> data;
	std::vector <double> weights;
	double dt2, dx2, dy2, csq, omega, tau, alpha2, cur_time;
	double sigma_0;
public:
	lb_method_2DQ5() {
		//normalize
		cur_time = 0;

		alpha2 = 1;
		dt2 = 1;
		dx2 = 1;
		dy2 = 1;
		sigma_0 = 0.04;

		std::cout << "dt=" << dt << std::endl;
		std::cout << "dx=" << dx << std::endl;
		std::cout << "alpha=" << alpha << std::endl;
		tau = alpha * 2 * dt / dx / dx + 0.5;
		std::cout << "tau=" << tau << std::endl;
		omega = dt2 / tau;
		std::cout << "omega=" << omega << std::endl;

		//set weights
		weights.resize(5);
		weights[0] = 0;
		weights[1] = 0.25;
		weights[2] = 0.25;
		weights[3] = 0.25;
		weights[4] = 0.25;

		//set size of arrs
		data.resize(n);
		for (int i = 0; i < n; ++i) {
			data[i].resize(m);
			for (int j = 0; j < m; ++j) {
				data[i][j].f.resize(5);
				data[i][j].feq.resize(5);
				data[i][j].q.resize(5);
			}
		}

		//Set arrs x/y
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

		//Start distr
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				data[i][j].rho = sin(data[i][j].x) * sin(data[i][j].y);
				//data[i][j].rho = exp(-((data[i][j].x - 0.5) * (data[i][j].x - 0.5) + (data[i][j].y - 0.5) * (data[i][j].y - 0.5)) / (2 * sigma_0 * sigma_0));
				for (int k = 0; k < 5; ++k) {
					data[i][j].f[k] = data[i][j].rho * weights[k];
				}
			}
		}
	}


	void collision() {
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				for (int k = 0; k < 5; ++k)
					data[i][j].q[k] = weights[k] * right_part(data[i][j].x, data[i][j].y, cur_time) * (1.0 - 1.0 / (2.0 * tau));

				data[i][j].rho = 0;
				for (int k = 0; k < 5; ++k)
					data[i][j].rho += data[i][j].f[k];

				data[i][j].rho += dt / 2.0 * right_part(data[i][j].x, data[i][j].y, cur_time);

				for (int k = 0; k < 5; ++k) {
					data[i][j].feq[k] = weights[k] * data[i][j].rho;
					data[i][j].f[k] = (1 - omega) * data[i][j].f[k] + omega * data[i][j].feq[k] + dt * data[i][j].q[k];;
				}
			}
		}
	}

	void streaming() {
		for (int j = m - 1; j > 0; --j) {
			for (int i = 0; i < n; ++i) {
				data[i][j].f[2] = data[i][j - 1].f[2];
			}
		}

		for (int j = m - 1; j >= 0; --j) {
			for (int i = n - 1; i > 0; --i) {
				data[i][j].f[1] = data[i - 1][j].f[1];
			}
		}

		for (int j = 0; j < m - 1; ++j) {
			for (int i = n - 1; i >= 0; --i) {
				data[i][j].f[4] = data[i][j + 1].f[4];
			}
		}

		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n - 1; ++i) {
				data[i][j].f[3] = data[i + 1][j].f[3];
			}
		}
	}

	void bound_cond() {
		for (int j = 0; j < m; ++j) {
			data[0][j].f[3] = data[0][j].f[1];
			data[n - 1][j].f[1] = data[n - 1][j].f[3];
		}

		for (int i = 0; i < n; ++i) {
			data[i][0].f[4] = data[i][0].f[2];
			data[i][m - 1].f[2] = data[i][m - 1].f[4];
		}

	}

	void Write_To_TecPlot_Cons(int step, double time)
	{
		std::string filename = "out\\out_" + std::to_string(step) + "_Q4.dat";

		std::ofstream myStream;
		int k1 = n;
		int k2 = m;

		myStream.fixed;
		myStream.precision(10);

		myStream.open(filename, 'w');
		myStream << "TITLE=\"OUT\"" << std::endl;
		myStream << "VARIABLES=\"X\",\"Y\",\"T\"" << std::endl;
		myStream << "ZONE T=\"D2Q5\", I = " << k1 << ", J = " << k2 << ", ZONETYPE=\"ORDERED\", DATAPACKING=\"BLOCK\"" << std::endl;

		myStream << "STRANDID=3" << std::endl;
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
	}

	double calc_l2_norm() {
		double res = 0.0;
		for (int i = 1; i < n - 1; ++i)
			for (int j = 1; j < m - 1; ++j)
				res += abs((data[i][j].rho - data[i][j].ex) * (data[i][j].rho - data[i][j].ex));

		return sqrt(res * dx * dy);
	}

	double calc_c_norm() {
		double res = 0.0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				if (abs(data[i][j].rho - data[i][j].ex) > res)
					res = abs(data[i][j].rho - data[i][j].ex);

		return res;
	}

	void calculate() {
		double sigma_D;
		for (int kk = 0; kk < time_steps + 1; ++kk) {
			cur_time += dt;

			collision();
			streaming();
			bound_cond();

			sigma_D = sqrt(2 * kk * dt);
			for (int j = 0; j < m; ++j) {
				for (int i = 0; i < n; ++i) {
					data[i][j].ex = sin(data[i][j].x) * sin(data[i][j].y) * exp(-cur_time);
					//data[i][j].ex = data[i][j].x * (pi - data[i][j].x) * sin(data[i][j].y) * exp(-cur_time);
					//data[i][j].ex = (sigma_0 * sigma_0) / (sigma_0 * sigma_0 + sigma_D * sigma_D) * exp(-((data[i][j].x - 0.5) * (data[i][j].x - 0.5) + (data[i][j].y - 0.5) * (data[i][j].y - 0.5)) / (2 * (sigma_0 * sigma_0 + sigma_D * sigma_D)));
				}
			}

			//if ((kk % 10) == 0)
				//Write_To_TecPlot_Cons(kk, kk * dt);
		}

		std::cout << "c_norm=" << calc_c_norm() << std::endl;
		std::cout << "l2_norm=" << calc_l2_norm() << std::endl;
	}
};

int main(int argc, char** argv) {
	std::cout << std::fixed << std::setprecision(10);
	dx = len / (n - 1);
	dy = len / (m - 1);
	dt = dx * dx / 4 / alpha;
	time_steps = total_time / dt;
	std::cout << "time_steps=" << time_steps << std::endl;

	lb_method solution_lb = lb_method();
	solution_lb.calculate();

	lb_method_2DQ5 solution = lb_method_2DQ5();
	solution.calculate();

	return 0;
}