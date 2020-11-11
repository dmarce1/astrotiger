#include <astrotiger/polytrope.hpp>
#include <astrotiger/options.hpp>

#include <cmath>

double lane_emden(double r, double dr0, double n) {

	const auto N = int(r / dr0) + 1;
	const auto dr = r / N;
	double theta = 1.0;
	double theta_dot = 0.0;

	const auto dtheta_dr = [](double theta, double theta_dot, double r) {
		return theta_dot;
	};

	const auto dtheta_dot_dr = [n](double theta, double theta_dot, double r) {
		double term1 = (r == 0.0 ? 0.0 : (NDIM - 1) * theta_dot / r);
		return -(term1 + std::pow(theta, n));
	};

	for (int i = 0; i < N; i++) {
		double r = i * dr;
		const double k1 = dtheta_dr(theta, theta_dot, r);
		const double l1 = dtheta_dot_dr(theta, theta_dot, r);
		if (theta + k1 * dr < 0.0) {
			theta = 0.0;
			break;
		}
		const double k2 = dtheta_dr(theta + k1 * dr, theta_dot + l1 * dr, r + dr);
		const double l2 = dtheta_dot_dr(theta + k1 * dr, theta_dot + l1 * dr, r + dr);
		theta += 0.5 * (k1 + k1) * dr;
		theta_dot += 0.5 * (l2 + l2) * dr;
		if (theta < 0.0) {
			theta = 0.0;
			break;
		}
//		printf("%e %e %e\n", r, theta, theta_dot);
	}
	return theta;

}
