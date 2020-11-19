#include <astrotiger/cosmos.hpp>
#include <astrotiger/options.hpp>

#include <cmath>

struct cosmos {
	double t, tau, a, adot;
	template<class A>
	void serialize(A &&arc, unsigned) {
		arc & t;
		arc & tau;
		arc & a;
		arc & adot;
	}
	cosmos();
	void advance_to_time(double, int = 10000);
	void advance_to_scalefactor(double);
};

static double H0;

cosmos::cosmos() {
	a = 1.0;
	t = 0.0;
	tau = 0.0;
	adot = H0;
}

double dadt(double adot) {
	return adot;
}

double dadotdt(double a, double adot) {
	return H0 * H0 * (-0.5 * opts.omega_m / std::pow(a,NDIM-1) + a * (1.0 - opts.omega_m));
}

double dtaudt(double a) {
	return 1.0 / a;
}

void cosmos::advance_to_time(double t2, int Nstep) {
	double t1 = t;
	double dt = (t2 - t1) / Nstep;

	for (int i = 0; i < Nstep; i++) {
		const auto da1 = dadt(adot) * dt;
		const auto dadot1 = dadotdt(a, adot) * dt;
		const auto da2 = dadt(adot + dadot1 * 0.5) * dt;
		const auto dadot2 = dadotdt(a + da1 * 0.5, adot + dadot1 * 0.5) * dt;
		const auto da3 = dadt(adot + dadot2 * 0.5) * dt;
		const auto dadot3 = dadotdt(a + da2 * 0.5, adot + dadot2 * 0.5) * dt;
		const auto da4 = dadt(adot + dadot3) * dt;
		const auto dadot4 = dadotdt(a + da3, adot + dadot3) * dt;
		const auto dtau1 = dtaudt(a) * dt;
		const auto dtau2 = dtaudt(a + da1 * 0.5) * dt;
		const auto dtau3 = dtaudt(a + da2 * 0.5) * dt;
		const auto dtau4 = dtaudt(a + da3) * dt;
		a += (da1 + 2.0 * da2 + 2.0 * da3 + da4) / 6.0;
		adot += (dadot1 + 2.0 * dadot2 + 2.0 * dadot3 + dadot4) / 6.0;
		tau += (dtau1 + 2.0 * dtau2 + 2.0 * dtau3 + dtau4) / 6.0;
	}
	t = t2;

}


static cosmos C;

void cosmos_init() {
	C = cosmos();
}

void cosmos_reset_time() {
	C.t = 0.0;
}

void cosmos::advance_to_scalefactor(double a2) {
	double a1 = a2;
	while (std::abs(a / a2 - 1.0) > 1.0e-6) {
		double dt = (a2 - a) / adot / 1000.0;
		advance_to_time(t + dt, 100);
	}

}

double cosmos_set_z(double z, double h0) {
	H0 = h0;
	cosmos_init();
	C.advance_to_scalefactor(1.0 / (z + 1.0));
	return -C.t;
}

double cosmos_a() {
	return C.a;
}

double cosmos_adot() {
	return C.adot;
}

void cosmos_advance(double t) {
	C.advance_to_time(t - opts.tmax);

}
