/*
 * fileio.cpp
 *
 *  Created on: Nov 17, 2020
 *      Author: dmarce1
 */

#include <astrotiger/fileio.hpp>
#include <astrotiger/options.hpp>
#include <astrotiger/cosmos.hpp>

// This header structure was copied from N-GenIC

struct io_header {
	std::uint32_t npart[6]; /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
	double mass[6]; /*!< mass[1] gives the particle mass */
	double time; /*!< time (=cosmological scale factor) of snapshot */
	double redshift; /*!< redshift of snapshot */
	std::int32_t flag_sfr; /*!< flags whether star formation is used (not available in L-Gadget2) */
	std::int32_t flag_feedback; /*!< flags whether feedback from star formation is included */
	std::uint32_t npartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
	 the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder. */
	std::int32_t flag_cooling; /*!< flags whether radiative cooling is included */
	std::int32_t num_files; /*!< determines the number of files that are used for a snapshot */
	double BoxSize; /*!< Simulation box size (in code units) */
	double Omega0; /*!< matter density */
	double OmegaLambda; /*!< vacuum energy density */
	double HubbleParam; /*!< little 'h' */
	std::int32_t flag_stellarage; /*!< flags whether the age of newly formed stars is recorded and saved */
	std::int32_t flag_metals; /*!< flags whether metal enrichment is included */
	std::int32_t hashtabsize; /*!< gives the size of the hashtable belonging to this snapshot file */
	char fill[84]; /*!< fills to 256 Bytes */
};

#define FREAD_ASSERT(b) { \
	if( (b) != 1) { \
		printf( "Unable to read file. source %s line %i", __FILE__, __LINE__); \
	} \
}

static std::vector<particle> parts;

const std::vector<particle>& fileio_get_particles() {
	return parts;
}

void fileio_init_read() {
	std::int32_t dummy;
	std::string filename = "ics";
	FILE *fp = fopen(filename.c_str(), "rb");
	if (!fp) {
		printf("Unable to load %s\n", filename.c_str());
		abort();
	}
	io_header header;
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	FREAD_ASSERT(fread(&header, sizeof(header), 1, fp));
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	const std::uint64_t total_parts = std::uint64_t(header.npartTotal[1]) + (std::uint64_t(header.npartTotal[2]) << std::uint64_t(32));
	opts.omega_m = header.Omega0;
	opts.omega_b = 0.17 * header.Omega0;
	opts.part_mass = header.mass[1];
	opts.m_tot = opts.part_mass * total_parts;
	const auto Gcgs = 6.672e-8;
	const auto Hcgs = 3.2407789e-18;
	const auto mtot = opts.part_mass * total_parts;
	opts.code_to_g = 1.98e44;
	opts.code_to_cm_per_s = 3e10;
	opts.code_to_cm = std::pow((8.0 * M_PI * Gcgs * mtot * opts.code_to_g) / (3.0 * opts.omega_m * Hcgs * Hcgs), 1.0 / 3.0);
	opts.code_to_s = opts.code_to_cm / opts.code_to_cm_per_s;
	opts.H0 = Hcgs * opts.code_to_s;
	opts.G = Gcgs / pow(opts.code_to_cm, 3) * opts.code_to_g * pow(opts.code_to_s, 2);
	opts.z = header.redshift;
	opts.tmax = cosmos_set_z(opts.z, opts.H0);
	printf("Reading %li particles\n", total_parts);
	printf("code_to_cm =    %e\n", opts.code_to_cm);
	printf("code_to_g  =    %e\n", opts.code_to_g);
	printf("code_to_s  =    %e\n", opts.code_to_s);
	printf("tmax  =        %e\n", opts.tmax);
	printf("G          =    %e\n", opts.G);
	printf("H0         =    %e\n", opts.H0);
	printf("Z =             %e\n", header.redshift);
	printf("particle mass = %e\n", header.mass[1]);
	printf("mtot =          %e\n", mtot);
	printf("Omega_m =       %e\n", header.Omega0);
	printf("Omega_lambda =  %e\n", header.OmegaLambda);
	printf("Hubble Param =  %e\n", header.HubbleParam);
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	for (int i = 0; i < header.npart[1]; i++) {
		float x, y, z;
		FREAD_ASSERT(fread(&x, sizeof(float), 1, fp));
		FREAD_ASSERT(fread(&y, sizeof(float), 1, fp));
		FREAD_ASSERT(fread(&z, sizeof(float), 1, fp));
		if (x > 1.0 || x < 0.0) {
			printf("Particle x out of range %e!\n", x);
			abort();
		}
		if (y > 1.0 || y < 0.0) {
			printf("Particle y out of range %e!\n", y);
			abort();
		}
		if (z > 1.0 || z < 0.0) {
			printf("Particle z out of range %e!\n", z);
			abort();
		}
		if (x == 1.0) {
			x = 0.0;
		}
		if (y == 1.0) {
			y = 0.0;
		}
		if (z == 1.0) {
			z = 0.0;
		}
		particle part;
		if ( NDIM > 0) {
			part.x[0] = x;
		}
		if ( NDIM > 1) {
			part.x[1] = y;
		}
		if (NDIM > 3) {
			part.x[2] = z;
		}
		part.rung = 0;
		part.m = opts.part_mass;
//		part.out = 0;
		parts.push_back(part);
	}
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	const auto c0 = std::pow(1.0 / (1.0 + header.redshift), 0.5);
	for (auto &part : parts) {
		float vx, vy, vz;
		FREAD_ASSERT(fread(&vx, sizeof(float), 1, fp));
		FREAD_ASSERT(fread(&vy, sizeof(float), 1, fp));
		FREAD_ASSERT(fread(&vz, sizeof(float), 1, fp));
		if ( NDIM > 0) {
			part.v[0] = vx * c0;
		}
		if ( NDIM > 1) {
			part.v[1] = vy * c0;
		}
		if ( NDIM > 2) {
			part.v[2] = vz * c0;
		}
	}
	FREAD_ASSERT(fread(&dummy, sizeof(dummy), 1, fp));
	fclose(fp);
}
