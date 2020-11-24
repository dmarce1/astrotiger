#pragma once

#define NDIM 2



#define SILO_CHECK(b) \
	if( b != 0 ) { \
		printf( "SILO call failed in %s on line %i\n", __FILE__, __LINE__); \
	}

#define hydro_i 0
