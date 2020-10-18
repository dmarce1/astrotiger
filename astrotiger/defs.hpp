#pragma once

#define NDIM 2

#if NDIM == 3
#define MULTI_FOR_BEGIN( i, box ) \
		{	\
			vect<index_type> i = box.min; \
			for( i[0] = box.min[0]; i[0] < box.max[0]; i[0]++ ) { \
				for( i[1] = box.min[1]; i[1] < box.max[1]; i[1]++ ) { \
					for( i[2] = box.min[2]; i[2] < box.max[2]; i[2]++ ) {
#define MULTI_FOR_END }}}}
#else
#if NDIM==2
#define MULTI_FOR_BEGIN( i, box ) \
		{ \
			vect<index_type> i = box.min; \
			for( i[0] = box.min[0]; i[0] < box.max[0]; i[0]++ ) { \
				for( i[1] = box.min[1]; i[1] < box.max[1]; i[1]++ ) {
#define MULTI_FOR_END }}}
#else
#define MULTI_FOR_BEGIN( i, box ) \
{ \
		vect<index_type> i = box.min; \
		for( i[0] = box.min[0]; i[0] < box.max[0]; i[0]++ ) {
#define MULTI_FOR_END }}
#endif
#endif


#define hydro_i 0
