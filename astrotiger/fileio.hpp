/*
 * fileio.hpp
 *
 *  Created on: Nov 17, 2020
 *      Author: dmarce1
 */

#ifndef ASTROTIGER_FILEIO_HPP_
#define ASTROTIGER_FILEIO_HPP_

#include <astrotiger/particles.hpp>
#include <vector>

void fileio_init_read();
const std::vector<particle>& fileio_get_particles();

#endif /* ASTROTIGER_FILEIO_HPP_ */
