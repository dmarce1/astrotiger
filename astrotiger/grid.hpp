#pragma once

#include <astrotiger/multi_array.hpp>

class grid {

public:
	virtual void resize(double dx, multi_range box_) = 0;
};
