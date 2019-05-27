#ifndef TransferOutput_hpp
#define TransferOutput_hpp

#include "AbstractOutput.hpp"

class TransferOutput : public AbstractOutput {
public:
	TransferOutput(int writePeriod) {
		this->writePeriod = writePeriod;
	};

	void writeParaview(Mesh *mesh, int step);
};

#endif
