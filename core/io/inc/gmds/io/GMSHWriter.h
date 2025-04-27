#pragma once

#include <gmds/io/IWriter.h>
#include <gmds/utils/CommonTypes.h>
#include "GMDSIo_export.h"

#include <map>
#include <sstream>

namespace gmds {

    class GMDSIo_API GMSHWriter : public IWriter {
    public:

        /** @brief Constructor
         *
         * @param AMeshService an implementation of an io service to write data
         * 					   into a mesh
         */
        explicit GMSHWriter(IMeshIOService *AMeshService);

	private:

        void initialize(const std::string &AFileName) override;
        void writeNodes() override;
        void writeFaces() override;
        void finalize() override;
    };
}