/***************************************************************************
*   Copyright 2011 - 2012 Jeffrey Kelling <j.kelling@hzdr.de>
*                  Helmholtz-Zentrum Dresden-Rossendorf
*                  Institute of Ion Beam Physics and Materials Research
*
*	This file is part of CudaKpz.
*
*   CudaKPZ is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   CudaKPZ is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with CudaKPZ.  If not, see <http://www.gnu.org/licenses/>.
***************************************************************************/

#ifndef KMC_Z_LIB_HELPER_H
#define KMC_Z_LIB_HELPER_H

#include <iostream>

/*! \brief Wrapper for zlib.
 * \ingroup DGio
 */
namespace zWrapper
{
	/*! Chunk size used with deflate, optimum value stated in zlib doc (256kB). */
	const int CHUNK = 1<<18;
	/*! Same value as Z_DEFAULT_COMPRESSION, but we don't want to include zlib.h here.
	 */
	const int DEFAULT_COMPRESSION = -1;

	/*! Write compressed data (using deflate) into stream.
	 * \param data pointer to the data to write
	 * \param bytes number of bytes to write
	 * \param o stream to write to
	 * \param level zlib compression level
	 */
	bool writeCompressed (const char* data, int bytes, std::ostream& o, int level = DEFAULT_COMPRESSION);
	/*! Read compressed data (using deflate) from stream.
	 * \param data pointer to where to pot the read data
	 * \param bytes number of bytes to retrieve (uncompressed)
	 * \param i stream to read from
	 */
	bool readCompressed (char* data, int bytes, std::istream& i);

};

#endif
