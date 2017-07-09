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

#include "zlibWrapper.h"

#include <zlib.h>
#include <cstdlib>

bool zWrapper::writeCompressed (const char* data, int bytes, std::ostream& o, int level)
{
	z_stream strm;
	char out[CHUNK];

	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.next_in = (Bytef*)data;
	strm.avail_in = bytes;
	if(deflateInit(&strm, level) != Z_OK)
	{
		std::cerr << "Deflate init failed.\n";
		return false;
	}

	do {
		strm.next_out = (Bytef*)out;
		strm.avail_out = CHUNK;

		if(deflate(&strm, Z_FINISH) == Z_STREAM_ERROR)
		{
			std::cerr << "Fatal runtime error in inflate: Z_STREAM_ERROR.\n";
			exit(1);
		}

		const int have = CHUNK - strm.avail_out;
		o.write(out, have);
		
	} while(strm.avail_out == 0);

	deflateEnd(&strm);
	return true;
}

bool zWrapper::readCompressed (char* data, int bytes, std::istream& i)
{
	z_stream strm;
	char in[CHUNK];

	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.next_in = Z_NULL;
	strm.avail_in = 0;
	strm.next_out = (Bytef*)data;
	strm.avail_out = bytes;
	if(inflateInit(&strm) != Z_OK)
	{
		std::cerr << "Inflate init failed.\n";
		return false;
	}

	int ret;
	do {
		i.read(in, CHUNK);
		strm.avail_in = i.gcount();
		if(strm.avail_in == 0)
			break;
		strm.next_in = (Bytef*)in;

		int lastOut;
		do {
			lastOut = strm.avail_out;

			ret = inflate(&strm, Z_NO_FLUSH);
			switch(ret)
			{
				case Z_STREAM_ERROR:
					std::cerr << "Fatal runtime error in inflate: Z_STREAM_ERROR.\n";
					exit(1);
				case Z_NEED_DICT:
					std::cerr << "Runtime error in inflate: Z_NEED_DICT.\n";
					inflateEnd(&strm);
					return false;
				case Z_DATA_ERROR:
					std::cerr << "Runtime error in inflate: Z_DATA_ERROR.\n";
					inflateEnd(&strm);
					return false;
				case Z_MEM_ERROR:
					std::cerr << "Runtime error in inflate: Z_MEM_ERROR.\n";
					inflateEnd(&strm);
					return false;
			}

		} while(lastOut != strm.avail_out);

	} while(ret != Z_STREAM_END);

	inflateEnd(&strm);
	return ret == Z_STREAM_END;
}
