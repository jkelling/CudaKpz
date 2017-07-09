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

#ifndef FILEHELPER_H
#define FILEHELPER_H

#include <string>

/*! \{
 * \defgroup DGioffflags Filetypes
 * \ingroup DGio
 */
/*! \brief Filetype: Invalid file.
 * \ingroup DGioffflags
 */
const int FT_INVALID = 0;
/*! \brief Filetype: .xyz file.
 * \ingroup DGioffflags
 */
const int FT_XYZ = 1;
/*! \brief Filetype: .bit file.
 * \ingroup DGioffflags
 */
const int FT_BIT = 2;
//!\}

/*! Determine file type by analyzing the head.
 * \param filename filename
 * \return typeflag, one of FT_INVALID, FT_XYZ and FT_BIT
 * \ingroup DGio
 */
int getFiletype(const char* filename);
/*! Generate a filename that does not already exist by adding an index to the
 * name.
 * \param a filename, this argument if modified to return the result
 * \param suffix Suffix to be added after the filename and index. If 0
 * everything after the last occurrence of '.' in the name will be used as
 * suffix.
 * \ingroup DGio
 */
void collisionFreeFilename(std::string& a, const char* suffix = 0);

#endif
