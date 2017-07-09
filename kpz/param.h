/***************************************************************************
*   Copyright 2011 - 2015 Jeffrey Kelling <j.kelling@hzdr.de>
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

#pragma once

// #define KMC_SWITCH_BENCHMARK
#define KPZ_ASYNC_MEASUREMENT
// available options: SLCG64 TinyMT32
#ifndef KPZ_SWITCH_RNG // allow Makefile to override this setting
// #define KPZ_SWITCH_RNG SLCG64
#define KPZ_SWITCH_RNG TinyMT32
#endif
/* To use PRAND set PRAND_PATH to aesb checkout root when calling make. */

#ifndef __OPENCL_VERSION__
/* Block layer DD, does only affect SchedulerDT for now.
 * Scheduler is fixed to KPZ_INNER_DT. */
#define KPZ_INNER_DD KPZ_INNER_DBSH

// the next to do not apply if KPZ_INNER_DC is enabled
// curcial to get valid statistics with DB at block level
// ~ x.5
// #define WAIT_IF_DEAD // switch no longer available
//#define OPTIMIZED_WAIT // no measureable effect

// do not switch off, the other case is not implmented
#define STORE_DISORDER_IN_REG

/*! Cache random word when update only requires single random bits for pq decision.
 *
 * Implemented in Kpz::SchedulerSCA<GpuRng>::caller_mcsPQSyncRng
 */
// #define USE_CACHED_RND_OPT
/*! Variant of USE_CACHED_RND_OPT, no effect USE_CACHED_RND_OPT is not in effect. */
// #define USE_DYNAMIC_RND_CACHE
#endif
