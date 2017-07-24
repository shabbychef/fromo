
/*
 
  This file is part of fromo.
  
  fromo is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  fromo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.
  
  You should have received a copy of the GNU Lesser General Public License
  along with fromo.  If not, see <http://www.gnu.org/licenses/>.

	Welford's object.

  Created: 2017.07.23
  Copyright: Steven E. Pav, 2017
  Author: Steven E. Pav <steven@gilgamath.com>
  Comments: Steven E. Pav
*/

#include "common.h"

#undef HAS_WTS
// na rm block//FOLDUP
#undef NA_RM
// check weight block//FOLDUP
#undef CHECK_WT
// normalize block//FOLDUP
#undef NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#define NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#undef NORMALIZE
//UNFOLD
#define CHECK_WT
// normalize block//FOLDUP
#undef NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#define NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#undef NORMALIZE
//UNFOLD
#undef CHECK_WT
//UNFOLD
#define NA_RM
// check weight block//FOLDUP
#undef CHECK_WT
// normalize block//FOLDUP
#undef NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#define NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#undef NORMALIZE
//UNFOLD
#define CHECK_WT
// normalize block//FOLDUP
#undef NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#define NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#undef NORMALIZE
//UNFOLD
#undef CHECK_WT
//UNFOLD
#undef NA_RM
//UNFOLD
#define HAS_WTS
// na rm block//FOLDUP
#undef NA_RM
// check weight block//FOLDUP
#undef CHECK_WT
// normalize block//FOLDUP
#undef NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#define NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#undef NORMALIZE
//UNFOLD
#define CHECK_WT
// normalize block//FOLDUP
#undef NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#define NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#undef NORMALIZE
//UNFOLD
#undef CHECK_WT
//UNFOLD
#define NA_RM
// check weight block//FOLDUP
#undef CHECK_WT
// normalize block//FOLDUP
#undef NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#define NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#undef NORMALIZE
//UNFOLD
#define CHECK_WT
// normalize block//FOLDUP
#undef NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#define NORMALIZE

#define MAX_ORDER ORDER_ONE
#include "welfords.hpp"
#define MAX_ORDER ORDER_TWO
#include "welfords.hpp"
#define MAX_ORDER ORDER_BEYOND
#include "welfords.hpp"

#undef NORMALIZE
//UNFOLD
#undef CHECK_WT
//UNFOLD
#undef NA_RM
//UNFOLD
#undef HAS_WTS

//for vim modeline: (do not edit)
// vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
