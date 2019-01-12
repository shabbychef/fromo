
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

	Kahan's summation object.

  Created: 2017.07.23
  Copyright: Steven E. Pav, 2017
  Author: Steven E. Pav <steven@gilgamath.com>
  Comments: Steven E. Pav
*/

#ifndef __DEF_KAHAN__
#define __DEF_KAHAN__

#include "common.h"

#define KAHAN_ADD(_sumx_,_err_,_newx_,_nxtv_,_tmpv_) \
    _tmpv_    = (_newx_) - (_err_);                    \
    _nxtv_    = (_sumx_) + _tmpv_;                     \
    (_err_)   = (_nxtv_ - (_sumx_)) - _tmpv_;          \
    (_sumx_)  = _nxtv_;

#define KAHAN_SUB(_sumx_,_err_,_newx_,_nxtv_,_tmpv_) \
    _tmpv_    = -(_newx_) - (_err_);                   \
    _nxtv_    = (_sumx_) + _tmpv_;                     \
    (_err_)   = (_nxtv_ - (_sumx_)) - _tmpv_;          \
    (_sumx_)  = _nxtv_;

// Kahan compensated summation object.//FOLDUP
template<class T>
class Kahan {
    public:
        T m_val;
    private:
        T m_errs;
    public:
        inline Kahan() : m_val(0), m_errs(0) {}

        inline Kahan(const T &t): m_val(t), m_errs(0) {}

        inline T as() const { return m_val; }

        inline Kahan& operator += (const T& rhs) {
            return add(rhs);
        }
        inline Kahan& operator -= (const T& rhs) {
            return add(-rhs);
        }
        inline Kahan& operator += (const Kahan<T>& rhs) {
            return join(rhs);
        }
        inline Kahan& operator -= (const Kahan<T>& rhs) {
            return unjoin(rhs);
        }
        inline Kahan& operator= (T rhs) {
            m_val = rhs;
            m_errs = 0;
            return *this;
        }
        // these are prefix;
        inline Kahan& operator++ () {
            return add(T(1));
        }
        inline Kahan& operator-- () {
            return add(T(-1));
        }
    private:
        inline Kahan& add(const T& rhs) {
            T tmpv, nxtv;
            KAHAN_ADD(m_val,m_errs,rhs,nxtv,tmpv)
            return *this;
        }
        inline Kahan& join(const Kahan<T>& rhs) {
            T tmpv, nxtv;
            tmpv = rhs.m_val - m_errs - rhs.m_errs;
            nxtv = m_val + tmpv;
            m_errs = (nxtv - m_val) - tmpv;
            m_val = nxtv;
            return *this;
        }
        inline Kahan& unjoin(const Kahan<T>& rhs) {
            T tmpv, nxtv;
            tmpv = -rhs.m_val - m_errs + rhs.m_errs;
            nxtv = m_val + tmpv;
            m_errs = (nxtv - m_val) - tmpv;
            m_val = nxtv;
            return *this;
        }
};
// specialization to int which do not require special treatment
template<>
class Kahan<int> {
    public:
        int m_val;
    public:
        inline Kahan() : m_val(0) {}

        inline Kahan(const int &t): m_val(t) {}

        inline int as() const { return m_val; }

        inline Kahan& operator += (const int& rhs) {
            return add(rhs);
        }
        inline Kahan& operator -= (const int& rhs) {
            return add(-rhs);
        }
        inline Kahan& operator += (const Kahan<int>& rhs) {
            return add(rhs.m_val);
        }
        inline Kahan& operator -= (const Kahan<int>& rhs) {
            return add(-rhs.m_val);
        }
        inline Kahan& operator= (int rhs) {
            m_val = rhs;
            return *this;
        }
        // these are prefix;
        inline Kahan& operator ++ () {
            m_val ++;
            return *this;
        }
        inline Kahan& operator -- () {
            m_val --;
            return *this;
        }
    private:
        inline Kahan& add(const int& rhs) {
            m_val += rhs;
            return *this;
        }
};
//UNFOLD

#endif /* __DEF_KAHAN__ */

//for vim modeline: (do not edit)
// vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=//%s:tags=.c_tags;:syn=cpp:ft=cpp:mps+=<\:>:ai:si:cin:nu:fo=croql:cino=p0t0c5(0:
