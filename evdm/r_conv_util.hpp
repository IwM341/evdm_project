#pragma once
#include "utils/type_tr.hpp"
#include "utils/math_external.hpp"
#include "measure.hpp"

namespace evdm {
    template <typename T>
    std::pair<T, T> _line_parabole_solver(T x0, T x1, T y0, T y1, T x_p_0, T C_p) {
        auto TanL = (y1 - y0) / (x1 - x0);
        auto T2 = TanL * TanL;
        auto _sqrt_s = std::sqrt(C_p * (C_p - 4 * TanL * y0 - 4 * T2 * x_p_0 + 4 * T2 * x0));
        auto _b_s = 2 * T2 * x0 - 2 * TanL * y0 + C_p;
        auto delta_sqrt_div_2t2 =
            2 * T2 * x0 * x0 - 4 * TanL * x0 * y0 + 2 * C_p * x_p_0 + 2 * y0 * y0;
        if (TanL == 0) {
            return { x_p_0,std::numeric_limits<T>::infinity() };
        }
        if (!(_sqrt_s >= 0)) {
            return { std::numeric_limits<T>::infinity(),
                    std::numeric_limits<T>::infinity() };
        }
        else {
            if (_b_s > 0) {
                return { delta_sqrt_div_2t2 / (_sqrt_s + _b_s),
                    (_b_s + _sqrt_s) / (2 * T2) };
            }
            else if (_b_s < 0) {
                return { (_b_s - _sqrt_s) / (2 * T2),
                    delta_sqrt_div_2t2 / (_b_s - _sqrt_s) };
            }
            else {
                return { -_sqrt_s / (2 * T2) ,_sqrt_s / (2 * T2) };
                return { -_sqrt_s / (2 * T2) ,_sqrt_s / (2 * T2) };
            }
        }
    }


    template <typename EL_Functype, typename T>
    inline bin_dedl_t<T> _bin_cut(EL_Functype&& LEf, bin_dedl_t<T> const& m_bin, same_t<T> phi, same_t<T> r) {
        auto bin1 = m_bin;
        auto& [e0, e1] = std::get<0>(bin1);
        auto l0 = std::get<1>(bin1).left;
        auto [em0, em1] =
            _line_parabole_solver(e0, e1, (T)LEf(-e0) * l0, (T)LEf(-e1) * l0, -phi, r * r);
        em0 = std::max(em0, -(T)phi);
        e0 = std::max((em0 < e1 ? em0 : e1), e0);
        e1 = std::min((em1 > e0 ? em1 : e0), e1);
        return bin1;
    };
};