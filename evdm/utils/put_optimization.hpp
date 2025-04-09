#pragma once
#include <algorithm>
#include <type_traits>
#include <vector>
#include <tuple>
#include <grob/multigrid.hpp>


template <typename HistoType,typename IterType>
auto put_values(IterType ValuesBegin, IterType ValuesEnd, HistoType& H) {
	typedef std::decay_t<decltype((*ValuesBegin)[0])> T;
	typedef std::pair<T,T> ret_type;
	
	T evap = 0;
	T sum = 0;

	auto const& Grid = H.Grid;
	auto const &gridE = Grid.grid();

	if (ValuesBegin == ValuesEnd) {
		return (T)0;
	}
	size_t i_E_0 = gridE.pos((*ValuesBegin)[0]);
	size_t GEsize = gridE.size();
	std::sort(ValuesBegin, ValuesEnd, [](auto const& x, auto const& y) {return x[0] < y[0]; });


	
	auto it_start = ValuesBegin;
	auto EndVec = ValuesEnd;
	while (it_start != ValuesEnd) {
		for (; i_E_0 < GEsize; ++i_E_0) {
			if (gridE[i_E_0].right > (*it_start)[0]) {
				break;
			}
		}
		if (i_E_0 == GEsize) {
			break;
		}
		auto Ei = gridE[i_E_0].right;
		auto it_next = it_start;
		for (; (it_next != EndVec) && (*it_next)[0] < Ei; ++it_next) {}
		
		//Calc for Ls:
		auto const&  gridL = Grid.inner(i_E_0);
		size_t i_L_0 = gridL.pos((*it_start)[1]);
		size_t GLsize = gridL.size();
		std::sort(it_start, it_next, [](auto const& x, auto const& y) {return x[1] < y[1]; });
		while (it_start != it_next) {
			for (; i_L_0 < GLsize-1; ++i_L_0) {
				if (gridL[i_L_0].right > (*it_start)[1]) {
					break;
				}
			}
			
			for (; i_L_0 < GLsize;++i_L_0) {
				T bin_value = 0;
				T bin_value2 = 0;
				auto Ltmp = gridL[i_L_0].right;
				
				for (; it_start < it_next && (*it_start)[1] < Ltmp; ++it_start) {
					bin_value += (*it_start)[2];
				}
				sum += bin_value;
				H[{i_E_0, i_L_0}] += bin_value;
			}
		}
	}
	
	for (auto it = it_start; it != EndVec; ++it) {
		auto v = (*it)[2];
		sum += v;
		evap += v;
	}
	return evap;

}


template <typename HistoType, typename IterType>
auto put_values2(IterType ValuesBegin, IterType ValuesEnd, HistoType& H) {
	typedef std::decay_t<decltype((*ValuesBegin)[0])> T;
	typedef std::pair<T, T> ret_type;

	//T evap = 0;
	T sum = 0;
	T sum2 = 0;
	auto const& Grid = H.Grid;
	auto const& gridE = Grid.grid();

	if (ValuesBegin == ValuesEnd) {
		return ret_type(0, 0);
	}
	size_t i_E_0 = gridE.pos((*ValuesBegin)[0]);
	size_t GEsize = gridE.size();
	std::sort(ValuesBegin, ValuesEnd, [](auto const& x, auto const& y) {return x[0] < y[0]; });



	auto it_start = ValuesBegin;
	auto EndVec = ValuesEnd;
	while (it_start != ValuesEnd) {
		for (; i_E_0 < GEsize; ++i_E_0) {
			if (gridE[i_E_0].right > (*it_start)[0]) {
				break;
			}
		}
		if (i_E_0 == GEsize) {
			break;
		}
		auto Ei = gridE[i_E_0].right;
		auto it_next = it_start;
		for (; (it_next != EndVec) && (*it_next)[0] < Ei; ++it_next) {}

		//Calc for Ls:
		auto const& gridL = Grid.inner(i_E_0);
		size_t i_L_0 = gridL.pos((*it_start)[1]);
		size_t GLsize = gridL.size();
		std::sort(it_start, it_next, [](auto const& x, auto const& y) {return x[1] < y[1]; });
		while (it_start != it_next) {
			for (; i_L_0 < GLsize - 1; ++i_L_0) {
				if (gridL[i_L_0].right > (*it_start)[1]) {
					break;
				}
			}

			for (; i_L_0 < GLsize; ++i_L_0) {
				T bin_value = 0;
				T bin_value2 = 0;
				auto Ltmp = gridL[i_L_0].right;

				for (; it_start < it_next && (*it_start)[1] < Ltmp; ++it_start) {
					auto ds = (*it_start)[2];
					bin_value += ds;
					bin_value2 += (*it_start)[3];
				}
				sum += bin_value;
				sum2 += bin_value2;
				H[{i_E_0, i_L_0}] += bin_value;
			}
		}
	}
	//Calc Evaporation
	/*
	for (auto it = it_start; it != EndVec; ++it) {
		auto v = (*it)[2];
		sum += v;
		sum2 += v * v;
		evap += v;
	}
	Evap += evap;
	*/
	return ret_type(sum, sum2);
}