#ifndef PRINT_GRID_HPP
#define PRINT_GRID_HPP
#include <iostream>
#include <Eigen/Dense>
#include "grid_variants.hpp"
#include "triangulation.hpp"
#include "measure.hpp"
namespace evdm{
    
    struct tri_sizes {
        size_t pts;///num of points
        size_t trs;///num of triangles
    };
    
    template <typename T, GridEL_type _grid_type>
    tri_sizes GridEl_TriSizes(
        GridEL<T, _grid_type> const& _el_grid
    ) {
        auto const& grid = _el_grid.inner(0);
        return {
            1 + grid.grid().size() + grid.size(),
            2* grid.size() - grid.inner().back().size()
        };
    }

    namespace tuple_operations{
        
        template <typename _Tuple1,typename _Tuple2, size_t...Is,typename Operator_t>
        inline void tuple_assign(_Tuple1 & X,_Tuple2 const & Y,Operator_t && assignator,std::index_sequence<Is...>){
            (assignator(std::get<Is>(X),std::get<Is>(Y)),...);
        }

        #define DECLARE_TUPLE_ASSIGN_OPERATION(op)\
        template <typename _Tuple1,typename _Tuple2>\
        inline decltype(auto)  operator op(_Tuple1 && X,_Tuple2 const & Y){\
            static_assert(std::tuple_size_v<std::decay_t<_Tuple1>> == \
                std::tuple_size_v<std::decay_t<_Tuple1>>,\
                "assigning tuples requires equal size");\
            tuple_assign(X,Y,[](auto & x,auto const & y){x op y;},\
                std::make_index_sequence<std::tuple_size_v<std::decay_t<_Tuple1>>>{});\
            return std::forward<_Tuple1>(X);\
        }
        DECLARE_TUPLE_ASSIGN_OPERATION(+=)
        DECLARE_TUPLE_ASSIGN_OPERATION(-=)
        DECLARE_TUPLE_ASSIGN_OPERATION(*=)
        DECLARE_TUPLE_ASSIGN_OPERATION(/=)
    };
    template < 
        typename HistoValuesType,typename LE_Functype,
        typename ValueArray,typename VertexArray,typename TriangleArray
    >
    void DirstributionPrinting(
        HistoValuesType const & Ni,
        size_t ptype,
        LE_Functype && LE,
        ValueArray && Values,
        VertexArray && Vertexes,
        TriangleArray && Triangles,
        const tri_sizes m_trs
    ) {
        auto const& grid = Ni.Grid.inner(0);
        
        typedef double T;
        
        Eigen::VectorX<T> Products(m_trs.pts);
        Products.setZero();
        Eigen::SparseMatrix<T> GrammMatrix(m_trs.pts,m_trs.pts);
        GrammMatrix.reserve(9*m_trs.pts);
        GrammMatrix.setZero();
        //the first triangle is absolutly the same as first rectangele;
        size_t triangle_index = 0;
        size_t vertex_index = 0;

        auto PointEL = [&](auto _e,auto _l){
            return Point2<T>{_e,LE(-_e)*_l};
        };
        /// ****************** e - energy, ln - normlized momentum
        auto getV = [&](size_t i)->Point2<T> {
            return Point2<T>{Vertexes[2*i],Vertexes[2*i+1]};
        };
        auto getT = [&](size_t i)->std::tuple<size_t,size_t,size_t>{
            return {Triangles[3*i],Triangles[3*i+1],Triangles[3*i+2]};
        };
        auto get3Prods = [&](size_t i0,size_t i1,size_t i2){
            return std::tie(Products[i0],Products[i1],Products[i2]);
        };

        auto append_vertex = [&](auto e,auto ln)->size_t{
            if (vertex_index >= m_trs.pts) {
                throw std::out_of_range("vertex_index out of range");
            }
            Vertexes[2*vertex_index] = e;
            Vertexes[2*vertex_index + 1] = (T)ln*LE(-e);
            return vertex_index++;
        };
        auto append_triangle = [&](size_t i0,size_t i1,size_t i2){
            if (triangle_index >= m_trs.trs) {
                throw std::out_of_range("triangle_index out of range");
            }
            Triangles[3*triangle_index] = i0; 
            Triangles[3*triangle_index+1] = i1; 
            Triangles[3*triangle_index+2] = i2;
            return triangle_index++;
        };
        auto L_limits = [&](size_t i_e,auto l0,auto l1,auto l2)->std::pair<size_t,size_t>{
            typedef decltype(l0) m_T;
            auto [l_min,l_max] = (l1 >= l0) ? 
                std::make_pair(std::min( (m_T) l0,(m_T) l2 ),std::max( (m_T) l1,(m_T) l2 )) :
                std::make_pair(std::min( (m_T) l1,(m_T) l2 ),std::max( (m_T) l0,(m_T) l2 ));
            
            return {grid.inner(i_e).pos(l_min),grid.inner(i_e).pos(l_max)+1};
        };

        auto mu_EL = get_bin_mes(measure_dEdL{},LE);

        auto _e_start = grid.grid()[0].left;
        auto _e_start_1 = grid.grid()[0].right;
        size_t i_start = append_vertex(_e_start,0);
        size_t i_right = append_vertex(_e_start_1,0);
        size_t i_prev_e = i_right;

        ///PART 0 FIRST LEFTEST BIN
        for(size_t i=0;i<grid.inner(0).size();++i){
            size_t index = append_vertex(_e_start_1,grid.inner(0)[i].right);
            auto Hvalue = Ni[grob::make_MI(ptype,0ULL,i)];
            append_triangle(i_start,i_prev_e,index);

            Products[i_start] += Hvalue/3;
            Products[i_prev_e] += Hvalue/3;
            Products[index] += Hvalue/3;
            i_prev_e = index;
        }
        ///PART 1 ITERATING OVER E
        size_t start_index_prev_e = i_right;
        using namespace tuple_operations;
        for(size_t i=1;i<grid.grid().size();++i){
            auto [e_0,e_1] = grid.grid()[i];
            size_t l_prev = grid.inner(i-1).size();
            size_t l_tmp = grid.inner(i).size();
            size_t l_min = std::min(l_prev,l_tmp);
            size_t l_max = std::max(l_prev,l_tmp);

            bool left_single_point = l_tmp >= l_prev;

            size_t start_index_tmp_e = append_vertex(e_1,0);
            for(size_t j_l =0;j_l < l_min;++j_l){
                auto [l_0_d,l_0_u] = grid.inner(i-1)[j_l];
                auto [l_1_d,l_1_u] = grid.inner(i)[j_l];

                auto P_0_d = PointEL(e_0,l_0_d); //start_index_prev_e + j_l
                size_t I0d = start_index_prev_e + j_l;
                auto P_0_u = PointEL(e_0,l_0_u); //start_index_prev_e + j_l + 1
                size_t I0u = start_index_prev_e + j_l+1;
                auto P_1_d = PointEL(e_1,l_1_d); //start_index_tmp_e + j_l
                size_t I1d = start_index_tmp_e + j_l;
                auto P_1_u = PointEL(e_1,l_1_u); //start_index_tmp_e + j_l + 1
                size_t I1u = start_index_tmp_e + j_l+1;
                size_t _check = append_vertex(e_1,l_1_u);
                append_triangle(I0d,I0u,I1d);
                append_triangle(I0u,I1d,I1u);

                std::array<Point2<T>,3> Tr1{P_0_d,P_0_u,P_1_d};
                std::array<Point2<T>,3> Tr2{P_0_u,P_1_d,P_1_u};
                auto S1 = TriangleArea(P_0_d,P_0_u,P_1_d);
                auto S2 = TriangleArea(P_0_u,P_1_d,P_1_u);
                auto [j_h_min_tr1,j_h_max_tr1] = L_limits(i,l_0_d,l_0_u,l_1_d);
                for(size_t j_h = j_h_min_tr1;j_h < j_h_max_tr1;++j_h ){
                    auto N_h = Ni[grob::make_MI(ptype,i,j_h)];
                    auto Hvalue = N_h/mu_EL(grid[grob::make_MI(i,j_h)]);
                    
                    auto [l_h_d,l_h_u] = grid.inner(i)[j_h];
                    std::array<Point2<T>,4> Bin{PointEL(e_0,l_h_d),PointEL(e_0,l_h_u),
                                                PointEL(e_1,l_h_u),PointEL(e_1,l_h_d)};
                    auto ISecArea1 = IntersectionPolygon(Tr1,Bin);
                    auto ISecArea2 = IntersectionPolygon(Tr2,Bin);
                    auto _integrals1 = TriangleIntegrals(Tr1,ISecArea1,Hvalue);
                    auto _integrals2 = TriangleIntegrals(Tr2,ISecArea2,Hvalue);
                    get3Prods(I0d,I0u,I1d) += _integrals1;
                    get3Prods(I0u,I1d,I1u) += _integrals2;
                }
            }
            size_t i_singl = (left_single_point ? i-1 : i);
            size_t i_multiple = (left_single_point ? i : i-1);
            auto e_single =  (left_single_point ? e_0 : e_1);
            auto e_multiple =  (left_single_point ? e_1 : e_0);
            
            size_t index_single = (left_single_point ? start_index_prev_e : start_index_tmp_e);
            size_t index_multiple = (left_single_point ? start_index_tmp_e : start_index_prev_e);

            
            auto l_s = grid.inner(i_singl)[l_min-1].right;
            auto Ps = PointEL(e_single,l_s);
            size_t Is = index_single + l_min;

            for (size_t j_l =l_min;j_l < l_max;++j_l){
                auto [l_d,l_u] = grid.inner(i_multiple)[j_l];
                auto Pd = PointEL(e_multiple,l_d);
                size_t Imd = index_multiple+j_l;
                auto Pu = PointEL(e_multiple,l_u);
                size_t Imu = index_multiple+j_l + 1;

                //_chck should be Imu
                if(left_single_point){
                    size_t _check = append_vertex(e_multiple,l_u);
                }
                append_triangle(Is,Imd,Imu);

                std::array<Point2<T>,3> Tr{Ps,Pd,Pu};

                auto [j_h_min_tr1,j_h_max_tr1] = L_limits(i,l_s,l_d,l_u);
                for(size_t j_h = j_h_min_tr1;j_h < j_h_max_tr1;++j_h ){
                    auto Hvalue = Ni[grob::make_MI(ptype,i,j_h)]/mu_EL(grid[grob::make_MI(i,j_h)]);
                    
                    auto [l_h_d,l_h_u] = grid.inner(i)[j_h];
                    std::array<Point2<T>,4> Bin{PointEL(e_0,l_h_d),PointEL(e_0,l_h_u),
                                                PointEL(e_1,l_h_u),PointEL(e_1,l_h_d)};
                    auto ISecArea = IntersectionPolygon(Tr,Bin);
                    get3Prods(Is,Imd,Imu) += TriangleIntegrals(Tr,ISecArea,Hvalue);
                }
            }
            start_index_prev_e = start_index_tmp_e;
        }

        std::cout << "triangle_index: " << triangle_index << "/" <<m_trs.trs << std::endl;
        std::cout << "vertex_index: " << vertex_index << "/" << m_trs.pts <<std::endl;
        ///Fill Gramm Matrix:
        for(size_t i=0;i<m_trs.trs;++i){
            auto [i0,i1,i2] = getT(i);
            Point2<T> P0 = getV(i0);
            Point2<T> P1 = getV(i1);
            Point2<T> P2 = getV(i2);
            T S = TriangleArea(P0,P1,P2);
            GrammMatrix.coeffRef(i0,i0) += S/6;
            GrammMatrix.coeffRef(i1,i1) += S/6;
            GrammMatrix.coeffRef(i2,i2) += S/6;

            GrammMatrix.coeffRef(i0,i1) += S/12;
            GrammMatrix.coeffRef(i0,i2) += S/12;
            GrammMatrix.coeffRef(i1,i2) += S/12;

            GrammMatrix.coeffRef(i1,i0) += S/12;
            GrammMatrix.coeffRef(i2,i0) += S/12;
            GrammMatrix.coeffRef(i2,i1) += S/12;
        }
        
        
        Eigen::ConjugateGradient<decltype(GrammMatrix), Eigen::Upper> solver;
        solver.compute(GrammMatrix);
        if(solver.info()!=Eigen::Success) {
            throw std::runtime_error("error trying to compute matrix");
        }
        Eigen::VectorX<T> VertSolution;
        VertSolution = solver.solve(Products);
        if(solver.info()!=Eigen::Success) {
            //std::cout << "Products:\n" << Products << std::endl;
            //std::cout << "GrammMatrix:\n" << GrammMatrix << std::endl;
            throw std::runtime_error("error trying to solve matrix");
        }
        for(int i=0;i<VertSolution.size();++i){
            Values[i] = VertSolution[i];
        }
    }

    template <typename T, GridEL_type _grid_type>
    tri_sizes GridEl_TriSizes_1order(
        GridEL<T, _grid_type> const& _el_grid
    ) {
        auto const& grid = _el_grid.inner(0);
        return {
            grid.size()*4,
            grid.size()*2
        };
    }

    template <
        typename HistoValuesType, typename LE_Functype,
        typename ValueArray, typename XArray_t, typename YArray_t, 
        typename TriangleArray,typename measure_t = measure_dEdL
    >
    void DirstributionPrinting_1order(
        HistoValuesType const& Ni,
        size_t ptype,
        LE_Functype&& LE,
        ValueArray&& Values,
        XArray_t && XArray,
        YArray_t&& YArray,
        TriangleArray&& Triangles,
        const tri_sizes m_trs,
        measure_t m_measure_inst = {}
    ) {
        using T = std::decay_t<decltype(Ni.Values[0])>;

        auto const& grid = Ni.Grid.inner(0);
        size_t _ptype_stride = Ni.Grid.grid().size();
        auto & m_Vals = Ni.Values;
        
        size_t bin_index = 0;
        for (auto m_bin : grid) {
            auto [e0, e1] = std::get<0>(m_bin);
            auto [l0, l1] = std::get<1>(m_bin);
            auto Lm0 = LE(-e0);
            auto Lm1 = LE(-e1);
            auto L00 = l0 * Lm0;
            auto L01 = l1 * Lm0;
            auto L10 = l0 * Lm1;
            auto L11 = l1 * Lm1;

            size_t I0 = 4 * bin_index;

            XArray[I0] = e0;
            XArray[I0+1] = e0;
            XArray[I0+2] = e1;
            XArray[I0+3] = e1;

            YArray[I0] = L00;
            YArray[I0 + 1] = L01;
            YArray[I0 + 2] = L10;
            YArray[I0 + 3] = L11;
            for (size_t i = 0; i < 3; ++i) {
                Triangles[6 * bin_index + i] = I0 + i;
            }
            for (size_t i = 0; i < 3; ++i) {
                Triangles[6 * bin_index+3 + i] = I0 + i + 1;
            }
            T m_dens = m_Vals[_ptype_stride * ptype + bin_index]/
                mes_down(measure_t{}, m_bin, Lm0, Lm1);
            for (size_t i = 0; i < 4; ++i) {
                Values[4 * bin_index+i] = m_dens;
            }
            ++bin_index;
        }
    }

};
#endif//PRINT_GRID_HPP