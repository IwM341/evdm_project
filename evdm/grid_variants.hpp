#ifndef GRID_VARIANTS_HPP
#define GRID_VARIANTS_HPP
#include <grob/multigrid.hpp>
#include <grob/grid_gen.hpp>
#include <variant>
#include "utils/variant_tools.hpp"
#include <concepts>
namespace evdm{

    enum class GridEL_type {
        GridCUU, GridCVV
    };

    template <GridEL_type _value>
    struct GridEL_type_t {
        constexpr static GridEL_type value = _value;
    };
    
    template <typename mGridEL_type_t>
    concept GridEL_type_C = requires{
        {mGridEL_type_t::value}->std::same_as<GridEL_type>;
    };

    

    template <typename T,GridEL_type>
    struct Grid_types;

    template <typename T>
    struct Grid_types<T, GridEL_type::GridCUU> {
        typedef grob::Point<grob::Rect<T>, grob::Rect<T>> bin_t;
        typedef grob::GridUniformHisto<T> GridE;
        typedef grob::GridUniformHisto<T> GridL;
        typedef grob::GridRangeLight<size_t> Range;
        typedef grob::MultiGrid<GridE, std::vector<GridL>> InnerEL_t;
        typedef grob::ConstValueVector<InnerEL_t> MGContainer_t;
        typedef grob::MultiGrid< Range, MGContainer_t> GridEL_t;

        /// @brief greates GridCUU grid
        /// @tparam NL_E_Lambda 
        /// @param ptypes number of m-states
        /// @param Emin  minimum enegy
        /// @param Ne number of bins in E grid
        /// @param Nl_E_func lambda : number of bins in L grid, depending on E
        /// @return
        template <typename NL_E_Lambda>
        inline static GridEL_t construct(size_t ptypes, T Emin, size_t Ne, NL_E_Lambda&& Nl_E_func) {
            auto E_Grid = GridE(Emin, 0, std::max((size_t)2, Ne + 1));
            return grob::mesh_grids(grob::GridRangeLight<size_t>(ptypes),
                grob::make_grid_f(E_Grid,
                    [&](size_t i) {
                        T t_e = 1 - E_Grid[i].center() / Emin;
                        return GridL(0.0, 1.0,
                            std::max((size_t)2, (size_t)Nl_E_func(t_e) + 1)
                        );
                    }
                )
            );
        }

        static GridEL_t Refine(GridEL_t const& G,size_t Er,size_t Lr) {
            if (Er == 0 || Lr == 0) {
                throw std::runtime_error("refine parameter Er or Lr is zero");
            }
            Range R1 = G.grid();
            InnerEL_t const & grid = G.inner(0);
            auto const& mgrid_e = grid.grid().unhisto();
            GridE Egrid_new = GridE(
                mgrid_e.front(), 
                mgrid_e.back(), 
                (mgrid_e.size() - 1) * Er + 1
            );
            const size_t Nl_old = grid.inner().size();
            std::vector<GridL> GridL_s(Nl_old *Er);
            for (size_t i = 0; i < Nl_old; ++i) {
                auto const& mgrid_l = grid.inner(i).unhisto();
                for (size_t j = 0; j < Er; ++j) {
                    GridL_s[i*Er + j] = GridL(
                        mgrid_l.front(),
                        mgrid_l.back(),
                        (mgrid_l.size() - 1) * Lr + 1
                    );
                }
            }
            return GridEL_t(R1,
                MGContainer_t(
                    InnerEL_t(
                        std::move(Egrid_new), std::move(GridL_s)
                    ), R1.size()
                )
            );
        }

        inline const char* type() {
            return "GridCUU";
        }
    };

    template <typename T>
    struct Grid_types<T, GridEL_type::GridCVV> {
        typedef grob::Point<grob::Rect<T>, grob::Rect<T>> bin_t;
        typedef grob::GridVectorHisto<T> GridE;
        typedef grob::GridVectorHisto<T> GridL;
        typedef grob::GridRangeLight<size_t> Range;
        typedef grob::MultiGrid<GridE, std::vector<GridL>> InnerEL_t;
        typedef grob::ConstValueVector<InnerEL_t> MGContainer_t;
        typedef grob::MultiGrid< Range, MGContainer_t> GridEL_t;

        /// @brief  greates GridCVV grid
        /// @param ptypes number of m-states
        /// @param Emin minimum enegy
        /// @param Ne number of bins in E grid
        /// @param E_density_view lambda : function of t ( where t in [0,1]), 
        /// which indicates density of buins in E grid
        /// @param Nl_E_func lambda : number of bins in L grid, depending on E 
        /// @param L_density_view lambda : function of (E,t), where E - energy, 
        /// t - density param, which indicates density of buins in L grid
        /// @return 
        template <typename E_density_view_t, typename NL_E_Lambda, typename L_density_view_t>
        static GridEL_t construct(size_t ptypes, T Emin, size_t Ne, E_density_view_t&& E_density_view,
            NL_E_Lambda&& Nl_E_func, L_density_view_t&& L_density_view) {

            //std::cout << "Emin = " << Emin << std::endl;
            auto E_Grid = GridE(
                grob::density_vector_nt<T>(Emin, 0,
                    E_density_view,
                    std::max((size_t)2, Ne + 1))
            );

            size_t i = E_Grid.pos(-0.5); // try to find and shift nearest  to 0.5 point
            if (std::abs(E_Grid.unhisto()[i] + 0.5) > std::abs(E_Grid.unhisto()[i + 1] + 0.5)) {
                i = i + 1;
            }
            E_Grid.unhisto()[i] = -0.5;
            //std::cout << "Egird = " << E_Grid << std::endl;
            return grob::mesh_grids(grob::GridRangeLight<size_t>(ptypes),
                grob::make_grid_f(E_Grid,
                    [&](size_t i) {
                        auto _E = 1 - E_Grid[i].center() / Emin;
                        //std::cout << "_E = " << _E << ", Emin = " << Emin << std::endl;
                        return GridL(
                            grob::density_vector_nt<T>(0.0f, 1.0f,
                                [&](auto l) {
                                    return L_density_view(_E, l);
                                },
                                std::max((size_t)2, (size_t)Nl_E_func(_E) + 1)
                            )
                        );
                    }
                )
            );
        }

        static GridEL_t Refine(GridEL_t const& G, size_t Er, size_t Lr) {
            if (Er == 0 || Lr == 0) {
                throw std::runtime_error("rifine parameter Er or Lr is zero");
            }
            Range R1 = G.grid();
            InnerEL_t const& grid = G.inner(0);
            auto const& mgrid_e = grid.grid().unhisto();
            
            std::vector<T> Epoints( (mgrid_e.size() - 1)* Er + 1);
            std::vector<GridL> Lpoints;
            Lpoints.reserve(Epoints.size());

            Epoints.back() = mgrid_e.back();
            for (size_t i = 0; i < mgrid_e.size() - 1; ++i) {
                T h = (mgrid_e[i + 1] - mgrid_e[i]) / Er;
                std::vector<T> const & ml_g = grid.inner(i);
                for (size_t j = 0; j < Er; ++j) {
                    Epoints[Er * i + j] = mgrid_e[i] + h * j;

                    std::vector<T> Lg((ml_g.size() - 1) * Lr + 1);
                    Lg.back() = ml_g.back();
                    for (size_t k = 0; k < ml_g.size() - 1;++k) {
                        T h_l = (ml_g[k + 1] - ml_g[k]) / Er;
                        for (size_t m = 0; m < Lr; ++m) {
                            Lg[k * Lr + m] = ml_g[k] + h_l * m;
                        }
                    }
                    Lpoints.push_back(std::move(Lg));
                }
            }

            return GridEL_t(R1,
                MGContainer_t(
                    InnerEL_t(
                        std::move(Epoints), std::move(Lpoints)
                    ),R1.size()
                )
            );
        }

        inline const char* type() {
            return "GridCVV";
        }
    };

    template <typename T, GridEL_type GridEL_type_ >
    struct GridEL : public Grid_types<T, GridEL_type_>::GridEL_t 
    {
        using Helper = Grid_types<T, GridEL_type_>;
        using GridELbase = typename Grid_types<T, GridEL_type_>::GridEL_t;
        using GridELbase::GridELbase;
        using GridE = typename Helper::GridE;
        using GridL = typename Helper::GridL;

        GridEL(GridELbase m_grid) :GridELbase(std::move(m_grid)) {}
        using Constructor = Grid_types<T, GridEL_type_>;
        
        inline size_t size1()const {
            return this->inner(0).size();
        }
        inline size_t ptypes()const {
            return this->grid().size();
        }
        inline size_t size_e()const {
            return this->inner(0).grid().size();
        }
        GridEL refine(size_t Ne, size_t Nl) const{
            return Helper::Refine(*this, Ne, Nl);
        }
    };

    struct printing_sizes {
        size_t verticals;
        size_t horisontals;
        static constexpr size_t bottom = 1;
        constexpr inline size_t size() const{
            return verticals + horisontals + bottom;
        }
    };
    template <typename T, GridEL_type GridEL_type_>
    printing_sizes GridEL_printing_size(GridEL<T, GridEL_type_> const & _el_grid)
    {
        auto const& grid = _el_grid.inner(0);
        return printing_sizes{
            grid.grid().size() + 1,
            grid.size()
        };
    }

    template <typename T, GridEL_type GridEL_type_,typename NormFunctype,typename Array_t>
    void GridEL_printing_fill(
        GridEL<T, GridEL_type_> const& _el_grid,
        NormFunctype && LE_or_1_func,
        Array_t&& X0, Array_t&& X1,
        Array_t&& Y0, Array_t&& Y1,
        printing_sizes _sizes)
    {
        auto const & _grid = _el_grid.inner(0);
        
        size_t max_size = _sizes.size();
        size_t array_index = 0;
        auto append_value = [&](T e0, T e1, T l0, T l1) {
            if (array_index < max_size) {
                X0[array_index] = e0;
                X1[array_index] = e1;
                Y0[array_index] = LE_or_1_func(e0) * l0;
                Y1[array_index] = LE_or_1_func(e1) * l1;
                array_index++;
            }
        };

        append_value(_grid.grid().front().left, _grid.grid().back().right, 0, 0);
        append_value(_grid.grid().front().left, _grid.grid().front().left, 0, 1);

        for (size_t i = 0; i < _grid.grid().size(); ++i) {
            auto _e0 = _grid.grid()[i].left;
            auto _e1 = _grid.grid()[i].right;
            append_value(_e1, _e1, 0, 1);
            for (size_t j = 0; j < _grid.inner(i).size(); ++j) {
                auto l1 = _grid.inner(i)[j].right;
                append_value(_e0, _e1, l1, l1);
            }
        }
    }
};

/*
    template <typename T>
    struct GridEL{
        typedef grob::Point<grob::Rect<T>,grob::Rect<T>> bin_t;

        typedef grob::GridRangeLight<size_t> Range;
        typedef grob::GridUniformHisto<T> GridU;
        typedef grob::GridVectorHisto<T> GridV;    

        typedef grob::MultiGrid< Range,
                    grob::ConstValueVector<
                            grob::MultiGrid<GridU,std::vector<GridU>>
                        >
                    > GridCUU;
        typedef grob::MultiGrid< Range,
                    grob::ConstValueVector<
                            grob::MultiGrid<GridV,std::vector<GridV>>
                        >
                    > GridCVV;
        
        typedef std::variant<GridCUU,GridCVV> GridType;
        static const char * get_type_name(size_t id) {
            switch(id){
                case 0:
                    return "GridCUU";
                case 1:
                    return "GridCVV";
                default:
                    return "ErrorType";
            }
        }
        GridType Grid;
        template <typename Grid_t>
        GridEL(Grid_t &&Grid):Grid(std::forward<Grid_t>(Grid)){}
        size_t size()const{
            return std::visit([](auto const & grd){return grd.size();},Grid);
        }
        size_t size_e()const {
            return std::visit([](auto const& grd) {return grd.inner(0).size(); }, Grid);
        }
        size_t size1()const {
            return std::visit([](auto const& grd) {return grd.inner(0).size(); }, Grid);
        }
        size_t ptypes()const {
            return std::visit([](auto const& grd) {return grd.grid().size(); }, Grid);
        }
        /// @brief greates GridCUU grid
        /// @tparam NL_E_Lambda 
        /// @param ptypes number of m-states
        /// @param Emin  minimum enegy
        /// @param Ne number of bins in E grid
        /// @param Nl_E_func lambda : number of bins in L grid, depending on E
        /// @return 
        template <typename NL_E_Lambda>
        static GridCUU grid_var_CUU(size_t ptypes,T Emin,size_t Ne,NL_E_Lambda && Nl_E_func){
            
            auto E_Grid = GridU(Emin,0,std::max((size_t)2,Ne+1));
            return grob::mesh_grids(grob::GridRangeLight<size_t>(ptypes),
                 grob::make_grid_f(E_Grid,
                    [&](size_t i){ 
                        T t_e = 1 - E_Grid[i].center() /Emin;
                        return GridU(0.0,1.0,
                            std::max((size_t)2,(size_t)Nl_E_func(t_e)+1)
                        );
                    }
                )
            );
        }

        ***
        /// @brief greates GridCUU grid
        /// @tparam NL_E_Lambda 
        /// @param ptypes number of m-states
        /// @param Emin  minimum enegy
        /// @param Ne number of bins in E grid
        /// @param Nl_E_func lambda : number of bins in L grid, depending on E
        /// @return 
        template <typename Ne_array_t,typename NL_E_Lambda_array_t>
        static GridUU grid_var_UU(size_t ptypes,double Emin,Ne_array_t const & Ne,NL_E_Lambda_array_t const & Nl_E_func){
            
            return grob::make_grid(grob::GridRangeLight<size_t>(ptypes),
                [&](size_t i){
                    auto E_Grid = grob::GridUniformHisto<double>(Emin,0,std::max(Ne[i]+1,2));
                    return grob::make_grid_f(E_Grid,
                        [&](auto const & rect_E){
                            return grob::GridUniformHisto<double>(0.0,1.0,
                                std::max(2,Nl_E_func[i](rect_E.center())+1)
                            );
                        }
                    );
                }
            );
        }*//*

        /// @brief  greates GridCVV grid
        /// @param ptypes number of m-states
        /// @param Emin minimum enegy
        /// @param Ne number of bins in E grid
        /// @param E_density_view lambda : function of t ( where t in [0,1]), 
        /// which indicates density of buins in E grid
        /// @param Nl_E_func lambda : number of bins in L grid, depending on E 
        /// @param L_density_view lambda : function of (E,t), where E - energy, 
        /// t - density param, which indicates density of buins in L grid
        /// @return 
        template <typename E_density_view_t,typename NL_E_Lambda,typename L_density_view_t>
        static GridCVV grid_var_CVV(size_t ptypes,T Emin,size_t Ne,E_density_view_t && E_density_view,
                                    NL_E_Lambda && Nl_E_func,L_density_view_t && L_density_view){
            
            auto E_Grid = GridV(
                grob::density_vector_nt<T>(Emin, 0, 
                    [&](auto _e) {return E_density_view(1-_e/Emin) ; },
                    std::max((size_t)2, Ne + 1))
                );

            size_t i = E_Grid.pos(-0.5); // try to find and shift nearest  to 0.5 point
            if(std::abs(E_Grid.unhisto()[i] + 0.5) > std::abs(E_Grid.unhisto()[i+1] + 0.5) ){
                i = i+1;
            }
            E_Grid.unhisto()[i] = -0.5;

            return grob::mesh_grids(grob::GridRangeLight<size_t>(ptypes),
                grob::make_grid_f(E_Grid,
                    [&](size_t i){
                        auto _E = 1 - E_Grid[i].center() /Emin;
                        return GridV(
                            grob::density_vector_nt<T>(0.0f,1.0f,
                                [&](auto l){
                                    return L_density_view(_E,l);
                                },
                                std::max((size_t)2,(size_t)Nl_E_func(_E)+1)
                            )
                        );
                    }
                )
            );
        }



        template <typename Serializator>
        auto Serialize(Serializator && S)const{
            auto m_type = stools::Serialize(Grid.index(),S);

            return S.MakeDict(
                std::array<std::string,3>{"type_id","type_name","grid"},
                std::array<decltype(m_type),3>{
                    m_type,
                    stools::Serialize(get_type_name(Grid.index()),S),
                    std::visit([&S](auto const & _grid){return _grid.Serialize(S);},Grid)
                }
            );
        }

        template <typename Object_t,typename DeSerializator_t>
        static auto DeSerialize(Object_t && Obj,DeSerializator_t && DS){
            size_t type_id = stools::DeSerialize<size_t>(DS.GetProperty(Obj,"type_id"),DS);
            return make_variant(type_id,
                [&](auto m_type){
                    return stools::DeSerialize<typename decltype(m_type)::type>
                        (DS.GetProperty(Obj,"grid"),DS);
                }
            );
        }




    };*/


#endif//GRID_VARIANTS_HPP