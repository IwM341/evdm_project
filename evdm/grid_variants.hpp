#ifndef GRID_VARIANTS_HPP
#define GRID_VARIANTS_HPP
#include <grob/multigrid.hpp>
#include <grob/grid_gen.hpp>
#include <variant>
#include "utils/variant_tools.hpp"
namespace evdm{

    enum class GridEL_variants {
        GridCUU, GridCVV
    };


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

        /*
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
        }*/

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


        /*
        /// @brief  greates GridCVV grid
        /// @param ptypes number of m-states
        /// @param Emin minimum enegy
        /// @param Ne number of bins in E grid
        /// @param E_density_view lambda : function of t ( where t in [0,1]), 
        /// which indicates density of bins in E grid
        /// @param Nl_E_func lambda : number of bins in L grid, depending on E 
        /// @param L_density_view lambda : function of (E,t), where E - energy, 
        /// t - density param, which indicates density of buins in L grid
        /// @return 
        template <typename Ne_array_t,typename E_density_view_array_t,typename NL_E_array_t ,typename L_density_array_view_t>
        static GridCVV grid_var_VV(size_t ptypes,double Emin,Ne_array_t const & Ne,E_density_view_array_t const& E_density_view,
                                    NL_E_array_t const& Nl_E_func,L_density_array_view_t const& L_density_view){
            
            return grob::make_grid(grob::GridRangeLight<size_t>(ptypes),
                [](size_t i){
                    auto E_Grid = grob::GridVectorHisto<double>(
                        grob::density_vector_nt<double>(Emin,0,E_density_view[i],std::max(2,Ne[i]+1))
                    );
                    return grob::make_grid_f(E_Grid,
                        [d](auto const & rect_E){
                            return grob::GridVectorHisto<double>(
                                grob::density_vector_nt<double>(0.0,1.0,
                                    std::max(2,Nl_E_func[i](rect_E.center())+1),
                                    [&rect_E](auto l){
                                        return L_density_view[i](rect_E.center(),l);
                                    }
                                )
                            );
                        }
                    );
                }
            );
        }*/

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




    };



};


#endif//GRID_VARIANTS_HPP