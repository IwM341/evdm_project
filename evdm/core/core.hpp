#ifndef CORE_HPP
#define CORE_HPP
#include "body_potential.hpp"
#include "grid_variants.hpp"
#include "histo_norm.hpp"
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <memory>
#include "dynamics/dynamics.hpp"
#include "r_conversion.hpp"
#include "traj_pool.hpp"

namespace evdm{

    template <typename T>
    using BodyModel = std::shared_ptr<Body<T>>;

    struct forward_value{
        template <typename T>
        static decltype(auto) forward(T && x){
            return std::forward<T>(x);
        } 
    };
    struct forward_shared{
        template <typename T>
        static decltype(auto) forward(T && x){
            return std::make_shared<std::decay_t<T>>(std::forward<T>(x));
        } 
    };

    template <typename T,typename forwarder = forward_value,
                typename IterableType,typename U>
    auto load_body_model(IterableType const &_array,size_t Rpoints,U Velocity){
        std::vector<T> Rho_values(Rpoints);
        for(size_t i=0;i<Rho_values.size();++i){
            Rho_values[i] = _array[i];
        }
        typename Body<T>::RFunc1_t Rho(
            typename Body<T>::GridR((T)0,(T)1,Rpoints),
            std::move(Rho_values)
        );
        return forwarder::forward(
                    Body<T>::FromRho(Velocity,std::move(Rho))
                );
    }

    template <typename forwarder = forward_value, typename T,typename U>
    auto load_body_model(std::vector<T> X,U Velocity) {
        std::vector<T> Rho_values(std::move(X));
        size_t size = Rho_values.size();
        
        typename Body<T>::RFunc1_t Rho(
            typename Body<T>::GridR((T)0, (T)1, size),
            std::move(Rho_values)
        );
       
        
        return forwarder::forward(
            Body<T>::FromRho(Velocity,std::move(Rho))
        );
    }

    template <typename T,typename forwarder = forward_value,
                typename Callable_t, typename U>
    auto make_body_model(Callable_t &&F,size_t Rpoints, U Velocity){
        std::vector<T> Rho_values(Rpoints);
        T h_inv = (T)1.0 / (Rpoints - 1);
        for(size_t i=0;i<Rho_values.size();++i){
            Rho_values[i] = F(i* h_inv);
        }
        typename Body<T>::RFunc1_t Rho(
            typename Body<T>::GridR((T)0,(T)1,Rpoints),
            std::move(Rho_values)
        );
        return forwarder::forward(
                    Body<T>::FromRho(Velocity,std::move(Rho))
                );
    }



    template <typename T>
    struct EL_Func{
        BodyModel<T> body;
        std::shared_ptr<typename Body<T>::LE_func_t> LE;

        static T get_vtype() {
            static_assert("EL_Func::get_vtype() is not callable");
        }
        EL_Func(BodyModel<T> const& body, decltype(LE) const& LE):
            body(body),LE(LE){}
        EL_Func(BodyModel<T> const& body,size_t Ne_bins):
            body(body),
            LE(
                std::make_shared<
                    typename Body<T>::LE_func_t
                >(body->get_le(Ne_bins))
            ){}
    };
    
    struct TrajPoolInitParams_t {
        bool is_static;
        float T_error;
        size_t traj_bins;
        size_t n_e_max;
        size_t n_l_max;
        TrajPoolInitParams_t(
            bool is_static = true,
            float T_error = 0.05,
            size_t traj_bins = 100,
            size_t n_e_max = 2,
            size_t n_l_max = 2
        ) :
            is_static(is_static),
            T_error(T_error),
            traj_bins(traj_bins),
            n_e_max(n_e_max),
            n_l_max(n_l_max) {}
    };

    template <typename BodyModel_vt,typename GridEL_vt, GridEL_type _grid_type>
    struct EL_Grid {
        using BM_vt = BodyModel_vt;
        using GEL_vt = GridEL_vt;
        using GridEL_t = GridEL<GridEL_vt, _grid_type>;
        using Body_t = Body<BM_vt>;



        static constexpr BodyModel_vt get_bm_vtype() {
            static_assert("EL_Grid::get_bm_vtype() is not callable");
        }
        static constexpr GEL_vt get_grid_vtype() {
            static_assert("EL_Grid::get_grid_vtype() is not callable");
        }
        static constexpr GridEL_type get_grid_type() {
            return _grid_type;
        }

        BodyModel<BM_vt> body; ///shared ptr to Body Model
        std::shared_ptr<GridEL_t> Grid; /// shared ptr to grid variants
        std::shared_ptr<typename Body<BM_vt>::LE_func_t > _LE; /// shared ptr to L(E) struct function
        std::shared_ptr<std::vector<bin_traj_pool_t<GridEL_vt>>> _TrajPools;

        template <typename GridType>
        EL_Grid(
            BodyModel<BM_vt> const& body, GridType&& EL_Grid,
            size_t NE_e_size = 100, TrajPoolInitParams_t TrajPoolInitParams = TrajPoolInitParams_t()) :
            body(body),
            Grid(std::make_shared<GridEL_t>(std::forward<GridType>(EL_Grid))),
            _LE(std::make_shared<typename Body<BM_vt>::LE_func_t>(body->get_le(NE_e_size)))
        {
            auto const& mEL_Grid = Grid->inner(0);
            //_TrajPools = std::make_shared<>
            std::vector<bin_traj_pool_t<GridEL_vt>> mTP;
            mTP.reserve(mEL_Grid.size());
            for (auto m_bin : mEL_Grid) {
                mTP.push_back(
                    GetTrajPool<GridEL_vt>(
                        TrajPoolInitParams.is_static,
                        m_bin, *body, *_LE,
                        TrajPoolInitParams.T_error,
                        TrajPoolInitParams.traj_bins,
                        TrajPoolInitParams.n_e_max,
                        TrajPoolInitParams.n_l_max
                    )
                );
            }
            _TrajPools = std::make_shared<
                std::vector<bin_traj_pool_t<GridEL_vt>>
            >(std::move(mTP));
        }

        inline EL_Func<BM_vt> getLE()const {
            return EL_Func<BM_vt>(body, _LE);
        }
        inline auto const & getBodyLE()const {
            return *_LE;
        }
        auto LE()const {
            return _LE->i_lm();
        }
        auto const & getLE_inner_grid() const {
            return Grid->inner(0);
        }

        using TrajPools_t = std::vector<bin_traj_pool_t<GridEL_vt>>;
        TrajPools_t const& TrajPools() const {
            return *_TrajPools;
        }

    };


    template <GridEL_type grid_type,typename BM_vt,typename GridType>
    auto make_core_GridEL(BodyModel<BM_vt> const& body,GridType && in_EL_Grid, size_t NE_e_size = 100){
        using EL_t = decltype(in_EL_Grid.inner(0).grid()[0].left);
        //in_EL_Grid.inner(0) -- is EL grid
        //in_EL_Grid.inner(0).grid() -- is E grid
        //in_EL_Grid.inner(0).grid()[0] -- is Rect<T>
        //Rect<T>.left is T
        // we need only T
        return EL_Grid<BM_vt,EL_t, grid_type>(body,std::forward<GridType>( in_EL_Grid));
    }


    template <typename T,typename Body_vt,typename GridEL_vt, GridEL_type grid_type>
    struct Distribution{
        using GEL_vt = GridEL_vt;
        using GridEL_t = EL_Grid<Body_vt,GEL_vt, grid_type>;
        
        GridEL_t Grid;
        typedef Eigen::VectorX<T> _Vec_t;
        typedef Eigen::VectorBlock<Eigen::VectorX<T>> Vec_t;
        //using const_Vec_t = Eigen::VectorX<T>::ConstSegmentReturnType;
        typedef Eigen::VectorBlock<const Eigen::VectorX<T>> const_Vec_t;
    private: 
        std::shared_ptr<Eigen::VectorX<T>> _Values;
        Vec_t _ValuesView;
        const_Vec_t _const_ValuesView;
        size_t padding;

        static Eigen::VectorX<T> process_vector(Eigen::VectorX<T> X,size_t _size,size_t padding) {
            size_t delta = _size % padding;
            if (X.size() >= _size && X.size() % padding == 0) {
                return std::move(X);
            }
            else {
                Eigen::VectorX<T> _X(delta ? _size + padding - delta : _size);
                _X.setZero();
                _X.segment(0, X.size()).noalias() = X;
                return std::move(_X);
            }
        }

    public:

        
        void check_ptype(int ptype) const{
            if (!(ptype < grid().grid().size())) {
                throw std::runtime_error("ptype should be less than max ptypes");
            }
        }
        static constexpr Body_vt get_body_vtype() {
            static_assert("Distribution::get_body_vtype() is not callable");
        };
        static constexpr T get_vtype() 
        {
            static_assert("Distribution::get_vtype() is not callable");
        }
        static constexpr GridEL_vt get_grid_vtype() {
            static_assert("Distribution::get_grid_vtype() is not callable");
        }
        static constexpr GridEL_type get_grid_type() {
            return grid_type;
        }
        
        inline _Vec_t& raw_vector() {
            return *_Values;
        }
        inline _Vec_t const& raw_vector()const {
            return *_Values;
        }

        inline auto const &grid()const{
            return *(Grid.Grid);
        }
        inline Body<Body_vt> const& body() const {
            return *Grid.body;
        }
        inline auto body_ptr() const {
            return Grid.body;
        }

        inline auto const& TrajPools() const {
            return Grid.TrajPools();
        }

        inline Vec_t values() {
            return _ValuesView;
        } 
        inline const_Vec_t values() const {
            return _const_ValuesView;
        }
        inline auto as_histo() {
            return grob::make_histo_view(grid(), values());
        }
        inline auto as_histo() const {
            return grob::make_histo_view(grid(), values());
        }

        Distribution(GridEL_t const &Grid, Eigen::VectorX<T> Values = {},size_t padding = 128) : 
            Grid(Grid),_Values(std::make_shared< Eigen::VectorX<T>>(process_vector(std::move(Values), grid().size(), padding))),
            _ValuesView(_Values->segment(0,grid().size())),
            _const_ValuesView( static_cast<const Eigen::VectorX<T> & >(*_Values).segment(0, grid().size()))
        {}
        
        template <typename U>
        Distribution<U, Body_vt, GridEL_vt, grid_type> as_type() const {
            return Distribution<U, Body_vt, GridEL_vt, grid_type>(
                Grid,_Values->template cast<U>()
                );
        }

        T count(int ptype)const {
            if (ptype < 0) {
                return _const_ValuesView.sum();
            }
            else {
                check_ptype(ptype);
                size_t N_size = grid().inner(0).size();
                //Eigen::VectorXf Values;
                return _const_ValuesView.segment(ptype * N_size,N_size).sum();
            }
        }
        void CopyBroadcast(size_t ptype_src, size_t ptype_dst,T weight = 1) {
            check_ptype(ptype_src);
            check_ptype(ptype_dst);
            size_t N_size = grid().inner(0).size();
            _ValuesView.segment(ptype_dst * N_size, N_size) =
                _ValuesView.segment(ptype_src * N_size, N_size);
        }
        void SummBroadcast(size_t ptype_src, size_t ptype_dst, T weight = 1) {
            check_ptype(ptype_src);
            check_ptype(ptype_dst);
            size_t N_size = grid().inner(0).size();
            _ValuesView.segment(ptype_dst * N_size, N_size) +=
                _ValuesView.segment(ptype_src * N_size, N_size)* weight;
        }
        Eigen::VectorBlock<Eigen::VectorX<T>> block(int ptype) {
            size_t N_size = grid().inner(0).size();
            if (ptype >= 0) {
                check_ptype(ptype);
                return _Values->segment(ptype * N_size, N_size);
            }
            else {
                return _ValuesView;
            }
            
        }
        const Eigen::VectorBlock<const Eigen::VectorX<T>> block(int ptype) const{
            size_t N_size = grid().inner(0).size();
            if (ptype >= 0) {
                check_ptype(ptype);
                return static_cast<const Eigen::VectorX<T> &>(*_Values).
                    segment(ptype * N_size, N_size);
            }
            else {
                return _const_ValuesView;
            }
        }

        template <typename Gen_t,typename RGrid_t>
        auto r_dens(
            std::vector<size_t> ptypes, Gen_t && G, RGrid_t && RGrid,
            size_t Nmk_per_bin = 1000,
            progress_omp_function<> m_prog_bar_func = progress_omp_function<>()
        ) const {
            for (auto ptype : ptypes) {
                check_ptype(ptype);
            }
            return convert_r_density(
                as_histo(),TrajPools(), 
                ptypes, Grid.LE(), body().Phi,body()._DD_F_C(1),
                G, std::forward<RGrid_t>(RGrid), Nmk_per_bin, m_prog_bar_func
            );
        }
    };
    template <
        typename T,typename Init_t, 
        typename Body_vt, typename GridEL_vt, 
        GridEL_type grid_type
    >
    auto make_Distribution(
        EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid, Init_t&& init) 
    {
        Distribution<T, Body_vt, GridEL_vt, grid_type> Dstr(Grid);
        auto LE = Dstr.Grid.LE();
        if constexpr (!std::is_same_v<bool, std::decay_t<Init_t>>) {
            auto H1 = Dstr.as_histo();
            for (auto [bin, val] : H1) {
                auto type = std::get<0>(bin);
                auto [e0, e1] = std::get<0>(bin.tail());
                auto [l0, l1] = std::get<1>(bin.tail());
                auto Lm0 = LE(-e0);
                auto Lm1 = LE(-e1);
                val = evdm::integrateAB2([&](auto e) {
                    auto Lm = LE(-e);
                    return evdm::integrateAB2([&](auto l) {
                            return init(e, l);
                    }, Lm * l0, Lm * l1, 1);
                }, e0, e1, 1);
            }
        } 
        return Dstr;
    }
    
    template <typename V1,typename BT1,typename T1, GridEL_type GT1,
              typename V2,typename BT2,typename T2, GridEL_type GT2,
                typename Deg_t,typename BinMeasure_t>
    auto distrib_norm(
        Distribution<V1, BT1,T1, GT1> const& D1,
        Distribution<V2, BT2,T2, GT2> const& D2,
        Deg_t deg,
        BinMeasure_t mu
    ) {
        auto m_bin_norm = get_bin_mes(mu, D1.Grid.LE());
        return histo_norm(
            D1.as_histo(),
            D2.as_histo(),
            deg, m_bin_norm
        );
    }

    template <
        typename T,
        typename Body_vt, typename GridEL_vt,
        GridEL_type grid_type
    >
    auto make_Distribution_data(
        EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid, 
        T * _data,size_t _size,size_t padding
    ) {
        size_t _size_grid = Grid.Grid->size();
        size_t _new_size = _size_grid + padding - _size_grid % padding;
        Eigen::VectorX<T> X(_new_size);
        X.setZero();
        size_t _copy_size = std::min(_size, _new_size);
        Eigen::Map<Eigen::VectorXf, Eigen::Unaligned> Data_Map(_data, _copy_size);
        X.segment(0, _copy_size) = Data_Map;
        return Distribution<T, Body_vt, GridEL_vt,grid_type>(
            Grid,
            std::move(X)
        );
    }

    template <typename T>
    using Matrix_t = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

    template <typename T, typename Body_vt, 
        typename GridEL_vt, GridEL_type grid_type>
    struct GridMatrix{
        using GEL_vt = GridEL_vt;
        using GridEL_t = EL_Grid<Body_vt, GEL_vt, grid_type>;

        GridEL_t Grid;

        typedef Matrix_t<T> _Mat_t;
        typedef Eigen::Block<Matrix_t<T>> Mat_t;
        //using const_Vec_t = Eigen::VectorX<T>::ConstSegmentReturnType;
        typedef Eigen::Block<const Matrix_t<T>> const_Mat_t;
    private:
        std::shared_ptr<_Mat_t> _Values;
        Mat_t _ValuesView;
        const_Mat_t _const_ValuesView;
        size_t padding;


    static _Mat_t process_matrix(
            _Mat_t X, size_t _N, 
            size_t padding) {
        size_t delta = _N % padding;
        if (X.cols() >= _N && X.size() % padding == 0) {
            return std::move(X);
        }
        else {
            size_t _N_true = delta ? _N + padding - delta : _N;
            _Mat_t _X(_N_true, _N_true);
            _X.setZero();
            _X.block(0,0, X.cols(), X.rows()).noalias() = X;
            return std::move(_X);
        }
    }

    public:


        static constexpr Body_vt get_body_vtype() {
            static_assert("Distribution::get_body_vtype() is not callable");
        };
        static constexpr T get_vtype()
        {
            static_assert("Distribution::get_vtype() is not callable");
        }
        static constexpr GridEL_vt get_grid_vtype() {
            static_assert("Distribution::get_grid_vtype() is not callable");
        }
        static constexpr GridEL_type get_grid_type() {
            return grid_type;
        }


        void check_ptype(int ptype) const {
            if (!(ptype < grid().grid().size())) {
                throw std::runtime_error("ptype should be less than max ptypes");
            }
        }


        inline Body<Body_vt> const& body() const {
            return *Grid.body;
        }
        inline auto body_ptr() const {
            return Grid.body;
        }
        inline auto const& grid()const {
            return *(Grid.Grid);
        }
        inline Mat_t values() {
            return _ValuesView;
        }
        inline const_Mat_t values() const {
            return _const_ValuesView;
        }
        inline _Mat_t & raw_matrix() {
            return *_Values;
        }
        inline _Mat_t const& raw_matrix() const{
            return *_Values;
        }

        template <typename Index>
        inline auto as_histo(Index const & _I) {
            size_t _shift = grid().LinearIndex(_I);
            return grob::make_histo_view(grid(), _ValuesView.col(_shift));
        }
        
        template <typename Index>
        inline auto as_histo(Index const& _I) const{
            size_t _shift = grid().LinearIndex(_I);
            return grob::make_histo_view(grid(), _ValuesView.col(_shift));
        }

        inline auto as_histo_i(size_t i) {
            return grob::make_histo_view(grid(), _ValuesView.col(i));
        }
        inline auto as_histo_i(size_t i) const {
            return grob::make_histo_view(grid(), _ValuesView.col(i));
        }

        inline auto as_histo_bank(size_t ptype_in, size_t ptype_out) {
            check_ptype(ptype_in);
            check_ptype(ptype_out);
            auto const& _grid = grid().inner(0); //geom grid without ptypes
            size_t inner_size = _grid.size(); //size of geom grid
            size_t ptypes = grid().grid().size(); //number of ptypes
            return grob::make_grid_object_ref(
                _grid,
                grob::as_container(
                    [&_grid, m_block = block(ptype_in, ptype_out)](size_t i)mutable {
                        auto m_col = m_block.col(i);
                        return grob::make_histo_view(_grid,
                            grob::vector_view(m_col.data(), m_col.size())
                        );
                    }, _grid.size()
                )
            );
        }

        GridMatrix(GridEL_t const& Grid,
            Matrix_t<T> Values = {},size_t padding = 128) :
            Grid(Grid), _Values(
                std::make_shared< Matrix_t<T>>(
                    process_matrix(std::move(Values), Grid.Grid->size(), padding)
                )), 
            _ValuesView(_Values->block(0, 0, grid().size(), grid().size())),
            _const_ValuesView(
                static_cast<_Mat_t const &>(*_Values)
                    .block(0, 0, grid().size(), grid().size())
            ) {}

        template <typename U>
        GridMatrix<U, Body_vt, GridEL_vt, grid_type> as_type() const {
            return GridMatrix<U, Body_vt, GridEL_vt, grid_type>(
                Grid, _Values->template cast<U>()
            );
        }

        void CopyBroadcast(
            size_t ptype_in_src, size_t ptype_out_src, 
            size_t ptype_in_dst, size_t ptype_out_dst , T weight = 1) {

            check_ptype(ptype_in_src);
            check_ptype(ptype_out_src);
            check_ptype(ptype_in_dst);
            check_ptype(ptype_out_dst);

            size_t N_size = grid().inner(0).size();
            
            _ValuesView.block(ptype_out_dst * N_size, ptype_in_src* N_size,N_size, N_size) =
                _ValuesView.block(ptype_out_src * N_size, ptype_in_src * N_size, N_size, N_size)
                *weight;
        }
        void SummBroadcast(
            size_t ptype_in_src, size_t ptype_out_src,
            size_t ptype_in_dst, size_t ptype_out_dst, T weight = 1) {

            check_ptype(ptype_in_src);
            check_ptype(ptype_out_src);
            check_ptype(ptype_in_dst);
            check_ptype(ptype_out_dst);

            size_t N_size = grid().inner(0).size();

            _ValuesView.block(ptype_out_dst * N_size, ptype_in_src * N_size, N_size, N_size) =
                _ValuesView.block(ptype_out_src * N_size, ptype_in_src * N_size, N_size, N_size)
                * weight;
        }

        inline Eigen::Block<Matrix_t<T>> block(int ptype_in, int ptype_out){

            

            size_t N_size = grid().inner(0).size();
            size_t N_size_full = grid().size();
            if (ptype_in < 0 || ptype_out < 0) {
                return _ValuesView;
            }
            else {
                check_ptype(ptype_in);
                check_ptype(ptype_out);
                return raw_matrix().block(
                    ptype_out * N_size, ptype_in * N_size,
                    N_size, N_size
                );
            }
        }
        inline Eigen::Block<const Matrix_t<T>> block(int ptype_in, int ptype_out)const{

            
            
            size_t N_size = grid().inner(0).size();
            size_t N_size_full = grid().size();
            if (ptype_in < 0 || ptype_out < 0) {
                return _const_ValuesView;
            }
            else {
                check_ptype(ptype_in);
                check_ptype(ptype_out);
                return raw_matrix().block(
                    ptype_out * N_size, ptype_in * N_size,
                    N_size, N_size
                );
            }
        }
    };

    template <typename T, typename Body_vt, typename GridEL_vt, GridEL_type grid_type>
    auto make_Matrix(EL_Grid<Body_vt, GridEL_vt, grid_type> const& Grid) {
        GridMatrix<T, Body_vt, GridEL_vt, grid_type> Dstr(Grid);
        return Dstr;
    }
};
#endif//CORE_HPP