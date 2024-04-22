#pragma once

template <typename T>
using BodyModel = std::shared_ptr<Body<T>>;

struct forward_value {
    template <typename T>
    static decltype(auto) forward(T&& x) {
        return std::forward<T>(x);
    }
};
struct forward_shared {
    template <typename T>
    static decltype(auto) forward(T&& x) {
        return std::make_shared<std::decay_t<T>>(std::forward<T>(x));
    }
};

template <typename T, typename forwarder = forward_value,
    typename IterableType, typename U>
auto load_body_model(IterableType const& _array, size_t Rpoints, U Velocity) {
    std::vector<T> Rho_values(Rpoints);
    for (size_t i = 0; i < Rho_values.size(); ++i) {
        Rho_values[i] = _array[i];
    }
    typename Body<T>::RFunc1_t Rho(
        typename Body<T>::GridR((T)0, (T)1, Rpoints),
        std::move(Rho_values)
    );
    return forwarder::forward(
        Body<T>::FromRho(Velocity, std::move(Rho))
    );
}

template <typename forwarder = forward_value, typename T, typename U>
auto load_body_model(std::vector<T> X, U Velocity) {
    std::vector<T> Rho_values(std::move(X));
    size_t size = Rho_values.size();

    typename Body<T>::RFunc1_t Rho(
        typename Body<T>::GridR((T)0, (T)1, size),
        std::move(Rho_values)
    );


    return forwarder::forward(
        Body<T>::FromRho(Velocity, std::move(Rho))
    );
}

template <typename T, typename forwarder = forward_value,
    typename Callable_t, typename U>
auto make_body_model(Callable_t&& F, size_t Rpoints, U Velocity) {
    std::vector<T> Rho_values(Rpoints);
    T h_inv = (T)1.0 / (Rpoints - 1);
    for (size_t i = 0; i < Rho_values.size(); ++i) {
        Rho_values[i] = F(i * h_inv);
    }
    typename Body<T>::RFunc1_t Rho(
        typename Body<T>::GridR((T)0, (T)1, Rpoints),
        std::move(Rho_values)
    );
    return forwarder::forward(
        Body<T>::FromRho(Velocity, std::move(Rho))
    );
}



template <typename T>
struct EL_Func {
    BodyModel<T> body;
    std::shared_ptr<typename Body<T>::LE_func_t> LE;

    static T get_vtype() {
        static_assert("EL_Func::get_vtype() is not callable");
    }
    EL_Func(BodyModel<T> const& body, decltype(LE) const& LE) :
        body(body), LE(LE) {}
    EL_Func(BodyModel<T> const& body, size_t Ne_bins) :
        body(body),
        LE(
            std::make_shared<
            typename Body<T>::LE_func_t
            >(body->get_le(Ne_bins))
        ) {}
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