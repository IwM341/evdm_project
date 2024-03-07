#ifndef VARIANT_TOOLS_HPP
#define VARIANT_TOOLS_HPP
#include <variant>
#include <tuple>

#ifndef MAX_RETURN_ALIGN
#define MAX_RETURN_ALIGN 4*sizeof(max_align_t)
#endif
namespace evdm {

    template <typename T>
    struct self_type{
        typedef T type;
    };
    namespace __visit_detail{

        template <typename T>
        struct aligned_wrapper{
            static constexpr size_t _size = sizeof(T);
            static constexpr size_t _full_size = alignof(T) + sizeof(T);
            char data[_full_size];
            size_t _shift;
            template <typename...Args>
            aligned_wrapper(Args&&...args){
                _shift = alignof(T) - (size_t)(data)%alignof(T);
                new (data + _shift) T(std::forward<Args>(args)...);
            }
            T const & self() const &{
                return *reinterpret_cast<T*>(data+_shift);
            }
            T & self()&{
                return *reinterpret_cast<T*>(data+_shift);
            }
            T && self() &&{
                return std::move(*reinterpret_cast<T*>(data+_shift));
            }

            operator T const &() const & {
                return self();
            }   
            operator T &() &{
                return self();
            }
            operator T &&() && {
                return self();
            }

            ~aligned_wrapper(){
                reinterpret_cast<T*>(data+_shift)->~T();
            }
        };
        template <typename VariantOrWrapperType>
        struct get_variant_type{
            template <typename VariantType>
            inline static decltype(auto) id(VariantType && V){
                return std::forward<VariantType>(V);
            }
        };

        
        template <typename VariantType>
        struct get_variant_type<aligned_wrapper<VariantType>>{
            template <typename WrapperType>
            inline static decltype(auto) id(WrapperType && V){
                return std::forward<WrapperType>(V).self();
            }
        };

        template <typename variant_type,typename InitializerLambda,size_t...I>
        constexpr inline auto make_variant_impl(size_t i,InitializerLambda && Init,std::index_sequence<I...>){
            constexpr size_t variantSize = sizeof...(I);
            typedef typename std::conditional<MAX_RETURN_ALIGN < alignof(variant_type),
                aligned_wrapper<variant_type>,
                variant_type
            >::type wrap_var_t;
            std::array<wrap_var_t (*)(InitializerLambda &&),variantSize> Funcs {
                [](InitializerLambda && _Init){
                    typedef typename std::variant_alternative<I,variant_type>::type type;
                    return wrap_var_t (_Init(self_type<type>{}));
                }...
            };
            return Funcs[i](std::forward<InitializerLambda>(Init));
        }

        template <typename variant_type, typename _TpForward,size_t...I>
        variant_type make_variant_alt_impl(size_t i,_TpForward&& _Tuple, std::index_sequence<I...>) {
            constexpr size_t variantSize = sizeof...(I);
            std::array<variant_type(*)(_TpForward&& AlterTp), variantSize> Funcs{
                [](_TpForward&& AlterTp) {
                    return variant_type(std::get<I>(std::forward<_TpForward>(AlterTp)));
                }...
            };
            return Funcs[i](std::forward<_TpForward>(_Tuple));
        }

    };
    
    /// @brief 
    /// @tparam variant_type std::variant of type
    /// @tparam InitializerLambda
    /// @param index index of std::variant
    /// @param Init lambda, which have template parametr of type self_type, indicating the return type
    /// @return std::varint = Init(self_type<type>{})
    template <typename variant_type,typename InitializerLambda>
    auto make_variant(size_t index,InitializerLambda && Init){
        return __visit_detail::make_variant_impl<variant_type>(index,Init,
        std::make_index_sequence<std::variant_size<variant_type>::value>{});
    }

    /// <summary>
    /// make variant from alternatives with particular
    /// </summary>
    /// <param name="index">index of alternateive</param>
    /// <param name="...Alternatives">all alternamtives</param>
    /// <returns></returns>
    template <
        typename...Alt_ts,
        typename variant_type = std::variant<std::decay_t<Alt_ts>...>
    >
    variant_type make_variant_alt(size_t index, Alt_ts &&...Alternatives) {
        if (index >= sizeof...(Alternatives)) {
            throw std::out_of_range("in make_variant_alt index out of range");
        }
        return __visit_detail::make_variant_alt_impl<variant_type>(
            index, std::forward_as_tuple(std::forward< Alt_ts>(Alternatives)...),
            std::make_index_sequence<sizeof...(Alternatives)>{}
        );
    }

    /// <summary>
    /// converts choosen alternative to type
    /// </summary>
    /// <tparam name="return_type">return type</tparam>
    /// <param name="...Alternatives">all alternamtives</param>
    /// <returns></returns>
    template <
        typename return_type,
        typename...Alt_ts
    >
    return_type make_choose_alt(size_t index, Alt_ts &&...Alternatives) {
        if (index >= sizeof...(Alternatives)) {
            throw std::out_of_range("in make_variant_alt index out of range");
        }
        return __visit_detail::make_variant_alt_impl<return_type>(
            index, std::forward_as_tuple(std::forward< Alt_ts>(Alternatives)...),
            std::make_index_sequence<sizeof...(Alternatives)>{}
        );
    }



    template <typename VariantWrapper>
    inline constexpr decltype(auto) vmove(
        VariantWrapper && V
    ) {
        return __visit_detail::get_variant_type<
            typename std::decay<VariantWrapper>::type
        >::id(std::move(V));
    }

    template <typename VariantWrapper>
    inline constexpr decltype(auto) vforward(
        std::remove_reference_t<VariantWrapper>&& V
    ){
        return __visit_detail::get_variant_type<
            typename std::decay<VariantWrapper>::type
        >::id(std::forward<VariantWrapper>(V));
    }

    template <typename VariantWrapper>
    inline constexpr decltype(auto) vforward(
        std::remove_reference_t<VariantWrapper>& V
    ) {
        static_assert(!std::is_lvalue_reference<VariantWrapper>::value, "template argument"
            "substituting VariantWrapper should be an lvalue reference type");
        return __visit_detail::get_variant_type<
            typename std::decay<VariantWrapper>::type
        >::id(std::forward<VariantWrapper>(V));
    }


    /// @brief makes x = y
    /// @param x std variant<Args1...>
    /// @param y std variant<Args2...>, Args2... should be in Args1...
    template <typename Variant1,typename Variant2>
    inline void variant_assign(Variant1 & x,Variant2 && y){
        std::visit([&x]<class T>(T && value){
            x = std::forward<T>(value);
        },std::forward<Variant2>(y));
    }

    /// @brief makes x = y
    /// @tparam OutputType output type parametr
    /// @param y std variant<Args2...>, Args2... should be in Args1...
    template <typename OutputType,typename Variant2>
    inline constexpr auto variant_cast(Variant2 && y){
        typedef typename std::conditional<MAX_RETURN_ALIGN < alignof(OutputType),
            __visit_detail::aligned_wrapper<OutputType>,
            OutputType
        >::type wrap_var_t;
        return std::visit([]<class T>(T && value)->wrap_var_t{
            return wrap_var_t(std::forward<T>(value));
        },vforward<Variant2>(y));
    }
};

#endif//VARIANT_TOOLS_HPP`