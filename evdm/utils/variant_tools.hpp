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

    };
    
    /// @brief 
    /// @tparam variant_type std::variant of type
    /// @tparam InitializerLambda 
    /// @tparam ...Args 
    /// @param Init lambda, which have template parametr of type self_type, indicating the return type
    /// @param index index of std::variant
    /// @return std::varint = Init(self_type<type>{})
    template <typename variant_type,typename InitializerLambda>
    auto make_variant(size_t index,InitializerLambda && Init){
        return __visit_detail::make_variant_impl<variant_type>(index,Init,
        std::make_index_sequence<std::variant_size<variant_type>::value>{});
    }

    template <typename VariantWrapper>
    inline decltype(auto) vforward(VariantWrapper && V){
        return __visit_detail::get_variant_type<
            typename std::decay<VariantWrapper>::type
        >::id(std::forward<VariantWrapper>(V));
    }  


    /// @brief makes x = y
    /// @param x std variant<Args1...>
    /// @param y std variant<Args2...>, Args2... should be in Args1...
    template <typename Variant1,typename Variant2>
    inline void variant_assign(Variant1 & x,Variant2 && y){
        std::visit([&x](auto && value){
            x = std::forward<decltype(value)>(value);
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