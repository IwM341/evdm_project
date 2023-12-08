#ifndef MK_HPP
#define MK_HPP


namespace evdm{
    template <typename ResultType,typename DensityType = double>
    struct MCResult{
        ResultType Result;
        DensityType RemainDensity;
        MCResult(ResultType Result,double RemainDensity = 1.0):Result(Result),RemainDensity(RemainDensity){}
        MCResult():RemainDensity(0){}
    };


    template <typename T>
    struct MCIntegral{
        T Result;
        T Sigma;
        MCIntegral(T Result,T Sigma = 0):Result(Result),Sigma(Sigma){}

        operator T (){
            return Result;
        }
        operator std::pair<T,T> (){
            return std::pair<T,T>(Result,Sigma);
        }
    };

    template<typename GenType>
    auto MCIntegrate(GenType G,size_t N){
        typedef typename std::decay<decltype(G())>::type value_type;
        value_type sum = 0;
        value_type sum2 = 0;
        for (size_t i=0;i<N;++i){
            auto f = G();
            sum += f;
            sum2 += f*f;
        }
        sum = sum/N;
        sum2 = sum2/N;
        return MCIntegral<value_type>(sum,sqrt(sum2-sum*sum));
    }
};
#endif//MK_HPP