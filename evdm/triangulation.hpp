#ifndef TRIANGULATION_HPP
#define TRIANGULATION_HPP
#include <Eigen/Dense>

    namespace evdm{
    template <typename T>
    using Point2 =Eigen::Matrix<T,2,1>;

    template <typename T>
    using Mat2 =Eigen::Matrix<T,2,2>;
    
    /// @brief calculates triangle pasis functions into polygon. 
    /// @param MainTriVerts 3 vertices, clock counterwize
    /// @param Polygon array of  clock counterwize ordered points
    /// @return integral over intersection between polygon and triangle 
    /// of basis functions of polygon
    template <typename TriangleArray_t,typename PolyArray>
    inline auto IntersectionPolygon(
        TriangleArray_t MainTriVerts,
        PolyArray const & Polygon
    ){
        typedef std::decay_t<decltype(MainTriVerts[0][0])> T;
        std::vector<Point2<T>> Tmp(Polygon.size()+MainTriVerts.size());
        size_t Tmp_size = Polygon.size();
        for(size_t i=0;i<Polygon.size();++i){
            Tmp[i] = Polygon[i];
        }
        
        std::vector<Point2<T>> Result(Tmp.size());
        size_t Result_size = 0;

        T TypeEps = (std::is_same_v<T,double> ? 1e-11 : 1e-5);
        T MaxAbs =0;
        for(auto const & V:MainTriVerts){
            MaxAbs = std::max(MaxAbs,V.squaredNorm());
        }
        T MaxDelta = 0;
        for(size_t i=0;i<MainTriVerts.size()-1;++i){
            MaxDelta = std::max(
                MaxDelta,
                (MainTriVerts[i+1]-MainTriVerts[i]).squaredNorm());
        }

        auto Projector = [Norm = MaxDelta](const auto & p1, const auto & pr)->Point2<T>
        {
            return (p1.squaredNorm()*pr-p1.dot(pr)*p1)/Norm;
        };
        T eps = 2*MaxDelta*MaxAbs*TypeEps;
        auto InCondition = [eps](const auto & X,
                    const auto & P0,const auto & proj)->bool{
            return proj.dot(X-P0) + eps>=0;
        };
        auto Intersection = [](const auto & X1,const auto & V1,
                                const auto & X2,const auto & V2
                                )->Point2<T>
        {
            T _det = (V1[0]*V2[1]-V1[1]*V2[0]);
            auto dX = X2-X1;
            T t = (V2[1]*dX[0]-V2[0]*dX[1])/_det;
            return X1 + V1*t;
        };
        auto TN = MainTriVerts.size();
        for (size_t i=0;i<TN;++i){
            Result_size = 0;
            const auto &P0 = MainTriVerts[i];
            const auto &P1 = MainTriVerts[(i+1)%TN];
            const auto & P_r =  MainTriVerts[(i+2)%TN];
            Point2<T> p1 = P1-P0;
            Point2<T> pr = P_r-P0;
            auto prj = Projector(p1,pr);
            

            size_t _tmp_index = 0;
            bool last_cond = InCondition(Tmp[_tmp_index],P0,prj);
            if(last_cond){
                Result[Result_size++] = Tmp[_tmp_index];
                if(++_tmp_index == Tmp_size){
                    goto NextItem;
                }
                while(InCondition(Tmp[_tmp_index],P0,prj)){
                    Result[Result_size++] = Tmp[_tmp_index];
                    ++_tmp_index;
                    if(_tmp_index == Tmp_size){
                        goto NextItem;
                    }
                }
                Result[Result_size++] = Intersection(
                    Tmp[_tmp_index-1],
                    Tmp[_tmp_index]-Tmp[_tmp_index-1],
                    P0,p1);
                if(++_tmp_index == Tmp_size){
                    goto AddLastItem;
                } 
                
                while(!InCondition(Tmp[_tmp_index],P0,prj)){
                    ++_tmp_index;
                    if(_tmp_index == Tmp_size){
                        break;
                    }
                }
                AddLastItem:
                Result[Result_size++] = Intersection(
                    Tmp[_tmp_index-1],
                    Tmp[_tmp_index%Tmp_size]-Tmp[_tmp_index-1],
                    P0,p1); 
                while(_tmp_index < Tmp_size){
                    Result[Result_size++] = Tmp[_tmp_index];
                    ++_tmp_index;
                }
            } else {
                ++_tmp_index;
                while(!InCondition(Tmp[_tmp_index],P0,prj)){
                    ++_tmp_index;
                    if(_tmp_index == Tmp_size){
                        goto NextItem;
                    }
                }
                    
                Result[Result_size++] = Intersection(Tmp[_tmp_index-1],Tmp[_tmp_index]-Tmp[_tmp_index-1],P0,p1);
                Result[Result_size++] = Tmp[_tmp_index];
                if(++_tmp_index == Tmp_size){
                    goto NextItem;
                }
                while(InCondition(Tmp[_tmp_index],P0,prj)){
                    Result[Result_size++] = Tmp[_tmp_index];
                    ++_tmp_index;
                    if(_tmp_index == Tmp_size){
                        break;
                    }
                }
                Result[Result_size++] = Intersection(
                    Tmp[_tmp_index-1],
                    Tmp[_tmp_index%Tmp_size]-Tmp[_tmp_index-1],
                    P0,p1);
            }
            NextItem:
            std::swap(Tmp_size,Result_size);
            std::swap(Tmp,Result);
            Result_size = 0;
        }
        Tmp.resize(Tmp_size);
        return Tmp;
    }
    template <bool oriented = false,typename T>
    T TriangleArea(Point2<T> const & Reper,
                      Point2<T> const & P1,
                      Point2<T> const & P2){
        T sum = 0;
        const Point2<T> V1 = P1 - Reper; 
        const Point2<T> V2 = P2 - Reper;
        sum += (V1[0]*V2[1]-V1[1]*V2[0]);
        return (oriented ? sum/2 : std::abs(sum/2));
    }

    template <typename TriangleArray_t,typename PolyArray_t,typename WeightType = int>
    auto TriangleIntegrals(TriangleArray_t const & Tri,
                        PolyArray_t const & Poly,
                        WeightType _weight = 1){

        typedef std::decay_t<decltype(Tri[0][0])> T;
        const auto &P0 = Tri[0];
        const Point2<T> V1 = Tri[1] - P0;
        const Point2<T> V2 = Tri[2] - P0;

        Mat2<T> AlphaMat = 
            Mat2<T>({{V1[0],V2[0]},{V1[1],V2[1]}}).inverse();
        if(Poly.size() < 3){
            return std::tuple<T,T,T>(0,0,0);
        }
        const auto & PolyReper = Poly[0];

        T sum0 = 0;
        T sum1 = 0;
        T sum2 = 0;

        for(size_t i=1;i<Poly.size()-1;++i){
            T S = TriangleArea(PolyReper,Poly[i],Poly[i+1]);
            Point2<T> Pav = (PolyReper+Poly[i]+Poly[i+1])/3;
            Point2<T> al1al2 = AlphaMat*(Pav-P0);
            sum0 += S*(1 - al1al2[0] - al1al2[1]);
            sum1 += S*al1al2[0];
            sum2 += S*al1al2[1];
        }
        return std::tuple<T,T,T>(sum0*_weight,sum1*_weight,sum2*_weight);
    }


};
#endif//TRIANGULATION_HPP