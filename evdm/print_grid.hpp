#ifndef PRINT_GRID_HPP
#define PRINT_GRID_HPP
#include <iostream>
namespace evdm{
    /// @brief get number of line segments of grid, which could be plotted (line is [x0,x1],[y0,y1]) 
    /// @param Grid 2 dimensional LE multi grid
    /// @return number of line segments
    template <typename GridType>
    size_t grid_lines_num(GridType const & Grid){
        return Grid.Grid.size() +2 + Grid.size();
    }

    /// @brief filles arrays X0,X1,Y0,Y1 with lines of borders of grids
    /// @param Grid 
    /// @param LE_func dependence function L(E)
    /// @param X0 array of x coords of line segment begin
    /// @param X1 array of x coords of line segment end
    /// @param Y0 array of y coords of line segment begin
    /// @param Y1 array of y coords of line segment end
    template <typename GridType,typename LE_Functype,
            typename X0_array_t,typename X1_array_t,typename Y0_array_t,typename Y1_array_t  >
    void grid_lines_fill_arrays(GridType const & Grid,LE_Functype const & LE_func,X0_array_t & X0,X1_array_t & X1,
                                    Y0_array_t & Y0,Y1_array_t & Y1){
        
        size_t index = 0;
        // inserting bottom segment --- (Emin,0)->(0,0)
        // +1
        X0[0] = Grid.Grid.front();
        X1[0] = Grid.Grid.back();
        Y0[0] = 0;
        Y1[0] = 0;
        index ++;
        // inserting the leftest segment --- (Emin,0)->(Emin,L(Emin))
        // + 1
        X0[index] = Grid.Grid.front();
        X1[index] = Grid.Grid.front();
        Y0[index] = 0;
        Y1[index] = LE_func(Grid.Grid.front());
        index ++;

        // Ne vertical lines: (E_i,0)->(E_i,L(E_i))
        // + Ne
        for(auto [E_left,E_right] : Grid.Grid){
            X0[index] = E_right;
            X1[index] = E_right;
            Y0[index] = 0;
            Y1[index] = LE_func(E_right);
            index ++;
        }
        // horizontal lines: (E_i_l,E_i_r)->(L(E_i_l),L(E_i_r))
        // + Grid.size()
        for(auto [Elr,Llr] : Grid.Grid){
            X0[index] = Elr.left();
            X1[index] = Elr.right();
            Y0[index] = LE_func(left)*Llr.right();
            Y1[index] = LE_func(E_right)*Llr.right();
            index ++;
        }
    }
};
#endif//PRINT_GRID_HPP