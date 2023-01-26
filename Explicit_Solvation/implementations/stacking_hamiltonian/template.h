#ifndef _TEMPLATE_H_
#define _TEMPLATE_H_

template <typename T>
void print (T arr_type){

    for ( int i{0}; i<static_cast<int>(arr_type.size()); ++i){
        std::cout << arr_type[i] << " | ";
    }
    std::cout << std::endl; 
}



#endif 