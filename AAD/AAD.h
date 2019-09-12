//
//  AAD.h
//  AAD
//
//  Created by Troels Tang Karlsen on 04/09/2019.
//  Copyright © 2019 Troels Tang Karlsen. All rights reserved.
//

#ifndef AAD_h
#define AAD_h

#include "number.h"
#include "blockList.h"
#include "tape.h"

using namespace std;

// Put multiple numbers on tape at once
template<class IT>
inline void putOnTape(IT begin, IT end)
{
    for_each(begin, end, [](number& n){n.putOnTape();});
}

// Convert a range of numbers/doubles to doubles/numbers
template<class IT1, class IT2>
inline void convertCollection(IT1 srcBegin, IT1 srcEnd, IT2 destBegin)
{
    using destType = remove_reference_t<decltype(*destBegin)>;
    transform(srcBegin, srcEnd, destBegin,
              [](const auto& source){return destType(source);});
}



#endif /* AAD_h */
