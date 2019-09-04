//
//  AAD.cpp
//  
//
//  Created by Troels Tang Karlsen on 04/09/2019.
//

#include <stdio.h>

#include "number.h"

Tape globalTape;
thread_local Tape* number::tape = &globalTape;


