#ifndef wellStruct_H
#define wellStruct_H

#include "volFields.H"

namespace Foam
{

struct well
{
    word name;

    labelList cells;

    bool injector;
    scalar Fb_inj;
    scalar Cs_inj;

    bool bhpControl;

    scalar target;

    scalar rate;
    scalar bhp;
    
};

} // End namespace Foam

#endif
 