#include "crystal.h"
#include "global.h"

using namespace Vipster;

Mat Vipster::makeBravais(int ibrav, Vec axes, Vec angles)
{
    switch (ibrav) {
    case 1:
        // cubic
        return {{{{1, 0, 0}},
                 {{0, 1, 0}},
                 {{0, 0, 1}}}};
    case 2:
        // fcc
        return {{{{-0.5,  0.0,  0.5}},
                 {{ 0.0,  0.5,  0.5}},
                 {{-0.5,  0.5,  0.0}}}};
    case 3:
        // bcc
        return {{{{ 0.5,  0.5,  0.5}},
                 {{-0.5,  0.5,  0.5}},
                 {{-0.5, -0.5,  0.5}}}};
    case -3:
        // bcc, alternative
        return {{{{-0.5,  0.5,  0.5}},
                 {{ 0.5, -0.5,  0.5}},
                 {{ 0.5,  0.5, -0.5}}}};
    case 4:
        // hexagonal/trigonal
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        return {{{{ 1.0,          0.0, 0.0}},
                 {{-0.5, sqrt(3)*0.5, 0.0}},
                 {{ 0.0,          0.0, axes[2]}}}};
    case 5:
        // trigonal, 3fold axis z
        if(std::abs(angles[0]) >=1) throw Error("Step::makeCell: angles[0] is invalid");
        {
        auto tx = sqrt((1.-angles[0])/2.);
        auto ty = sqrt((1.-angles[0])/6.);
        auto tz = sqrt((1.+2*angles[0])/3.);
        return {{{{tx, -ty, tz}}, {{0, 2*ty, tz}},{{-tx, -ty, tz}}}};
        }
    case -5:
        // trigonal, 3fold axis <111>
        if(std::abs(angles[0]) >=1) throw Error("Step::makeCell: angles[0] is invalid");
        {
        auto ty = sqrt((1.-angles[0])/6.);
        auto tz = sqrt((1.+2*angles[0])/3.);
        auto u = (tz - 2*sqrt(2)*ty)/sqrt(3);
        auto v = (tz + sqrt(2)*ty)/sqrt(3);
        return {{{{u,v,v}}, {{v,u,v}}, {{v,v,u}}}};
        }
    case 6:
        // tetragonal
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        return {{{{1,0,0}}, {{0,1,0}}, {{0,0,axes[2]}}}};
    case 7:
        // bct
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        return {{{{ 0.5, -0.5, axes[2]*0.5}},
                 {{ 0.5,  0.5, axes[2]*0.5}},
                 {{-0.5, -0.5, axes[2]*0.5}}}};
    case 8:
        // orthorhombic
        if(float_comp(axes[1], 0.)) throw Error("Step::makeCell: axes[1] not provided");
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        return {{{{1,0,0}}, {{0,axes[1],0}}, {{0,0,axes[2]}}}};
    case 9:
        // bco
        if(float_comp(axes[1], 0.)) throw Error("Step::makeCell: axes[1] not provided");
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        return {{{{ 0.5, axes[1]*0.5, 0}},
                 {{-0.5, axes[1]*0.5, 0}},
                 {{ 0.0,         0.0, axes[2]}}}};
    case -9:
        // bco, alternative
        if(float_comp(axes[1], 0.)) throw Error("Step::makeCell: axes[1] not provided");
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        return {{{{ 0.5,-axes[1]*0.5, 0}},
                 {{ 0.5, axes[1]*0.5, 0}},
                 {{ 0.0,         0.0, axes[2]}}}};
    case 91:
        // bco, A type
        if(float_comp(axes[1], 0.)) throw Error("Step::makeCell: axes[1] not provided");
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        return {{{{ 0.5,         0.0,          0.0}},
                 {{ 0.0, axes[1]*0.5, -axes[2]*0.5}},
                 {{ 0.0, axes[1]*0.5,  axes[2]*0.5}}}};
    case 10:
        // fco
        if(float_comp(axes[1], 0.)) throw Error("Step::makeCell: axes[1] not provided");
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        return {{{{ 0.5,         0.0, axes[2]*0.5}},
                 {{ 0.5, axes[1]*0.5,         0.0}},
                 {{ 0.0, axes[1]*0.5, axes[2]*0.5}}}};
    case 11:
        // body-co
        if(float_comp(axes[1], 0.)) throw Error("Step::makeCell: axes[1] not provided");
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        return {{{{ 0.5,  axes[1]*0.5, axes[2]*0.5}},
                 {{-0.5,  axes[1]*0.5, axes[2]*0.5}},
                 {{-0.5, -axes[1]*0.5, axes[2]*0.5}}}};
    case 12:
        // monoclinic, unique axis z
        if(float_comp(axes[1], 0.)) throw Error("Step::makeCell: axes[1] not provided");
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        if(std::abs(angles[0]) >=1) throw Error("Step::makeCell: angles[0] is invalid");
        return {{{{ 1, 0, 0}},
                 {{ axes[1]*angles[0], axes[1]*sqrt(1.-pow(angles[1],2)), 0}},
                 {{ 0, 0, axes[2]}}}};
    case -12:
        // monoclinic, unique axis y
        if(float_comp(axes[1], 0.)) throw Error("Step::makeCell: axes[1] not provided");
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        if(std::abs(angles[1]) >=1) throw Error("Step::makeCell: angles[1] is invalid");
        return {{{{ 1, 0, 0}},
                 {{ 0, axes[1], 0}},
                 {{ axes[2]*angles[1], 0, axes[2]*sqrt(1.-pow(angles[1],2))}}}};
    case 13:
        // bcm, unique axis z
        if(float_comp(axes[1], 0.)) throw Error("Step::makeCell: axes[1] not provided");
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        if(std::abs(angles[0]) >=1) throw Error("Step::makeCell: angles[0] is invalid");
        return {{{{0.5, 0, -0.5*axes[2]}},
                 {{axes[1]*angles[1], 0, axes[1]*sqrt(1.-pow(angles[1],2))}},
                 {{0.5, 0, 0.5*axes[2]}}}};
    case -13:
        // bcm, unique axis y
        if(float_comp(axes[1], 0.)) throw Error("Step::makeCell: axes[1] not provided");
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        if(std::abs(angles[1]) >=1) throw Error("Step::makeCell: angles[1] is invalid");
        return {{{{0.5, -0.5*axes[1], 0}},
                 {{0.5, 0.5*axes[1], 0}},
                 {{axes[2]*angles[2], 0, axes[2]*sqrt(1.-pow(angles[2],2))}}}};
    case 14:
        // triclinic
        if(float_comp(axes[1], 0.)) throw Error("Step::makeCell: axes[1] not provided");
        if(float_comp(axes[2], 0.)) throw Error("Step::makeCell: axes[2] not provided");
        if(std::abs(angles[0]) >=1) throw Error("Step::makeCell: angles[0] is invalid");
        if(std::abs(angles[1]) >=1) throw Error("Step::makeCell: angles[1] is invalid");
        if(std::abs(angles[2]) >=1) throw Error("Step::makeCell: angles[2] is invalid");
        {
        auto singam = sqrt(1.-pow(angles[2],2));
        return {{{{ 1, 0, 0}},
                 {{axes[1]*angles[2], axes[2]*singam, 0}},
                 {{axes[2]*angles[1], axes[2]*(angles[0]-angles[1]*angles[2])/singam,
                   axes[2]*sqrt(1+2*angles[0]*angles[1]*angles[2]-
                                 pow(angles[0],2)-pow(angles[1],2)-
                                 pow(angles[2],2))/singam}}}};
        }
    default:
        throw Error("Step::makeCell unknown ibrav");
    }
}
