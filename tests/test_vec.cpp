#include "vec.h"
#include "configfile.h"
#include <ostream>

using namespace Vipster;

namespace std{
static ostream& operator<<(ostream& out, const Vec& v)
{
    out << "Vec{" << v[0] << ", " << v[1] << ", " << v[2] << "}";
    return out;
}

static ostream& operator<<(ostream& out, const Mat& m)
{
    out << "Mat{" << m[0] << ", " << m[1] << ", " << m[2] << "}";
    return out;
}
}

#include "catch.hpp"

TEST_CASE( "Vipster::Vec operators", "[vec]" ) {
    Vec v1{{1, 1, 1}};
    Vec v2{{1.5, 1.5, 1.5}};

    REQUIRE(v1 == Vec{{1,1,1}});
    REQUIRE(v1 != v2);

    SECTION( "Addition and subtraction" ) {
        Vec v3{{2.5, 2.5, 2.5}};
        REQUIRE( v1+0.5 == v2 );
        REQUIRE( 0.5+v1 == v2 );
        REQUIRE( v2-0.5 == v1 );
        v1 += 0.5;
        INFO( "In-place vec-float addition" );
        REQUIRE( v1 == v2 );
        v1 -= 0.5;
        INFO( "In-place vec-float subtraction" );
        REQUIRE( v1+0.5 == v2 );
        REQUIRE( v1+v2 == v3 );
        REQUIRE( v3-v2 == v1 );
        v1 += v2;
        INFO( "In-place vec-vec addition" );
        REQUIRE( v1 == v3 );
        v1 -= v2;
        INFO( "In-place vec-vec subtraction" );
        REQUIRE( v1+v2 == v3 );
    }

    SECTION( "Multiplication and division" ) {
        REQUIRE( v1*1.5 == v2 );
        REQUIRE( 1.5*v1 == v2 );
        REQUIRE( v2/1.5 == v1 );
        v1 *= 1.5;
        INFO( "In-place multiplication" );
        REQUIRE( v1 == v2 );
        v1 /= 1.5;
        INFO( "In-place division" );
        REQUIRE( v1*1.5 == v2 );
    }

    SECTION( "Linear algebra ops" ) {
        REQUIRE( float_comp(Vec_dot(v1,v2), 4.5f) );
        REQUIRE( float_comp(Vec_length(v1), sqrtf(3)) );
        REQUIRE( Vec_cross(v1, v2) == Vec{{0,0,0}} );
        Vec v3{{1,0,0}};
        REQUIRE( Vec_cross(v1, v3) == Vec{{0,1,-1}} );
        REQUIRE( Vec_cross(v1, v3) == -Vec_cross(v3,v1) );
    }
}

TEST_CASE( "Vipster::Mat operators", "[mat]" ) {
    Mat m1{{{{1,0,0}}, {{0,1,0}}, {{0,0,1}}}};

    SECTION( "Multiplication and divison" ) {
        Mat m2{{{{2,0,0}}, {{0,2,0}}, {{0,0,2}}}};
        REQUIRE( m1*2 == m2 );
        REQUIRE( m2/2 == m1 );
        m1 *= 2;
        REQUIRE( m1 == m2 );
        m1 /= 2;
        REQUIRE( m1*2 == m2 );
    }

    SECTION( "Linear algebra" ) {
        Vec v1{{1,2,3}};
        Mat m2{{{{1,2,3}}, {{0,1,0}}, {{0,0,1}}}};
        REQUIRE( Mat_trans(m2) == Mat{{{{1,0,0}},{{2,1,0}},{{3,0,1}}}} );
        REQUIRE( v1*m1 == v1 );
        REQUIRE( m1*v1 == v1 );
        REQUIRE( m2*v1 == Vec{{14,2,3}} );
        REQUIRE( v1*m2 == Vec{{1,4,6}} );
        REQUIRE( m1*m2 == m2 );
        REQUIRE( m2*m2 == Mat{{{{1,4,6}},{{0,1,0}},{{0,0,1}}}} );
        REQUIRE( Mat_det(m2) == 1.f );
        REQUIRE( Mat_inv(m2) == Mat{{{{1,-2,-3}},{{0,1,0}},{{0,0,1}}}} );
        REQUIRE( m2*Mat_inv(m2) == m1 );
        m2*=Mat_inv(m2);
        REQUIRE( m2 == m1 );
        Mat m3{{{{1,2,3}}, {{4,5,6}}, {{7,8,9}}}};
        REQUIRE( Mat_det(m3) == 0.f );
        REQUIRE_THROWS( Mat_inv(m3) );
    }
}
