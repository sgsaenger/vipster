#include <iostream>
#include <fstream>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <molecule.h>

namespace spirit = boost::spirit;
namespace qi = boost::spirit::qi;
using namespace qi;
namespace phx = boost::phoenix;

namespace std {
    ostream& operator<<(ostream& stream, const Vipster::Atom& at)
    {
        stream << "Atom: " << at.name << endl;
        stream << "at position: (" << at.coord[0] << ", " << at.coord[1] << ", " << at.coord[2] << ")" << endl;
        stream << "Charge: " << at.charge;
        stream << ", fixed: (" << at.fix[0] << ", " << at.fix[1] << ", " << at.fix[2];
        stream << "), hidden: " << at.hidden << endl;
        return stream;
    }
    ostream& operator<<(ostream& stream, const Vipster::Step& step)
    {
        stream << "Step with " << step.getNat() << " atoms" << endl;
        stream << "Comment: " << step.getComment() << endl;
        stream << "Atoms:" << endl;
        for(auto &a: step.getAtoms())
        {
            stream << a;
        }
        return stream;
    }
    ostream& operator<<(ostream& stream, const Vipster::Molecule& mol)
    {
        stream << "Molecule with " << mol.getNstep() << " steps" << endl;
        stream << "Name: " << mol.getName() << endl;
        stream << "Steps: " << endl;
        for(auto &s: mol.getSteps())
        {
            stream << s;
        }
        return stream;
    }
}
BOOST_FUSION_ADAPT_STRUCT(
        Vipster::Atom,
        (std::string, name)
        (float, coord[0])
        (float, coord[1])
        (float, coord[2])
)
BOOST_FUSION_ADAPT_ADT(
        Vipster::Step,
        (std::string, std::string, obj.getComment(), obj.setComment(val))
        (std::vector<Vipster::Atom>, std::vector<Vipster::Atom>, obj.getAtoms(), obj.newAtoms(val))
)

BOOST_FUSION_ADAPT_ADT(
        Vipster::Molecule,
        (std::vector<Vipster::Step>, std::vector<Vipster::Step>, obj.getSteps(), obj.newSteps(val))
)

template<typename Iterator>
struct xyz_parse_grammar
        : qi::grammar<Iterator, Vipster::Molecule(), blank_type>
{
    xyz_parse_grammar(): xyz_parse_grammar::base_type(mol)
    {
        mol = steps;
        steps = step % *eol;
        step %= omit[int_ >> eol]
                >> comment
                >> atoms
                >> eps[phx::bind(&Vipster::Step::setCellDim,_val,1,true,Vipster::AtomFmt::Angstrom)];
        comment = *(char_ - eol) >> eol;
        atoms = atom % eol;
        atom = name
               >> double_
               >> double_
               >> double_;
        name = +alnum;
    }
    rule<Iterator, Vipster::Molecule(), blank_type> mol;
    rule<Iterator, std::vector<Vipster::Step>(), blank_type> steps;
    rule<Iterator, Vipster::Step(), blank_type> step;
    rule<Iterator, std::string()> comment;
    rule<Iterator, std::vector<Vipster::Atom>(), blank_type> atoms;
    rule<Iterator, Vipster::Atom(),blank_type> atom;
    rule<Iterator, std::string()> name;
};

int main()
{
    Vipster::Molecule m{"AHA",0};
    std::string s = "4\n"
                    "ba rgl\n"
                    "C 0 0 0\n"
                    "H 1 2 3\n"
                    "O 1.2 3.4 5.6\n"
                    "C1 17  213 1.23\n"
                    "2\n"
                    "wurgls\n"
                    "C 1 2 3\n"
                    "F 3 5 6\n";
    std::ifstream f{"test.xyz"};
    typedef std::istreambuf_iterator<char> iter;
    spirit::multi_pass<iter> first = spirit::make_default_multi_pass(iter(f));
    xyz_parse_grammar<spirit::multi_pass<iter>> grammar;

    qi::phrase_parse(first, spirit::make_default_multi_pass(iter()), grammar, qi::blank, m);

    std::cout << m << std::endl;
    return 0;
}
