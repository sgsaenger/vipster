#include "ioplugins/xyz.h"
#define BOOST_SPIRIT_DEBUG
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_adt.hpp>


namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;

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

std::ostream& operator<< (std::ostream& stream, const Vipster::Step& s)
{
    stream << "Step with " << s.getNat() << " atoms." << std::endl;
    return stream;
}

template<typename Iterator>
struct xyz_parse_grammar
        : qi::grammar<Iterator, Vipster::Molecule(), qi::blank_type>
{
    xyz_parse_grammar(): xyz_parse_grammar::base_type(mol, "XYZ")
    {
        mol = steps > qi::eoi;
        mol.name("Molecule");
        steps = step % *qi::eol;
        steps.name("Steps");
        step %= qi::omit[qi::int_[qi::_a = qi::_1] > qi::eol]
                > comment
                > atoms(qi::_a)
                > qi::eps[phx::bind(&Vipster::Step::setCellDim,qi::_val,1,true,Vipster::AtomFmt::Angstrom)];
        step.name("Step");
        comment = *(qi::char_ - qi::eol);
        comment.name("Comment");
        atoms = qi::repeat(qi::_r1)[qi::eol > atom];
        atoms.name("Atoms");
        atom = name
               > qi::double_
               > qi::double_
               > qi::double_;
        atom.name("Atom");
        name = +qi::alnum;
        name.name("Element");
        qi::on_error<qi::fail>(
                mol,
                std::cout
                << phx::val("Error! Expecting ")
                << qi::_4
                << phx::val(" here: \"")
                << phx::construct<std::string>(qi::_3, qi::_2)
                << phx::val("\"")
                << std::endl
        );
//        BOOST_SPIRIT_DEBUG_NODE(mol);
//        BOOST_SPIRIT_DEBUG_NODE(steps);
//        BOOST_SPIRIT_DEBUG_NODE(step);
        BOOST_SPIRIT_DEBUG_NODE(atoms);
        BOOST_SPIRIT_DEBUG_NODE(atom);
        BOOST_SPIRIT_DEBUG_NODE(comment);
        BOOST_SPIRIT_DEBUG_NODE(name);
    }
    qi::rule<Iterator, Vipster::Molecule(), qi::blank_type> mol;
    qi::rule<Iterator, std::vector<Vipster::Step>(), qi::blank_type> steps;
    qi::rule<Iterator, Vipster::Step(), qi::locals<int>, qi::blank_type> step;
    qi::rule<Iterator, std::string()> comment;
    qi::rule<Iterator, std::vector<Vipster::Atom>(int), qi::blank_type> atoms;
    qi::rule<Iterator, Vipster::Atom(), qi::blank_type> atom;
    qi::rule<Iterator, std::string()> name;
};

Vipster::IOData xyz_file_parser(std::string fn, std::ifstream &file)
{
    Vipster::IOData d{
        Vipster::Molecule{fn,0},
        std::shared_ptr<Vipster::IO::BaseParam>()
    };

    typedef std::istreambuf_iterator<char> iter;
    boost::spirit::multi_pass<iter> first = boost::spirit::make_default_multi_pass(iter(file));

    xyz_parse_grammar<boost::spirit::multi_pass<iter>> grammar;

    qi::phrase_parse(first, boost::spirit::make_default_multi_pass(iter()), grammar, qi::blank, d.mol);

    return d;
}

const Vipster::IOPlugin Vipster::IO::XYZ =
{
    "xyz",
    "xyz",
    "xyz",
    &xyz_file_parser,
    nullptr
};
