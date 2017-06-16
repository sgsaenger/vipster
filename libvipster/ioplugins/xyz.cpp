#include "ioplugins/xyz.h"
#define BOOST_SPIRIT_DEBUG
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/fusion/adapted/std_array.hpp>


namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;

namespace boost { namespace spirit { namespace traits {
    template <typename T, size_t N>
        struct is_container<std::array<T, N>, void> : mpl::false_ { };
} } }


BOOST_FUSION_ADAPT_STRUCT(
        Vipster::Atom,
        name,
        coord
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
        name = +(qi::char_ - qi::space);
        name.name("Element");
        atom = name > qi::as<Vipster::Vec>()[qi::float_ > qi::float_ > qi::float_];
        atom.name("Atom");
        atoms = (atom % qi::eol)
                > qi::eps(qi::_r1 == phx::bind(&std::vector<Vipster::Atom>::size,qi::_val));
        atoms.name("Atoms");
        comment = *(qi::char_ - qi::eol) > qi::eol;
        comment.name("Comment");
        step %= qi::omit[qi::int_[qi::_a = qi::_1] > qi::eol]
                > comment
                > atoms(qi::_a);
        step.name("Step");
        steps = step % +qi::eol;
        steps.name("Steps");
        mol = steps;
        mol.name("Molecule");
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
    }
    qi::rule<Iterator, Vipster::Molecule(), qi::blank_type> mol;
    qi::rule<Iterator, std::vector<Vipster::Step>(), qi::blank_type> steps;
    qi::rule<Iterator, Vipster::Step(), qi::locals<int>, qi::blank_type> step;
    qi::rule<Iterator, std::string()> comment;
    qi::rule<Iterator, std::vector<Vipster::Atom>(int), qi::blank_type> atoms;
    qi::rule<Iterator, Vipster::Atom(), qi::blank_type> atom;
    qi::rule<Iterator, std::string()> name;
};

Vipster::IO::BaseData xyz_file_parser(std::string fn, std::ifstream &file)
{
    Vipster::IO::BaseData d;
    d.mol.setName(fn);

    typedef std::istreambuf_iterator<char> iter;
    boost::spirit::multi_pass<iter> first = boost::spirit::make_default_multi_pass(iter(file));

    xyz_parse_grammar<boost::spirit::multi_pass<iter>> grammar;

    qi::phrase_parse(first, boost::spirit::make_default_multi_pass(iter()), grammar, qi::blank, d.mol);
    d.mol.setFmtAll(Vipster::AtomFmt::Angstrom);

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
