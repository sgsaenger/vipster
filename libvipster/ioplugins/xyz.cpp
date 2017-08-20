#include "ioplugins/xyz.h"
#include <iomanip>
#define BOOST_SPIRIT_DEBUG
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/fusion/adapted/std_array.hpp>


namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;
using namespace Vipster;

namespace boost { namespace spirit { namespace traits {
    template <typename T, size_t N>
        struct is_container<std::array<T, N>, void> : mpl::false_ { };
} } }

BOOST_FUSION_ADAPT_STRUCT(
        Atom,
        name,
        coord
)

BOOST_FUSION_ADAPT_ADT(
        Step,
        (obj.getComment(), obj.setComment(val))
        (obj.getAtoms(), obj.newAtoms(val))
)

BOOST_FUSION_ADAPT_ADT(
        Molecule,
        (obj.getSteps(), obj.newSteps(val))
)

std::ostream& operator<< (std::ostream& stream, const Step& s)
{
    stream << "Step with " << s.getNat() << " atoms." << std::endl;
    return stream;
}

template<typename Iterator>
struct xyz_parse_grammar
        : qi::grammar<Iterator, Molecule(), qi::blank_type>
{
    xyz_parse_grammar(): xyz_parse_grammar::base_type(mol, "XYZ")
    {
        name = +(qi::char_ - qi::space);
        name.name("Element");
        atom = name >> qi::as<Vec>()[qi::float_ > qi::float_ > qi::float_];
        atom.name("Atom");
        atoms = (atom % qi::eol)
                > qi::eps(qi::_r1 == phx::bind(&std::vector<Atom>::size,qi::_val));
        atoms.name("Atoms");
        comment = *(qi::char_ - qi::eol) > qi::eol;
        comment.name("Comment");
        step %= qi::omit[qi::uint_[qi::_a = qi::_1] > qi::eol]
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
    qi::rule<Iterator, Molecule(), qi::blank_type> mol;
    qi::rule<Iterator, std::vector<Step>(), qi::blank_type> steps;
    qi::rule<Iterator, Step(), qi::locals<unsigned int>, qi::blank_type> step;
    qi::rule<Iterator, std::string()> comment;
    qi::rule<Iterator, std::vector<Atom>(unsigned int), qi::blank_type> atoms;
    qi::rule<Iterator, Atom(), qi::blank_type> atom;
    qi::rule<Iterator, std::string()> name;
};

std::shared_ptr<IO::BaseData> xyz_file_parser(std::string name, std::ifstream &file)
{
    auto d = std::make_shared<IO::BaseData>();
    d->mol.setName(name);

    typedef std::istreambuf_iterator<char> iter;
    boost::spirit::multi_pass<iter> first = boost::spirit::make_default_multi_pass(iter(file));

    xyz_parse_grammar<boost::spirit::multi_pass<iter>> grammar;

    qi::phrase_parse(first, boost::spirit::make_default_multi_pass(iter()), grammar, qi::blank, d->mol);
    d->mol.setFmtAll(AtomFmt::Angstrom);

    return d;
}

bool xyz_file_writer(const Molecule& m, std::ofstream &file, const IO::BaseParam*)
{
    const Step& s = m.getStep(0);
    file << s.getNat() << '\n';
    file << s.getComment() << '\n';
    file << std::fixed << std::setprecision(5);
    for(const Atom& at: s.getAtoms()){
        file << std::left << std::setw(3) << at.name << " "
             << std::right << std::setw(10) << at.coord[0] << " "
             << std::right << std::setw(10) << at.coord[1] << " "
             << std::right << std::setw(10) << at.coord[2] << '\n';
    }
    return true;
}

const IOPlugin IO::XYZ =
{
    "xyz",
    "xyz",
    "xyz",
    &xyz_file_parser,
    &xyz_file_writer
};
