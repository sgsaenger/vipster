#include "ioplugins/pwinput.h"
#define BOOST_SPIRIT_DEBUG
#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/phoenix.hpp>
//#include <boost/fusion/include/adapt_struct_named.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/fusion/adapted/std_array.hpp>
#include <boost/fusion/include/std_pair.hpp>

namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;
using namespace std;
using namespace Vipster;
using namespace qi;
using pwmap = std::map<std::string,std::string>;

namespace boost { namespace spirit { namespace traits {
    template <typename T, size_t N>
        struct is_container<std::array<T, N>, void> : mpl::false_ { };
} } }

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::IO::PWData,
        data,
        mol
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::PseEntry,
        m,
        PWPP
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::Atom,
        name,
        coord,
        fix
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::KPoints,
        active,
        mpg,
        discrete
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::KPoints::MPG,
        x,
        y,
        z,
        sx,
        sy,
        sz
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::KPoints::Discrete,
        properties,
        kpoints
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::DiscreteKPoint,
        pos,
        weight
)

BOOST_FUSION_ADAPT_ADT(
        Vipster::Molecule,
        (obj.getStep(0).getAtoms(), obj.getStep(0).newAtoms(val))
        (obj.getKPoints(),obj.setKPoints(val))
        (obj.getStep(0).getCellVec(), obj.getStep(0).setCellVec(val))
)

template<typename Iterator>
struct pwi_parse_grammar
        : grammar<Iterator, IO::PWData(), blank_type>
{
    pwi_parse_grammar(): pwi_parse_grammar::base_type(file)
    {
        //TODO: ignore comments
        //TODO: change skipper to remove eols
        //TODO: honor AtomFmt/CellFmt
        file = param > *eol > mol;
        file.name("PWScf input file");

        /* Namelists:
         * &NAME
         * key=value, key=value
         * key=value
         * !comment
         * /
         */
        nl_names.add
              ("control", &IO::PWParam::control)
              ("system", &IO::PWParam::system)
              ("electrons", &IO::PWParam::electrons)
              ("ions", &IO::PWParam::ions)
              ("cell", &IO::PWParam::cell)
        ;
        param = ('&' > no_case[nl_names[_a = _1]] > eol >
                 nl[phx::bind([](pwmap IO::PWParam::* v, IO::PWParam& p, pwmap m){
                        p.*v = m;}, _a, _val, _1)]
                 > eol > '/')
                % *eol;
        nl = nl_entry % (eol|',');
        nl_entry = key > qi::lit('=') > value;
        value = fstr | key;
        key = +(qi::graph - qi::char_("=/,'"));
        fstr = qi::char_('\'') >> qi::lexeme[*(qi::char_ - '\'')] > qi::char_('\'');

        /* Cards:
         * NAME {options}
         * #comment
         * !comment
         * <content>
         */
        mol.name("Molecule");
        //TODO: use permutation-op ^ when skipper is solved
        mol = species(_val) > *eol > positions > *eol > kpoints > *eol > cell;
        // Atomic species
        species.name("atomic_species");
        species = no_case["atomic_species"] > eol >
                  (pseentry(_r1) % eol) > eol;
        pseentry.name("PSE Entry");
        pseentry = (key > double_ > key)[phx::bind(
                    [](Molecule& mol, std::string k, float m, std::string p){
                      PseEntry& e = (*mol.pse)[k]; e.m = m; e.PWPP = p;
                     }, _r1, _1, _2, _3)];
        // Atomic positions
        positions.name("atomic_positions");
        positions = no_case["atomic_positions"] > omit[-no_case[atomfmt]] > eol >
                (atom % eol) > eol;
        atom.name("Atom");
        //TODO: parse fix-bools (0=true,1=false)
        atom = key > as<Vec>()[double_ > double_ > double_];
        atomfmt.add
            ("angstrom", AtomFmt::Angstrom)
            ("alat", AtomFmt::Alat)
            ("bohr", AtomFmt::Bohr)
            ("crystal", AtomFmt::Crystal)
        ;
        // K-Points
        kpoints.name("k_points");
        kpoints = no_case["k_points"] > (
                  no_case["gamma"] |
                  (no_case["automatic"] > attr(KPointFmt::MPG) > kmpg) |
                  (attr(KPointFmt::Discrete) > attr(KPoints::MPG{}) > kdisc)
                  ) > eol;
        kmpg.name("MPG");
        kmpg = eol > int_ > int_ > int_ > float_ > float_ > float_;
        kdisc.name("Discrete KPoints");
        kdisc = (no_case[kdiscfmt]|attr(KPoints::Discrete::none)) > eol > omit[int_] > eol > (disckpoint % eol);
        kdiscfmt.add
            ("tpiba", KPoints::Discrete::none)
            ("tpiba_b", KPoints::Discrete::band)
            ("tpiba_c", KPoints::Discrete::contour)
            ("crystal", KPoints::Discrete::crystal)
            ("crystal_b", (KPoints::Discrete::Properties)(
                 KPoints::Discrete::crystal|KPoints::Discrete::band))
            ("crystal_c", (KPoints::Discrete::Properties)(
                 KPoints::Discrete::crystal|KPoints::Discrete::contour))
        ;
        disckpoint.name("Discrete KPoint");
        disckpoint = as<Vec>()[float_ > float_ > float_] > float_;
        // Cell
        cell.name("Cell parameters");
        cell = no_case["cell_parameters"] > -omit[no_case[atomfmt]] > eol >
               as<Vec>()[float_ > float_ > float_] > eol >
               as<Vec>()[float_ > float_ > float_] > eol >
               as<Vec>()[float_ > float_ > float_];

        //DEBUG
        on_error<fail>(
                file,
                cout
                << phx::val("Error! Expecting ")
                << _4
                << phx::val(" here: \"")
                << phx::construct<std::string>(_3, _2)
                << phx::val("\"")
                << endl
        );
    }

    rule<Iterator, IO::PWData(), blank_type> file;
    //PWParam
    rule<Iterator, IO::PWParam(), locals<pwmap IO::PWParam::*>, blank_type> param;
    rule<Iterator, pwmap(), blank_type> nl;
    symbols<char, pwmap IO::PWParam::*> nl_names;
    rule<Iterator, pair<std::string,std::string>(), blank_type> nl_entry;
    rule<Iterator, std::string()> key;
    rule<Iterator, std::string()> value;
    rule<Iterator, std::string()> fstr;
    //Molecule
    rule<Iterator, Molecule(), blank_type> mol;
    //
    rule<Iterator, unused_type(Molecule), blank_type> species;
    rule<Iterator, unused_type(Molecule), locals<PseEntry*>, blank_type> pseentry;
//    rule<Iterator, pair<std::string,PseEntry>(), blank_type> pseentry;
    //
    rule<Iterator, vector<Atom>, blank_type> positions;
    rule<Iterator, Atom(), blank_type> atom;
    symbols<char, AtomFmt> atomfmt;
    //
    rule<Iterator, KPoints(), blank_type> kpoints;
    rule<Iterator, KPoints::MPG(), blank_type> kmpg;
    rule<Iterator, KPoints::Discrete(), blank_type> kdisc;
    symbols<char, KPoints::Discrete::Properties> kdiscfmt;
    rule<Iterator, DiscreteKPoint(), blank_type> disckpoint;
    //
    rule<Iterator, Mat, blank_type> cell;
};

Vipster::IO::BaseData pwi_file_parser(std::string fn, std::ifstream &file)
{
    // set up data
    Vipster::IO::PWData d;
    d.mol.setName(fn);
    d.mol.newStep();

    // set up iterators
    typedef std::istreambuf_iterator<char> iter;
    boost::spirit::multi_pass<iter> first = boost::spirit::make_default_multi_pass(iter(file));
    boost::spirit::multi_pass<iter> last = boost::spirit::make_default_multi_pass(iter());

    pwi_parse_grammar<boost::spirit::multi_pass<iter>> grammar;

    phrase_parse(first, last, grammar, blank, d);

    return d;
}

const Vipster::IOPlugin Vipster::IO::PWInput =
{
    "pwi",
    "pwi",
    "pwi",
    &pwi_file_parser,
    nullptr
};
