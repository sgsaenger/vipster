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

template<typename Iterator>
struct pwi_parse_grammar
        : grammar<Iterator, IO::PWData(), blank_type>
{
    pwi_parse_grammar(): pwi_parse_grammar::base_type(file)
    {
        //TODO: ignore comments
        //TODO: change skipper to remove eols
        //TODO: honor CellFmt
        file = param > *eol > mol;

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
        param       = ('&' > no_case[nl_names[_a = _1]] > eol
                      > nl[phx::bind(
                        [](pwmap IO::PWParam::* v, IO::PWParam& p, pwmap m)
                        { p.*v = m; }, _a, _val, _1)]
                      > eol > '/') % *eol;
        nl          = nl_entry % (eol|',');
        nl_entry    = key > qi::lit('=') > value;
        value       = fstr | key;
        key         = +(qi::graph - qi::char_("=/,'"));
        fstr        = qi::char_('\'')
                      >> qi::lexeme[*(qi::char_ - '\'')]
                      > qi::char_('\'');

        /* Cards:
         * NAME {options}
         * #comment
         * !comment
         * <content>
         */
        mol         = species(phx::ref(_val))
                      ^ positions(phx::ref(_val))
                      ^ kpoints(phx::ref(_val))
                      ^ cell(phx::ref(_val));
        // Atomic species
        species     = no_case["atomic_species"] > eol
                      > (pseentry(phx::ref(_r1)) % eol) > *eol;
        pseentry    = (key >> float_ >> key)
                      [phx::bind(
                       [](Molecule& mol, std::string k, float m, std::string p){
                         PseEntry& e = (*mol.pse)[k]; e.m = m; e.PWPP = p;
                        }, _r1, _1, _2, _3)];
        // Atomic positions
        //TODO: parse fix-bools (0=true,1=false)
        positions   = (no_case["atomic_positions"] >
                      no_case[atomfmt] > eol >
                      (atom % eol) > *eol)
                      [phx::bind(
                       [](Molecule& mol, AtomFmt f, vector<Atom> at){
                          Step& s = mol.getStep(0); s.setFmt(f); s.newAtoms(at);
                       }, _r1, _1, _2)];
        atom        = key >> as<Vec>()[float_ > float_ > float_];
        atomfmt.add
            ("angstrom", AtomFmt::Angstrom)
            ("alat", AtomFmt::Alat)
            ("bohr", AtomFmt::Bohr)
            ("crystal", AtomFmt::Crystal)
        ;
        // K-Points
        kpoints     = no_case["k_points"] > (
                      no_case["gamma"] |
                      (no_case["automatic"] > attr(KPointFmt::MPG) > kmpg) |
                      (attr(KPointFmt::Discrete) > attr(KPoints::MPG{}) > kdisc)
                      ) > *eol;
        kmpg        = eol > int_ > int_ > int_ > float_ > float_ > float_;
        kdisc       = (no_case[kdiscfmt]|attr(KPoints::Discrete::none))
                      > eol > omit[int_]
                      > eol > (disckpoint % eol);
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
        disckpoint  = as<Vec>()[float_ > float_ > float_] > float_;
        // Cell
        cell        = as<Mat>()[no_case["cell_parameters"]
                      > -omit[no_case[atomfmt]] > eol
                      > as<Vec>()[float_ > float_ > float_] > eol
                      > as<Vec>()[float_ > float_ > float_] > eol
                      > as<Vec>()[float_ > float_ > float_] > *eol]
                      [phx::bind([](Molecule& m, Mat v){
                        m.getStep(0).setCellVec(v);
                        },_r1,_1)];

        //DEBUG
        file.name("PWScf input file");
        mol.name("Molecule");
        species.name("atomic_species");
        pseentry.name("PSE Entry");
        positions.name("atomic_positions");
        atom.name("Atom");
        kpoints.name("k_points");
        kmpg.name("MPG");
        kdisc.name("Discrete KPoints");
        disckpoint.name("Discrete KPoint");
        cell.name("Cell parameters");
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
        BOOST_SPIRIT_DEBUG_NODES((species)(pseentry)(positions));
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
    rule<Iterator, unused_type(Molecule&), blank_type> species;
    rule<Iterator, unused_type(Molecule&), locals<PseEntry*>, blank_type> pseentry;
    //
    rule<Iterator, unused_type(Molecule&), blank_type> positions;
    rule<Iterator, Atom(), blank_type> atom;
    symbols<char, AtomFmt> atomfmt;
    //
    rule<Iterator, unused_type(Molecule&), blank_type> kpoints;
    rule<Iterator, KPoints::MPG(), blank_type> kmpg;
    rule<Iterator, KPoints::Discrete(), blank_type> kdisc;
    symbols<char, KPoints::Discrete::Properties> kdiscfmt;
    rule<Iterator, DiscreteKPoint(), blank_type> disckpoint;
    //
    rule<Iterator, unused_type(Molecule&), blank_type> cell;
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
    auto cdm = d.data.system.find("celldm(1)");
    if (cdm != d.data.system.end()){
        d.mol.getStep(0).setCellDim(stof(cdm->second));
        d.data.system.erase(cdm);
    }
    d.data.system.erase("nat");
    d.data.system.erase("ntyp");

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
