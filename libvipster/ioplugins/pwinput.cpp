#include "ioplugins/pwinput.h"
#define BOOST_SPIRIT_DEBUG
#define BOOST_SPIRIT_USE_PHOENIX_V3
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/fusion/include/adapt_struct_named.hpp>
#include <boost/fusion/include/adapt_adt.hpp>
#include <boost/fusion/include/std_pair.hpp>

namespace qi = boost::spirit::qi;
namespace phx = boost::phoenix;
using namespace std;
using namespace Vipster;

std::ostream& operator<< (ostream& stream, const map<string, string>m)
{
    for(auto& p:m){
        stream << p.first << " has value " << p.second << endl;
    }
    return stream;
}

std::ostream& operator<< (ostream& stream, const Vipster::IO::PWParam& p)
{
    stream << "control" << endl << p.control << endl;
    stream << "system" << endl << p.system << endl;
    stream << "electrons" << endl << p.electrons << endl;
    stream << "ions" << endl << p.ions << endl;
    stream << "cell" << endl << p.cell << endl;
    return stream;
}

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::Atom,
        (std::string, name)
        (float, coord[0])
        (float, coord[1])
        (float, coord[2])
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::DiscreteKPoint,
        (float, pos[0])
        (float, pos[1])
        (float, pos[2])
        (float, weight)
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
        (Vipster::DiscreteKPoint::Properties, properties)
        (std::vector<Vipster::DiscreteKPoint>, kpoints)
)


BOOST_FUSION_ADAPT_STRUCT(
        Vipster::KPoints,
        active,
        mpg,
        discrete
)

BOOST_FUSION_ADAPT_ADT(
        Vipster::Step,
        (std::vector<Vipster::Atom>, std::vector<Vipster::Atom>, obj.getAtoms(), obj.newAtoms(val))
)

BOOST_FUSION_ADAPT_ADT(
        Vipster::Molecule,
        (obj.getSteps(), obj.newSteps(val))
        (obj.getKPoints(), obj.setKPoints(val))
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::PseEntry,
        m,
        PWPP
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::IO::PWParam,
        control,
        system,
        electrons,
        ions,
        cell
)

BOOST_FUSION_ADAPT_STRUCT(
        Vipster::IO::PWData,
        data,
        mol
)

template<typename Iterator>
struct pwi_parse_grammar
        : qi::grammar<Iterator, IO::PWData(), qi::locals<int,int,int>, qi::blank_type>
{
    pwi_parse_grammar(): pwi_parse_grammar::base_type(file)
    {
        key = *(qi::graph - qi::char_("=/,'"));
        key.name("Key (or unquoted value)");
        fstr = qi::char_('\'') >> qi::lexeme[*(qi::char_ - '\'')] > qi::char_('\'');
        fstr.name("Quoted value");
        value = fstr | key;
        value.name("Value");
        nl_entry = key >> qi::lit('=') >> value;
        nl_entry.name("Namelist KV-pair");
        namelist = qi::omit[qi::no_case[qi::string(qi::_r1)]] >> *qi::eol
                   >> nl_entry % (*(qi::eol >> -(qi::lit('!') >> *(qi::char_ - qi::eol) > qi::eol)) | ',')
                   >> *qi::eol >> '/';
        namelist.name("Namelist");
        param = namelist(string{"&control"}) > *qi::eol >
                namelist(string{"&system"}) > *qi::eol >
                namelist(string{"&electrons"}) > *qi::eol >
                -(namelist(string{"&ions"}) > *qi::eol) >
                -(namelist(string{"&cell"}));
        param.name("PW parameter set");
        //TODO: hier weitermachen
        //IDEE: semantic action nach param -> nat,ntyp,ibrav an mol uebergeben?!
        //IDEE2: ^-operator um reihenfolge beliebig zu machen?
        name = +(qi::char_ - qi::space);
        pseentry = qi::double_ > name;
        pseentry.name("PSEEntry");
        psepair = name > pseentry;
        psepair.name("PSE KV-pair");
        species = qi::no_case[qi::lit("atomic_species")] > qi::eol
                > (psepair % qi::eol)
                > qi::eps(qi::_r1 == phx::bind(&vector<pair<string,PseEntry>>::size,qi::_val));
        species.name("Atomic species");
        atom = name > qi::double_ > qi::double_ > qi::double_;
        atom.name("Atom");
        format = qi::no_case[qi::lit("crystal")][phx::ref(fmt) = AtomFmt::Crystal]
               | qi::no_case[qi::lit("alat")][phx::ref(fmt) = AtomFmt::Alat]
               | qi::no_case[qi::lit("bohr")][phx::ref(fmt) = AtomFmt::Bohr]
               | qi::no_case[qi::lit("angstrom")][phx::ref(fmt) = AtomFmt::Angstrom];
        format.name("Atom format");
        atoms = qi::no_case[qi::lit("atomic_positions")] > -format > qi::eol
                > (atom % qi::eol)
                > qi::eps(qi::_r1 == phx::bind(&vector<Atom>::size,qi::_val));
        atoms.name("Atomic positions");
        positions = qi::as<Step>()[atoms(qi::_r1)];
        kpointmpg = qi::eol > qi::int_ > qi::int_ > qi::int_ > qi::double_ > qi::double_ > qi::double_;
        disckpoint = qi::double_ > qi::double_ > qi::double_ > qi::double_;
        disckpoints = (disckpoint % qi::eol);
        kpointdisc = ((qi::no_case[qi::lit("crystal_b")] > qi::attr(DiscreteKPoint::crystal|DiscreteKPoint::band))
                     |(qi::no_case[qi::lit("crystal")] > qi::attr(DiscreteKPoint::crystal))
                     |(qi::no_case[qi::lit("tpiba_b")] > qi::attr(DiscreteKPoint::band))
                     |(-qi::no_case[qi::lit("tpiba")] > qi::attr(0))
                     ) > qi::eol > qi::omit[qi::int_] > qi::eol > (disckpoint % qi::eol);
        kpoints = qi::no_case[qi::lit("k_points")] >
                   (qi::no_case[qi::lit("gamma")]
                   |(qi::no_case[qi::lit("automatic")] > qi::attr(KPointFmt::MPG) > kpointmpg)
                   |(qi::attr(KPointFmt::Discrete) > qi::attr(KPoints::MPG{}) > kpointdisc));
        mol.name("Molecule");
        mol = qi::omit[species(qi::_r1)] > *qi::eol >
              positions(qi::_r2) > *qi::eol >
              kpoints > *qi::eol;// >
//              cell(qi::_r3) > *qi::eol >
//              occupations > *qi::eol >
//              constraints > *qi::eol >
//              forces > *qi::eol;
        file %= param >
                qi::eps[qi::_a = phx::bind([](const IO::PWData& p){return stoi(p.data.system.at("ntyp"));}, qi::_val)] >
                qi::eps[qi::_b = phx::bind([](const IO::PWData& p){return stoi(p.data.system.at("nat"));}, qi::_val)] >
                qi::eps[qi::_c = phx::bind([](const IO::PWData& p){return stoi(p.data.system.at("ibrav"));}, qi::_val)] >
                *qi::eol > mol(qi::_a,qi::_b,qi::_c);
        file.name("PWScf input file");
        qi::on_error<qi::fail>(
                file,
                std::cout
                << phx::val("Error! Expecting ")
                << qi::_4
                << phx::val(" here: \"")
                << phx::construct<std::string>(qi::_3, qi::_2)
                << phx::val("\"")
                << std::endl
        );
//        BOOST_SPIRIT_DEBUG_NODE(file);
//        BOOST_SPIRIT_DEBUG_NODE(param);
//        BOOST_SPIRIT_DEBUG_NODE(namelist);
//        BOOST_SPIRIT_DEBUG_NODE(nl_entry);
//        BOOST_SPIRIT_DEBUG_NODE(value);
//        BOOST_SPIRIT_DEBUG_NODE(key);
//        BOOST_SPIRIT_DEBUG_NODE(fstr);
//        BOOST_SPIRIT_DEBUG_NODE(mol);
//        BOOST_SPIRIT_DEBUG_NODE(name);
//        BOOST_SPIRIT_DEBUG_NODE(species);
//        BOOST_SPIRIT_DEBUG_NODE(psepair);
//        BOOST_SPIRIT_DEBUG_NODE(pseentry);
//        BOOST_SPIRIT_DEBUG_NODE(atoms);
//        BOOST_SPIRIT_DEBUG_NODE(atom);
//        BOOST_SPIRIT_DEBUG_NODE(format);
    }
    AtomFmt fmt = AtomFmt::Alat;
    Mat m={{{1,0,0},{0,1,0},{0,0,1}}};
    qi::rule<Iterator, IO::PWData(), qi::locals<int,int,int>, qi::blank_type> file;
//Parameter-set, composed of F90-namelists
    qi::rule<Iterator, IO::PWParam(), qi::blank_type> param;
    qi::rule<Iterator, map<string, string>(string), qi::blank_type> namelist;
    qi::rule<Iterator, pair<string, string>(), qi::blank_type> nl_entry;
    qi::rule<Iterator, string()> value;
    qi::rule<Iterator, string()> key;
    qi::rule<Iterator, string()> fstr;
//Molecule
    qi::rule<Iterator, Molecule(int,int,int), qi::blank_type> mol;
    qi::rule<Iterator, string()> name;
    qi::rule<Iterator, vector<pair<string,PseEntry>>(int), qi::blank_type> species;
    qi::rule<Iterator, pair<string,PseEntry>(), qi::blank_type> psepair;
    qi::rule<Iterator, PseEntry(), qi::blank_type> pseentry;
    qi::rule<Iterator, vector<Step>(int), qi::blank_type> positions;
    qi::rule<Iterator, vector<Atom>(int), qi::blank_type> atoms;
    qi::rule<Iterator, qi::blank_type> format;
    qi::rule<Iterator, Atom(), qi::blank_type> atom;
    qi::rule<Iterator, KPoints(), qi::blank_type> kpoints;
    qi::rule<Iterator, KPoints::MPG(), qi::blank_type> kpointmpg;
    qi::rule<Iterator, KPoints::Discrete(), qi::blank_type> kpointdisc;
    qi::rule<Iterator, DiscreteKPoint(), qi::blank_type> disckpoint;
//    qi::rule<Iterator, qi::unused_type(int), qi::blank_type> cell;
//    qi::rule<Iterator, qi::unused_type(), qi::blank_type> occupations;
//    qi::rule<Iterator, qi::unused_type(), qi::blank_type> constraints;
//    qi::rule<Iterator, qi::unused_type(), qi::blank_type> forces;
};

Vipster::IO::BaseData pwi_file_parser(std::string fn, std::ifstream &file)
{
    Vipster::IO::PWData d;
    d.mol.setName(fn);
//    Vipster::IO::PWParam* p = dynamic_cast<Vipster::IO::PWParam*>(&d.data);

    typedef std::istreambuf_iterator<char> iter;
    boost::spirit::multi_pass<iter> first = boost::spirit::make_default_multi_pass(iter(file));

    pwi_parse_grammar<boost::spirit::multi_pass<iter>> grammar;

    qi::phrase_parse(first, boost::spirit::make_default_multi_pass(iter()), grammar, qi::blank, d);

//    cout << *p << endl;
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
